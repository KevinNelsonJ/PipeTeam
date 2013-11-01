#include <fstream>
#include <iostream>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>

namespace Poission
{
    using namespace dealii;

    template<int dim> class PoissonEqn {
        public:
        PoissonEqn();
        void run();
        
        private:
        void make_grid();
        void setup_system();
        void assemble_system();
        void solve();
        void output_results() const;
        
        Triangulation<dim> triangulation;
        FE_Q<dim>          fe;
        DoFHandler<dim>    dof_handler;
        
        SparsityPattern      sparsity_pattern;
        SparseMatrix<double> system_matrix;
        Vector<double>       solution;
        Vector<double>       system_rhs;
    };
    
    template<int dim> PoissonEqn<dim>::PoissonEqn() :
        fe(1),                     // use piecwise polynomial elements of order 1
        dof_handler(triangulation) // knows dimensionality because of triangulation 
        {}
        
    template<int dim> void PoissonEqn<dim>::make_grid() {

        // makes a simple hypercube grid to solve the system on
        GridGenerator::hyper_cube(triangulation, -1, 1);
        triangulation.refine_global(5);
        
        std::cout << "Finished make_grid()"
                  << std::endl;
        std::cout << "Number of active cells: "
                  << triangulation.n_active_cells()
                  << std::endl;
        std::cout << "Total number of cells: "
                  << triangulation.n_cells()
                  << std::endl
                  << std::endl;
    }
    
    template<int dim> void PoissonEqn<dim>::setup_system() {
        // determine number of dofs necessary based on mesh and order of test functions
        // and then number them
        dof_handler.distribute_dofs(fe);
        
        // generate our sparsity pattern
        CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
        DoFTools::make_sparsity_pattern(dof_handler, c_sparsity);
        sparsity_pattern.copy_from(c_sparsity);
    
        // reinitialize matricies and vectors with the new sparsity pattern
        system_matrix.reinit(sparsity_pattern);
        solution.reinit(dof_handler.n_dofs());
        system_rhs.reinit(dof_handler.n_dofs());  
        
        std::cout << "Finished setup_system()"
                  << std::endl;
        std::cout << "Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl
                  << std::endl; 
    }
    
    template<int dim> void PoissonEqn<dim>::assemble_system() {

        // use quadrature scheme with 2 quadrature points
        // FEValues can provide the quadrature points for the scheme on the correct test functions
        QGauss<dim>   quadrature_formula(2); 
        FEValues<dim> fe_values(fe, quadrature_formula, 
                                update_values | update_gradients | update_JxW_values);
        
        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        const unsigned int n_q_points = quadrature_formula.size();
        
        // We're going to compute each cell's contribution to the system independantly, then add
        // back to the overall system matrix.  Those contributions are stored here.
        FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double>     cell_rhs(dofs_per_cell);
        
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                       endc = dof_handler.end();
        // Iterate over the active cells                                    
        for(; cell!=endc; ++cell) {
        
            // Reinitialize cell-independant data back to 0
            fe_values.reinit(cell);
            cell_matrix = 0;
            cell_rhs = 0;
            
            // compute the dot product of the gradient of the trial and test functions on this cell
            // this is the mass matrix
            for(unsigned int i=0; i<dofs_per_cell; ++i) {
                for(unsigned int j=0; j<dofs_per_cell; ++j) {
                    for(unsigned int q_point=0; q_point<n_q_points; ++q_point) {
                        cell_matrix(i,j) += (fe_values.shape_grad(i, q_point) * 
                                             fe_values.shape_grad(j, q_point) *
                                             fe_values.JxW(q_point));
                    }
                }
            }
            
            // compute the local forcing function                                 
            for(unsigned int i=0; i<dofs_per_cell; ++i) {
                for(unsigned int q_point=0; q_point<n_q_points; ++q_point) {
                    cell_rhs(i) += (fe_values.shape_value(i, q_point) *
                                    1 *
                                    fe_values.JxW(q_point));
                }
            }
            
            cell->get_dof_indices(local_dof_indices);
            
            // add the contribution of this cell back to the system mass matrix
            for(unsigned int i=0; i<dofs_per_cell; ++i) {
                for(unsigned int j=0; j<dofs_per_cell; ++j) {
                    system_matrix.add(local_dof_indices[i],
                                      local_dof_indices[j],
                                      cell_matrix(i,j));
                }
            }
            
            // add the contribution of this cell back to the forcing function
            for(unsigned int i=0; i<dofs_per_cell; ++i) {
                system_rhs(local_dof_indices[i]) += cell_rhs(i);
            }
        }
        
        // apply our dirchlet boundary values (of zero) to the system
        std::map<types::global_dof_index,double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler, 0, ZeroFunction<dim>(), boundary_values);
        MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs);
        
        std::cout << "Finished assemble_system()"
                  << std::endl
                  << std::endl;
    }
    
    template<int dim> void PoissonEqn<dim>::solve() {
        // solve the system for the heat energy at each point on the mesh
        SolverControl solver_control(1000, 1e-12);
        SolverCG<>    solver(solver_control);
        solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
        
        std::cout << "Finished solve()"
                  << std::endl
                  << std::endl;
    }
    
    template<int dim> void PoissonEqn<dim>::output_results() const {
        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(solution, "solution");
        data_out.build_patches();
        
        const std::string filename = "data/solution-" +
                                     Utilities::int_to_string(dim,1) + 
                                     "D.vtu";
        
        std::ofstream output(filename.c_str());
        data_out.write_vtu(output);
        
        std::cout << "Finished output_results()"
                  << std::endl
                  << std::endl;
    }
    
    template<int dim> void PoissonEqn<dim>::run() {
        make_grid();
        setup_system();
        assemble_system();
        solve();
        output_results();
    }
}

int main() {
    try {
        
        using namespace dealii;
        using namespace Poission;
        
        deallog.depth_console(2);
        PoissonEqn<2> poisson_problem;
        poisson_problem.run();
        return 0;
        
    } catch(std::exception &exc) {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    } catch(...) {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }
}
        
        
    
