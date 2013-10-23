#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/base/utilities.h>

namespace Wave
{

    using namespace dealii;
    
    template<int dim> class WaveEquation {
        public:
        WaveEquation();
        void run();
        
        private:
        void setup_system();
        void solve_u();
        void solve_v();
        void output_results() const;
        
        Triangulation<dim> triangulation;
        FE_Q<dim>          fe;
        DoFHandler<dim>    dof_handler;
        ConstraintMatrix   constraints;

        SparsityPattern        sparsity_pattern;
        SparseMatrix<double> mass_matrix;
        SparseMatrix<double> laplace_matrix;
        SparseMatrix<double> matrix_u;
        SparseMatrix<double> matrix_v;
        
        Vector<double>         solution_u, solution_v;
        Vector<double>         old_solution_u, old_solution_v;
        Vector<double>         system_rhs;
        
        double       time, time_step;
        unsigned int timestep_number;
        const double theta;
    };
    
    template<int dim> class InitialValuesU : public Function<dim> {
        public:
        InitialValuesU() : Function<dim>() {}
        
        virtual double value(const Point<dim> &p, const unsigned int component = 0) const {
            Assert(component == 0, ExcInternalError());
            return 0;
        }
    };
    
    template<int dim> class InitialValuesV : public Function<dim> {
        public:
        InitialValuesV() : Function<dim>() {}
        
        virtual double value(const Point<dim> &p, const unsigned int component = 0) const {
            Assert(component == 0, ExcInternalError());
            return 0;
        }
    };
    
    template<int dim> class RightHandSide : public Function<dim> {
        public:
        RightHandSide() : Function<dim>() {}
        virtual double value(const Point<dim> &p, const unsigned int component = 0) const {
            Assert(component == 0, ExcInternalError());
            return 0;
        }
    }
    
    template<int dim> class BoundaryValuesU : Function<dim> {
        public:
        BoundaryValuesU() : Function<dim>() {}
        virtual double value(const Point<dim> &p, const unsigned int component = 0) const {
            Assert(component == 0, ExcInternalError());
            
            if((this->get_time() <= 0.5) && 
               (p[0] < 0) && 
               (p[1] <  1/3.0) &&
               (p[1] > -1/3.0)) {
            
                // This sinewave is the boundary conditions on the system
                return std::sin(this->get_time() * 4 * numbers::PI);
            } else {
                return 0;
            }
        }
    }
    
    template<int dim> class BoundaryValuesU : Function<dim> {
        public:
        BoundaryValuesU() : Function<dim>() {}
        virtual double value(const Point<dim> &p, const unsigned int component = 0) const {
            Assert(component == 0, ExcInternalError());
            
            if((this->get_time() <= 0.5) && 
               (p[0] < 0) && 
               (p[1] <  1/3.0) &&
               (p[1] > -1/3.0)) {
            
                // This is the time derivitave of the boundary conditions on U
                return std::cos(this->get_time() * 4 * numbers::PI) * 4 * numbers::PI;
            } else {
                return 0;
            }
        }
    }
    
    //Constructor definition
    template<int dim> WaveEquation<dim>::WaveEquation() :
        fe(1),
        dof_handler(triangulation),
        time_step(1/64.0),
        theta(0.5) {}
        
    template<int dim> void WaveEquation<dim>::setup_system() {
        GridGenerator::hyper_cube(triangulation, -1, 1);
        triangulation.refine_local(7);
        
        sparsity_pattern.reinit(dof_handler.n_dofs(), 
                                dof_handler.n_dofs(), 
                                dof_handler.max_couplings_between_dofs());
        DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
        sparsity_pattern.compress();
        
        mass_matrix.reinit(sparsity_pattern);
        laplace_matrix.reinit(sparsity_pattern);
        matrix_u.reinit(sparsity_pattern);
        matrix_v.reinit(sparsity_pattern);
        
        MatrixCreator::create_mass_matrix(dof_handler, QGauss<dim>(3), mass_matrix);
        MatrixCreator::create_laplace_matrix(dof_handler, QGauss<dim>(3), laplace_matrix);
        
        solution_u.reinit(dof_handler.n_dofs());
        solution_v.reinit(dof_handler.n_dofs());
        old_solution_u.reinit(dof_handler.n_dofs());
        old_solution_v.reinit(dof_handler.n_dofs());
        system_rhs.reinit(dof_handler.n_dofs());
        
        // QUESTION: when do we normally need the constraints matrix?
        
        constraints.close();
    }
    
    // QUESTION: Why is resolution based on system_rhs?
    // QUESTION: How do we tell when a preconditioner is a good thing to use?
    // QUESTION: What exactly does the preconditoner precondition?
    
    template<int dim> void WaveEquation<dim>::solve_u() {
        SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
        SolverCG<>    cg(solver_control);
        cg.solve(matrix_u, solution_u, system_rhs, PreconditionIdentity());
    
    template<int dim> void WaveEquation<dim>::solve_v() {
        SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
        SolverCG<>    cg(solver_control);
        cg.solve(matrix_v, solution_v, system_rhs, PreconditionIdentity());
    }
    
    template<int dim> void WaveEquation<dim>::output_results() const {
        DataOut<dim> data_out;
        
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(solution_u, "U");
        data_out.add_data_vector(solution_v, "V");
        data_out.build_patches();
        
        const std::string filename = "data/solution-" + 
                                     Utilities::int_to_string(timestep_number, 3) + 
                                     ".vtu";
                                     
        std::ofstream output(filename.c_str());
        data_out.write_vtu(output);
    }
    
    template<int dim> void WaveEquation<dim>::run() {
        setup_system();
        
        VectorTools::project(dof_handler, constraints, QGauss<dim>(3), InitialValues<dim>(), old_solution_u);    
        VectorTools::project(dof_handler, constraints, QGauss<dim>(3), InitialValues<dim>(), old_solution_v);
        
        Vector<double> tmp(solution_u.size());
        Vector<double> forcing_terms(solution_u.size());
        
        for (timestep_number=1, time=time_step; 
             time <= 5;
             time+=time_step, ++timestep_number) {
             
            mass_matrix.vmult(system_rhs, old_solution_u);
            mass_matrix.vmult(tmp, old_solution_v);
            system_rhs.add(time_step, tmp);
            
            laplace_matrix.vmult(tmp, old_solution_u);
            system_rhs.add(-theta *(1-theta) * time_step * time_step, tmp);
            
            RightHandSide<dim> rhs_function;
            rhs_function.set_time(time);
            VectorTools::craete_right_hand_side(dof_handler, QGauss<dim>(2), rhs_function, tmp);
            
            forcing_terms = tmp;
            forcing_terms *= theta * time_step;
            
            rhs_function.set_time(time-time_step);
            VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(2), rhs_function, tmp)
            
            forcing_terms.add ((1-theta) * time_step, tmp);

            system_rhs.add (theta * time_step, forcing_terms);
            
             
    
        
        
        
        
}

int main() {

    try {
    
        using namespace dealii;
        using namespace Wave;
        
        deallog.depth_console(0);
        
        Wave::InitialValuesU<2> a = Wave::InitialValuesU<2>();
        Wave::InitialValuesV<2> b = Wave::InitialValuesV<2>();
    
    } catch (std::exception &exc) {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    } catch (...) {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
      }
    return 0;
}
    
    
    
    
        
