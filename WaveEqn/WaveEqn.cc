// If VERBOSE is defined, the system will print out a lot more information about
// what it's doing.
#define VERBOSE


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
  
  // Class Defn for keeping track of all of the data for solving the wave equation
  template<int dim> class WaveEquation {
    public:
    WaveEquation();
    void run();
    
    private:
    void generate_mesh();
    void setup_system();
    void solve_u();
    void solve_v();
    void output_results() const;
    
    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;
    ConstraintMatrix   constraints;

    SparsityPattern      sparsity_pattern;
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
    
  // Class Defn for providing initial values to the U vector (position of membrane)
  //   Currently, initial values are set to 0.
  template<int dim> class InitialValuesU : public Function<dim> {
    public:
    InitialValuesU() : Function<dim>() {}
      
    virtual double value(const Point<dim> &p, const unsigned int component = 0) const {
      Assert(component == 0, ExcInternalError());
      return 0;
    }
  };
    
  // Class Defn for providing initial values to the V vector (velocity of membrane)
  //   Currently, initial values are set to 0.
  template<int dim> class InitialValuesV : public Function<dim> {
    public:
    InitialValuesV() : Function<dim>() {}
    
    virtual double value(const Point<dim> &p, const unsigned int component = 0) const {
      Assert(component == 0, ExcInternalError());
      return 0;
    }
  };
    
  // Class Defn for providing the forcing function vector
  //   This would provide energy inside the membrane
  //   Currently, the forcing function is set to zero.
  template<int dim> class RightHandSide : public Function<dim> {
    public:
    RightHandSide() : Function<dim>() {}
    virtual double value(const Point<dim> &p, const unsigned int component = 0) const {
      Assert(component == 0, ExcInternalError());
      return 0;
    }
  };
  
  // Class Defn for providing boundary values on U (the position of the membrane)
  //   The membrane is clamped to 0 on all sides.  
  //   For half a second, 2/3 of the right hand side of the membrane is snapped in a sine wave
  //   pattern.  Then it returns to 0 like the rest of the conditions.
  template<int dim> class BoundaryValuesU : public Function<dim> {
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
  };
    
  // Class Defn for providing boundary values on V (the velocity of the membrane)
  //   The membrane is clamped to 0 on all sides.
  //   The initial boundary conditions for the velocity are the time derivative of the
  //   boundary conditions for the U vector.
  template<int dim> class BoundaryValuesV : public Function<dim> {
    public:
    BoundaryValuesV() : Function<dim>() {}
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
  };
  
  // Constructor Defn for wave equation class
  //   We initialize the test function order to 1, give the DoF Handler the location of the mesh,
  //   set the time step size, and choose theta such that the Crank-Nichols time discretization
  //   method is used.
  template<int dim> WaveEquation<dim>::WaveEquation() :
    fe(1),
    dof_handler(triangulation),
    time_step(1/64.0),
    theta(0.5) {}
  
  // Function Defn for generating the mesh
  //   builds a hypercube (rectangle in 2d) on which the wave propegates    
  template<int dim> void WaveEquation<dim>::generate_mesh() {
    GridGenerator::hyper_cube(triangulation, -1, 1);
    triangulation.refine_global(7);

    dof_handler.distribute_dofs(fe);
    
    #ifdef VERBOSE
    std::cout << "Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl << std::endl;
    #endif
  }
    
  // Function Defn for building matricies and initializing system
  //   Builds a sparsity pattern which is used on both the mass and laplace matrix
  //   Creates the mass and laplace matrix which are do not vary with time
  //   Reinitializes the solution vectors with the correct number of entries    
  template<int dim> void WaveEquation<dim>::setup_system() {
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
    
    // Guarenteed no hanging nodes because we used a uniformly refined mesh
    // so we don't need any entries in this matrix, so just close this matrix
    constraints.close();
  }
    
  // Function Defn for solving the U vector
  //   QUESTION: why is the solver resolution based on the l2norm of the rhs?
  template<int dim> void WaveEquation<dim>::solve_u() {
    SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
    SolverCG<>    cg(solver_control);
    cg.solve(matrix_u, solution_u, system_rhs, PreconditionIdentity());
    
    #ifdef VERBOSE
    std::cout << "   u-equation: " << solver_control.last_step()
              << " CG iterations."
              << std::endl;
    #endif
  }
  
  // Function Defn for solving the V vector
  //   QUESTION: why is the solver resolution based on the l2norm of the rhs?
  template<int dim> void WaveEquation<dim>::solve_v() {
    SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
    SolverCG<>    cg(solver_control);
    cg.solve(matrix_v, solution_v, system_rhs, PreconditionIdentity());
    
    #ifdef VERBOSE
    std::cout << "   v-equation: " << solver_control.last_step()
              << " CG iterations."
              << std::endl;  
    #endif  
  }
    
  // Function Defn for outputting the results of a single time step
  //   We write the output as a collection of .vtu files with lexiocographically ordered names
  //   by including the timestep value.  ParaView, the program I use to analyze the data
  //   is smart enough to open all the files as a group for animation.
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
  
  // Function Defn for running the system
  //   Computes the linear equation representing the U solution and solves for the U vector
  //   Computes the linear equation representing the V solution and solves for the V vector
  //   Outputs the solution
  //   Iterates time by 1 time step
  template<int dim> void WaveEquation<dim>::run() {
    generate_mesh();
    setup_system();
    
    VectorTools::project(dof_handler, constraints, QGauss<dim>(3), InitialValuesU<dim>(), old_solution_u);    
    VectorTools::project(dof_handler, constraints, QGauss<dim>(3), InitialValuesV<dim>(), old_solution_v);
    
    Vector<double> tmp(solution_u.size());
    Vector<double> forcing_terms(solution_u.size());
    
    for (timestep_number=1, time=time_step; 
       time <= 5;
       time+=time_step, ++timestep_number) {
      
      #ifdef VERBOSE
      std::cout << "Time step " << timestep_number
                << " at t=" << time
                << std::endl;
      #endif
      
      // Compute the vector that represents the right hand side of the equation
      // used to solve the U vector 
      mass_matrix.vmult(system_rhs, old_solution_u);
      mass_matrix.vmult(tmp, old_solution_v);
      system_rhs.add(time_step, tmp);
      
      laplace_matrix.vmult(tmp, old_solution_u);
      system_rhs.add(-theta *(1-theta) * time_step * time_step, tmp);
      
      RightHandSide<dim> rhs_function;
      rhs_function.set_time(time);
      VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(2), rhs_function, tmp);
      
      forcing_terms = tmp;
      forcing_terms *= theta * time_step;
      
      rhs_function.set_time(time-time_step);
      VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(2), rhs_function, tmp);
      
      forcing_terms.add ((1-theta) * time_step, tmp);

      system_rhs.add (theta * time_step, forcing_terms);
      
      // Compute the matrix used in solving the U vector and apply the boundary values 
      // to be applied to the matrix
      // NOTE: this section is put in it's own scope so that temporary values calculated
      //       cannot accidently be reused
      // QUESTION: unsure how we know we have to apply bounary values to this particular matrix
      //           and why they aren't put in the weak form like usual?
      {
        BoundaryValuesU<dim> boundary_values_u_function;
        boundary_values_u_function.set_time(time);
        std::map<types::global_dof_index,double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler, 
                                                 0, 
                                                 boundary_values_u_function, 
                                                 boundary_values);
                                                 
        // Compute the matrix used in solving the U vector
        matrix_u.copy_from(mass_matrix);
        matrix_u.add(theta*theta*time_step*time_step, laplace_matrix);
        MatrixTools::apply_boundary_values(boundary_values, matrix_u, solution_u, system_rhs);
      }
      solve_u();
      
      // Compute the vector that represents the right hand side of the equation
      // used to solve the V vector      
      laplace_matrix.vmult(system_rhs, solution_u);
      system_rhs *= -theta * time_step;
      
      mass_matrix.vmult(tmp, old_solution_v);
      system_rhs += tmp;
      
      laplace_matrix.vmult(tmp, old_solution_u);
      system_rhs.add(-time_step*(1-theta), tmp);
      
      system_rhs += forcing_terms;
      
      // Compute the matrix used in solving the V vector and apply the boundary values 
      // to be applied to the matrix
      // NOTE: this section is put in it's own scope so that temporary values calculated
      //       cannot accidently be reused
      {
        BoundaryValuesV<dim> boundary_values_v_function;
        boundary_values_v_function.set_time(time);
        
        std::map<types::global_dof_index,double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler, 
                                                  0, 
                                                  boundary_values_v_function, 
                                                  boundary_values);
        matrix_v.copy_from(mass_matrix);
        MatrixTools::apply_boundary_values(boundary_values, matrix_v, solution_v, system_rhs);
      }
      solve_v();  
      
      output_results();
      
      old_solution_u = solution_u;
      old_solution_v = solution_v;
      
      #ifdef VERBOSE
      std::cout << "   Total energy: "
                << (mass_matrix.matrix_norm_square (solution_v) +
                    laplace_matrix.matrix_norm_square (solution_u)) / 2
                << std::endl << std::endl;
      #endif
    }
  }
}

int main() {

    try {
    
        using namespace dealii;
        using namespace Wave;
        
        int depth = 0;
        #ifdef VERBOSE
        depth = 2;
        #endif
        
        deallog.depth_console(depth);
        
        WaveEquation<2> wave_equation_solver;
        wave_equation_solver.run();
    
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
    
    
    
    
        
