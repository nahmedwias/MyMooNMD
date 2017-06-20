/**
 *
 * @name   NSE2D_axialsymmetric.h
 * @brief  System class for an axisymmetric NSE2D system.
 * It has some subtle, yet severe differences from the "standard"
 * NSE2D problem, this is why I decided to put it into its own class.
 * Will lead to dupe code to some extent, I'm afraid.
 *
 * @author     Clemens Bartsch
 * @date       16.06.17
 */

#ifndef USER_PROJECTS_INC_NSE2D_AXISYMMETRIC_H_
#define USER_PROJECTS_INC_NSE2D_AXISYMMETRIC_H_

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Residuals.h>
#include <ParameterDatabase.h>
#include <Solver.h>
#include <Example_NSE2D.h>
#include <PostProcessing2D.h>
#include <MainUtilities.h> // FixedSizeQueue

#include <utility>
#include <array>


class NSE2D_axisymmetric
{
    /// @brief store a complete system on a particular grid.
    struct System_per_grid
    {
      /** @brief Finite Element space for the axial velocity */
      TFESpace2D velocity_space_z;

      /** @brief Finite Element space for the radial velocity */
      TFESpace2D velocity_space_r;

      /** @brief Finite Element space for the pressure */
      TFESpace2D pressure_space;

      /** The system matrix. */
      BlockFEMatrix matrix;

      /** @brief The right hand side vector */
      BlockVector rhs;
      /** @brief The solution vector. */
      BlockVector solution;
      /** @brief Finite Element function for axial velocity component. */
      TFEFunction2D u_z;
      /** @brief Finite Element function for radial velocity component. */
      TFEFunction2D u_r;
      /** @brief Finite Element function for pressure */
      TFEFunction2D p;

      /** @brief constructor */
      System_per_grid(const Example_NSE2D& example, TCollection& coll,
                      std:: pair <int,int> velocity_pressure_orders);

      //! These things are not to be copied or moved.
      System_per_grid(const System_per_grid&) = delete;
      System_per_grid(System_per_grid&&) = delete;
      System_per_grid& operator=(const System_per_grid&) = delete;
      System_per_grid& operator=(System_per_grid&&) = delete;
      ~System_per_grid() = default;
    };

    /// @brief a complete system on each grid.
    std::deque<System_per_grid> systems;

    /** @brief Definition of the used example */
    const Example_NSE2D example;

    /// @brief The database object. It contains control parameters.
    ParameterDatabase db;

    /** @brief output handling (vtk and case files) */
    PostProcessing2D outputWriter;

    /** @brief a solver object which will solve the linear system*/
    Solver<BlockFEMatrix, BlockVector> solver;

    //! @brief Store the current defect vector.
    BlockVector defect;

    ///@brief The norms of residuals from up to 10 previous iterations
    FixedSizeQueue<10, Residuals> oldResiduals;

    /** @brief store the initial residual so that the nonlinear iteration can
     *         be stopped as soon as a desired reduction is achieved
     */
    double initial_residual;

    /// Store the initial rhs norm, so that the nonlinear iteration can be
    /// stopped relative to the initial right hand side.
    double initial_rhs_norm;

    /** @brief Errors to be accesed from outside the class
     * The array is filled during the function call NSE2D::output()
     * Currently, the errors store the L2 and H1 errors of the velocity
     * and pressure
     */
    std::array<double, int(4)> errors;

    /** @brief set parameters in database
     *
     * This functions checks if the parameters in the database are meaningful
     * and resets them otherwise. The hope is that after calling this function
     * this class is fully functional.
     *
     * If some parameters are set to unsupported values, an error occurs and
     * throws an exception.
     */
    void set_parameters();

    /** @brief set the velocity and pressure orders
     *
     * This function sets the corresponding velocity and
     * pressure orders. The pressure order is set if it is
     * not specified by the readin file. Default is -4711
     */
    void get_velocity_pressure_orders(std::pair <int,int>
                   &velocity_pressure_orders);

    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;

  public:
    /// COnstrcutor, taking a grid or grid hiearchie, a controlling database and an example.
    NSE2D_axisymmetric(std::list<TCollection*> grids, const ParameterDatabase& param_db,
          const Example_NSE2D _example);

    /**
     * Assemble those parts which do not contain nonlinearities
     * i.e. a Stokes problem. When solving a Navier Stokes problem
     * then this must be called once before entering the nonlinear loop.
     */
    void assemble_linear_terms();

    /**
     * Assemble the nonlinear term. Need not be used when this is
     * a Stokes problem, or once per nonlinear iteration if this
     * is a Navier Stokes problem.
     */
    void assemble_nonlinear_term();

    /** @brief solve the system */
    void solve();

    /**
     * @brief Compute the defect Ax-b and store its norm in NSE2D::oldResiduals
     *
     * where A is the current matrix, x is the current solution and b is the
     * right hand side. Call this function after assembling the nonlinear
     * matrix with the current solution.
     */
    void compute_norms_of_residuals();

    /** @brief check if one of the stopping criteria is fulfilled
     *
     * either converged, maximun number of iterations reached, or slow
     * convergence
     *
     * @param iteration_counter current iterate
     */
    bool stop_iteration(unsigned int iteration_counter);

    /**
     * @brief measure errors and write pictures
     *
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk files are created.
     *
     * @param i suffix for output file name, -1 means no suffix
     */
    void output(int i = -1);

    // getters and setters
    const BlockFEMatrix & get_matrix() const
    { return systems.front().matrix; }
    BlockFEMatrix & get_matrix()
    { return systems.front().matrix; }

    const BlockVector & get_rhs() const
    { return systems.front().rhs; }
    BlockVector & get_rhs()
    { return systems.front().rhs; }

    const TFEFunction2D & get_axial_velocity() const
    { return systems.front().u_r; }
    TFEFunction2D & get_axial_velocity()
    { return systems.front().u_r; }

    const TFEFunction2D & get_radial_velocity() const
    { return systems.front().u_r; }
    TFEFunction2D & get_radial_velocity()
    { return systems.front().u_r; }

    const TFEFunction2D & get_pressure() const
    { return this->systems.front().p; }
    TFEFunction2D & get_pressure()
    { return this->systems.front().p; }

    const TFESpace2D & get_axial_velocity_space() const
    { return this->systems.front().velocity_space_z; }
    const TFESpace2D & get_radial_velocity_space() const
    { return this->systems.front().velocity_space_r; }

    const TFESpace2D & get_pressure_space() const
    { return this->systems.front().pressure_space; }

    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    BlockVector & get_solution()
    { return this->systems.front().solution; }

    const LoopInfo& get_it_solver_info()
    {return solver.get_solver_loop_info();}

    unsigned int get_size() const
    { return this->systems.front().solution.length(); }

    const Example_NSE2D & get_example() const
    { return example; }

    /// @brief get the current residuals
    /// @details updated in NSE2D::computeNormsOfResiduals which in turn is
    /// called from NSE2D::stopIt
    const Residuals& getResiduals() const;
    /// @brief get the current impuls residual
    /// @details updated in NSE2D::computeNormsOfResiduals which in turn is
    /// called from NSE2D::stopIt
    double getImpulsResidual() const;
    /// @brief get the current mass residual
    /// @details updated in NSE2D::computeNormsOfResiduals which in turn is
    /// called from NSE2D::stopIt
    double getMassResidual() const;
    /// @brief get the current residual
    /// @details updated in NSE2D::computeNormsOfResiduals which in turn is
    /// called from NSE2D::stopIt
    double getFullResidual() const;
    /// @brief return the computed errors
    /// @details updated in NSE2D::stopIt
    std::array<double, int(4)> get_errors() const;
};

#endif /* USER_PROJECTS_INC_NSE2D_AXISYMMETRIC_H_ */
