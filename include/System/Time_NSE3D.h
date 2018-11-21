/** ***************************************************************************
 *
 * @name       Time_NSE3D
 * @brief      store everything needed to solve a time-dependent
 *             Navier-Stokes flow problem in 3D
 *
 *             Store matrix, right hand side, FE spaces, FE functions and
 *             the solution vector of a NSE problem in 3D. When a multigrid
 *             solver is used then all of this is stored on multiple grids.
 *
 * @author     Najib Alia
 * @date       2016/04/16
 * @history    2016/04/16 
 * @history    2016/06/16 (multigrid enabled) Naveed
 * @history    2016/10/26 (multigrid enabled)
 *
 ******************************************************************************/

#ifndef INCLUDE_SYSTEM_TIME_NSE3D_H_
#define INCLUDE_SYSTEM_TIME_NSE3D_H_

#include <FEVectFunct3D.h>
#include <FEFunction3D.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Residuals.h>
#include <FESpace3D.h>
#include <Example_TimeNSE3D.h>
#include <TimeDiscretizations.h>
#include "LocalAssembling.h"
#include <ParameterDatabase.h>
#include <Solver.h>
#include <DataWriter.h>
#include <MainUtilities.h> // FixedSizeQueu

#include <vector>
#include <deque>
#include <list>
#include <utility>
#include <array>

class Time_NSE3D
{

  enum class Matrix{Type14, Type1, Type2, Type3, Type4};

  protected:

    /** @brief store a complete system on a particular grid
     *
     * This combines a matrix, rhs, solution, spaces and functions
     * needed to describe one time-dependent Navier Stokes flow problem in 3D.
     * In MPI case the parallel infrastructure is stored, too.
     */
    struct System_per_grid
    {
      /** @brief Finite Element space for the velocity */
      TFESpace3D velocitySpace_;
      /** @brief Finite Element space for the pressure */
      TFESpace3D pressureSpace_;

      /** @brief the system matrix (depends strongly on
       *         TDatabase::ParamDB->NSTYPE)
       *  [ A11  A12  A13  B1T ]
       *  [ A21  A22  A23  B2T ]
       *  [ A31  A32  A33  B3T ]
       *  [ B1   B2   B3   C   ]
       */
      BlockFEMatrix matrix_;
      /** @brief Mass matrix
       *  [ M11  0    0    0 ]
       *  [ 0    M22  0    0 ]
       *  [ 0    0    M33  0 ]
       *  [ 0    0    0    0 ]
       */
      BlockFEMatrix massMatrix_;
      /** @brief the right hand side vector */
      BlockVector rhs_;
      /** @brief solution vector with two components. */
      BlockVector solution_;
      /** @brief Finite Element function for velocity */
      TFEVectFunct3D u_;
      /** @brief Finite Element function for pressure */
      TFEFunction3D p_;

      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] maxSubDomainPerDof Only in MPI case! The maximal number of
       * processes which share a single dof. The value is calculated in the
       * main program and handed down to the FESpaces. Getting rid of this
       * construction is a TODO .
       */
      System_per_grid(const Example_TimeNSE3D& ex, TCollection& coll,
                      std::pair<int, int> order, Time_NSE3D::Matrix type);

      // System_per_grid is not supposed to be copied or moved
      // until underlying classes realize the rule of zero.

      //! Delete copy constructor. No copies allowed.
      System_per_grid( const System_per_grid& ) = delete;

      //! Delete move constructor. No moves allowed.
      System_per_grid( System_per_grid&& ) = delete;

      //! Delete copy assignment operator. No copies allowed.
      System_per_grid& operator=( const System_per_grid& ) = delete;

      //! Default move assignment operator. No moves allowed.
      System_per_grid& operator=( System_per_grid&& ) = delete;

      //! Default destructor. Does most likely cause memory leaks.
      ~System_per_grid() = default;
    };

    /** @brief a local parameter database which controls this class
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * default ParMooN parameters). Note that this usually does not include
     * other parameters such as solver parameters. Those are only in the
     * Solver object.
     */
    ParameterDatabase db_;
    
    /** @brief output object */
    DataWriter3D outputWriter;

    /** @brief A complete system on each involved grid.
     *
     * Note that the size of this double ended queue is at least one and
     * larger only when a multgrid preconditioner is used.
     */
    std::deque<System_per_grid> systems_;

    /** @brief Definition of the used example. */
    const Example_TimeNSE3D& example_;

    /** @brief a solver object which will solve the linear system
     *
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver_;

    //! @brief An array to store the current defect.
    BlockVector defect_;
    
    ///@brief The norm of residuals from up to 10 previous iterations
    FixedSizeQueue<10, Residuals> old_residual_;

    /*! @brief The initial residual. Stored so that the nonlinear iteration can
    !        be stopped as soon as a desired reduction is achieved */
    double initial_residual_;
    
    /** @brief Errors, held in ready to be accessed from outside the class
     * The array is filled during the function call TNSE3D::output()
     * Currently, the errors store the L2 and H1-semi errors of the velocity
     * (errors.at(0) is L2 and errors.at(1) is H1-semi)
     * and the pressure (errors.at(2) is L2 and errors.at(3) is H1-semi).
     * From 5 to 8, errors corrected with time tau*0.5
     */
    std::array<double, int(12)> errors_;

    /** @brief right hand side vector from previous time step (on finest mesh)*/
    BlockVector old_rhs_;

    /** @brief solution vectors from previous time steps (on finest mesh)*/
    BlockVector solution_m1;
    BlockVector solution_m2;

    /** @brief constructs a solution vector extrapolated from previous steps
     * Currently, it is used for IMEX-Scheme: 2u(t-1)-u(t-2). */
    BlockVector extrapolated_solution_;

    /** @brief old time step length used to scale the pressure blocks */
    double oldtau_;
    
    /// @brief time stepping scheme object to access everything
    TimeDiscretization time_stepping_scheme;
    
    /// @brief is the mass matrix and right hand side is solution
    /// dependent: for example the SUPG method
    bool is_rhs_and_mass_matrix_nonlinear;
    /// @brief system right hand side which passes to the solver
    BlockVector rhs_from_time_disc;
    
    /// this is for setting the discretization type globally, since the DISCTYPE 
    /// is removed fromt the global database but it's not simple/possible to use
    /// different space discretization in the LocalAssembling2D
    int space_disc_global;

    /** @brief check parameters in database
    *
    * This functions checks if the parameters in the database are meaningful.
    * If some parameters are set to unsupported values, an error occurs and
    * throws an exception.
    */
    void check_and_set_parameters();
    
    /** @brief set the velocity and pressure orders
     *
     * This function sets the corresponding velocity and
     * pressure orders. The pressure order is set if it is
     * not specified by the reading file. Default is -4711
     *
     * Tried to stay with the function GetVelocityPressureSpace3D()
     * in the MainUtilities file.
     */
    void get_velocity_pressure_orders(std::pair <int,int>
                   &velocity_pressure_orders);
    
    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;

 public:

    /** @brief Standard constructor of an NSE3D problem.
     *
     * @note The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and
     * functions as well as matrices, solution and right hand side vectors are
     * initialized.
     *
     * @param domain
     * @param example The example to perform
     */
    Time_NSE3D(const TDomain& domain, const ParameterDatabase& param_db, 
               const Example_TimeNSE3D& example);
    
// ======================================================================
    /** @brief This returns the number of the current time step.
     * This counter is set at 0 before the time loop and is incremented at each
     * time step (but not at each sub-step) in the main program.
     * It can be useful to give info to the members of the class. It is for example
     * used in IMEX scheme to detect when we passed 2 time steps, so that we
     * are guaranteed to have saved both old_solution_ and old_solution2_ correctly.    */
    int current_step_;

// ======================================================================
    /** @brief Assemble all the matrices and rhs before the time iterations
    *
    * This includes the assembling of the Stiffness matrix, the Mass matrix,
    * the additional matrix K in case of SUPG stabilization, and rhs.
    * This assembling occurs just once, before entering any loop. It assembles
    * linear terms only.
    */
    void assemble_initial_time_m();
    void assemble_initial_time();

    /** @brief Assemble the rhs only
    * 1. Assembling the right hand side only
    * 2. Scaling of the B-Blocks due to time stepping
    * This function will prepare the right hand side during the time
    * discretization but should be outside the nonlinear loop.
    */
    void assemble_rhs();
    void assemble_rhs_m();

    /** @brief Assemble the nonlinear terms
     * Assemble the nonlinear terms. Need not be used when this is
     * a Stokes problem, or once per nonlinear iteration if this
     * is a Navier Stokes problem.
     * The matrix blocks to which the nonlinear term contribute are reset
     * to zero and then completely reassembled, including the linear and
     * nonlinear terms.
     */
    void assemble_nonlinear_term();

    /** @brief Assemble the whole system matrix which will be passed
     * to the solvers.
     */
    void assemble_system();

    /** @brief Solve the current linear system. Nonlinear
     * loop is outside of this class.
    */
    void solve();

    /** @brief check if one of the stopping criteria is fulfilled
     *
     * either converged, maximum number of iterations reached, or slow
     * convergence
     *
     * @param iteration_counter current iterate
     *
     * @note For the sequential case, this is the copy-paste NSE2D
     * (with exception of slightly different compute_residual methods).
     */
    bool stop_it(unsigned int iteration_counter);

    /** @brief Compute the defect Ax-b, and the residuals and store it all.
     * This method is also the one displaying the residuals.
     * Updates defect and old_residuals.
     * A is the current matrix, x is the current solution and b is the
     * right hand side. Call this function after assembling the nonlinear
     * matrix with the current solution.
     */
    void compute_residuals();

    /** @brief Measure errors and draw a nice VTK picture,
     * if requested to do so.
     */
    void output(int m);

    /** @brief Construct the extrapolated solution.
     * At the moment, only IMEX is implemented. */
    void construct_extrapolated_solution();

    /** */
    bool imex_scheme(bool print_info);

/* ******************************************************************************/
    // Declaration of special member functions - delete all but destructor.
    // This problem class will be used to request the whole process of
    // putting up and solving a single flow problem. It is, due to
    // its procedural nature, not supposed to be copied.

    //! Delete copy constructor. No copies allowed.
    Time_NSE3D(const Time_NSE3D&) = delete;

    //! Delete move constructor. No moves allowed.
    Time_NSE3D(Time_NSE3D&&) = delete;

    //! Delete copy assignment operator. No copies allowed.
    Time_NSE3D& operator=(const Time_NSE3D&) = delete;

    //! Default move assignment operator. No moves allowed.
    Time_NSE3D& operator=(Time_NSE3D&&) = delete;

    //! Default destructor. Does most likely cause memory leaks.
    ~Time_NSE3D() = default;

/* *******************************************************************************/

    //// getters

    /// Get the velocity space.
    const TFESpace3D     & get_velocity_space() const
    { return this->systems_.front().velocitySpace_; }

    /// Get the pressure space.
    const TFESpace3D     & get_pressure_space() const
    { return this->systems_.front().pressureSpace_; }

    const TFEVectFunct3D& get_velocity() const
    { return this->systems_.front().u_; }

    TFEVectFunct3D& get_velocity()
    { return this->systems_.front().u_; }

    TFEFunction3D *get_velocity_component(int i);

    const TFEFunction3D& get_pressure() const
    { return this->systems_.front().p_; }

    TFEFunction3D& get_pressure()
    { return this->systems_.front().p_; }

    /// Get number of degrees of freedom.
    const int get_size() const
    { return this->systems_.front().solution_.length(); }
    
    /// @brief get the space discretization type
    int get_space_disc_global() {return space_disc_global;}
    
    /// Get the stored database.
    const ParameterDatabase & get_db() const
    { return db_; }

    /// @brief Get the current residuals  (updated in compute_residuals)
    const Residuals& get_residuals() const;
    /// @brief get the current impulse residual (updated in compute_residuals)
    double get_impulse_residual() const;
    /// @brief get the current mass residual (updated in compute_residuals)
    double get_mass_residual() const;
    /// @brief get the current residual (updated in compute_residuals)
    double get_full_residual() const;

    /** @brief return the computed errors (computed in output()) */
    std::array<double, int(6)> get_errors() const;
    
    /// @brief get the time stepping object
    TimeDiscretization& get_time_stepping_scheme()
    {return time_stepping_scheme;}
    const TimeDiscretization& get_time_stepping_scheme() const
    {return time_stepping_scheme;}
    
private:
  
  /// @brief this routines wraps up the call to Assemble3D
  void call_assembling_routine(Time_NSE3D::System_per_grid& s,
                               LocalAssembling_type type);
  /// @brief set the matrices and right hand side depending on the
  /// assemling routines, nstypes and the methods
  void set_matrices_rhs(Time_NSE3D::System_per_grid& s,
                        LocalAssembling_type type,
                        std::vector<TSquareMatrix3D*> &sqMat,
                        std::vector<TMatrix3D*> &reMat,
                        std::vector<double*> &rhs_array);
  /// @brief set the spaces depending on disc types
  void set_arrays(Time_NSE3D::System_per_grid& s,
                  std::vector<const TFESpace3D*> &spaces,
                  std::vector< const TFESpace3D* >& spaces_rhs,
                  std::vector< TFEFunction3D*> &functions);  
  
  /// @brief
  void restrict_function();
};




#endif /* INCLUDE_SYSTEM_TIME_NSE3D_H_ */
