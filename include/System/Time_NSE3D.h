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
#include <Example_NSE3D.h>

#include <NSE_MultiGrid.h>

#include <MainUtilities.h> // FixedSizeQueu

#include <vector>
#include <deque>
#include <list>
#include <utility>
#include <array>

#ifdef _MPI
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

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

#ifdef _MPI
      /** @brief A parallel FE mapper storing parallel information
       * for the velocity degrees of freedom.*/
      TParFEMapper3D parMapperVelocity_;

      /** @brief A parallel FE mapper storing parallel information
       * for the pressure degrees of freedom.*/
      TParFEMapper3D parMapperPressure_;

      /** @brief A parallel FE communicator taking care for the MPI
       * communication between the velocity dofs on this grid. */
      TParFECommunicator3D parCommVelocity_;

      /** @brief A parallel FE communicator taking care for the MPI
       * communication between the pressure dofs on this grid. */
      TParFECommunicator3D parCommPressure_;
#endif

      /** TODO Implement the constructor used when ifdef_MPI
       * @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] maxSubDomainPerDof Only in MPI case! The maximal number of
       * processes which share a single dof. The value is calculated in the
       * main program and handed down to the FESpaces. Getting rid of this
       * construction is a TODO .
       */
//#ifdef _MPI
//      System_per_grid(const Example_NSE3D& example,
//                    TCollection& coll, std::pair<int, int> order, Time_NSE3D::Matrix type,
//                    int maxSubDomainPerDof);
//#else
      System_per_grid(const Example_NSE3D& example, TCollection& coll,
                      std::pair<int, int> order, Time_NSE3D::Matrix type);
//#endif

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

    /** @brief A complete system on each involved grid.
     *
     * Note that the size of this double ended queue is at least one and
     * larger only when a multgrid preconditioner is used.
     */
    std::deque<System_per_grid> systems_;

    /** @brief Definition of the used example. */
    const Example_NSE3D& example_;

    /** @brief A multigrid object which is set to nullptr in case it is not
     *         needed.
     */
    std::shared_ptr<TNSE_MultiGrid> multigrid_;
    /// This sorry thing is needed for multigrid with NSTypes 1 or 3, where
    /// transposed blocks are not stored explicitely...sad but true.
    std::vector<std::shared_ptr<TStructure>> transposed_B_structures_;

    //! @brief An array to store the current defect.
    BlockVector defect_;
    
    ///@brief The norm of residual from the previous iteration
    double old_residual_;

    /*! @brief The initial residual. Stored so that the nonlinear iteration can
    !        be stopped as soon as a desired reduction is achieved */
    double initial_residual_;
    
    /** @brief Errors, held in ready to be accessed from outside the class
     * The array is filled during the function call NSE3D::output()
     * Currently, the errors store the L2 and H1-semi errors of the velocity
     * (errors.at(0) is L2 and errors.at(1) is H1-semi)
     * and the pressure (errors.at(2) is L2 and errors.at(3) is H1-semi).
     */
    std::array<double, int(4)> errors_;

    /** @brief right hand side vector from previous time step (on finest mesh)*/
    BlockVector old_rhs_;

    BlockVector old_solution_;

    /** old time step length used to scale the pressure blocks*/
    double oldtau_;

    /** TODO Implement this method
     * @brief set the velocity and pressure orders
     *
     * This function sets the corresponding velocity and
     * pressure orders. The pressure order is set if it is
     * not specified by the readin file. Default is -4711
     *
     * Tried to stay with the function GetVelocityPressureSpace3D()
     * in the MainUtilities file.
     */
    void get_velocity_pressure_orders(std::pair <int,int>
                   &velocity_pressure_orders);
 public:

    /** @brief constructor
    * This constructor calls the other constructor creating an Example_NSE3D
    * object.
    */
   Time_NSE3D(const TDomain& domain, int reference_id = -4711);

    /** TODO Implement the constructor used when ifdef_MPI
     * @brief Standard constructor of an NSE3D problem.
     *
     * @note The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and
     * functions as well as matrices, solution and right hand side vectors are
     * initialized.
     *
     * @param domain The domain this problem lives on.
     * @param example The example to perform
     */
//#ifdef _MPI
//    Time_NSE3D(const TDomain& domain, const Example_NSE3D& example,
//          int maxSubDomainPerDof);
//#else
    Time_NSE3D(const TDomain& domain, const Example_NSE3D& example,
               int reference_id = -4711);
//#endif
    
    
// ======================================================================
     /**
//     * TODO Implement this method. What's the difference with "set_parameters()"
//     * in TNSE2D class?
//     * @brief Check whether the program will be working with the
//     * current input parameters.
//     *
//     * ParMooN is work in progress, and so is this class. This method
//     * checks the parameters stored in the database and stops execution
//     * of the program, if some of these do not match.
//     * The method is a little makeshift and the cases caught here are various,
//     * but basically it is intended to stop execution of cases with
//     * parameter combinations which are not implemented for NSE3D or
//     * are currently known to be problematic.
//     *
//     * This is not yet a guarantee for a functioning program, but is
//     * intended to be, someday. Eventually this method and the like
//     * will be moved to TDatabase.
//     */
//    static void check_parameters();

   /** @brief Assemble all the matrices and rhs before the time iterations
    *
    * This includes the assembling of the Stiffness matrix, the Mass matrix,
    * the additional matrix K in case of SUPG stabilization, and rhs.
    * This assembling occurs just once, before entering any loop. It assembles
    * linear terms only.
    */
    void assemble_initial_time();

//    * @brief Assemble the rhs only
//    * 1. Assembling the right hand side only
//    * 2. Scaling of the B-Blocks due to time stepping
//    * This function will prepare the right hand side during the time
//    * discretization but should be outside the nonlinear loop.
//    */
    void assemble_rhs();

//    /** TODO Implement this method.
//     * @brief Assemble the whole system matrix which will be passed
//     * to the solvers.
//     */
//    void assemble_system();
//
//    /** TODO Implement this method.
//     * @brief Descale matrices
//     * This function will descale all A-blocks which were scaled
//     * during the function call time Time_NSE3D::assemble_system().
//     */
//    void descale_matrices();
//
//    /** TODO Implement this method.
//     * @brief Assemble the nonlinear terms
//     * Assemble the nonlinear terms. Need not be used when this is
//     * a Stokes problem, or once per nonlinear iteration if this
//     * is a Navier Stokes problem.
//     * The matrix blocks to which the nonlinear term contribute are reset
//     * to zero and then completely reassembled, including the linear and
//     * nonlinear terms.
//     */
//    void assemble_nonlinear_term();
//
//    /** TODO Implement this method.
//    *Solve the current linear system. Nonlinear loop is outside of this class.
//    */
//    void solve();
//
//    /** TODO Implement this method.
//     *  @brief check if one of the stopping criteria is fulfilled
//     *
//     * either converged, maximum number of iterations reached, or slow
//     * convergence
//     *
//     * @param iteration_counter current iterate
//     *
//     * @note For the sequential case, this is the copy-paste NSE2D
//     * (with exception of slightly different compute_residual methods).
//     */
//    bool stop_it(unsigned int iteration_counter);
//
//    /** TODO Implement this method.
//     * ! Measure errors and draw a nice VTK picture, if requested to do so.
//    ! @param i suffix for output file name, -1 means no suffix. */
//    void output(int i = -1);
//
//    /** TODO Implement this method.
//     * @brief Compute the defect Ax-b, and the residuals and store it all.
//     *
//     * Updates defect and old_residuals.
//     * A is the current matrix, x is the current solution and b is the
//     * right hand side. Call this function after assembling the nonlinear
//     * matrix with the current solution.
//     */
//    void compute_residuals();
//
///*******************************************************************************/
//    // Declaration of special member functions - delete all but destructor.
//    // This problem class will be used to request the whole process of
//    // putting up and solving a single flow problem. It is, due to
//    // its procedural nature, not supposed to be copied.
//
//    //! Delete copy constructor. No copies allowed.
//    Time_NSE3D(const Time_NSE3D&) = delete;
//
//    //! Delete move constructor. No moves allowed.
//    Time_NSE3D(Time_NSE3D&&) = delete;
//
//    //! Delete copy assignment operator. No copies allowed.
//    Time_NSE3D& operator=(const Time_NSE3D&) = delete;
//
//    //! Default move assignment operator. No moves allowed.
//    Time_NSE3D& operator=(Time_NSE3D&&) = delete;
//
//    //! Default destructor. Does most likely cause memory leaks.
//    ~Time_NSE3D() = default;
//
///*******************************************************************************/
//    /**
//     * @brief initialize multigrid levels for different NSTYPE's
//     *
//     * @param: level
//     * @param: grid to be added
//     */
//    TNSE_MGLevel* mg_levels(int level, System_per_grid& s);
//
//    /**
//     * TODO Implement this method.
//     * @brief multigrid solver
//     * preparing the stuff used to call multigrid solver
//     */
//    void mg_solver();
//
///********************************************************************************/
//// getters
//   const BlockVector    & get_rhs()      const
//    { return this->systems_.front().rhs_; }
//
//   const TFEVectFunct3D & get_velocity() const
//    { return this->systems_.front().u_; }
//
//   const TFEFunction3D  & get_pressure() const
//     { return this->systems_.front().p_; }
   const TFESpace3D     & get_velocity_space() const
     { return this->systems_.front().velocitySpace_; }
   const TFESpace3D     & get_pressure_space() const
     { return this->systems_.front().pressureSpace_; }
//   const BlockVector    & get_solution() const
//     { return this->systems_.front().solution_; }
//
//   const int get_size() const
//     { return this->systems_.front().solution_.length(); }
//
//   const Example_NSE2D  & get_example()  const
//     { return example_; }
//
////  // try not to use the below methods which are not const
////  TFEVectFunct3D & get_velocity()
////    { return this->systems_.front().u_; }
////  TFEFunction3D *get_velocity_component(int i);
//
////  /// @brief Get the current residuals  (updated in compute_residuals)
////  const Residuals& get_residuals() const;
////  /// @brief get the current impuls residual (updated in compute_residuals)
////  double get_impuls_residual() const;
////  /// @brief get the current mass residual (updated in compute_residuals)
////  double get_mass_residual() const;
////  /// @brief get the current residual (updated in compute_residuals)
////  double get_full_residual() const;
//
//  /** TODO Implement this method.
//   * @brief return the computed errors (computed in output())
//   */
//  std::array<double, int(4)> get_errors() const;
};




#endif /* INCLUDE_SYSTEM_TIME_NSE3D_H_ */
