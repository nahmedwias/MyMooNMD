/** ***************************************************************************
 *
 * @name       NSE3D
 * @brief      store everything needed to solve a (stationary/instationary)
 *             Navier-Stokes flow problem
 *
 *             Store matrix, right hand side, FE spaces, FE functions and
 *             the solution vector of a NSE problem in 3D. When a multigrid
 *             solver is used then all of this is stored on multiple grids.
 *
 * @author     Clemens Bartsch, Naveed Ahmed
 * @date       2015/12/07 
 * @history    2016/02/21
 *
 ******************************************************************************/

#ifndef INCLUDE_SYSTEM_NSE3D_H_
#define INCLUDE_SYSTEM_NSE3D_H_

#include <FEVectFunct3D.h>
#include <FEFunction3D.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Residuals.h>
#include <FESpace3D.h>
#include <Example_NSE3D.h>

#include <ParameterDatabase.h>
#include <Solver.h>
#include <DataWriter.h>

#include <MainUtilities.h> // FixedSizeQueu

#include <vector>
#include <deque>
#include <list>
#include <utility>
#include <array>

class NSE3D
{

  enum class Matrix{Type14, Type1, Type2, Type3, Type4};

  protected:

    /** @brief store a complete system on a particular grid
     *
     * This combines a matrix, rhs, solution, spaces and functions
     * needed to describe one Navier Stokes flow problem in 3D.
     * In MPI case the parallel infrastructure is stored, too.
     */
    struct System_per_grid
    {
      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] maxSubDomainPerDof Only in MPI case! The maximal number of
       * processes which share a single dof. The value is calculated in the
       * main program and handed down to the FESpaces. Getting rid of this
       * construction is a TODO .
       */
      System_per_grid(const Example_NSE3D& example, TCollection& coll,
                      std::pair<int, int> order, NSE3D::Matrix type);

      /** @brief Finite Element space for the velocity */
      std::shared_ptr<TFESpace3D> velocitySpace_;
      /** @brief Finite Element space for the pressure */
      std::shared_ptr<TFESpace3D> pressureSpace_;

      /** @brief the system matrix (depends strongly on
       *         TDatabase::ParamDB->NSTYPE)
       *  [ A11  A12  A13  B1T ]
       *  [ A21  A22  A23  B2T ]
       *  [ A31  A32  A33  B3T ]
       *  [ B1   B2   B3   C   ]
       */
      BlockFEMatrix matrix_;
      /** @brief the right hand side vector */
      BlockVector rhs_;
      /** @brief solution vector with two components. */
      BlockVector solution_;
      /** @brief Finite Element function for velocity */
      TFEVectFunct3D u_;
      /** @brief Finite Element function for pressure */
      TFEFunction3D p_;

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
     * systems_.at(0) holds the finest grid.
     */
    std::deque<System_per_grid> systems_;

    /** @brief Definition of the used example. */
    const Example_NSE3D example_;
    
    /** @brief a local parameter database which controls this class
     * 
     * The database given to the constructor will be merged into this one. Only 
     * parameters which are of interest to this class are stored (and the 
     * default ParMooN parameters). Note that this usually does not include 
     * other parameters such as solver parameters. Those are only in the 
     * Solver object.
     */
    ParameterDatabase db;
    
    /** @brief output object */
    DataWriter3D outputWriter;
    
    /** @brief a solver object which will solve the linear system
     * 
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;

    //! @brief An array to store the current defect.
    BlockVector defect_;
    
    ///@brief The norms of residuals from up to 10 previous iterations
    FixedSizeQueue<10, Residuals> old_residuals_;

    //! @brief The initial residual. Stored so that the nonlinear iteration can
    //!        be stopped as soon as a desired reduction is achieved
    double initial_residual_;
    
    /** @brief Errors, held in ready to be accessed from outside the class
     * The array is filled during the function call NSE3D::output()
     * Currently, the errors store the L2 and H1-semi errors of the velocity
     * (errors.at(0) is L2 and errors.at(1) is H1-semi), the L2 error in the 
     * divergence (errors.at(2)), and the pressure (errors.at(3) is L2 and 
     * errors.at(4) is H1-semi).
     */
    std::array<double, int(5)> errors_;

    /** @brief set the velocity and pressure orders
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

    /** @brief Standard constructor of an NSE3D problem.
     *
     * @note The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and
     * functions as well as matrices, solution and right hand side vectors are
     * initialized.
     *
     * @param domain the computational domain to get the grid(s)
     * @param param_db parameters controlling this class
     */
    NSE3D(const TDomain& domain, const ParameterDatabase& param_db);

    /** @brief Standard constructor of an NSE3D problem.
     *
     * @note The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and
     * functions as well as matrices, solution and right hand side vectors are
     * initialized.
     *
     * @param domain the computational domain to get the grid(s)
     * @param param_db parameters controlling this class
     * @param example The example to perform
     */
    NSE3D(const TDomain& domain, const ParameterDatabase& param_db,
          const Example_NSE3D& example);


    /** @brief return a database with all parameters necessary for 
     * Navier--Stokes.
     */
    static ParameterDatabase default_NSE_database();
    
    /**
     * @brief Check whether the program will be working with the
     * current input parameters.
     *
     * ParMooN is work in progress, and so is this class. This method
     * checks the parameters stored in the database and stops execution
     * of the program, if some of these do not match.
     * The method is a little makeshift and the cases caught here are various,
     * but basically it is intended to stop execution of cases with
     * parameter combinations which are not implemented for NSE3D or
     * are currently known to be problematic.
     *
     * This is not yet a guarantee for a functioning program, but is
     * intended to be, someday. Eventually this method and the like
     * will be moved to TDatabase.
     */
    void check_parameters();

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
    void assemble_non_linear_term();
    
    //! Solve the current linear system. Nonlinear loop is outside of this class.
    void solve();

    /**
     * @brief Compute the defect Ax-b, and the residuals and store it all.
     *
     * Updates defect and old_residuals.
     * A is the current matrix, x is the current solution and b is the
     * right hand side. Call this function after assembling the nonlinear
     * matrix with the current solution.
     */
    void compute_residuals();

    /** @brief check if one of the stopping criteria is fulfilled
     * 
     * either converged, maximun number of iterations reached, or slow 
     * convergence
     * 
     * @param iteration_counter current iterate
     *
     * @note For the sequential case, this is the copy-paste NSE2D
     * (with exception of slightly different compute_residual methods).
     */
    bool stop_it(unsigned int iteration_counter);

    /**! Measure errors and draw a nice VTK picture, if requested to do so.
    ! @param i suffix for output file name, -1 means no suffix. */
    void output(int i = -1);

/*******************************************************************************/
    // Declaration of special member functions - delete all but destructor.
    // This problem class will be used to request the whole process of
    // putting up and solving a single flow problem. It is, due to
    // its procedural nature, not supposed to be copied.

    //! Delete copy constructor. No copies allowed.
    NSE3D(const NSE3D&) = delete;

    //! Delete move constructor. No moves allowed.
    NSE3D(NSE3D&&) = delete;

    //! Delete copy assignment operator. No copies allowed.
    NSE3D& operator=(const NSE3D&) = delete;

    //! Default move assignment operator. No moves allowed.
    NSE3D& operator=(NSE3D&&) = delete;

    //! Default destructor. Does most likely cause memory leaks.
    ~NSE3D() = default;
    
/********************************************************************************/
// getters

    const BlockFEMatrix & get_matrix() const
    { return this->systems_.front().matrix_; }
    
    BlockFEMatrix & get_matrix()
    { return this->systems_.front().matrix_; }
    
    const BlockVector & get_rhs() const
    { return this->systems_.front().rhs_; }
    
    BlockVector & get_rhs()
    { return this->systems_.front().rhs_; }
    
    const TFEVectFunct3D& get_velocity() const
    { return this->systems_.front().u_; }

    TFEVectFunct3D& get_velocity()
    { return this->systems_.front().u_; }

    TFEFunction3D *get_velocity_component(int i);

    const TFEFunction3D& get_pressure() const
    { return this->systems_.front().p_; }

    TFEFunction3D& get_pressure()
    { return this->systems_.front().p_; }

    const TFESpace3D& get_velocity_space() const
    { return *this->systems_.front().velocitySpace_; }

    const TFESpace3D& get_pressure_space() const
    { return *this->systems_.front().pressureSpace_; }
    
    const int get_size(){return this->systems_.front().solution_.length();}

    const ParameterDatabase & get_db() const
    { return db; }
    
    /// @brief Get the current residuals  (updated in compute_residuals)
    const Residuals& get_residuals() const;
    /// @brief get the current impuls residual (updated in compute_residuals)
    double get_impuls_residual() const;
    /// @brief get the current mass residual (updated in compute_residuals)
    double get_mass_residual() const;
    /// @brief get the current residual (updated in compute_residuals)
    double get_full_residual() const;
    
    /// @brief return the computed errors (computed in output())
    std::array<double, int(5)> get_errors() const;

 private:
    void output_problem_size_info() const;

};




#endif /* INCLUDE_SYSTEM_NSE3D_H_ */
