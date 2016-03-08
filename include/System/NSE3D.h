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

#include <NSE_MultiGrid.h>
#include <FESpace3D.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <FEVectFunct3D.h>
#include <FEFunction3D.h>
#include <Example_NSE3D.h>

#include <vector>
#include <deque>
#include <list>
#include <utility>

#ifdef _MPI
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

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
    struct SystemPerGrid
    {
      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] maxSubDomainPerDof Only in MPI case! The maximal number of
       * processes which share a single dof. The value is calculated in the
       * main program and handed down to the FESpaces. Getting rid of this
       * construction is a TODO .
       */
#ifdef _MPI
      SystemPerGrid(const Example_NSE3D& example,
                    TCollection& coll, std::pair<int, int> order, NSE3D::Matrix type, 
                    int maxSubDomainPerDof);
#else
      SystemPerGrid(const Example_NSE3D& example, TCollection& coll, std::pair<int, int> order, 
                    NSE3D::Matrix type);
#endif

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


      // SystemPerGrid is not supposed to be copied or moved
      // until underlying classes realize the rule of zero.

      //! Delete copy constructor. No copies allowed.
      SystemPerGrid( const SystemPerGrid& ) = delete;

      //! Delete move constructor. No moves allowed.
      SystemPerGrid( SystemPerGrid&& ) = delete;

      //! Delete copy assignment operator. No copies allowed.
      SystemPerGrid& operator=( const SystemPerGrid& ) = delete;

      //! Default move assignment operator. No moves allowed.
      SystemPerGrid& operator=( SystemPerGrid&& ) = delete;

      //! Default destructor. Does most likely cause memory leaks.
      ~SystemPerGrid() = default;
    };

    /** @brief A complete system on each involved grid.
     *
     * Note that the size of this double ended queue is at least one and
     * larger only when a multgrid preconditioner is used.
     */
    std::deque<SystemPerGrid> systems_;

    /** @brief Definition of the used example. */
    const Example_NSE3D& example_;

    /** @brief A multigrid object which is set to nullptr in case it is not
     *         needed.
     */
    std::shared_ptr<TNSE_MultiGrid> multigrid_;
    
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
    
    /** @brief an array to store defect, so that we don't have to reallocate
     *         so often
     */
    BlockVector defect;
    /***/
    double oldResidual;
    
    /** @brief store the initial residual so that the nonlinear iteration can 
     *         be stopped as soon as a desired reduction is achieved
     */
    double initial_residual;
    
    /** @brief
     */
    std::array<double, int(4)> errors;
 public:

    /** @brief The standard constructor, can be used for multigrid and non-multigrid.
     *
     * It is user responsibility to call this constructor with correct TCollections.
     * Here is a "safe" way to do it.
     *
     * When not using multigrid, this must be of length 1,
     * containing a reference to a collection gained e.g. by calling
     * domain.GetCollection(It_Finest, 0) in the main program.
     *
     * In multigrid case, the collections vector must contain TCollections gained
     * by calling domain.GetCollection(It_Finest, 0) in the main program repeatedly,
     * directly after the corresponding refinement step (and after Domain_Crop in MPI
     * case).
     *
     * @param[in] collections A hierarchy of cell collections used for multigrid solve,
     * or just one collection in non-multigrid case. Ordered by fineness of the grid -
     * finest collection first!
     *
     * @param[in] example The example which is to be calculated.
     *
     * @param[in] maxSubDomainPerDof Only in MPI case! the maximal number of
     * processes which share a single dof. The value is calculated in the
     * main program and handed down to the FESpaces. Getting rid of this
     * construction is a TODO .
     */
#ifdef _MPI
    NSE3D(std::list<TCollection* > collections, const Example_NSE3D& example, int maxSubDomainPerDof);
#else

    NSE3D(std::list<TCollection* > collections, const Example_NSE3D& _example);
#endif

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
    static void check_parameters();

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
    
    /** @brief check if one of the stopping criteria is fulfilled
     * 
     * either converged, maximun number of iterations reached, or slow 
     * convergence
     * 
     * @param iteration_counter current iterate
     */
    bool stop_it(unsigned int iteration_counter);

    //! Solve the current linear system. Nonlinear loop will be outside of this class.
    void solve();

    //! Measure errors and draw a nice VTK picture, if requested to do so.
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

/*******************************************************************************/
    /**
     * @brief initialize multigrid levels for different NSTYPE's
     * 
     * @param: level
     * @param: grid to be added
     */
    TNSE_MGLevel* mg_levels(int level, SystemPerGrid& s);
    
    /**
     * @brief multigrid solver
     * preparing the stuff used to call multigrid solver
     */
    void mg_solver();
    
/********************************************************************************/
// getters
   const TFEVectFunct3D & get_velocity() const
    { return this->systems_.front().u_; }
    TFEVectFunct3D & get_velocity()
    { return this->systems_.front().u_; }    
    TFEFunction3D *get_velocity_component(int i)
    { 
      if(i==0) 
        return this->systems_.front().u_.GetComponent(0);
      else if(i==1) 
        return this->systems_.front().u_.GetComponent(1);
      else 
        return this->systems_.front().u_.GetComponent(2);
    }
    TFEFunction3D & get_pressure()
    { return this->systems_.front().p_; }
    const TFESpace3D & get_velocity_space() const
    { return this->systems_.front().velocitySpace_; }
    const TFESpace3D & get_pressure_space() const
    { return this->systems_.front().pressureSpace_; }
    
    const int get_size(){return this->systems_.front().solution_.length();}
    
    /// @brief return the computed errors 
    std::array<double, int(4)> get_errors();

};




#endif /* INCLUDE_SYSTEM_NSE3D_H_ */
