#ifndef __SYSTEM_NAVIERSTOKES_H__
#define __SYSTEM_NAVIERSTOKES_H__

#include "BlockFEMatrix.h"
#ifdef __2D__
#include "Example_NSE2D.h"
#include "FEFunction2D.h"
#include "FEVectFunct2D.h"
#include "BoundaryAssembling2D.h"
#else
#include "Example_NSE3D.h"
#include "FEFunction3D.h"
#include "FEVectFunct3D.h"
#include "BoundaryAssembling3D.h"
#endif
#include "BlockVector.h"
#include "ParameterDatabase.h"
#include "Solver.h"
#include "templateNames.h"
#include "DataWriter.h"
#include "Residuals.h"
#include "MainUtilities.h" // FixedSizeQueue

#include <deque>
#include <array>

template<int d>
class NavierStokes
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FEVectFunct = typename Template_names<d>::FEVectFunct;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_NSE = typename Template_names<d>::Example_NSE;
    
    /** @brief Standard constructor of an NSE3D problem.
     *
     * @note The domain must have been refined a couple of times already if you
     * want to use multigrid. On the finest level the finite element spaces and
     * functions as well as matrices, solution and right hand side vectors are
     * initialized.
     *
     * @param domain the computational domain to get the grid(s)
     * @param param_db parameters controlling this class
     */
    NavierStokes(const TDomain& domain, const ParameterDatabase& param_db);
    
    /** @brief Standard constructor of an NSE3D problem.
     *
     * @note The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and
     * functions as well as matrices, solution and right hand side vectors are
     * initialized.
     *
     * @param domain the computational domain to get the grid(s)
     * @param param_db parameters controlling this class
     * @param example The example to use
     */
    NavierStokes(const TDomain& domain, const ParameterDatabase& param_db,
                 Example_NSE example);
    
    /** @brief return a database with all parameters necessary for 
     * Navier--Stokes.
     */
    static ParameterDatabase default_nse_database();
    
    /** @brief Assemble those parts which do not contain nonlinearities
     * i.e. a Stokes or Brinkman problem. When solving a Navier-Stokes problem,
     * then this must be called once before entering the nonlinear loop.
     */
    void assemble_linear_terms();

    /** @brief Assemble terms on boundary edges (2D) or faces (3D).
     */
    void assemble_boundary_terms();

    /** @brief Assemble the nonlinear term. Need not be used when this is
     * a Stokes or Brinkman problem. Call this once per nonlinear iteration if
     * this is a Navier-Stokes problem.
     */
    void assemble_nonlinear_term();
    
    /** @brief Solve the current linear system.
     *
     * Nonlinear loop is outside of this class. 
     */
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

    
    // Declaration of special member functions - delete all but destructor.
    // This problem class will be used to request the whole process of
    // putting up and solving a single flow problem. It is, due to
    // its procedural nature, not supposed to be copied.
    
    // copy constructor is protected
    //! Delete move constructor.
    NavierStokes(NavierStokes&&) = delete;

    //! Delete copy assignment operator.
    NavierStokes& operator=(const NavierStokes&) = delete;

    //! Delete move assignment operator.
    NavierStokes& operator=(NavierStokes&&) = delete;

    //! Destructor.
    ~NavierStokes() = default;

    
    // getters and setters
    const BlockFEMatrix & get_matrix() const
    { return this->systems.front().matrix; }
    BlockFEMatrix & get_matrix()
    { return this->systems.front().matrix; }
    
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    BlockVector & get_rhs()
    { return this->systems.front().rhs; }
    
    const FEVectFunct & get_velocity() const
    { return this->systems.front().u; }
    FEVectFunct & get_velocity()
    { return this->systems.front().u; }
    
    // try not to use this as it is not const
    FEFunction *get_velocity_component(int i);
    
    const FEFunction & get_pressure() const
    { return this->systems.front().p; }
    FEFunction & get_pressure()
    { return this->systems.front().p; }

    const FEFunction & get_exact_pressure() const
    { return this->systems.front().p_exact; }
    FEFunction & get_exact_pressure()
    { return this->systems.front().p_exact; }
    
    const FESpace & get_velocity_space() const
    { return *this->systems.front().velocity_space.get(); }
    const FESpace & get_pressure_space() const
    { return *this->systems.front().pressure_space.get(); }
    
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    BlockVector & get_solution()
    { return this->systems.front().solution; }
    
    const LoopInfo& get_it_solver_info()
    {return solver.get_solver_loop_info();}
    
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    
    const Example_NSE & get_example() const
    { return example; }
    
    void add_to_output(const FEVectFunct* fe_vector_fct)
    { outputWriter.add_fe_vector_function(fe_vector_fct); }
    void add_to_output(const FEFunction* fe_fct)
    { outputWriter.add_fe_function(fe_fct); }
    
    const ParameterDatabase & get_db() const
    { return db; }
    
    /// @brief get the current residuals 
    /// @details updated in NavierStokes<d>::computeNormsOfResiduals which in 
    /// turn is called from NavierStokes<d>::stopIt.
    const Residuals& get_residuals() const;
    /// @brief get the current impuls residual
    /// @details updated in NavierStokes<d>::computeNormsOfResiduals which in
    /// turn is called from NavierStokes<d>::stopIt.
    double get_impuls_residual() const;
    /// @brief get the current mass residual
    /// @details updated in NavierStokes<d>::computeNormsOfResiduals which in
    /// turn is called from NavierStokes<d>::stopIt.
    double get_mass_residual() const;
    /// @brief get the current residual
    /// @details updated in NavierStokes<d>::computeNormsOfResiduals which in
    /// turn is called from NavierStokes<d>::stopIt.
    double get_full_residual() const;
    /// @brief reset the residuals.
    /// use this if you want to use this object again for a second nonlinear 
    /// iteration.
    void reset_residuals();
    
    /// @brief return the computed errors (computed in output())
    std::array<double, 6> get_errors() const;
    
  protected:
    
    /// @brief default copy constructor (useful in derived classes)
    NavierStokes(const NavierStokes&) = default;
    
    /** @brief store a complete system on a particular grid
     *
     * This combines a matrix, rhs, solution, spaces and functions
     * needed to describe one stationary Navier-Stokes flow problem.
     */
    struct System_per_grid
    {
      /** @brief Finite Element space for the velocity */
      std::shared_ptr<FESpace> velocity_space;
      /** @brief Finite Element space for the pressure */
      std::shared_ptr<FESpace> pressure_space;

      /** @brief the system matrix (depends strongly on
       *         TDatabase::ParamDB->NSTYPE)
       *  [ A11  A12  A13  B1T ]
       *  [ A21  A22  A23  B2T ]
       *  [ A31  A32  A33  B3T ]
       *  [ B1   B2   B3   C   ]
       */
      BlockFEMatrix matrix;
      /** @brief the right hand side vector */
      BlockVector rhs;
      /** @brief solution vector with two components. */
      BlockVector solution;
      /** @brief Finite Element function for velocity */
      FEVectFunct u;
      /** @brief Finite Element function for pressure */
      FEFunction p;

      /** @brief Finite Element function for exact pressure (if available)*/
      FEFunction p_exact;
      double *values_exact_p;
      
      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] order the order for the velocity and pressure space
       */
      System_per_grid(const Example_NSE& example, TCollection& coll,
                      std::pair<int, int> order);

      // System_per_grid is not supposed to be copied or moved
      // until underlying classes realize the rule of zero.

      //! copy constructor.
      System_per_grid(const System_per_grid&);

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
    std::deque<System_per_grid> systems;

    /** @brief Definition of the used example. */
    const Example_NSE example;
    
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
    DataWriter<d> outputWriter;
    
    /** @brief a solver object which will solve the linear system
     * 
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;

    //! @brief An array to store the current defect.
    BlockVector defect;
    
    ///@brief The norms of residuals from up to 10 previous iterations
    FixedSizeQueue<10, Residuals> old_residuals;

    /** @brief store the initial residual so that the nonlinear iteration can 
     * be stopped as soon as a desired reduction is achieved.
     */
    double initial_residual;
    
    /** @brief Store the initial rhs norm, so that the nonlinear iteration can
     * be stopped relative to the initial right hand side.
     */
    double initial_rhs_norm;
    
    /** @brief Errors, held in ready to be accessed from outside the class.
     * The array is filled during the function call NavierStokes<d>::output()
     * Currently, the errors store the L2 and H1-semi errors of the velocity
     * (errors.at(0) is L2 and errors.at(2) is H1-semi), the L2 error in the 
     * divergence (errors.at(1)), and the pressure (errors.at(3) is L2 and 
     * errors.at(4) is H1-semi). Finally, at errors[5] there is a natural norm
     * for pressure stabilized methods (such as PSPG or GLS).
     */
    std::array<double, 6> errors;

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
    
    /**
     * @brief Check whether the program will be working with the
     * current input parameters.
     *
     * ParMooN is work in progress, and so is this class. This method
     * checks the parameters stored in the database and stops execution
     * of the program, if some of these do not match.
     * The method is a little makeshift and the cases caught here are various,
     * but basically it is intended to stop execution of cases with
     * parameter combinations which are not implemented or
     * are currently known to be problematic.
     *
     * This is not yet a guarantee for a functioning program, but is
     * intended to be, someday. Eventually this method and the like
     * will be moved to TDatabase.
     */
    void check_parameters();
    
    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;
    
    /** @brief modify computed pressure due to the used nonlinear form
     * 
     * Using the rotational or emac form of the nonlinear term leads to a 
     * pressure which includes an additional term 0.5*|u|^2. This modification 
     * is substracted, so that the resulting pressure is comparable to the ones
     * from other nonlinear forms.
     */
    void adjust_pressure();
};


#endif // __SYSTEM_NAVIERSTOKES_H__
