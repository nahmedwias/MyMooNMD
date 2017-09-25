#ifndef TIME_NSE3D_MERGED_H
#define TIME_NSE3D_MERGED_H


/** ************************************************************************
 *
 * @name         Time_NSE3D_Merged
 * @brief        store everything needed to solve a time dependent convection
 *               diffusion reaction problem
 *               Stores matrix, right hand side, FE spaces, FE functions
 *               and the solution vector
 *
 * @author       Naveed Ahmed
 * @History      30.08.2017
**************************************************************************/

#include <FEVectFunct3D.h>

#include <BlockFEMatrix.h>
#include <BlockVector.h>

#include <FESpace3D.h>
#include <Example_TimeNSE3D.h>

#include <Multigrid.h>
#include <Solver.h>

#include <MainUtilities.h>

#include <ParameterDatabase.h>
#include <LocalAssembling3D.h>
#include <TimeDiscretization.h>
#include <Residuals.h>

class Time_NSE3D_Merged
{
  enum class Matrix{Type14, Type1, Type2, Type3, Type4};
protected:
  /** @brief store a complete system on a paticular grid.
     *
     * This combines a matrix, rhs, solution, spaces and functions
     * needed to describe a Time Navier-Stokes problem in 2D
     */
    struct System_per_grid
    {
      /** @brief Finite element space for velocity*/
      TFESpace3D velocity_space;
      /** @brief Finite element space for pressure*/
      TFESpace3D pressure_space;

      /** @brief system matrix
       *  [ A11  A12  B1 ]
       *  [ A21  A22  B2 ]
       *  [ B3   B4   C  ]
      */
      BlockFEMatrix matrix;
      /** @brief mass matrix: this will be the standard mass matrix
       * for the standard Galerkin scheme. However, for the SUPG/RBVMS
       * schems, this includes all the terms which are related to the
       * time derivatives in the fully discrete scheme Eq: (45)
       * Ahmed, Rebollo, John and Rubino (2015)
       *  [ M11  M12  0 ]
       *  [ M21  M22  0 ]
       *  [ 0     0   0 ]
       * This additionaly created due to the struture of the
       * residual based Variational Multiscale Method:
       */
      BlockFEMatrix mass_matrix;
      /** @brief right hand side vector*/
      BlockVector rhs;
      /** @brief solution vector*/
      BlockVector solution;
      /** @brief Finite element function for velocity*/
      TFEVectFunct3D u;
      /** @brief Finite element function for pressure*/
      TFEFunction3D p;
      /** @brief
       * old solution for the computation of the residual
       * that passes as an FEFunction to the local assembling
       * routines
       */
      BlockVector solution_m1;
      BlockVector solution_m2;
      TFEVectFunct3D u_m1;
      TFEFunction3D p_old;
      TFEVectFunct3D u_m2;

      BlockVector combined_old_sols;
      TFEVectFunct3D comb_old_u;

      BlockVector extrapolate_sol;
      TFEVectFunct3D extrapolate_u;

      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] maxSubDomainPerDof Only in MPI case! The maximal number of
       * processes which share a single dof. The value is calculated in the
       * main program and handed down to the FESpaces. Getting rid of this
       * construction is a TODO .
       */
#ifdef _MPI
      System_per_grid(const Example_TimeNSE3D& example, TCollection& coll, std::pair< int, int > order, 
                      Time_NSE3D_Merged::Matrix type, int maxSubDomainPerDof);
#else
      System_per_grid(const Example_TimeNSE3D& ex, TCollection& coll,
                      std::pair<int, int> order, Time_NSE3D_Merged::Matrix type);
#endif
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

    /** @brief A complete system on each involved grid.
     *
     * Note that the size of this double ended queue is at least one and
     * larger only when a multgrid preconditioner is used. In the multigrid
     * case, the systems_ are ordered from front to back 
     * (finest to coarses)
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

    /** @brief solution vector from previous time step (on finest mesh)*/
    BlockVector old_solution_;
    
    BlockVector old_rhs_w_old_sol;
    /** @brief constructs a solution vector extrapolated from previous steps
     * Currently, it is used for IMEX-Scheme: 2u(t-1)-u(t-2). */
    BlockVector extrapolated_solution_;

    /** @brief old time step length used to scale the pressure blocks */
    double oldtau_;
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
    
    /** @brief projection space used for VMS method*/
    std::shared_ptr<TFESpace3D> projection_space_;
    
    /** @brief finite element function for vms projection*/
    // can we rename it to large scales?? also check BlockVector!! currently just vector
    std::vector<double> vms_small_resolved_scales; 
    std::shared_ptr<TFEVectFunct3D> vms_small_resolved_scales_fefct;
    /** matrices for turbulence model*/
    std::array<std::shared_ptr<FEMatrix>, int(7)> matrices_for_turb_mod;
    
    // piecewise constant space containing the labels of the local projection space
    std::shared_ptr<TFESpace3D> label_for_local_projection_space_;
    // vector used for indicating local projection space 
    std::vector<double> label_for_local_projection;
    // finite element function for local projection space
    std::shared_ptr<TFEFunction3D> label_for_local_projection_fefct;
    // prepare spaces and matrices for the vms projection
    void set_matrices_vms(TCollection *coll_);
    
    ///
    BlockVector rhs_from_time_disc;
public:
  #ifdef _MPI
  Time_NSE3D_Merged(std::list<TCollection* > collections_, const ParameterDatabase& param_db, 
               const Example_TimeNSE3D& example, int maxSubDomainPerDof);
#else
    Time_NSE3D_Merged(std::list<TCollection* > collections_, const ParameterDatabase& param_db, 
               const Example_TimeNSE3D& example);
#endif
    
    /** @brief check parameters in database
    *
    * This functions checks if the parameters in the database are meaningful.
    * If some parameters are set to unsupported values, an error occurs and
    * throws an exception.
    */
    void check_and_setparameters();
    
    /** @brief 
     */
    TimeDiscretization time_stepping_scheme;
    
    /** @brief Assemble all the matrices and rhs before the time iterations
    *
    * This includes the assembling of the Stiffness matrix, the Mass matrix,
    * the additional matrix K in case of SUPG stabilization, and rhs.
    * This assembling occurs just once, before entering any loop. It assembles
    * linear terms only.
    */
    void assemble_initial_time();
    
    /** @brief Assemble the rhs only
     * TODO: update according to the 2D case
    * 1. Assembling the right hand side and matrices 
    * 2. Scaling of the B-Blocks due to time stepping
    * This function will prepare the right hand side during the time
    * discretization but should be outside the nonlinear loop.
    */
    void assemble_matrices_rhs(unsigned int it_counter);
    
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
    
    /** @brief Solve the current linear system. Nonlinear
     * loop is outside of this class.
    */
    void solve();

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
    void output(int m, int &image);
    
    /** */
    bool imex_scheme(bool print_info);
    
    /* *******************************************************************************/

    //// getters

    /// Get the velocity space.
    const TFESpace3D     & get_velocity_space() const
    { return this->systems_.front().velocity_space; }

    /// Get the pressure space.
    const TFESpace3D     & get_pressure_space() const
    { return this->systems_.front().pressure_space; }

    const TFEVectFunct3D& get_velocity() const
    { return this->systems_.front().u; }

    TFEVectFunct3D& get_velocity()
    { return this->systems_.front().u; }

    const TFEFunction3D& get_pressure() const
    { return this->systems_.front().p; }

    TFEFunction3D& get_pressure()
    { return this->systems_.front().p; }

    /// Get number of degrees of freedom.
    const int get_size() const
    { return this->systems_.front().solution.length(); }

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
    
    private:
  /// this routines wraps up the call to Assemble3D  
  void call_assembling_routine(Time_NSE3D_Merged::System_per_grid& s, 
                               LocalAssembling3D_type la_type);
  
  /// @brief set the spaces depending on disc types
  void set_arrays(Time_NSE3D_Merged::System_per_grid& s, 
        std::vector<const TFESpace3D*> &spaces, std::vector< const TFESpace3D* >& spaces_rhs,
        std::vector< TFEFunction3D*> &functions, std::vector<BoundCondFunct3D*> &bc, 
                               std::vector<BoundValueFunct3D*> &bv);
  
  /**
   * This will initialize the matrices depending on DISCTYPE and 
   * NSTYPE.
   * @param s system on a particular grid
   * @param type local assembling type
   * @param sqMatrices square matrices for different DISCTYPE and different 
   * NSTYPE
   * @param reMatrirces rectangular matrices 
   */
  void prepare_matrices_rhs(Time_NSE3D_Merged::System_per_grid& s, 
         LocalAssembling3D_type la_type, std::vector< TSquareMatrix3D* >& sqMatrices, 
         std::vector< TMatrix3D* >& rectMatrices, std::vector< double* >& rhs_array);
  /** 
   * @brief restrict the function to on every grid
  * nonliear assembling requires an approximate velocity
  * on every grid
  */
  void restrict_function();
};

#endif // TIME_NSE3D_MERGED_H
