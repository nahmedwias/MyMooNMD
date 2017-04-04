/** ************************************************************************
 * 
 * @name         Time_NSE2D
 * @brief        store everything needed to solve a time dependent convection 
 *               diffusion reaction problem
 *               Stores matrix, right hand side, FE spaces, FE functions 
 *               and the solution vector
 *               
 * @author       Naveed Ahmed
 * @History      16.09.2015
**************************************************************************/

#ifndef __TIME_NSE2D__
#define __TIME_NSE2D__

#include <FEVectFunct2D.h>

#include <BlockFEMatrix.h>
#include <BlockVector.h>

#include <FESpace2D.h>
#include <Example_TimeNSE2D.h>
#include <LocalAssembling2D.h>

#include <Multigrid.h>
#include <Solver.h>

#include <MainUtilities.h>

#include <ParameterDatabase.h>
#include <PostProcessing2D.h>

#include <vector>
#include <deque>
#include <utility>
#include <array>

class Time_NSE2D
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
      TFESpace2D velocity_space;
      /** @brief Finite element space for pressure*/
      TFESpace2D pressure_space;
      
      /** @brief system matrix
       *  [ A11  A12  B1 ]
       *  [ A21  A22  B2 ]
       *  [ B3   B4   C  ]
      */
      BlockFEMatrix matrix;
      /** @brief Mass matrix
       *  [ M11  0    0 ]
       *  [ 0    M22  0 ]
       *  [ 0    0    0 ]
       */
      BlockFEMatrix Mass_Matrix;
      /** @brief right hand side vector*/
      BlockVector rhs;
      /** @brief solution vector*/
      BlockVector solution;
      /** @brief Finite element function for velocity*/
      TFEVectFunct2D u;
      /** @brief Finite element function for pressure*/
      TFEFunction2D p;

      /** @brief constructor*/
      System_per_grid(const Example_TimeNSE2D& example, TCollection& coll, 
                      std::pair<int,int> order, Time_NSE2D::Matrix type);
    };
    
    /* FEFunction equal to property fields, construction
     * based on the above BlockVectors
     */
    TFEFunction2D* rho_fefunct = nullptr;
    TFEFunction2D* mu_fefunct = nullptr;

    /** @brief a local parameter database which controls this class
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * default ParMooN parameters). Note that this usually does not include
     * other parameters such as solver parameters. Those are only in the
     * Solver object.
     */
    ParameterDatabase db;

    /** @brief class for output handling */
    PostProcessing2D outputWriter;
    
    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    Example_TimeNSE2D example;
    
    /** @brief a solver object which will solve the linear system
     * 
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;
    
    /** @brief an array to store defect, so that we don't have to reallocate
     *         so often
     */
    BlockVector defect;

    /**
    * @brief store the square root of the residual from previous iteration
    */
    double oldResidual;

    /** @brief store the initial residual so that the nonlinear iteration can
    *         be stopped as soon as a desired reduction is achieved
    */
    double initial_residual;

    /** @brief store errors  */
    std::vector<double> errors;

    /** @brief right hand side vector from previous time step (on finest mesh)*/
    BlockVector old_rhs;
      
    BlockVector old_solution;
    
    /** old time step length used to scale the pressure blocks*/
    double oldtau;
      
    /** @brief parameter to store the space_discretization_type
     * used in TNSE2D. It can be used throughout the methods,
     * especially in the calls to LocalAssembling object.
     * Its main purpose is to convert the "space_discretization_type"
     * into an "int", so that it matches the usual DISCTYPE: GALERKIN,
     * SUPG, SMAGORINSKY, VMS_PROJECTION,...etc
     * The default is 1 (=GALERKIN) but the real value is set
     * in the call to set_parameters().
     */
    int disctype = 1;

    /*
     * @brief this detects if the rho and mu coefficients are given
     *  as fe_functions in the local assembling objects, instead of
     *  taking eps=1/reynolds_number from LinCoeffs from example
     */
    bool with_variable_fluid_properties = 0;

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

    /** @brief get velocity and pressure space*/
    void get_velocity_pressure_orders(std::pair <int,int> &velocity_pressure_orders);
    
    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;

    /** Below all the members needed for Projection-Based VMS */
    /** @brief finite element function for vms projection*/
    // can we rename it to large scales?? also check BlockVector!! currently just vector
    std::shared_ptr<TFESpace2D> projection_space_;
    std::vector<double> vms_small_resolved_scales;
    std::shared_ptr<TFEVectFunct2D> vms_small_resolved_scales_fefct;
    // piecewise constant space containing the labels of the local projection space
    std::shared_ptr<TFESpace2D> label_for_local_projection_space_;
    std::vector<double> label_for_local_projection;
    std::shared_ptr<TFEFunction2D> label_for_local_projection_fefct;
    /** matrices for turbulence model
     *
     * (  A11  A12  B1T  G11tilde   G12tilde      0       )   ( u1  )
     * (  A21  A22  B2T     0       G22tilde   G24tilde   )   ( u2  )
     * (  B1   B2    0      0          0          0       ) . (  p  )
     * (  G11   0    0      M          0          0       )   ( g11 )
     * (  G21  G22   0      0         M/2         0       )   ( g12 )
     * (   0   G42   0      0          0          M       )   ( g22 )
     *
     *  It remains 5 matrices due to the following relations:
     *  G21 = G42/2 ;  G12tilde = G24tilde/2 ;
     *  G22 = G11/2 ;  G22tilde = G11tilde/2 ; and the fifth matrix is M
     *
     * */
    std::array<std::shared_ptr<FEMatrix>, int(5)> matrices_for_turb_mod;

  public:

    /** @brief constructor
     * This constructor calls the other constructor creating an Example_TimeNSE2D
     * object. 
     */
    Time_NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
               int reference_id = -4711);

    /** @brief constructor 
     * 
     * The domain must have been refined a couple of times already. On the 
     * finest level the finite element spaces and functions as well as 
     * matrices, solution and right hand side vectors are initialized. 
     * 
     * The reference_id can be used if only the cells with the give reference_id
     * should be used. The default implies all cells.
     */
    Time_NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
               const Example_TimeNSE2D& ex, int reference_id = -4711);
    
    /** @brief Assemble all the matrices and rhs before the time iterations
     * 
     * This includes the assembling of: Stiff_matrix, Mass_Matrix, 
     * (additional matrixK in case of SUPG stabilization), rhs
     */
    void assemble_initial_time();
    
    /** @brief 
     * 1. assembling the right hand side
     * 2. scaling of the B-blocks due to time stepping
     * this function will prepare the right hand side during the time 
     * discretization but should be outside the nonlinear loop
     */
    void assemble_rhs();
    
    /** @brief Assemble the system matrix
     * This function will prepare the system which will be 
     * used for solvers
     */
    void assemble_system();
    
    /** descale matrices
     * This function will descale all A-blocks which were scaled
     * during the function call Time_NSE2D::assemble_system():
     */
    void deScaleMatrices();

    /** @brief assemble nonlinear term
     * 
     * The matrix blocks to which the nonlinear term contributes are reset to 
     * zero and then completely reassembled, including the linear and nonlinear
     * terms. If this->assemble() has been called before, the matrix is now set
     * up correctly. 
     */
    void assemble_nonlinear_term(unsigned int it_counter = 0);
    
    /** @brief solve the system */
    void solve();

    /** @brief check if one of the stopping criteria is fulfilled
     * 
     * either converged, maximun number of iterations reached, or slow 
     * convergence
     * 
     * @param it_counter current iterate
     */
    bool stopIte(unsigned int it_counter);
    
    /** @brief
     * compute errors and write solution
     */
    void output(int m);
    
    // getters and setters
    /*const BlockMatrixNSE2D & get_matrix() const
    { return this->systems.front().matrix; }*/
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }

    const TFEVectFunct2D & get_velocity() const
    { return this->systems.front().u; }
    // try not to use this as it is not const
    TFEFunction2D *get_velocity_component(int i)
    { return (i==0) ? this->systems.front().u.GetComponent(0)
                    : this->systems.front().u.GetComponent(1); }
    const TFEFunction2D & get_pressure() const
    { return this->systems.front().p; }
    TFEFunction2D & get_pressure()
    { return this->systems.front().p; }
    const TFESpace2D & get_velocity_space() const
    { return this->systems.front().velocity_space; }
    const TFESpace2D & get_pressure_space() const
    { return this->systems.front().pressure_space; }
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    const Example_TimeNSE2D & get_example() const
    { return example; }
    const ParameterDatabase & get_db() const
    { return db; }
    /// @brief return the computed errors at each discre time point
    std::array<double, int(6)> get_errors() const;







    /** ***************BELOW THIS LINE, USER SPECIFIC CODE ********/
  public:

    /*
     * @brief this sets the boolean with_variable_fluid_properties
     */
    void set_bool_variable_properties(bool rho_mu_variable)
    { this->with_variable_fluid_properties = rho_mu_variable; }

    /*
     * @brief this sets the rho and mu fefunctions coming from VOF
     */
    void set_rho_mu_fefunct(TFEFunction2D* rho, TFEFunction2D* mu)
    { this->rho_fefunct = rho; this->mu_fefunct = mu; }

    /**@brief this is a call to assemble2Dslipbc(), which sets up
     * everything needed for "Slip with friction, Penetration
     * with resistance" Boundary conditions in RHS and matrices
     */
    void apply_slip_penetration_bc(bool change_A_offdiagonal_blocks = false,
                                   bool change_B_Mass_blocks = false);

    void assemble_massmatrix_withfields(TFEFunction2D* rho_field=nullptr);

// ======================================================================
    // The following members are used for IMEX
/** @brief This returns the number of the current time step.
     * This counter is set at 0 before the time loop and is incremented at each
     * time step (but not at each sub-step) in the main program.
     * It can be useful to give info to the members of the class. It is for example
     * used in IMEX scheme to detect when we passed 2 time steps, so that we
     * are guaranteed to have saved both old_solution_ and old_solution2_ correctly.    */
    int current_step_;

    /** @brief Construct the extrapolated solution.
     * At the moment, only IMEX is implemented. */
    void construct_extrapolated_solution();

    bool imex_scheme(bool print_info);

    /** @brief constructs a solution vector extrapolated from previous steps
     * Currently, it is used for IMEX-Scheme: 2u(t-1)-u(t-2). */
  protected:
    BlockVector extrapolated_solution_;
// ======================================================================


  private:
    /* @brief Set the fe_functions to prepare
     * the local assembling object. It is basically to improve
     * the readability of the code.   */
    void prepare_fefunc_for_localassembling(Time_NSE2D::System_per_grid& s,
                                            std::vector<TFEFunction2D*> &fe_functions);

    /* @brief Set the spaces and matrices to prepare
     * the call to Assemble. It is basically to improve
     * the readability of the code.   */
    void prepare_spaces_and_matrices_for_assemble(Time_NSE2D::System_per_grid& s,
                                                  LocalAssembling2D_type type,
                                                  std::vector<const TFESpace2D*> &spaces,
                                                  std::vector<const TFESpace2D*> &spaces_rhs,
                                                  std::vector<TSquareMatrix2D*> &sqMat,
                                                  std::vector<TMatrix2D*> &rectMat,
                                                  std::vector<double*> &rhs);

};

#endif // __TIME_NSE2D__
