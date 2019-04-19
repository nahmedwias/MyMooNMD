#ifndef VOLUMEOFFLUID_H
#define VOLUMEOFFLUID_H

#include "vector"
#include "../include/General/templateNames.h"
#include "../include/Matrix/BlockFEMatrix.h"
#include "../include/Matrix/BlockVector.h"
#include "../../include/General/ParameterDatabase.h"
#include "../../include/Geometry/Domain.h"
#include "Example_TnseTcd2D.h"
#include <deque>
#include <array>
#include "MainUtilities.h" // FixedSizeQueue
#include "Residuals.h"
#include "DataWriter.h"
#include "Solver.h"
#include "TimeDiscretizations.h"
#include "LocalAssembling.h"
#include "LinesEval.h"


#ifdef __2D__
#include "../include/FE/FEFunction2D.h"
#include "../include/FE/FEVectFunct2D.h"
#else
#include "../include/FE/FEFunction3D.h"
#include "../include/FE/FEVectFunct3D.h"
#endif
/**
 * @todo write docs
 */

template <int d>
class MultiPhase
{
public:
  using FEFunction = typename Template_names<d>::FEFunction;
  using FEVectFunct = typename Template_names<d>::FEVectFunct;
  using FESpace = typename Template_names<d>::FESpace;
  
public:
  /** @brief see the other constructor */
  MultiPhase(const TDomain& domain, const ParameterDatabase& param_db);
  /** @brief constructor
     *
     * The domain must have been refined a couple of times already. On the
     * finest level the finite element spaces and functions as well as
     * matrices, solution and right hand side vectors are initialized.
     * 
     * @param domain the computational domain to get the grid(s)
     * @param param_db parameters controlling this class
     * @param example The example to use
     */
  MultiPhase(const TDomain& domain, const ParameterDatabase& param_db,
                     const Example_TnseTcd2D& ex);
  
  static ParameterDatabase default_coupled_database();
  /// @brief get the current residuals 
  /// @details updated in MultiPhase<d>::compute_residuals() which in 
  /// turn is called from MultiPhase<d>::stop_it().
  const Residuals& get_residuals() const;
  /// @brief get the current impuls residual
  /// @details updated in MultiPhase<d>::compute_residuals() which in 
  /// turn is called from MultiPhase<d>::stop_it().
  double get_impuls_residual() const;
  /// @brief get the current mass residual
  /// @details updated in MultiPhase<d>::compute_residuals() which in 
  /// turn is called from MultiPhase<d>::stop_it().
  double get_mass_residual() const;
  /// @brief get the current residual
  /// @details updated in MultiPhase<d>::compute_residuals() which in 
  /// turn is called from MultiPhase<d>::stop_it().
  double get_full_residual() const;
  /// @brief reset the residuals.
  /// use this if you want to use this object again for a second nonlinear 
  /// iteration.
  void reset_residuals();

  /** @brief Assemble all the matrices before the time iterations
   * 
   * This includes the assembling of: Stiff_matrix, Mass_Matrix, 
   * (additional matrixK in case of SUPG stabilization), rhs
   */
  void assemble_initial_time();
  
  
  /** @brief Assemble those parts which do not contain nonlinearities

   * like a Stokes or Brinkman problem. When solving a Navier-Stokes problem,

   * then this must be called once before entering the nonlinear loop.

   */
  void assemble_linear_terms();
  
  /** @brief assemble the matrices
   * this function will assemble the stiffness matrix and rhs
   * In addition the system matrix and the rhs which passes to the solver 
   * are also prepared within the function
   */
  void assemble();
  
  /** @brief This will assemble the right-hand side, will prepare the 
   * right-hand side using the class TimeDiscretization object. 
   * The pressure blocks are scaled properly and the nonlinear matrices
   * are assembled again. The system matrix is assemble using the 
   * object of class TimeDiscretization and modify the matrices if needed.
   */
  void assemble_matrices_rhs(unsigned int it_counter);
 
  
  /** @brief check if the semi-implicit scheme is used */
  bool imex_scheme();
  
  /**
   * @brief Compute the defect Ax-b and store its norm.
   *
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
   * @param it_counter current iterate
   */
  bool stop_it(unsigned int it_counter);

  /** @brief solve the system
   */
  void solve();
  
  /** @brief measure errors and write solution
   * 
   */
  void output();

   // getters and setters
  
   const FEVectFunct & get_velocity() const
  { return this->systems.front().u; }
  FEVectFunct & get_velocity()
  { return this->systems.front().u; }
  
    const FEFunction & get_pressure() const
  { return this->systems.front().p; }
  FEFunction & get_pressure()
  { return this->systems.front().p; }

private:
  struct System_per_grid
  {
    /** @brief Finite element space for velocity*/
    std::shared_ptr<FESpace> velocity_space;
    /** @brief Finite element space for pressure*/
    std::shared_ptr<FESpace> pressure_space;
    /** @brief Finite element space for pressure*/
    std::shared_ptr<FESpace> temp_space;
    
    /** @brief system matrix                     |    [ A11  A12  A13  B1T ]
       *                       [ A11  A12  B1 ]    |    [ A21  A22  A23  B2T ]
       *                       [ A21  A22  B2 ]    |    [ A31  A32  A33  B3T ]
       *                       [ B3   B4   C  ]    |    [ B1   B2   B3   C   ]
      */
    BlockFEMatrix matrix_NS;
    /***/
    BlockFEMatrix matrix_temp;
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
    BlockFEMatrix mass_matrixNS;
    
    /** @brief Stiffness Matrix for Convection Diffusion*/
    BlockFEMatrix stiffness_matrixtemp;
    /****/
    BlockFEMatrix mass_matrixtemp;
    
    /** @brief right hand side vector*/
    BlockVector rhs_NS;
    /** @brief solution vector*/
    BlockVector solution_NS;
    /** @brief Finite element function for velocity*/
    FEVectFunct u;
    /** @brief Finite element function for pressure*/
    FEFunction p;
    
    /** @brief right hand side vector for Convection Diffusion*/
    BlockVector rhs_temp;
    
   /** @brief solution vector for Convection Diffusion*/
    BlockVector solution_temp;
    
    /** @brief Finite element function for temprature*/
    FEFunction temp;
    /** @brief
      * old solution for the computation of the residual
      * that passes as an FEFunction to the local assembling
      * routines
      * 
      * @todo why don't we store a vector of size 
      * time_stepping_scheme.n_old_solutions() here?
    */
    BlockVector solution_temp_m1;
    FEFunction temp_u_m1;
    BlockVector solution_temp_m2;
    FEFunction temp_u_m2;
    
    
    /** @brief
      * old solution for the computation of the residual
      * that passes as an FEFunction to the local assembling
      * routines
      * 
      * @todo why don't we store a vector of size 
      * time_stepping_scheme.n_old_solutions() here?
      */
    BlockVector solution_NS_m1;
    FEFunction u_m1;
    FEFunction p_m1;
    BlockVector solution_NS_m2;
    FEFunction u_m2;
    FEFunction p_m2;      
    
    BlockVector time_avg_sol;
    FEVectFunct u_time_avg;
    FEFunction  p_time_avg;

    BlockVector combined_old_sols_NS;
    FEVectFunct comb_old_u;

    BlockVector extrapolate_sol_NS;
    FEVectFunct extrapolate_u;     
    /** @brief constructor*/
    System_per_grid(const Example_TnseTcd2D& example, TCollection& coll,
                      std::tuple<int,int,int> order);    
    
      /** @brief copy constructor*/
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
  
  /** @brief a local parameter database which controls this class
   * 
   *
   * The database given to the constructor will be merged into this one. Only
   * parameters which are of interest to this class are stored (and the
   * default ParMooN parameters). Note that this usually does not include
   * other parameters such as solver parameters. Those are only in the
   * Solver object.
   */
  ParameterDatabase db;
   
  /** @brief a complete system on each grid
  *
  * Note that the size of this deque is at least one and larger only in case
  * of multigrid.
  */
  std::deque<System_per_grid> systems;
  
   
  /** @brief Definition of the used example */
  Example_TnseTcd2D example;
    
  /** @brief an array to store defect, so that we don't have to reallocate
   *         so often
   */
  BlockVector defect;

  ///@brief The norms of residuals from up to 10 previous iterations
  FixedSizeQueue<10, Residuals> old_residuals;
   /** @brief store the initial residual so that the nonlinear iteration can
   *         be stopped as soon as a desired reduction is achieved
   */
  double initial_residual;
  

  /** @brief a solver object which will solve the linear system
   * 
   * Storing it means that for a direct solver we also store the factorization
   * which is usually not necessary.
   */
  Solver<BlockFEMatrix, BlockVector> solver;
  

  /// this is for setting the discretization type globally, since the DISCTYPE 
  /// is removed fromt the global database but it's not simple/possible to use
  /// different space discretization in the LocalAssembling2D
  int space_disc_global;
  
 /** @brief old right hand side vectior 
  * this will be used to save the right hand side from the 
  * previous time step that will be used for different 
  * time stepping schemes
  */
  BlockVector old_rhs_temp;

  BlockVector old_rhs_NS;
  
  
  static constexpr int n_errors = 5;
 /** @brief store the errors to compute accumulated error norms */
  std::array<double, n_errors> errors_temp;
    
 /** @brief store the errors to compute accumulated error norms */
  std::array<double, n_errors> errors_NS;
    
 /** @brief output object */
  DataWriter<d> outputWriter;
    
  /// @brief time stepping scheme object to access everything
  TimeDiscretization time_stepping_scheme;
  
  /// @brief is the mass matrix and right hand side is solution
  /// dependent: for example the SUPG method
  bool is_rhs_and_mass_matrix_nonlinear;
  /// @brief system right hand side which passes to the solver
  BlockVector rhs_temp_from_time_disc;
  
    /// @brief system right hand side which passes to the solver
  BlockVector rhs_NS_from_time_disc;
  
    
  /** @brief set parameters in database
   *
   * This functions checks if the parameters in the database are meaningful
   * and resets them otherwise. The hope is that after calling this function
   * this class is fully functional.
   *
   * If some parameters are set to unsupported values, an error occurs and
   * throws an exception.
   */
  void check_and_set_parameters();
  
  /** @brief write some information (number of cells, dofs, ...) */
  void output_problem_size_info() const;
  
     
  /** @brief get velocity and pressure space*/
  void get_velocity_pressure_orders(
    std::tuple<int, int, int> &velocity_pressure_orders);
    
  /**  @brief set parameters in database
     * Scramble together the parameters which Assemble2D/Assemble3D needs and
     * call it. Is only put here to keep the code of assemble() slender.
     * assemble() should take care of the right choice of the LocalAssembling
     * object and whether it fits the system matrix' block structure
     * @param s The sytem where rhs, mass and stiffness matrices are to be assembled.
     * @param la_stiff The local assembling object of choice.
     * @param assemble_both If true, both stiffness (+rhs) and mass matrix are
     * assembled, if false only stiffness matrix and rhs.
     */
  void call_assembling_routine_ConvDiff(MultiPhase::System_per_grid& s,
                                 LocalAssembling<d>& la_stiff,
                                 bool assemble_both);  
  
  
  /// @brief this routines wraps up the call to Assemble2D
  void call_assembling_routine_NS(MultiPhase<d>::System_per_grid& s, 
                                  LocalAssembling_type type);
   
  /// @brief set the spaces depending on disc types
  void set_arrays_NS(MultiPhase<d>::System_per_grid& s,
                  std::vector<const FESpace*> &spaces,
                  std::vector<const FESpace*>& spaces_rhs_NS,
                  std::vector<FEFunction*>& functions);


  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD; 
  
   /// @brief set the matrices and right hand side depending on the
  /// assemling routines, for navier stokes part
  void set_matrices_rhs_NS(MultiPhase<d>::System_per_grid& s,
                        LocalAssembling_type type,
                        std::vector<SquareMatrixD*> &sqMat,
                        std::vector<MatrixD*> &reMat,
                        std::vector<double*> &rhs_array_NS);
  
 /// @brief Special case for the SUPG and RBVMS method. The right-hand side is 
 /// nonlinear and need a complete re-assembling before passing to the solver 
 /// for solving.
  void assemble_rhs_nonlinear();
  

 /// @brief restrict the function to every (coarse) grid
 /// nonlinear assembling requires an approximate velocity on every grid
  void restrict_function(); 
 
 /// update matrices for local projection stabilization
  void update_matrices_lps(MultiPhase<d>::System_per_grid& s);
  
    /** @brief modify matrices according to the Slip type boundary
     * conditions
     * If the mass matrix (also BT's are independent of solution )
     * is independent of solution then it only needs to modify only
     * once:
     * NOTE: mass matrix and BT's are solution dependent for
     * residual-VMS, and SUPG case.
     * Nonlinear matrices needs to be modify within each
     * time step and also within non-linear iteration
     */
  void modify_slip_bc(bool BT_Mass = false, bool slip_A_nl = false); 

 /**
     * Apply an algebraic flux correction scheme to the assembled matrix.
     * Should be called within the assemble routine, after the assembling
     * of pure mass and stiffness matrix and right hand side
     * has been performed with the INTERNAL_FULL_MATRIX_STRUCTURE switch on.
     *
     *
     * Which afc algorithm is performed is determined by switching over
     * the algebraic_flux_correction parameter of the database. (So far only
     * fem-fct-cn is enabled.)
     */
  void do_algebraic_flux_correction();
    
    /** @brief LineEval object to store the lines where evalation has to be done
     */
  LinesEval<d> Lines;
    
  
      /** @brief calculate the time averaging of the solution in time_avg_sol */
  void time_averaging();  

};

#endif // VOLUMEOFFLUID_H
