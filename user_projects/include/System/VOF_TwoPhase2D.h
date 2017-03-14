/** ************************************************************************
 * 
 * @name         VOF_TwoPhase2D
 * @brief        store everything needed to solve a two-phase problem
 *               modelled with Volume Of Fluid
 *               
 * @author       Najib Alia
 * @History      27.02.2017
**************************************************************************/

#ifndef __VOF_TwoPhase2D__
#define __VOF_TwoPhase2D__

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <FESpace2D.h>
#include <FEFunction2D.h>
#include <Time_NSE2D.h>
#include <Time_CD2D.h>
#include <Example_TimeNSE2D.h>
#include <Example_TimeCD2D.h>

#include <MainUtilities.h>
#include <ParameterDatabase.h>
#include <vector>
#include <deque>
#include <utility>
#include <array>

class VOF_TwoPhase2D
{
  public:
    Example_TimeNSE2D example_tnse2d_;
    Time_NSE2D tnse2d_;
    Example_TimeCD2D example_tcd2d_;
    Time_CD2D phaseconvection2d_;

    /* example number of vof=example of tnse=tcd */
    int example_number_;
    /* rhol = constant density of liquid phase
     * default value is 1 */
    double rhol_ = 1;
    /* mul = constant dyn. visco of liquid phase
     * default value is 1 */
    double mul_  = 1;
    /* rhog = constant density of gas phase
     * default value is 0 */
    double rhog_ = 0;
    /* mug = constant dyn. visco of gas phase
     * default value is 0 */
    double mug_  = 0;

    /* Boolean parameters which activates and
     * controls different features of this class
     */
    bool tnse_variable_fluid_ = false;
    bool solve_convection_    = false;
    bool nse2cd_coupling_     = false;
    bool cd2nse_coupling_     = false;

    /* Vector equal to property fields at the nodes */
    BlockVector rho_vector_;
    BlockVector mu_vector_;
    BlockVector unity_vector_; // same structure and equals to 1

    /* FEFunction equal to property fields, construction
     * based on the above BlockVectors
     */
    TFEFunction2D rho_fefunction_;
    TFEFunction2D mu_fefunction_;

  public:
    /** @brief constructor*/
    VOF_TwoPhase2D(const TDomain& domain,
                   const ParameterDatabase& param_db_tnse,
                   const ParameterDatabase& param_db_tcd
                   );

    /*************************************************************/
   /**
    * Special member functions mostly deleted
    * ...needs to be optimized
    */
   //! Delete copy constructor.
    VOF_TwoPhase2D(const VOF_TwoPhase2D&) = delete;

   //! Delete move constructor.
    VOF_TwoPhase2D(VOF_TwoPhase2D&&) = delete;

   //! Delete copy assignment operator.
    VOF_TwoPhase2D& operator=(const VOF_TwoPhase2D&) = delete;

   //! Delete move assignment operator.
    VOF_TwoPhase2D& operator=(VOF_TwoPhase2D&&) = delete;

   //! Default destructor. Most likely causes memory leaks.
   ~VOF_TwoPhase2D() = default;

   /*************************************************************/

   /* Check that the input parameters are consistent,
    * and correct/throw errors if not the case
    */
   void manage_example_parameters();

   /* Update the BlockVectors rho and mu with the phase fraction
    * vector, via the equation: rho = rhol.phi + rhog.(1-phi)
    * idem for mu
    */
   void update_field_vectors();

   /* Write the vectors in a file for output */
   void output_vectors(std::string filename_phi,
                       std::string filename_rho,
                       std::string filename_mu);

   /* Print some info, mostly useful after the constructor */
   void output_initial_info();

};

#endif // __VOF_TwoPhase2D__

//
//    /** @brief store a complete system on a paticular grid.
//     *
//     * This combines a matrix, rhs, solution, spaces and functions
//     * needed to describe a Time Navier-Stokes problem in 2D
//     */
//    struct System_per_grid
//    {
//      /** @brief Finite element space for velocity*/
//      TFESpace2D velocity_space;
//      /** @brief Finite element space for pressure*/
//      TFESpace2D pressure_space;
//
//      /** @brief system matrix
//       *  [ A11  A12  B1 ]
//       *  [ A21  A22  B2 ]
//       *  [ B3   B4   C  ]
//      */
//      BlockFEMatrix matrix;
//      /** @brief Mass matrix
//       *  [ M11  0    0 ]
//       *  [ 0    M22  0 ]
//       *  [ 0    0    0 ]
//       */
//      BlockFEMatrix Mass_Matrix;
//      /** @brief right hand side vector*/
//      BlockVector rhs;
//      /** @brief solution vector*/
//      BlockVector solution;
//      /** @brief Finite element function for velocity*/
//      TFEVectFunct2D u;
//      /** @brief Finite element function for pressure*/
//      TFEFunction2D p;
//
//      /** @brief constructor*/
//      System_per_grid(const Example_TimeNSE2D& example, TCollection& coll,
//                      std::pair<int,int> order, Time_NSE2D::Matrix type);
//    };
//
//    /** @brief a local parameter database which controls this class
//     *
//     * The database given to the constructor will be merged into this one. Only
//     * parameters which are of interest to this class are stored (and the
//     * default ParMooN parameters). Note that this usually does not include
//     * other parameters such as solver parameters. Those are only in the
//     * Solver object.
//     */
//    ParameterDatabase db;
//
//    /** @brief class for output handling */
//    PostProcessing2D outputWriter;
//
//    /** @brief a complete system on each grid
//     *
//     * Note that the size of this deque is at least one and larger only in case
//     * of multigrid.
//     */
//    std::deque<System_per_grid> systems;
//
//    /** @brief Definition of the used example */
//    Example_TimeNSE2D example;
//
//    /** @brief a solver object which will solve the linear system
//     *
//     * Storing it means that for a direct solver we also store the factorization
//     * which is usually not necessary.
//     */
//    Solver<BlockFEMatrix, BlockVector> solver;
//
//    /** @brief an array to store defect, so that we don't have to reallocate
//     *         so often
//     */
//    BlockVector defect;
//
//    /**
//    * @brief store the square root of the residual from previous iteration
//    */
//    double oldResidual;
//
//    /** @brief store the initial residual so that the nonlinear iteration can
//    *         be stopped as soon as a desired reduction is achieved
//    */
//    double initial_residual;
//
//    /** @brief store errors  */
//    std::vector<double> errors;
//
//    /** @brief right hand side vector from previous time step (on finest mesh)*/
//    BlockVector old_rhs;
//
//    BlockVector old_solution;
//
//    /** old time step length used to scale the pressure blocks*/
//    double oldtau;
//
//    /** @brief set parameters in database
//    *
//    * This functions checks if the parameters in the database are meaningful
//    * and resets them otherwise. The hope is that after calling this function
//    * this class is fully functional.
//    *
//    * If some parameters are set to unsupported values, an error occurs and
//    * throws an exception.
//    */
//    void set_parameters();
//
//    /** @brief get velocity and pressure space*/
//    void get_velocity_pressure_orders(std::pair <int,int> &velocity_pressure_orders);
//
//    /** @brief write some information (number of cells, dofs, ...) */
//    void output_problem_size_info() const;
//
//  public:
//
//    /** @brief constructor
//     * This constructor calls the other constructor creating an Example_TimeNSE2D
//     * object.
//     */
//    Time_NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
//               int reference_id = -4711);
//
//    /** @brief constructor
//     *
//     * The domain must have been refined a couple of times already. On the
//     * finest level the finite element spaces and functions as well as
//     * matrices, solution and right hand side vectors are initialized.
//     *
//     * The reference_id can be used if only the cells with the give reference_id
//     * should be used. The default implies all cells.
//     */
//    Time_NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
//               const Example_TimeNSE2D& ex, int reference_id = -4711);
//
//    /** @brief Assemble all the matrices and rhs before the time iterations
//     *
//     * This includes the assembling of: Stiff_matrix, Mass_Matrix,
//     * (additional matrixK in case of SUPG stabilization), rhs
//     */
//    void assemble_initial_time();
//
//    /** @brief
//     * 1. assembling the right hand side
//     * 2. scaling of the B-blocks due to time stepping
//     * this function will prepare the right hand side during the time
//     * discretization but should be outside the nonlinear loop
//     */
//    void assemble_rhs();
//
//    /** @brief Assemble the right hand side only
//     *  @param la: LocalAssembling2D object
//     *  @param rhs: rhs vector to be assembled
//    */
//    void AssembleRhs();
//
//    /** @brief Assemble the system matrix
//     * This function will prepare the system which will be
//     * used for solvers
//     */
//    void assemble_system();
//
//    /** descale matrices
//     * This function will descale all A-blocks which were scaled
//     * during the function call Time_NSE2D::assemble_system():
//     */
//    void deScaleMatrices();
//
//    /** @brief assemble nonlinear term
//     *
//     * The matrix blocks to which the nonlinear term contributes are reset to
//     * zero and then completely reassembled, including the linear and nonlinear
//     * terms. If this->assemble() has been called before, the matrix is now set
//     * up correctly.
//     */
//    void assemble_nonlinear_term();
//
//    /** @brief solve the system */
//    void solve();
//
//    /** @brief check if one of the stopping criteria is fulfilled
//     *
//     * either converged, maximun number of iterations reached, or slow
//     * convergence
//     *
//     * @param it_counter current iterate
//     */
//    bool stopIte(unsigned int it_counter);
//
//    /** @brief
//     * compute errors and write solution
//     */
//    void output(int m);
//
//    // getters and setters
//    /*const BlockMatrixNSE2D & get_matrix() const
//    { return this->systems.front().matrix; }*/
//    const BlockVector & get_rhs() const
//    { return this->systems.front().rhs; }
//
//    const TFEVectFunct2D & get_velocity() const
//    { return this->systems.front().u; }
//    // try not to use this as it is not const
//    TFEFunction2D *get_velocity_component(int i)
//    { return (i==0) ? this->systems.front().u.GetComponent(0)
//                    : this->systems.front().u.GetComponent(1); }
//    const TFEFunction2D & get_pressure() const
//    { return this->systems.front().p; }
//    const TFESpace2D & get_velocity_space() const
//    { return this->systems.front().velocity_space; }
//    const TFESpace2D & get_pressure_space() const
//    { return this->systems.front().pressure_space; }
//    const BlockVector & get_solution() const
//    { return this->systems.front().solution; }
//    unsigned int get_size() const
//    { return this->systems.front().solution.length(); }
//    const Example_TimeNSE2D & get_example() const
//    { return example; }
//    const ParameterDatabase & get_db() const
//    { return db; }
//    /// @brief return the computed errors at each discre time point
//    std::array<double, int(6)> get_errors();
//
//
//
//
//
//
//
//    /** ***************BELOW THIS LINE, USER SPECIFIC CODE ********/
//  public:
//
//    /** @brief assemble matrix,
//     *
//     * This assembles everything which is not related to the nonlinear term.
//     * I.e. it assembles a Stokes matrix.
//     */
//    void assemble_initial_time_withfields(TFEFunction2D* rho_field=nullptr,
//                                          TFEFunction2D* mu_field=nullptr);
//
//    /** @brief assemble nonlinear term
//     *
//     * The matrix blocks to which the nonlinear term contributes are reset to
//     * zero and then completely reassembled, including the linear and nonlinear
//     * terms. If this->assemble() has been called before, the matrix is now set
//     * up correctly.
//     */
//    void assemble_nonlinear_term_withfields(TFEFunction2D* rho_field=nullptr,
//                                            TFEFunction2D* mu_field=nullptr);
//
//
//    void assemble_rhs_withfields(TFEFunction2D* rho_field=nullptr,
//                                 TFEFunction2D* mu_field=nullptr);
//
//    void assemble_massmatrix_withfields(TFEFunction2D* rho_field=nullptr);
//
//// ======================================================================
//    // The following members are used for IMEX
///** @brief This returns the number of the current time step.
//     * This counter is set at 0 before the time loop and is incremented at each
//     * time step (but not at each sub-step) in the main program.
//     * It can be useful to give info to the members of the class. It is for example
//     * used in IMEX scheme to detect when we passed 2 time steps, so that we
//     * are guaranteed to have saved both old_solution_ and old_solution2_ correctly.    */
//    int current_step_;
//
//    /** @brief Construct the extrapolated solution.
//     * At the moment, only IMEX is implemented. */
//    void construct_extrapolated_solution();
//
//    bool imex_scheme(bool print_info);
//
//    /** @brief constructs a solution vector extrapolated from previous steps
//     * Currently, it is used for IMEX-Scheme: 2u(t-1)-u(t-2). */
//  protected:
//    BlockVector extrapolated_solution_;
//// ======================================================================
//
//  protected: // these members are used when interpolating
//    // rho and mu field into same space as velocity_space
//    std::vector<double> entries_rho_scalar_field;
//    std::vector<double> entries_mu_scalar_field;
////    TFEFunction2D interpolated_rho_field;
////    TFEFunction2D interpolated_mu_field;
//
