/** ************************************************************************* 
 * @name         Time_NSE2D_BDF
 * @brief        store everything needed to solve a time dependent convection 
 *               diffusion reaction problem
 *               Stores matrix, right hand side, FE spaces, FE functions 
 *               and the solution vector
 *               
 * @author       Naveed Ahmed
 * @History      16.09.2015
 **************************************************************************/

#ifndef TIME_NSE2D_BDF_H
#define TIME_NSE2D_BDF_H


#include <FEVectFunct2D.h>

#include <BlockFEMatrix.h>
#include <BlockVector.h>

#include <FESpace2D.h>
#include <Example_TimeNSE2D.h>

#include <Multigrid.h>
#include <Solver.h>

#include <MainUtilities.h>

#include <ParameterDatabase.h>
#include <PostProcessing2D.h>
#include <Residuals.h>
#include <LocalAssembling2D.h>

#include <vector>
#include <deque>
#include <utility>
#include <array>


class Time_NSE2D_BDF
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
      
      /** @brief MatrixK includes all the terms which are related to the 
       * time derivatives in the fully discrete scheme Eq: (45)
       * Ahmed, Rebollo, John and Rubino (2015)
       * Weighted Mass matrix
       *  [ M11  M12  0 ]
       *  [ M21  M22  0 ]
       *  [ 0     0   0 ]
       * This additionaly created due to the struture of the 
       * residual based Variational Multiscale Method: In all 
       * other methods, for example SUPG or Galerkin, the Mass_Matrix
       * is used.
       */
      BlockFEMatrix MatrixK;
      
      /** @brief 
       * old solution for the computation of the residual
       * that passes as an FEFunction to the local assembling
       * routines
       */
      BlockVector solution_m1;
      BlockVector solution_m2;
      /** @brief Finite element function for old velocities*/
      TFEVectFunct2D u_m1;
      TFEVectFunct2D u_m2;
      TFEFunction2D p_old;
      /** @brief linear combination of the old solutions
       * this is needed for the assembling of the pressure 
       * right-hand side when the SUPG or residual based
       * vms methods are considered
       */      
      BlockVector combined_old_sols;
      /// the corresponding FE function
      TFEVectFunct2D comb_old_u;
      
      /// In the case of implicit explicit schemes
      /// one needs extrapolation of the solution
      BlockVector extrapolate_sol;
      /// the corresponding fe function
      TFEVectFunct2D extrapolate_u;

      /** @brief constructor*/
      System_per_grid(const Example_TimeNSE2D& example, TCollection& coll, 
                      std::pair<int,int> order, Time_NSE2D_BDF::Matrix type);
    };
    
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

    ///@brief The norms of residuals from up to 10 previous iterations
    FixedSizeQueue<10, Residuals> oldResiduals;
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
public:
    Time_NSE2D_BDF(const TDomain& domain, const ParameterDatabase& param_db,
                       const Example_TimeNSE2D& ex, int reference_id=0);
    
    void assemble_initial_time();
    
    /** @brief 
     * 1. assembling the right hand side
     * 2. scaling of the B-blocks due to time stepping
     * this function will prepare the right hand side during the time 
     * discretization but should be outside the nonlinear loop
     */
    void assemble_rhs(bool rhs_assemble = true);
    
    /** @brief Assemble the system matrix
     * This function will prepare the system which will be 
     * used for solvers
     */
    void assemble_system();
    
    /** @brief assemble nonlinear term
     * 
     * The matrix blocks to which the nonlinear term contributes are reset to 
     * zero and then completely reassembled, including the linear and nonlinear
     * terms. If this->assemble() has been called before, the matrix is now set
     * up correctly. 
     */
    void assemble_nonlinear_term();
    
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
    
    ParameterDatabase pre_step_time_db;
    /** descale matrices
     * This function will descale all A-blocks which were scaled
     * during the function call Time_NSE2D::assemble_system():
     */
    void deScaleMatrices();
    
    /** @brief
     * compute errors and write solution
     */
    void output(int m);
    /** @brief this returns the number of the current time step 
     * In the main programe:: initialize this parameter to 0 before
     * the time iterations
     * This will also serves for the semi-implicit schemes.
     * For example, in the BDF2 scheme where in the first step
     * backward Euler time stepping is used to get the solution
     * at the second time step to be sure that one have at least 
     * 2 initial solutions
     */
    int current_step_;
    /** @brief check if the semi-implicit scheme is used
     */
    bool imex_scheme(bool print_info);

    
    void modify_slip_bc(bool modify_all);
    // getters and setters
    /*const BlockMatrixNSE2D & get_matrix() const
    { return this->systems.front().matrix; }*/
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }

    const TFEVectFunct2D & get_velocity() const
    { return this->systems.front().u; }
    
    const TFEVectFunct2D & get_velocity_old() const 
    { return this->systems.front().u_m1;}
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
    std::array<double, int(6)> get_errors();

    
private:
  /// this routines wraps up the call to Assemble2D  
  void call_assembling_routine(Time_NSE2D_BDF::System_per_grid& s, LocalAssembling2D_type type);
  /// set the matrices and right hand side depending on the 
  /// assemling routines, nstypes and the methods 
  void set_matrices_rhs(Time_NSE2D_BDF::System_per_grid& s, LocalAssembling2D_type type, 
                    std::vector<TSquareMatrix2D*> &sqMat, std::vector<TMatrix2D*> &reMat, 
                    std::vector<double*> &rhs);
  /// set the spaces depending on disc types
  void set_arrays(Time_NSE2D_BDF::System_per_grid& s, std::vector<const TFESpace2D*> &spaces, 
                  std::vector< const TFESpace2D* >& spaces_rhs, 
                  std::vector< TFEFunction2D*> &functions);
};

#endif // TIME_NSE2D_BDF_H
