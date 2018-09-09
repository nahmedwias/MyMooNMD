#ifndef TIME_NSE2D_MERGE_H
#define TIME_NSE2D_MERGE_H

/** ************************************************************************
 *
 * @name         Time_NSE2D
 * @brief        store everything needed to solve a time dependent convection
 *               diffusion reaction problem
 *               Stores matrix, right hand side, FE spaces, FE functions
 *               and the solution vector
 *
 * @author       Naveed Ahmed
 * @History      20.02.2017
**************************************************************************/


#include <FEVectFunct2D.h>

#include <BlockFEMatrix.h>
#include <BlockVector.h>

#include <FESpace2D.h>
#include <Example_TimeNSE2D.h>

#include <Multigrid.h>
#include <Solver.h>

#include <MainUtilities.h>

#include <ParameterDatabase.h>
#include <DataWriter.h>
#include <Residuals.h>
#include "LocalAssembling.h"
#include <TimeDiscretizations.h>


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
      TFEVectFunct2D u;
      /** @brief Finite element function for pressure*/
      TFEFunction2D p;
      /** @brief
       * old solution for the computation of the residual
       * that passes as an FEFunction to the local assembling
       * routines
       */
      BlockVector solution_m1;
      ///@brief Finite element functions for velocity and pressure
      TFEVectFunct2D u_m1;
      TFEFunction2D p_m1;
      ///@brief solution from two previous time steps
      ///this is used for the BDF2 scheme
      BlockVector solution_m2;
      ///@brief the corresponding finite element functions for 
      ///velocity and pressure
      TFEVectFunct2D u_m2;
      TFEFunction2D p_m2;
      ///@brief linear combination of old solutions
      ///this particularly is used for the BDF2 and SUPG/RBVMS
      ///when the equal-order elements are used
      BlockVector combined_old_sols;
      TFEVectFunct2D comb_old_u;
      
      /// @brief extrapolated solution for imex scheme
      BlockVector extrapolate_sol;
      TFEVectFunct2D extrapolate_u;
      
      /// @brief extrapolated pressure 
      TFEFunction2D extrapolate_p;
      /// @brief time derivative of the solution
      /// this is actually used in the computation of the residual
      /// for the residual based vms methods
      BlockVector time_deriv_sol;      
      TFEVectFunct2D u_td;

      /** @brief constructor*/
      System_per_grid(const Example_TimeNSE2D& example, TCollection& coll,
                      std::pair<int,int> order, Time_NSE2D::Matrix type);
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
    DataWriter2D outputWriter;

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
    
    /// this parameter is used to write the solution after the number of time 
    /// steps: 
    int counter_for_time_scale;
public:
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


    //
    /** @brief Assemble all the matrices and rhs before the time iterations
     *
     * This includes the assembling of: Stiff_matrix, Mass_Matrix,
     * (additional matrixK in case of SUPG stabilization), rhs
     */
    void assemble_initial_time();

    /** @brief This will assemble the right-hand side, will prepare the 
     * right-hand side using the class TimeDiscretization object. 
     * The pressure blocks are scaled properly and the nonlinear matrices
     * are assembled again. The system matrix is assemble using the 
     * object of class TimeDiscretization and modify the matrices if needed.
     */
    void assemble_matrices_rhs(unsigned int it_counter);
    
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

    /**
     * @brief Compute the defect Ax-b and store its norm in NSE2D::oldResiduals
     *
     * where A is the current matrix, x is the current solution and b is the
     * right hand side. Call this function after assembling the nonlinear
     * matrix with the current solution.
     */
    void computeNormsOfResiduals();

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

    /** @brief this returns the number of the current time step
     * In the main programe:: initialize this parameter to 0 before
     * the time iterations
     * This will also serves for the semi-implicit schemes.
     * For example, in the BDF2 scheme where in the first step
     * backward Euler time stepping is used to get the solution
     * at the second time step to be sure that one have at least
     * 2 initial solutions
     */
    // int current_step_;
    /** @brief check if the semi-implicit scheme is used
     */
    bool imex_scheme(bool print_info);

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
      * some common declaration for postprocessing
      * this only needs on the finest grid
     */
    std::shared_ptr<TFESpace2D> vorticity_space;
    unsigned int n_vort_dofs;
    std::vector<double> vorticity;
    std::shared_ptr<TFEFunction2D> vorticity_funct;
    std::shared_ptr<TFEFunction2D> divergence;
    double zero_vorticity;
    const TFEFunction2D & get_vorticity_funct() const
    {return *vorticity_funct.get();}

    std::shared_ptr<TFESpace2D> stream_function_space;
    unsigned int n_psi;
    std::vector<double> psi;
    std::shared_ptr<TFEFunction2D> stream_function;

    void prepared_postprocessing(TCollection *coll);
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
    // @brief
    double getFullResidual() const;

    TimeDiscretization& get_time_stepping_scheme()
    {return time_stepping_scheme;}
    const TimeDiscretization& get_time_stepping_scheme() const
    {return time_stepping_scheme;}
    
    int get_space_disc_global() {return space_disc_global;}
    
private:
  /// @brief this routines wraps up the call to Assemble2D
  void call_assembling_routine(Time_NSE2D::System_per_grid& s,
                               LocalAssembling_type type);
  /// @brief set the matrices and right hand side depending on the
  /// assemling routines, nstypes and the methods
  void set_matrices_rhs(Time_NSE2D::System_per_grid& s,
                        LocalAssembling_type type,
                        std::vector<TSquareMatrix2D*> &sqMat,
                        std::vector<TMatrix2D*> &reMat,
                        std::vector<double*> &rhs);
  /// @brief set the spaces depending on disc types
  void set_arrays(Time_NSE2D::System_per_grid& s,
                  std::vector<const TFESpace2D*> &spaces,
                  std::vector< const TFESpace2D* >& spaces_rhs,
                  std::vector< TFEFunction2D*> &functions);  
  /// @brief restrict the function to on every grid
  /// nonliear assembling requires an approximate velocity
  /// on every grid
  void restrict_function();
  /// update matrices for local projection stabilization
  void update_matrices_lps(Time_NSE2D::System_per_grid& s);
  /// @brief Special case for the SUPG and RBVMS method. The right-hand side is 
  /// nonlinear and need a complete re-assembling before passing to the solver 
  /// for solving.
  void assemble_rhs_nonlinear();
};



#endif // TIME_NSE2D_MERGE_H
