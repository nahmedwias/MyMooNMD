/** ************************************************************************
 * 
 * @name         Time_CD2D
 * @brief        store everything needed to solve a time dependent convection 
 *               diffusion reaction problem
 *               Stores matrix, right hand side, FE spaces, FE functions 
 *               and the solution vector
 *               
 * @author       Naveed Ahmed, Ulrich Wilbrandt, Clemens Bartsch
 * @History      16.09.2015
**************************************************************************/

#ifndef __Time_CD2D__
#define __Time_CD2D__

#include <FEFunction2D.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Example_TimeCD2D.h>
#include <Multigrid.h>
#include <Domain.h>
#include <DataWriter.h>
#include <Solver.h>
#include "LocalAssembling.h"
#include "TimeDiscretizations.h"

#include <vector>
#include <deque>

class Time_CD2D
{
  protected:
    /** @brief store a complete system on a particular grid.
     * 
     * This combines a matrix, rhs, solution, spaces and functions 
     * needed to describe a Time CDR problem in 2D
     */
    struct System_per_grid
    {
      /** @brief Finite element space */
      TFESpace2D fe_space;
      /** @brief Stiffness Matrix */
      BlockFEMatrix stiff_matrix;
      /** @brief Mass matrix */
      BlockFEMatrix mass_matrix;
      /** @brief right hand side vector */
      BlockVector rhs;
      /** @brief solution vector */
      BlockVector solution;
      /** @brief Finite element function */
      TFEFunction2D fe_function;
      /** @brief
       * old solution for the computation of the residual
       * that passes as an FEFunction to the local assembling
       * routines
       */
      BlockVector solution_m1;
      TFEFunction2D u_m1;
      BlockVector solution_m2;
      TFEFunction2D u_m2;

      /** @brief constructor*/
      System_per_grid(const Example_TimeCD2D& example, TCollection& coll);

      /**
       * Special member functions mostly deleted,
       * for struct takes ownership of the bad
       * classes TFEFunction2D and TFESpace2D.
       */
      //! Delete copy constructor.
      System_per_grid(const System_per_grid&) = delete;

      //! Delete move constructor.
      System_per_grid(System_per_grid&&) = delete;

      //! Delete copy assignment operator.
      System_per_grid& operator=(const System_per_grid&) = delete;

      //! Delete move assignment operator.
      System_per_grid& operator=(System_per_grid&&) = delete;

      //! Default destructor. Most likely causes memory leaks.
      ~System_per_grid() = default;

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
    /** @brief a solver object which will solve the linear system
     * 
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;

    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger than that
     * only in case of multigrid (when it holds as many systems as there are
     * multigrid levels).
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    const Example_TimeCD2D example;
    
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
    
    /** @brief old right hand side vectior 
     * this will be used to save the right hand side from the 
     * previous time step that will be used for different 
     * time stepping schemes
     */
    BlockVector old_rhs;
    
    /** @brief store the errors to compute accumulated error norms */
    std::vector<double> errors;

    /// @brief class for handling (time dependent) output 
    DataWriter2D timeDependentOutput;
    
    /// @brief time stepping scheme object to access everything
    TimeDiscretization time_stepping_scheme;
    /// @brief system right hand side which passes to the solver
    BlockVector rhs_from_time_disc;
    
  public:
    /** @brief constructor
     * This constructor calls the other constructor creating an Example_CD2D
     * object. 
     */
    Time_CD2D(const TDomain& domain, const ParameterDatabase& param_db,
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
    Time_CD2D(const TDomain& domain, const ParameterDatabase& param_db,
    		const Example_TimeCD2D& ex, int reference_id = -4711);
    
    /** @brief Assemble all the matrices before the time iterations
     * 
     * This includes the assembling of: Stiff_matrix, Mass_Matrix, 
     * (additional matrixK in case of SUPG stabilization), rhs
     */
    void assemble_initial_time();
    
    /** @brief assemble the matrices
     * this function will assemble the stiffness matrix and rhs
     * In addition the system matrix and the rhs which passes to the solver 
     * are also prepared within the function
     */
     void assemble();

    /** @brief solve the system
     */
    void solve();
    
    /// @brief measure errors and write solution
    void output();
     // getters and setters
    const Example_TimeCD2D& get_example() const
    { return example; }
    const TFEFunction2D & get_function() const
    { return this->systems.front().fe_function; }
    TFEFunction2D & get_function()
    { return this->systems.front().fe_function; }
    const BlockFEMatrix & get_stiff_matrix() const
    { return this->systems.front().stiff_matrix; }
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    BlockVector & get_rhs()
    { return this->systems.front().rhs; }
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    const TFESpace2D & get_space() const
    { return this->systems.front().fe_space; }
    const ParameterDatabase & get_db() const
    { return db; }
    
    TimeDiscretization& get_time_stepping_scheme()
    {return time_stepping_scheme;}
    const TimeDiscretization& get_time_stepping_scheme() const
    {return time_stepping_scheme;}
    
    /**
    * @brief return the computed errors at each discre time point
    * 
    */
    std::array<double, int(3)> get_errors() const;


  private:
    /**
     * Apply an algebraic flux correction scheme to the assembled matrix.
     * Should be called within the assemble routine, after the assembling
     * of pure mass and stiffness matrix and right hand side
     * has been performed with the INTERNAL_FULL_MATRIX_STRUCTURE switch on.
     *
     * Which afc algorithm is performed is determined by switching over
     * the algebraic_flux_correction parameter of the database. (So far only
     * fem-fct-cn is enabled.)
     */
    void do_algebraic_flux_correction();

    /**
     * This wraps up the tedious call to Assemble2D and is only used
     * to avoid code duping. Is really not written very sophisticated,
     * use it with care.
     * @param block_mat should be one system's stiffness or mass matrix.
     * @param la A fittingly constructed LocalAssemble2D object which
     * is responsible for the assembling of the stiffness matrix and right hand
     * side.
     * @param assemble_both If true, both stiffness (+rhs) and mass matrix are
     * assembled, if false only stiffness matrix and rhs.
     */
    void call_assembling_routine(Time_CD2D::System_per_grid& system,
                                 LocalAssembling2D& la, bool assemble_both);
};

#endif
