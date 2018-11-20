/** ************************************************************************
 * 
 * @name         Time_CD3D
 * @brief        store everything needed to solve a time dependent convection 
 *               diffusion reaction problem
 *               Stores matrix, right hand side, FE spaces, FE functions 
 *               and the solution vector
 *               
 * @author       Naveed Ahmed
 * @History      13.05.2016
**************************************************************************/

#ifndef __Time_CD3D__
#define __Time_CD3D__

#include <FEFunction3D.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Example_TimeCD3D.h>
#include <Solver.h>
#include <DataWriter.h>
#include "LocalAssembling.h"
#include "TimeDiscretizations.h"

#include <vector>
#include <array>
#include <deque>

class Time_CD3D
{
  protected:
    /** @brief store a complete system on a particular grid.
     * 
     * This combines a matrix, rhs, solution, spaces and functions 
     * needed to describe a Time CDR problem in 3D
     */
    struct SystemPerGrid
    {
      /** @brief Finite element space */
      TFESpace3D feSpace_;
      /** @brief Stiffness Matrix */
      BlockFEMatrix stiffMatrix_;
      /** @brief Mass matrix */
      BlockFEMatrix massMatrix_;
      /** @brief right hand side vector */
      BlockVector rhs_;
      /** @brief solution vector */
      BlockVector solution_;

      /** @brief Finite element function */
      TFEFunction3D feFunction_;
      
      /** @brief
       * old solution for the computation of the residual
       * that passes as an FEFunction to the local assembling
       * routines
       */
      BlockVector solution_m1;
      TFEFunction3D u_m1;
      BlockVector solution_m2;
      TFEFunction3D u_m2;

      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] maxSubDomainPerDof Only in MPI case! The maximal number of
       * processes which share a single dof. The value is calculated in the
       * main program and handed down to the FESpaces. Getting rid of this
       * construction is a TODO .
       */
      SystemPerGrid(const Example_TimeCD3D& example, TCollection& coll);
      /**
       * Special member functions mostly deleted,
       * for struct takes ownership of the bad
       * classes TFEFunction3D and TFESpace3D.
       */
      //! Delete copy constructor.
      SystemPerGrid(const SystemPerGrid&) = delete;

      //! Delete move constructor.
      SystemPerGrid(SystemPerGrid&&) = delete;

      //! Delete copy assignment operator.
      SystemPerGrid& operator=(const SystemPerGrid&) = delete;

      //! Delete move assignment operator.
      SystemPerGrid& operator=(SystemPerGrid&&) = delete;

      //! Default destructor. Most likely causes memory leaks.
      ~SystemPerGrid() = default;

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
    std::deque<SystemPerGrid> systems_;
    
    
    /** @brief Definition of the used example */
    const Example_TimeCD3D example_;
    
    /** @brief old right hand side vectior 
     * this will be used to save the right hand side from the 
     * previous time step that will be used for different 
     * time stepping schemes
     */
    BlockVector old_rhs;
    
    /** @brief store the errors to compute accumulated error norms */
    std::array<double, 5> errors_;
    
    /** @brief output object */
    DataWriter3D outputWriter;
    
    /// @brief time stepping scheme object to access everything
    TimeDiscretization time_stepping_scheme;
    /// @brief system right hand side which passes to the solver
    BlockVector rhs_from_time_disc;

    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;

    /** @brief an internal integer to store the discretization type.
     * This integer corresponds to the old DISCTYPE, but it is taken
     * from the database of this class, and is not a global parameter
     * anymore. It is more convenient than using the long syntax
     * db["space_discretization_type"].is("blablabla"), and it is
     * used when calling the local assembling objects.
     * Setting this parameter is done in the method CheckParameters
     * which has now been renamed check_and_set_parameters.
     * This disctype is necessary to use correctly different disctype,
     * e.g. supg.
     */
    int disctype;

  public:
    /** @brief The standard constructor, can be used for multigrid and 
     * non-multigrid.
     *
     * @param[in] domain The computational domain providing the grids
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] example The example which is to be calculated.
     */
    Time_CD3D(const TDomain& domain, const ParameterDatabase &param_db,
	      const Example_TimeCD3D& _example);
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
    
    /** @brief measure errors and write solution
     * 
     */
    void output(int m);
    
    /**
     * @brief Check whether the program will be working with the
     * current input parameters.
     *
     * ParMooN is work in progress, and so is this class. This method
     * checks the parameters stored in the database and stops execution
     * of the program, if some of these do not match.
     * The method is a little makeshift and the cases caught here are various,
     * but basically it is intended to stop execution of cases with
     * parameter combinations which are not implemented for CD3D or
     * are currently known to be problematic.
     *
     * This is not yet a guarantee for a functioning program, but is
     * intended to be, someday. Eventually this method and the like
     * will be moved to TDatabase.
     */
    void check_and_set_parameters();
    
     // getters and setters
    const Example_TimeCD3D& get_example() const
    { return example_; }
    const TFEFunction3D & get_function() const
    { return this->systems_.front().feFunction_; }
    TFEFunction3D & get_function()
    { return this->systems_.front().feFunction_; }
    const BlockFEMatrix & get_stiff_matrix() const
    { return this->systems_.front().stiffMatrix_; }
    const BlockVector & get_rhs() const
    { return this->systems_.front().rhs_; }
    BlockVector & get_rhs()
    { return this->systems_.front().rhs_; }
    const BlockVector & get_solution() const
    { return this->systems_.front().solution_; }
    const TFESpace3D & get_space() const
    { return this->systems_.front().feSpace_; }
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
     *
     * Which afc algorithm is performed is determined by switching over
     * the algebraic_flux_correction parameter of the database. (So far only
     * fem-fct-cn is enabled.)
     */
    void do_algebraic_flux_correction();

    /**
     * Scramble together the parameters which Assemble3D needs and call it.
     * Is only put here to keep the code of assemble() slender.
     * assemble() should take care of the right choice of the LocalAssembling
     * object and whether it fits the system matrix' block structure
     * @param s The sytem where rhs, mass and stiffness matrices are to be assembled.
     * @param la_stiff The local assembling object of choice.
     * @param assemble_both If true, both stiffness (+rhs) and mass matrix are
     * assembled, if false only stiffness matrix and rhs.
     */
    void call_assembling_routine(Time_CD3D::SystemPerGrid& system,
                                 LocalAssembling3D& la_stiff, bool assemble_both);
};

#endif
