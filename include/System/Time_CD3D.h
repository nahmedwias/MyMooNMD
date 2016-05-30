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
#include <Example_CD3D.h>
#include <MultiGrid3D.h>
#include <Solver.h>

#include <vector>
#include <deque>

#ifdef _MPI
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

class LocalAssembling3D; //forward declaration

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
      /** Stores stiffness_matrix * solution of last time step -
       *  this is needed in (partly) explicit time stepping schemes.*/
      BlockVector old_Au;
      /** @brief Finite element function */
      TFEFunction3D feFunction_;
#ifdef _MPI
      /** @brief A parallel FE mapper storing parallel information for 
       *the one matrix of this grid.
       */
      TParFEMapper3D parMapper_;
      
      /** @brief A parallel FE communicator taking care for the 
       * MPI communication on this grid. 
       */
      TParFECommunicator3D parComm_;
#endif 
      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] maxSubDomainPerDof Only in MPI case! The maximal number of
       * processes which share a single dof. The value is calculated in the
       * main program and handed down to the FESpaces. Getting rid of this
       * construction is a TODO .
       */
#ifdef _MPI
      SystemPerGrid(const Example_CD3D& example, TCollection& coll, 
		      int maxSubDomainPerDof);
#else
      SystemPerGrid(const Example_CD3D& example, TCollection& coll);
#endif

      /**
       * Gives a non-const pointer to the one block which is stored
       * by matrix. FIXME Is terribly unsafe as it makes use
       * of a deprecated block matrix method and must be replaced soon.
       */
      [[deprecated]]TSquareMatrix3D* get_stiff_matrix_pointer();

      /**
       * Reset the stiffness matrix A to its 'pure' state before the
       * modifications due to a one-step/fractional-step theta scheme.
       * Sets A = 1/(tau*theta_1)*(A - mass)
       * This is for the case we want to reuse A in the next time step.
       *
       * @param tau The current time step length.
       * @param theta_1 The impliciteness parameter for the transport (e.g. 1 for bw Euler).
       */
      void descale_stiff_matrix(double tau, double theta_1);

      void update_old_Au();

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
    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger than that
     * only in case of multigrid (when it holds as many systems as there are
     * multigrid levels).
     */
    std::deque<SystemPerGrid> systems_;
    
    /** @brief Definition of the used example */
    const Example_CD3D example_;
    
    /** a multigrid object from new multigrid class;
     * it stays nullptr if not used
     */
    std::shared_ptr<Multigrid>  multigrid_;
            
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
    
    /** @brief store the errors to compute accumulated error norms */
    std::vector<double> errors_;
    
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
     *
     * @param[in] collections A hierarchy of cell collections used for multigrid solve,
     * or just one collection in non-multigrid case. Ordered by fineness of the grid -
     * finest collection first!
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] example The example which is to be calculated.
     * @param[in] maxSubDomainPerDof Only in MPI case! the maximal number of
     * processes which share a single dof. The value is calculated in the
     * main program and handed down to the FESpaces. Getting rid of this
     * construction is a TODO .
     */
#ifdef _MPI
    Time_CD3D(std::list<TCollection* >collections, const ParameterDatabase &param_db,
	      const Example_CD3D& _example, int maxSubDomainPerDof);
#else
    Time_CD3D(std::list<TCollection* >collections, const ParameterDatabase &param_db,
	      const Example_CD3D& _example);
#endif    
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
    
    /**
     * Descales the stiffness matrices from the modifications due to time
     * discretization (one step/fractional step theta).
     * This should be called after all solve() of the current time step
     * have been performed and the stiffness should be reused.
     *
     * Sets A := 1/(theta_1 * tau) * (A - M), where A is stiffness matrix and
     * M mass matrix, on all grids.
     *
     * As a side effect, it updates the "old_Au" value which will be needed
     * for the next time step.
     */
    void descale_stiffness();

    /** @brief solve the system
     */
    void solve();
    
    /** @brief measure errors and write solution
     * 
     */
    void output(int m, int& image);
    
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
    void checkParameters();
    
     // getters and setters
    const Example_CD3D& get_example() const
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
     * ALGEBRAIC_FLUX_CORRECTION (so far only 2: linear C-N FEM-FCT).
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