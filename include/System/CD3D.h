/** ***************************************************************************
 *
 * @name       CD3D
 * @brief      store everything needed to solve a convection-diffusion-reaction
 *             (cdr) problem
 *
 *             Store matrix, right hand side, FE spaces, FE functions and 
 *             the solution vector of a convection-diffusion problem. This 
 *             wraps up everything which is necessary to solve a convection 
 *             diffusion problem in 3D.
 *
 *             Here are some tasks for class CD3D
 *             \todo enable different discretizations/stabilizations
 *             \todo enable direct solver in MPI
 *             \todo test multigrid with all smoothers (so far: only SSOR)
 *
 * @author     Ulrich Wilbrandt, Clemens Bartsch
 * @date       09.06.15, 2015/10/20
 *
 ******************************************************************************/

#ifndef __CD3D_H__
#define __CD3D_H__

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <ParameterDatabase.h>
#include <Solver.h>

#include <PostProcessing3D.h>

#include <FEFunction3D.h>
#include <Example_CD3D.h>
#include <anderson.h>
#include <PostProcessing3D.h>

#include <vector>
#include <deque>
#include <array>
#include <list>

class LocalAssembling3D; //forward declaration
class Example_CD3D;


class CD3D
{
  protected:
    /** @brief store a complete system on a particular grid
     * 
     * This combines a matrix, rhs, solution, space and function needed to
     * describe one convection-diffusion-reaction problem in 3D.
     * In MPI case the parallel infrastructure is stored, too.
     */
    struct SystemPerGrid
    {
      /** @brief Finite Element space */
      std::shared_ptr<TFESpace3D> feSpace_;
      /** @brief the system matrix */
      BlockFEMatrix matrix_;
      /** @brief the right hand side vector */
      BlockVector rhs_;
      /** @brief solution vector with one component. */
      BlockVector solution_;
      /** @brief Finite Element function */
      TFEFunction3D feFunction_;
      /** @brief vector of weights for an AFC scheme */
      std::vector<double> afc_gamma;
      /** @brief entries of correction matrix for AFC schemes */
      std::vector<double> afc_matrix_D_entries;
      
      /** @brief constructor */
      SystemPerGrid( const Example_CD3D& example, TCollection& coll, 
                       int ansatz_order );


      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] maxSubDomainPerDof Only in MPI case! The maximal number of
       * processes which share a single dof. The value is calculated in the
       * main program and handed down to the FESpaces. Getting rid of this
       * construction is a TODO .
       */
#ifdef _MPI
      SystemPerGrid(const Example_CD3D& example,
                    TCollection& coll, int maxSubDomainPerDof);
#else
      SystemPerGrid( const Example_CD3D& example, TCollection& coll );
#endif

      // Declaration of special member functions - delete all but destructor.
      // TODO CB Apply rule of zero/five as soon as underlying classes do.

      /** Delete copy constructor. No copies allowed until
      TFESpace3D, BlockMatrixCD3D and TFEFunction3D obey rule of 0/5. */
      SystemPerGrid(const SystemPerGrid&) = delete;

      //! Delete move constructor. No moves allowed.
      SystemPerGrid(SystemPerGrid&&) = delete;

      //! Delete copy assignment operator. No copies allowed.
      SystemPerGrid& operator=(const SystemPerGrid&) = delete;

      //! Default move assignment operator. No moves allowed.
      SystemPerGrid& operator=(SystemPerGrid&&) = delete;

      //! Default destructor. Does most likely cause memory leaks.
      ~SystemPerGrid() = default;

    };
    
    /** @brief A complete system on each involved grid.
     * 
     * Note that the size of this double ended queue is at least one and
     * larger only when a multgrid preconditioner is used.
     */
    std::deque<SystemPerGrid> systems_;
    
    /** @brief Definition of the used example. */
    const Example_CD3D example_;
    
    /** @brief A multigrid object which is set to nullptr in case it is not
     *         needed.
     */
    std::shared_ptr<TMultiGrid3D> multigrid_;
    
    /** @brief a local parameter database which controls this class
     * 
     * The database given to the constructor will be merged into this one. Only 
     * parameters which are of interest to this class are stored (and the 
     * default ParMooN parameters). Note that this usually does not include 
     * other parameters such as solver parameters. Those are only in the 
     * CD3D::solver object.
     */
    ParameterDatabase db;
    
    /** @brief class for output handling */
    PostProcessing3D outputWriter;

    
    /** @brief a solver object which will solve the linear system
     * 
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;

    /** @brief Errors to be accesed from outside the class
     * The array is filled during the function call CD3D::output()
     * Currently, the errors store the global L2 and H1 errors.
     */
    std::array<double, int(2)> errors_;
    
    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;
    
     /** @brief an interpolant of the exact solution onto the fe space */
    std::vector<double> exact_interpolant;
    
    void call_assembling_routine(SystemPerGrid& s, LocalAssembling3D& la);


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
    CD3D(std::list<TCollection* > collections,
         const ParameterDatabase & param_db, const Example_CD3D& example,
         int maxSubDomainPerDof);
#else
    CD3D(std::list<TCollection* > collections, 
         const ParameterDatabase & param_db, const Example_CD3D& _example);
#endif
    
    /** @brief Assemble the system matrix resp. matrices in multigrid case.
     * 
     * Depending on 'this->db["space_discretization_type]' different (local)
     * assembling routines are supposed to be used - so far only standard
     * "galerkin" ("1") is supported.
     */
    void assemble();
    
    /** @brief Solve the system.
     *
     * Have a look at the code to see which solvers are already enabled.
     * More are to come.
     */
    void solve();
    
    /** 
     * @brief Measure errors and write nice pictures.
     * 
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk files are created.
     *
     * @param i Integer suffix for output file name. -1 means no suffix.
     */
    void output(int i = -1);
    
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

    /** @brief Get the system matrix on the currently finest grid.*/
    const BlockFEMatrix & getMatrix() const
    {
      return systems_.front().matrix_;
    }

    /** @brief Get the system matrix on the currently finest grid.*/
    BlockFEMatrix & getMatrix()
    {
      return systems_.front().matrix_;
    }

    /** @brief Get the right hand side on the currently finest grid.*/
    const BlockVector & getRhs() const
    {
      return systems_.front().rhs_;
    }

    /** @brief Get the right hand side on the currently finest grid.*/
    BlockVector & getRhs()
    {
      return systems_.front().rhs_;
    }


    /** @brief Get the fe function on the currently finest grid.*/
    const TFEFunction3D & getFunction() const
    {
      return systems_.front().feFunction_;
    }

    /** @brief Get the fe space on the currently finest grid.*/
    const TFESpace3D & getSpace() const
    {
      return *systems_.front().feSpace_;
    }

    /** @brief Get the solution vector on the currently finest grid.*/
    const BlockVector & getSolution() const
    {
      return systems_.front().solution_;
    }

    /** @brief Get the solution vector on the currently finest grid.*/
    BlockVector & getSolution()
    {
      return systems_.front().solution_;
    }

    /** @brief Get the used example.*/
    const Example_CD3D& getExample() const
    {
      return example_;
    }
    const ParameterDatabase & get_db() const
    {
      return db;
    }
    
     const int & get_global_disc_type() const
    {
      return global_space_type;
    }
    
    int & get_global_disc_type() 
    {
      return global_space_type;
    }

    /** ************************************************************************ */
    std::array< double, int(2) > get_errors() const
    {
      return errors_;
    }

    // Declaration of special member functions - delete all but destructor.
    // TODO CB Apply rule of zero/five as soon as underlying classes do.

    /** Delete copy constructor. This class may not be copied,
     * until its member classes fulfil the rule of 0/5. */
    CD3D(const CD3D&) = delete;

    /** Delete move constructor. This class may not be moved yet. */
    CD3D(CD3D&&) = delete;

    /** Delete copy assignment. This class may not be copied yet. */
    CD3D& operator=(const CD3D&) = delete;

    /** Delete move assignment. This class may not be moved yet. */
    CD3D& operator=(CD3D&&) = delete;
    
    int global_space_type;

    //! Default destructor. Does most likely cause memory leaks.
    ~CD3D() = default;
   
    
};

#endif // __CD3D_H__
