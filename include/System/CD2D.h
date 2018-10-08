/** ***************************************************************************
 *
 * @class       CD2D
 * @brief      store everything needed to solve a convection-diffusion-reaction
 *             (cdr) problem
 *
 * Store matrix, right hand side, FE spaces, FE functions and the solution 
 * vector of a convection-diffusion problem. This wraps up everything which is 
 * necessary to solve a convection diffusion problem in 2D.
 *
 * @author     Ulrich Wilbrandt
 * @date       06.09.13
 *
 ******************************************************************************/

#ifndef __CD2D_H__
#define __CD2D_H__

#include <Example_CD2D.h>
#include <FEFunction2D.h>
#include <Domain.h>
#include <deque>

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <ParameterDatabase.h>
#include <Solver.h>
#include <DataWriter.h>

class Multigrid;

class CD2D
{
  protected:
    
    /** @brief store a complete system on a particular grid
     * 
     * This combines a matrix, rhs, solution, space and function needed to 
     * describe one convection-diffusion-reaction problem in 2D.
     */
    struct System_per_grid
    {
      /** @brief Finite Element space */
      std::shared_ptr<const TFESpace2D> fe_space;

      /** @brief the system matrix */
      BlockFEMatrix matrix;
      /** @brief the right hand side vector */
      BlockVector rhs;
      /** @brief solution vector with one component. */
      BlockVector solution;
      /** @brief Finite Element function */
      TFEFunction2D fe_function;
      
      /** @brief constructor */
      System_per_grid( const Example_CD2D& example, TCollection& coll, 
                       int ansatz_order );

      // Special member functions. Disable copy/move, set destructor to default.
      // Will be changed only when the underlying
      // classes TFESpace2D and TFESpace2D follow rule of 0/5.

      //! Delete copy constructor.
      System_per_grid(const System_per_grid&) = delete;

      //! Delete move constructor.
      System_per_grid(System_per_grid&&) = delete;

      //! Delete copy assignment operator.
      System_per_grid& operator=(const System_per_grid&) = delete;

      //! Delete move assignment operator.
      System_per_grid& operator=(System_per_grid&&) = delete;

      //! Default destructor.
      ~System_per_grid() = default;
    };
    
    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;
    
    /** @brief Definition of the used example */
    const Example_CD2D example;
    
    /** @brief a local parameter database which constrols this class
     * 
     * The database given to the constructor will be merged into this one. Only 
     * parameters which are of interest to this class are stored (and the 
     * defualt ParMooN parameters). Note that this usually does not include 
     * other parameters such as solver parameters. Those are only in the 
     * CD2D::solver object.
     */
    ParameterDatabase db;

    /** @brief class for output handling */
    DataWriter2D outputWriter;

    /** @brief a solver object which will solve the linear system
     * 
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;
    
    /** @brief store the errors to access them from outside this class
     * 
     * This array is filled during a call to CD2D::output if the parameter
     * "output_compute_errors" is set to true. The exact solution is taken from
     * CD2D::example. If that example does not provide an exact solution,
     * typically it is set to be zero, so that this array contains the norms of
     * the solution instead of the error.
     * 
     * The errors are stored in the following order: 
     * 
     *  - L2 error
     *  - H1-semi
     *  - SD error (streamline diffusion, useful for SDFEM)
     *  - L_inf error
     */
    std::array<double, 4> errors;
    
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
    
  public:
    
    /** @brief constructor 
     * 
     * This constructor calls the other constructor creating an Example_CD2D
     * object for you. See there for more documentation.
     *
     * @param[in] domain The readily treated (refined/partitioned...) domain 
     *                   object. Must not go out of scope before CD2D does!
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] reference_id The cell reference id, of which cells to create
     *                         the TCollection.
     */
    CD2D(const TDomain& domain, const ParameterDatabase& param_db,
         int reference_id = -4711);
    
    /** @brief constructor 
     * 
     * All members are initialized, including systems, which has only one entry
     * usually.
     *
     * @param[in] domain The readily treated (refined/partitioned...) domain
     *                   object. Must not go out of scope before CD2D does!
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] example a description of the example to be used, this is 
     *                    copied into a local member.
     * @param[in] reference_id The cell reference id, of which cells to create
     *                         the TCollection.
     */
    CD2D(const TDomain& domain, const ParameterDatabase& param_db, 
         const Example_CD2D& example, int reference_id = -4711);
    
    /** @brief assemble matrix, 
     * 
     * depending on 'this->db["space_discretization_type]' different (local)
     * assembling routines are used. Also in case of multigrid the matrices
     * on all grids are assembled.
     */

    // Initialize/Declare CD2D database, called with the constructor (for the use outside CD2D.C)
    static ParameterDatabase get_default_CD2D_parameters();

    void assemble(TFEFunction2D* coefficient_function1 = nullptr, TFEFunction2D* coefficient_function2 = nullptr);
    
    /** @brief solve the system */
    void solve();
    
    /** 
     * @brief measure errors and write pictures 
     * 
     * The current errors will be printed out. If desired, further output, e.g.,
     * vtk or case files are created.
     * 
     * @param i suffix for vtk output file name, -1 means no suffix
     */
    void output(int i = -1);
    
    /// @name return computed errors
    ///
    /// You have to call CD2D::output for any of these to return a meaningful 
    /// value.
    //@{
    /// @brief return the computed L2 error.
    double get_L2_error() const;
    /// @brief return the computed H1-semi.
    double get_H1_semi_error() const;
    /// @brief return the streamline diffusion (SD) error.
    double get_SD_error() const;
    /// @brief return the maximum error over all quadrature points in all cells.
    double get_L_inf_error() const;
    //@}
    
    // getters and setters
    const BlockFEMatrix & get_matrix() const
    { return this->systems.front().matrix; }
    BlockFEMatrix & get_matrix()
    { return this->systems.front().matrix; }
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    BlockVector & get_rhs()
    { return this->systems.front().rhs; }
    const TFEFunction2D & get_function() const
    { return this->systems.front().fe_function; }
    
    const TFESpace2D & get_space() const
    { return *this->systems.front().fe_space; }
    
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    BlockVector & get_solution()
    { return this->systems.front().solution; }
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    const Example_CD2D& get_example() const
    { return example; }
    const ParameterDatabase & get_db() const
    { return db; }

    // Special member functions. Disable copy/move, set destructor to default.
    // Will be changed only when the underlying classes follow rule of 0/5.

    //! Delete copy constructor.
    CD2D(const CD2D&) = delete;

    //! Delete move constructor.
    CD2D(CD2D&&) = delete;

    //! Delete copy assignment operator.
    CD2D& operator=(const CD2D&) = delete;

    //! Delete move assignment operator.
    CD2D& operator=(CD2D&&) = delete;

    //! Destructor. Still leaks memory (esp. because of the multigrid objeect).
    ~CD2D();

  private:
    /**
     * Apply an algebraic flux correction scheme to the assembled matrix.
     * Should be called within the assemble routine, after the actual assembling
     * has been performed with the INTERNAL_FULL_MATRIX_STRUCTURE switch on.
     *
     * Which afc algorithm is performed is determined by switching over
     * ALGEBRAIC_FLUX_CORRECTION.
     */
    void do_algebraic_flux_correction();
};

#endif // __CD2D_H__
