#ifndef __SYSTEM_CONVECTIONDIFFUSION_H__
#define __SYSTEM_CONVECTIONDIFFUSION_H__

#include "BlockFEMatrix.h"
#ifdef __2D__
#include "Example_CD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_CD3D.h"
#include "FEFunction3D.h"
#endif
#include "BlockVector.h"
#include "ParameterDatabase.h"
#include "Solver.h"
#include "templateNames.h"
#include "DataWriter.h"

#include <deque>
#include <array>

template<int d>
class ConvectionDiffusion
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_CD = typename Template_names<d>::Example_CD;
  
  protected:
    
    /** @brief store a complete system on a particular grid
     * 
     * This combines a matrix, rhs, solution, space and function needed to 
     * describe one convection-diffusion-reaction problem in 2D.
     */
    struct SystemPerGrid
    {
      /** @brief Finite Element space */
      std::shared_ptr<const FESpace> fe_space;

      /** @brief the system matrix */
      BlockFEMatrix matrix;
      /** @brief the right hand side vector */
      BlockVector rhs;
      /** @brief solution vector with one component. */
      BlockVector solution;
      /** @brief Finite Element function */
      FEFunction fe_function;
      
      /** @brief constructor */
      SystemPerGrid(const Example_CD& example, TCollection& coll,
                    int ansatz_order);

      // Special member functions. Disable copy/move, set destructor to default.

      //! Delete copy constructor.
      SystemPerGrid(const SystemPerGrid&) = delete;

      //! Delete move constructor.
      SystemPerGrid(SystemPerGrid&&) = delete;

      //! Delete copy assignment operator.
      SystemPerGrid& operator=(const SystemPerGrid&) = delete;

      //! Delete move assignment operator.
      SystemPerGrid& operator=(SystemPerGrid&&) = delete;

      //! Default destructor.
      ~SystemPerGrid() = default;
    };
    
    /** @brief a complete system on each grid 
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<SystemPerGrid> systems;
    
    /** @brief Definition of the used example */
    const Example_CD example;
    
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
    DataWriter<d> outputWriter;

    /** @brief a solver object which will solve the linear system
     * 
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;
    
    static constexpr int n_errors = 4;
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
    std::array<double, n_errors> errors;
    
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
    
    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;
    
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
     */
    ConvectionDiffusion(const TDomain& domain,
                        const ParameterDatabase& param_db);
    
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
     */
    ConvectionDiffusion(const TDomain& domain,
                        const ParameterDatabase& param_db,
                        const Example_CD& example);
    
    /** @brief get a database with all parameters needed by this class, 
     *  initialized with default values.
     */
    static ParameterDatabase get_default_cd_database();

    /** @brief assemble matrix, 
     * 
     * depending on 'this->db["space_discretization_type]' different (local)
     * assembling routines are used. Also in case of multigrid the matrices
     * on all grids are assembled.
     */
    void assemble();
    
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
    const FEFunction & get_function() const
    { return this->systems.front().fe_function; }
    
    const FESpace & get_space() const
    { return *this->systems.front().fe_space; }
    
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    BlockVector & get_solution()
    { return this->systems.front().solution; }
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    const Example_CD& get_example() const
    { return example; }
    const ParameterDatabase & get_db() const
    { return db; }

    // Special member functions. Disable copy/move, set destructor to default.
    // Will be changed only when the underlying classes follow rule of 0/5.

    //! Delete copy constructor.
    ConvectionDiffusion(const ConvectionDiffusion&) = delete;

    //! Delete move constructor.
    ConvectionDiffusion(ConvectionDiffusion&&) = delete;

    //! Delete copy assignment operator.
    ConvectionDiffusion& operator=(const ConvectionDiffusion&) = delete;

    //! Delete move assignment operator.
    ConvectionDiffusion& operator=(ConvectionDiffusion&&) = delete;

    //! Destructor.
    ~ConvectionDiffusion() = default;

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

#endif // __SYSTEM_CONVECTIONDIFFUSION_H__
