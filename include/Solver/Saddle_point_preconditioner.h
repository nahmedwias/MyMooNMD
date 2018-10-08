#ifndef __SADDLE_POINT_PRECONDITIONER__
#define __SADDLE_POINT_PRECONDITIONER__

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <DirectSolver.h>
#include <Solver.h>
#include <Preconditioner.h>
#ifdef _MPI
#include <MumpsWrapper.h>
#endif

class ParameterDatabase;

/** @brief implement special preconditioners for saddle point problems
 * 
 * @Note Up to now this class is not very nicely written. This needs some 
 * refactoring.
 */
class Saddle_point_preconditioner : public Preconditioner<BlockVector>
{
  public:
    /// simple - semi-implicit method for pressure-linked equations
    /// lsc - least squares commutator
    /// bd_lsc - boundary corrected least squares commutator
    /// AL - augmented Lagrangian based preconditioner
	  /// mod_AL - modified augmented Lagrangian based preconditioner
    enum class type {simple, lsc, bd_lsc, AL, mod_AL};
    /// The database passed to the constuctor should be named this or have a 
    /// nested database with this name. Otherwise a default solver database is 
    /// used to initialize a solver object for the velocity subsystem.
    constexpr static char required_database_name[] = 
      "Saddle Point Preconditioner Database";
    
    /** @brief constructor for a given system matrix */
    explicit Saddle_point_preconditioner(const BlockFEMatrix & m,
                                         type t, const ParameterDatabase& db, const BlockVector& rhs = BlockVector());
    /** @brief don't use this constuctor. It is here only for compatability 
     * in the Solver class.
     * 
     * @warning Do not use this constructor. You need a BlockFEMatrix instead.
     */
    explicit Saddle_point_preconditioner(const BlockMatrix & m,
                                         type t, const ParameterDatabase& db);
    
    /** @brief destructor, delete all allocated memory */
    ~Saddle_point_preconditioner() = default;
    
    /** "brief update the members in this class after a change in the matrix
     * 
     * The idea is that during a nonlinear iteration only the velocity blocks 
     * change, not the blocks involving pressure, therefore much of what has been
     * computed in the constructor can be used without modification. This is true
     * for LSC and for SIMPLE. In particular is saves the costly computation of 
     * B*B^T.
     */
    void update();
    
    /** @brief applying the preconditioner */
    void apply(const BlockVector &z, BlockVector &r) const;
    
    /** @brief Method to apply the preconditioner in flexible gmres.
     *
     * So far i and j get ignored and simply apply(z,r) is called.
     *
     * @param i Number of current iteration since last restart in flexible gmres
     * @param j Number of current iteration in flexible gmres
     * @param z The right hand side of the preconditioning
     * @param r The obtained vector
     */
    void apply(int i, int j, const BlockVector &z, BlockVector &r) const
    {
      this->apply(z,r);
    }

    /* @brief get the augmentated system matrix in case of augmented Lagrangian preconditioner:
     * A_ij -> A_ij + gamma * B_i^T * W^{-1} * B_j */
    BlockMatrix& get_augmented_matrix()
    {
    	return this->augmented_matrix;
    }

    /* @brief get the augmentated rhs in case of augmented Lagrangian preconditioner:
     * (f_1,f_2,f_3,g) -> (f_1 + gamma * B_1^T * W^{-1} g, f_2 + gamma * B_2^T * W^{-1} g, g) */
    BlockVector get_augmented_blockvector(BlockVector right_hand_side);
  
  protected:
    // saddle point preconditioner (spp) type
    type spp_type;
    
    /** @brief the system matrix
     * 
     * This preconditioner should be an approximation to the inverse of M. We
     * do not handle any memory here. The matrix is not used in the `solve` 
     * method in this class. 
     */
    const BlockFEMatrix * M;
    
    /** @brief all velocity blocks as one TMatrix */
    BlockFEMatrix velocity_block;
   
  /** @brief pressure block (C) as one BlockFEMatrix copied from the system matrix*/
    BlockFEMatrix pressure_block;
   
    /** @brief pressure mass matrix (p_h,q_h) (for W in the augmented Lagrangian based preconditioner) */
    BlockFEMatrix pressure_mass;

    /** @brief the augmentation for the augmented Lagrangian based preconditioner: gamma * B^T * W^{-1} * B   */
    BlockMatrix augmented_matrix;

    /** rectangular matrix which has to be applied to the pressure block in rhs in case of augmented Lagrangian (see Solver.C) */
    BlockMatrix augmentation_matrix_for_rhs;

    /** @brief the block which represents the gradient of the pressure 
     * 
     * It is scaled by the 'inverse_diagonal' during the constructor in case of
     * this->type == Saddle_point_preconditioner::type::simple.
     */
    std::shared_ptr<TMatrix> gradient_block;
    /** @brief the block which represents the divergence of the pressure */
    std::shared_ptr<TMatrix> divergence_block;
    
    /** @brief storing a factorization for the 'velocity_block'.
     *
     *  In MPI case, either the velocity_solver will be used (if an
     *  iterative velocity solver is requested) or the MUMPS wrapper
     *  object velocity_mumps_wrapper, in the case the velocity problem
     *  should be solved with a direct solver. The respective other
     *  object will be left empty.
     */
    std::shared_ptr<Solver<BlockFEMatrix, BlockVector>> velocity_solver;
#ifdef _MPI
    std::shared_ptr<MumpsWrapper> velocity_mumps_wrapper;
#endif
    
    /** @brief the inverse of the diagonal of the system matrix */
    std::vector<double> inverse_diagonal;
   
    /** @brief the scaled (with this->gamma) inverse of the diagonal of the system matrix */
    std::vector<double> scaled_inverse_diagonal_of_W;
    std::vector<double> diagonal_of_W;

    /** @brief references to the fe spaces 
     * 
     * These are really needed only for the boundary corrected LSC. I think 
     * otherwise we would be fine with a BlockMatrix instead of a BlockFEMatrix
     * as well.
     */
    #ifdef __2D__
    const TFESpace2D * velocity_space;
    const TFESpace2D * pressure_space;
    #else
    const TFESpace3D * velocity_space;
    const TFESpace3D * pressure_space;
    #endif
    
    /** @brief damping factor, zero means no preconditioner is used 
     * 
     * Currently this is always 1.0, i.e., no damping.
     */
    double damping_factor;
   
/** @brief gamma: augmentation parameter; gamma = 0 means no augmentation is used */
    double gamma;

    /* LSC === */
    /** @brief the Poisson solver matrix (the part of the approximation of the
     * Schur complement)
     */
    std::shared_ptr<BlockMatrix> Poisson_solver_matrix;
    
#ifdef _SEQ
    /** @brief storing a factorization for the 'Poisson_solver_matrix' */
    std::shared_ptr<DirectSolver> Poisson_solver;
#endif
#ifdef _MPI
    /** @brief storing a factorization for the 'Poisson_solver_matrix' */
    std::shared_ptr<MumpsWrapper> Poisson_solver;
#endif
    
    /** @brief a vector to store some intermediate results
     * 
     * This is stored here only to avoid allocation and destruction during each
     * call to solve.
     */
    mutable BlockVector up_star;
    mutable BlockVector u_star;
    mutable BlockVector p_star;
    mutable BlockVector u_tmp;
    mutable BlockVector p_tmp;
    
    /* Boundary corrected LSC */
    
    /** The diagonal weighting matrix which contains the boundary correction. 
     * We store H^{-1} = D*D_Q^{-1}, where D is the boundary weighting and 
     * D_Q^{-1} the diagonal version of the velocity mass matrix.
     */
    std::vector<double> bdryCorrectionMatrix_;
    
    /** The Poisson solver matrix which contains the boundary correction. We
     * store B*H^{-1}B^T, where H^{-1} is just bdryCorrectionMatrix_.
     */
    std::shared_ptr<TMatrix> poissonMatrixBdry_;
    
    /** Stores a factorization for the extra Poisson solver matrix 
     * poissonMatrixBdry_.
     *
     * @note Seems this can only be done via pointer, for DirectSolver causes
     * trouble when using default constructor and setMatrix(..) later */
    std::shared_ptr<DirectSolver> poissonSolverBdry_;
    
    
    // methods
    
    /** @brief return an approximation to the Poisson solver matrix */
    std::shared_ptr<BlockMatrix> compute_Poisson_solver_matrix() const;
    
    /** @brief fill the member Saddle_point_preconditioner::inverse_diagonal.
     * 
     * This involves assembling of a mass matrix.
     */
    void fill_inverse_diagonal();
   

        /** @brief fill the member Saddle_point_preconditioner::scaled_inverse_diagonal_of_W and
         * Saddle_point_preconditioner::diagonal_of_W
       */
        void fill_AL_weight_W();

    /** @brief fill the member Saddle_point_preconditioner::pressure_mass.
       * This involves assembling of a mass matrix.
       */
    void fill_pressure_mass_matrix();

    void fill_augmented_matrix_and_rhs();


    /** @brief Sets up bdryCorrectionMatrix_. This is 2D specific and must be redone for 3D.
     *  @param m A const reference to the whole system matrix.
     */
    void computeBdryCorrectionMatrix(const BlockFEMatrix& m);
    
    /**
     * @brief Sets up the extra Poisson matrix (poissonMatrixBdry_) needed in bdry corrected LSC.
     *
     * Call only after Poisson_solver_matrix has been constructed, work with a copy of its structure
     */
    void computePoissonMatrixBdry();

    /**
     * @brief Call the velocity solver to solve for rhs and sol.
     * In MPI case, this decides whether velocity_solver or
     * velocity_mumps_wrapper must be used.
     */
    void solve_velocity(const BlockVector& rhs, BlockVector& sol) const;

};


#endif // __SADDLE_POINT_PRECONDITIONER__
