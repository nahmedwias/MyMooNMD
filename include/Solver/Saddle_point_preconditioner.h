#ifndef __SADDLE_POINT_PRECONDITIONER__
#define __SADDLE_POINT_PRECONDITIONER__

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <DirectSolver.h>
#include <Preconditioner.h>

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
    enum class type {simple, lsc, bd_lsc};
    
    /** @brief constructor for a given system matrix */
    explicit Saddle_point_preconditioner(const BlockFEMatrix & m,
                                         type t = type::lsc);
    /** @brief don't use this constuctor. It is here only for compatability 
     * in the Solver class.
     * 
     * @warning Do not use this constructor. You need a BlockFEMatrix instead.
     */
    explicit Saddle_point_preconditioner(const BlockMatrix & m,
                                         type t = type::lsc);
    
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
    void update(const BlockFEMatrix & m);
    
    /** @brief don't use this method. It is here only for compatability 
     * in the Solver class.
     * 
     * @warning Do not use this method. You need a BlockFEMatrix instead.
     */
    void update(const BlockMatrix & m);
    
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
  
  protected:
    // saddle point preconditioner (spp) type
    type spp_type;
    
    // solution strategy within the LSC preconditioner
    // 0: direct (umfpack)
    // 1: fgmres (with ssor preconditioning), then use SC_LIN_RED_FACTOR_SCALAR
    // in case of 1 you need a flexible solver, because the number of iterations
    // is in general not constant.
    unsigned int lsc_strategy;
    
    
    /** @brief the system matrix
     * 
     * This preconditioner should be an approximation to the inverse of M. We
     * do not handle any memory here. The matrix is not used in the `solve` 
     * method in this class. 
     */
    const BlockFEMatrix * M;
    
    /** @brief all velocity blocks as one TMatrix */
    BlockFEMatrix velocity_block;
    /** @brief the block which represents the gradient of the pressure 
     * 
     * It is scaled by the 'inverse_diagonal' during the constructor
     */
    std::shared_ptr<TMatrix> gradient_block;
    /** @brief the block which represents the divergence of the pressure */
    std::shared_ptr<TMatrix> divergence_block;
    /** @brief the Schur complement matrix (possibly an approximation to it) */
    std::shared_ptr<TMatrix> Schur_complement;
    
    /** @brief storing a factorization for the 'velocity_block' */
    std::shared_ptr<DirectSolver> velocity_solver;
    /** @brief storing a factorization for the 'Schur_complement' */
    std::shared_ptr<DirectSolver> Schur_solver;
    
    /** @brief the inverse of the diagonal of the system matrix */
    std::vector<double> inverse_diagonal;
    
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
    
    /** @brief damping factor, zero means no preconditioner is used */
    double damping_factor;
    
    /* LSC === */
    /** @brief the Poisson solver matrix (the part of the approximation of the
     * Schur complement)
     */
    std::shared_ptr<BlockMatrix> Poisson_solver_matrix;
    
    /** @brief storing a factorization for the 'Poisson_solver_matrix' */
    std::shared_ptr<DirectSolver> Poisson_solver;
    
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
    
    /** The Poisson solver matrix which contains the boundary correction. We store
     * B*H^{-1}B^T, where H^{-1} is just bdryCorrectionMatrix_.
     */
    std::shared_ptr<TMatrix> poissonMatrixBdry_;
    
    /** Stores a factorization for the extra Poisson solver matrix poissonMatrixBdry_.
     *
     * @note Seems this can only be done via pointer, for DirectSolver causes
     * trouble when using default constructor and setMatrix(..) later */
    std::shared_ptr<DirectSolver> poissonSolverBdry_;
    
    
    // methods
    
    /** @brief return an approximation to the Schur complement matrix */
    std::shared_ptr<TMatrix> compute_Schur_complement_approximation() const;
    
    /** @brief return an approximation to the Poisson solver matrix */
    std::shared_ptr<BlockMatrix> compute_Poisson_solver_matrix() const;
    
    /** @brief fill the member Saddle_point_preconditioner::inverse_diagonal.
     * 
     * This involves assembling of a mass matrix.
     */
    void fill_inverse_diagonal();
    
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
    
    /** @brief solve the system involving only the velocity block 
     * 
     * Use LSC_STRATEGY to change the way how this is solved
     */
    void solve_velocity_block(const double* rhs, double* solution) const;
};


#endif // __SADDLE_POINT_PRECONDITIONER__