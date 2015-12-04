#ifndef __BLOCKMATRIX__
#define __BLOCKMATRIX__

#include <Matrix.h>
#include <BlockPattern.h>
#include <memory>

class BlockVector;

/** ************************************************************************ 
*
* @class      BlockMatrix
* @brief      represent a matrix consisting of blocks which are sparse matrices
*
*             This is a purely algebraic object. And one which does frequently
*             appear in FEM. Identifying and exploiting a certain block structure
*             of matrices is a fundamental concept of solvers in CFD.
*             
*             \todo The blockmatrix does not follow the rule of 0/5 yet.
*             The unresolved problem is that it does own a vector of pointers
*             to its blocks, but these are only shallow copied in the copy
*             and move constructor. If we implement a copy constructor that
*             makes deep copies of the objects as TMatrix, we must make sure
*             in the derived classes, that their copy constructors do not call
*             the base class copy constructor (high overhead!) but instead perform
*             their own copy or move actions along with their knowledge about the actual
*             type ("FEMatrix") of the matrices.
*
* @author     Ulrich Wilbrandt
* @date       08.09.2015
*
****************************************************************************/
class BlockMatrix
{
  protected:
    /** @brief the block pattern, determining e.g. the number of blocks
     *
     *  A matrix with a fixed block pattern does not change it again,
     *  this is why we can store a shared_ptr to a const object here.
     */
    std::shared_ptr<const BlockPattern> block_pattern;

    /** @brief the blocks, each is a sparse matrix (consider vector<vector<>>)
     * 
     * We store pointers to sparse matrices because we allow constructors 
     * without setting the blocks. Then they are initially set to nullptr. 
     * Maybe those constructors are not really needed, then we could store 
     * objects rather than pointers in this vector.
     */
    std::vector<std::shared_ptr<TMatrix>> blocks;
    
    /** @brief number of active entries for each block
     * 
     * The size of this vector is the same as the size of the vector 'blocks'.
     * The i-th block in 'blocks' has 'actives[i]' active entrires. This enables
     * the methods like scale_active, and add_scaled_active.
     * 
     * If you use a constructor which does not create the matrices in 'blocks',
     * the default value in this vector will be
     * std::numeric_limits<unsigned int>::max().
     */
    std::vector<unsigned int> actives;
    
    /** @brief all blocks as one TMatrix
     * 
     * This object is only created upon request. If there is only one block then
     * it is combined_matrix = blocks[0] set.
     */
    std::shared_ptr<TMatrix> combined_matrix;
    
    /**
     * @brief find the index of the block where the (i,j)-th entry is at
     *
     * the indices i and j are replaced by the local indices of the block to
     * which (i,j) originally belonged to. This enables the following:
     *
     * unsigned int i,j; // set to some admissible value
     * unsigned int bI = this->blockOfIndex(i,j);
     * double v = (this->blocks[bT])(i,j);
     *
     * Then the value stored in 'v' is the value of this BlockMatrix at the 
     * indices (i,j) which were initially set.
     *
     * @param i row index
     * @param j column index
     */
    unsigned int block_of_index(unsigned int& i, unsigned int& j) const;
    
  public:
    /** @brief Creates an empty BlockMatrix without any properties.
     */
    BlockMatrix();
    
    /**
     * @brief Creates an empty BlockMatrix with nRows times nCols blocks.
     *
     * @param nRows - number of blocks per column
     * @param nCols - number of blocks per row
     */
    BlockMatrix(unsigned int nRows, unsigned int nCols);
    
    /**
     * @brief Creates a nRows times nCols BlockMatrix filled with blocks
     *
     * @param nRows - number of blocks per column
     * @param nCols - number of blocks per row
     * @param new_blocks - the blocks as a vector
     */
    BlockMatrix(unsigned int nRows, unsigned int nCols, 
                std::vector<std::shared_ptr<TMatrix>> new_blocks);
    
    /** @brief construct a BlockMatrix suitable for the given problem type
     * 
     * All blocks contain zeros only.
     * 
     * This constructor essentially creates a BlockPattern using the given
     * arguments, and then does the same as the constructor using only a 
     * BlockPattern.
     */
    BlockMatrix(const Problem_type, unsigned int space_dimension,
                bool mass_matrix = false);
    
    /** @brief construct a BlockMatrix using a given BlockPattern
     * 
     * All blocks contain zeros only.
     */
    BlockMatrix(std::shared_ptr<const BlockPattern>);

    /** @brief copy constructor */
    BlockMatrix(BlockMatrix&);
    
    /** @brief move constructor */
    BlockMatrix(BlockMatrix&&);
    
    /// @brief destructor, deleting every block
    ~BlockMatrix() noexcept;
    
    /** @brief Set all submatrices to zero
     * 
     * Possibly existing special matrices are not changed.
     */
    void reset();
    
    /** 
     * @brief adding a scaled matrix to this matrix
     * 
     * The summation is index-wise, i.e. M(i,j) += factor*A(i.j), where M is 
     * this matrix. 
     * 
     * Note that this only works if the sparsity structure is the same for this
     * matrix and A.
     * 
     * Possibly existing special matrices are not changed.
     */
    void add_scaled(const BlockMatrix &A, double factor = 1.0);
    
    /** 
     * @brief adding a scaled matrix to this matrix, but only the active entries
     * 
     * This does exactly the same as add_scaled, except that nonactive entries
     * are not changed.
     */
    void add_scaled_active(const BlockMatrix &A, double factor = 1.0);
    
    /** 
     * @brief scale this matrix
     * 
     * That means for each submatrix all entries are scaled.
     * 
     * Possibly existing special matrices are not changed.
     */
    void scale(double factor);
    
    /** 
     * @brief scale all active entries of this matrix
     * 
     * That means for each submatrix all active entries are scaled.
     * 
     * Possibly existing special matrices are not changed.
     */
    void scale_active(double factor);
    
    /** @brief compute y = Ax 
     *
     * write the matrix-vector product "Ax" to y if "A" is this matrix.
     *
     * @param x the BlockVector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication
     */
    void apply(const BlockVector & x, BlockVector & y) const;
    
    /** @brief compute y = y + a * Ax 
     *
     * add the matrix-vector product "Ax", scaled by "a", to y if "A" is this 
     * matrix.
     * 
     * This function can be used to compute the residual r = b - Ax, for example
     * BlockVector r(b);  // structural copy (no values copied yet)
     * r = b;             // copy values from b to r
     * A.apply_scaled_add(x, r, -1.0);
     * cout << "Norm of residual " << r.norm() << endl;
     *
     * @param x the BlockVector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication
     * @param a optional factor, default to 1.0
     */
    void apply_scaled_add(const BlockVector & x, BlockVector & y,
                          double a = 1.0) const;
    
    /** @brief return this BlockMatrix as one TMatrix
     * 
     * This returns a merged version of this martix. That means this matrix then
     * exists twice, as blocks and as a combined matrix. Changes to one of them
     * does not affect the other.
     * 
     * Usually this is used to pass this matrix to a solver.
     */
    std::shared_ptr<TMatrix> get_combined_matrix();
  
    /* getter/setter functions */
    std::shared_ptr<const BlockPattern> get_block_pattern() const
    { return block_pattern; }
    double & operator()(unsigned int i, unsigned int j);
    const double & operator()(unsigned int i, unsigned int j) const;
    unsigned int n_blocks() const { return block_pattern->n_blocks(); }
    unsigned int n_rows() const { return block_pattern->n_rows(); }
    unsigned int n_cols() const { return block_pattern->n_cols(); }
    /// @brief total number of rows (> n_block_rows)
    unsigned int n_total_rows() const; 
    /// @brief total number of columns(> n_block_columns)
    unsigned int n_total_cols() const;
     /// @brief total number of entries
    unsigned int n_total_entries() const;
    std::shared_ptr<const TMatrix> block(const unsigned int i) const;
    std::shared_ptr<TMatrix> block(const unsigned int i);
    /** @brief return the TMatrix located in the r-th block row and c-th block
     *         column
     */
    const TMatrix& block(const unsigned int r, const unsigned int c) const;
    
    /** @brief print some information on this BlockMatrix */
    void info(size_t verbose) const;
};



#endif //__BLOCKMATRIX__
