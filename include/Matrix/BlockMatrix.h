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
*             This is a purely algebraic object.
*             
* @author     Ulrich Wilbrandt
* @date       08.09.2015
*
****************************************************************************/
class BlockMatrix
{
  protected:
    /// @brief the block pattern, determining e.g. the number of blocks
    std::shared_ptr<BlockPattern> block_pattern;
    /** @brief the blocks, each is a sparse matrix (consider vector<vector<>>)
     * 
     * We store pointers to sparse matrices because we allow constructors 
     * without setting the blocks. Then they are initially set to nullptr. 
     * Maybe those constructors are not really needed, then we could store 
     * objects rather than pointers in this vector.
     */
    std::vector<std::shared_ptr<TMatrix>> blocks;
    
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
     * This constructor essentially creates a BlockPattern using the 
     * Problem_type and the finite element spaces, and then does the same as 
     * the constructor using only a BlockPattern.
     */
    BlockMatrix(const Problem_type, unsigned int space_dimension);
    
    /** @brief construct a BlockMatrix using a given BlockPattern
     * 
     * All blocks contain zeros only.
     */
    BlockMatrix(std::shared_ptr<BlockPattern>);
    
    /** @brief copy constructor */
    BlockMatrix(BlockMatrix&);
    
    /** @brief move constructor */
    BlockMatrix(BlockMatrix&&);
    
    /// @brief destructor, deleting every block
    ~BlockMatrix() noexcept;
    
    /** @brief Set all submatrices to zero
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
     */
    void add_scaled(const BlockMatrix &A, double factor = 1.0);
    
    /** 
     * @brief scale this matrix
     * 
     * That means for each submatrix all entries are scaled.
     */
    void scale(double factor);
    
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
