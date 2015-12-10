/*
 * ColoredBlockMatrix.h
 *
 *  Created on: Dec 10, 2015
 *      Author: bartsch
 */

#ifndef USER_PROJECTS_COLOREDBLOCKMATRIX_H_
#define USER_PROJECTS_COLOREDBLOCKMATRIX_H_

#include <Matrix.h>
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
*
* @author     Ulrich Wilbrandt, Clemens Bartsch
* @date       08.09.2015
*
****************************************************************************/
class BlockMatrix
{
  //! Store information of a certain grid cell in the block matrix.
  //! Will by default perform a shallow copy when copied, which is what we want.
  struct CellInfo
  {
      const size_t n_rows_;

      const size_t n_columns_;

      shared_ptr<TMatrix> block_in_cell_;

      size_t color_;

      bool is_transposed_;

      //! Default constructor, will be called implicitely.
      CellInfo();

      //! Constructo which sets the constant values correctly.
      CellInfo(size_t n_rows, size_t n_columns);

  };

  protected:

    //! The number of block rows. Equals the number of blocks/cells in each block column.
    size_t n_block_rows_;

    //! The number of block columns. Equals the number of blocks/cells in each block row.
    size_t n_block_columns_;

    /*!
     * A tableau of cell information. Has dimension
    *   n_block_rows x n_block_columns, and stores cell information.
    *   I.e. cell_info_grid[i][j] holds information on the cell in block row i,
    *   block column j.
    */
    std::vector<std::vector<CellInfo>> cell_info_grid_;

    /*! The number of colors currently existing in the block matrix,
     *  i.e. the number of physically stored TMatrices.
     */
    size_t n_colors_;

    /** @brief all blocks as one TMatrix
     *
     * This object is only created upon request. If there is only one block then
     * it is combined_matrix = blocks[0] set.
     *
     * Consider not storing but making up and returning on request.
     */
    std::shared_ptr<TMatrix> combined_matrix_;


  public:

    /**
     * @brief checks whether the coloring is correct - use in tests only
     *
     * The method checks, whether the coloring of the cells is correct, e.g.
     *
     *  - for each color from 0 to n_colors - 1 there is at least one cell
     *  in the matrix with that color
     *
     *  - no cell has color greater or equal n_colors
     *
     *  - all cells with the same color hold a pointer to the same matrix
     *
     *  - cells with different colors hold pointers to different matrices
     *
     *  - the colors are ordered in such a way that from left to right, top to bottom
     *    the colors of first appearing matrices are in ascending order.
     *
     *    Does nothing but throw if one of the rules stated above is broken.
     *    Is written with comprehensibility in mind, not performance. Use
     *    for testing purpose only.
     */
    void checkColoring();

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

    // Getter and setter

    //! @return The number of block rows.
    size_t get_n_block_rows() const {
      return n_block_rows_;
    }

    //! @return the number of block columns.
    size_t get_n_block_columns() const {
      return n_block_columns_;
    }


    // Special member functions.

    /** @brief copy constructor */
    BlockMatrix(BlockMatrix&);

    /** @brief move constructor */
    BlockMatrix(BlockMatrix&&);

    /// @brief destructor, deleting every block
    ~BlockMatrix() noexcept;

  public:
    /** @brief Creates an empty BlockMatrix without any properties.
     */
    BlockMatrix();

    /**
     * @brief Creates a BlockMatrix which is filled with fitting zero blocks.
     *
     * The number of block rows is the length of cell_row_numbers,
     * the number of block columns the length of cell_column_numbers.
     *
     * If e.g. cell_row_numbers[3] == 100, all cells in block row 3 will have
     * 100 matrix rows.
     *
     * @param cell_row_numbers holds the number of rows for all matrices in one block row
     * @param cell_column_numbers holds the number of columns for all matrices in one block column
     */
    BlockMatrix(std::vector<size_t> cell_row_numbers, std::vector<size_t> cell_column_numbers);

//    /**
//     * @brief Creates a nRows times nCols BlockMatrix filled with blocks
//     *
//     * @param nRows - number of blocks per column
//     * @param nCols - number of blocks per row
//     * @param new_blocks - the blocks as a vector
//     */
//    BlockMatrix(unsigned int nRows, unsigned int nCols,
//                std::vector<std::shared_ptr<TMatrix>> new_blocks);
//
//    /** @brief construct a BlockMatrix suitable for the given problem type
//     *
//     * All blocks contain zeros only.
//     *
//     * This constructor essentially creates a BlockPattern using the given
//     * arguments, and then does the same as the constructor using only a
//     * BlockPattern.
//     */
//    BlockMatrix(const Problem_type, unsigned int space_dimension,
//                bool mass_matrix = false);



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
     * @brief scale this matrix
     *
     * That means for each submatrix all entries are scaled.
     *
     * Possibly existing special matrices are not changed.
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

    double & operator()(unsigned int i, unsigned int j);
    const double & operator()(unsigned int i, unsigned int j) const;
    unsigned int n_blocks() const { return block_pattern->n_blocks(); }

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



#endif /* USER_PROJECTS_COLOREDBLOCKMATRIX_H_ */
