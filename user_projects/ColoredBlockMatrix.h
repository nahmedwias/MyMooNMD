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
#include <tuple>

class BlockVector;

/** ************************************************************************
 *
 * @class      ColoredBlockMatrix
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
class ColoredBlockMatrix
{
  public:

    /** @brief Creates an empty ColoredBlockMatrix without any properties.
     */
    ColoredBlockMatrix();

    /**
     * @brief Creates a ColoredBlockMatrix which is filled with fitting zero blocks.
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
    ColoredBlockMatrix(std::vector<size_t> cell_row_numbers, std::vector<size_t> cell_column_numbers);

    /**
     * @brief checks whether the coloring is correct - use in tests only
     *
     * The method checks, whether the coloring of the cells is correct, e.g.
     *
     *  - the set of assigned colors equals the set {0,...,n_colors_ -1}
     *
     *  - on the set of cells the equivalence relation "has the same color as"
     *    is equivalent to "holds a pointer to the same matrix as"
     *
     *  - the colors are ordered in such a way that from left to right, top to bottom
     *    the colors of first appearing matrices are in ascending order.
     *
     *    Does nothing but throw if one of the rules stated above is broken.
     *    Is written with comprehensibility in mind, not performance. Use
     *    for testing purpose only.
     *
     *    @note Please note that this method scales quadratically in the number of blocks
     *    as it employs checkEquivalenceOfRelations. So please do only use it for testing
     *    of other methods in this class.
     */
    void check_coloring() const;

    /*
     * Prints matrix coloring pattern and color_count_,
     * then runs the check. Thus one can visualize errors.
     */
    void print_and_check() const;

    /*
     * Print out a little picture of the current coloring pattern of the matrix.
     * Does not perform a check on consistency - use it for debugging!
     */
    void print_coloring_pattern() const;

    /*
     * Print out the current color count of the matrix.
     * Does not perform a check on consistency - use it for debugging!
     */
    void print_coloring_count() const;


    //! A datatype which stores the block_row the block_column
    //! and the mode where and how a cell should be modified.
    //! 'Mode' can have to values: true - 'transposed' and
    //! false - 'not transposed'
    typedef std::tuple<size_t, size_t, bool> grid_place_and_mode;

    /*
     * Add one matrix to several blocks at once.
     *
     * @param[in]  The matrix to be added.
     * @param[out] A vector of tuples, each of which encodes information on
     *             one block where summand is supposed to be added to,
     *             and the way - "true" is for transposed,
     *             "false" for not-transposed.
     *
     * @note If anyone has a better idea how to realize the input than with
     * a vector of tuples, that would be very welcome...
     */
    void add_matrix_to_blocks(
        const TMatrix& summand,
        std::vector<grid_place_and_mode> row_column_transpose_tuples);

    void add_scaled_matrix_to_blocks(
        const TMatrix& summand, double scaling_factor,
        std::vector<grid_place_and_mode> row_column_transpose_tuples);

    /*
     * Comment
     */
     void scale_blocks( double scaling_factor,
                        std::vector<grid_place_and_mode> row_column_transpose_tuples);

     /*
      * Replaces the blocks whose positions are given in grid_places by a copy of
      * new_block, as long as that fits into the given grid postions
      *
      * Expects conditions 1, 2 and 3 to hold and maintains them.
      *
      * @param[in] new_block The new block to be inserted.
      *
      * @param[in] row_column_transpose_tuples The places where the new block should to be inserted.
      */
     void replace_blocks(
         const TMatrix& new_block,
         std::vector<grid_place_and_mode> row_column_transpose_tuples);

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
      * Then the value stored in 'v' is the value of this ColoredBlockMatrix at the
      * indices (i,j) which were initially set.
      *
      * @param i row index
      * @param j column index
      */
     unsigned int block_of_index(unsigned int& i, unsigned int& j) const;

     // Getter and setter

     //! @return The number of different colors
     size_t get_n_colors() const
     {
       // color_count_ must be kept up-to-date
       return color_count_.size();
     }

     //! @return The number of block rows.
     size_t get_n_block_rows() const
     {
       return n_block_rows_;
     }

     //! @return the number of block columns.
     size_t get_n_block_columns() const
     {
       return n_block_columns_;
     }


     // Special member functions.

     /** @brief copy constructor */
     ColoredBlockMatrix(ColoredBlockMatrix&) = default;

     /** @brief move constructor */
     ColoredBlockMatrix(ColoredBlockMatrix&&) = default;

     /// @brief destructor, deleting every block
     ~ColoredBlockMatrix() = default;

     // Data members (and declaration of a special struct)


  protected:

     //! Store information of a certain grid cell in the block matrix.
     //! Will by default perform a shallow copy when copied, which is what we want.
     struct CellInfo
     {
       // todo Consider encapsulating n_rows_ and n_columns_

       // the number of rows a matrix stored here must possess
       size_t n_rows_;

       // the number of columns a matrix stored here must possess
       size_t n_columns_;

       //! The TMatrix currently associated with this cell.
       std::shared_ptr<TMatrix> block_;

       /*!
        * The color of the cell - std::numeric_limits<size_t>::max()
        * stands for 'uncolored'.
        */
       size_t color_;

       //! whether or not the stored block is viewed as transposed
       bool is_transposed_;

       //! A flag used in algorithms which change the coloring scheme
       // of the owning ColoredBlockMatrix
       enum class ReColoringFlag {SPLIT, MERGE, KEEP} re_color_;

       /*! Default constructor, will be called implicitely. Constructs
        *  a non-transposed, uncolored 0x0 object marked with ReColoringFlag:KEEP.
        */
       CellInfo();

       /*! Constructor which initializes n_rows and n_columns,
        *  Behaves like default parameter apart from that.
        */
       CellInfo(size_t n_rows, size_t n_columns);

     };

     //! The number of block rows. Equals the number of blocks/cells in each block column.
     size_t n_block_rows_;

     //! The number of block columns. Equals the number of blocks/cells in each block row.
     size_t n_block_columns_;

     /*!
      *  A tableau of cell information. Has dimension
      *   n_block_rows x n_block_columns, and stores cell information.
      *   I.e. cell_info_grid[i][j] holds information on the cell in block row i,
      *   block column j.
      */
     std::vector<std::vector<CellInfo>> cell_info_grid_;

     /*!
      * The size of the vector is the number of colors currently
      * existing in the block matrix, i.e. the number of physically
      * stored TMatrices.
      *
      * Storing these numbers enables quick checking whether a
      * modification of the blocks will require a modification
      * of the coloring scheme.
      *
      * color_count_[i] gives you the current number of cells
      * with color i in the block matrix
      */
     std::vector<size_t> color_count_;

     /** @brief all blocks as one TMatrix
      *
      * This object is only created upon request. If there is only one block then
      * it is combined_matrix = blocks[0] set.
      *
      * Consider not storing but making up and returning on request.
      */
     std::shared_ptr<TMatrix> combined_matrix_;


  private:

     /*
      * Checks an input vector of places and edits it:
      *   - sorts in alphanumeric order of grid places
      *   - removes duplicates
      *   - throws if the input contains duplicates of the type
      *     (a,b,true) , (a, b, false) - mode at place has to
      *     be unique!
      *
      * @param[in,out] The input vector of grid places and modes.
      *
      */
     static void check_and_edit_input(
         std::vector<grid_place_and_mode>& row_column_transpose_tuples);

     /*
      * Check whether the set of assigned colors equals the set
      * {0,...,n_colors_ -1}. If not so - just throw.
      */
     void check_color_count() const;

     /**
      * This method checks whether the coloring fulfils the internal ordering.
      * If not so it just throws.
      *
      * The method relies on conditions 1 and 2 being established.
      */
     void check_coloring_order() const;

     /*
      * Check whether on the set of cells the equivalence relation
      * "has the same color as" is equivalent to "holds a pointer
      * to the same matrix as". If not so - just throw.
      *
      * The result is independent of the check which
      * checkColorAssignment performs.
      *
      * @note That this method scales quadratically in the number of blocks,
      * so by n^4 for a nxn block matrix - so please do only use it if you
      * a) want to see your computer failing or
      * b) test other methods in this class.
      */
     void check_equivalence_of_relations() const;

     /*
      * Checks whether a matrix fits to the cells and mode specified in
      * row_column_transpose_tuples.
      * Throws if not so.
      *
      * @param[in] matrix The matrix whose fit into the grid shall be checked.
      *
      * @param[in] row_column_transpose_tuples The grid places and modes to check.
      * The method check_and_edit_input is called upon it to make sure it is
      * valid input.
      */
     void check_grid_fit(
         const TMatrix& matrix,
         std::vector<grid_place_and_mode>& row_column_transpose_tuples) const;

     /**
      * This method checks whether all cells have their re-coloring flags set to KEEP.
      * If not so it just throws.
      */
     void check_re_coloring_flags() const;

     /*!
      * Check if a given block fits into a given cell in the given transposed state.
      *
      */
     static bool does_block_fit_cell(
         const TMatrix& block, const CellInfo& cell, bool transposed );

     /*!
      * Check if for two cells the relations "has the same color as"
      * and "store a pointer to the same matrix" produce the same result.
      *
      * @param[in] first  the first CellInfo object to compare
      * @param[in] second the second CellInfo object to compare
      */
     static bool does_color_match_block(const CellInfo& first, const CellInfo& second);

     /**
      * Checks if a modification which affects a set of cells requires at least
      * one color class to be split in two.
      *
      * @return True if the modification requires to split at least one color in two.
      *
      * @param[out] color_to_split The lowest color which has to be split in two.
      * @param[in] row_column_tuples The index pairs of the cells affected by
      * the modification. Must not contain the same index pair twice! Need not be sorted.
      *
      * When returning true the parameters of the method should be handed over to
      * mark_for_color_split and then to split_color
      */
     bool does_modification_require_color_split(
         size_t& color_to_split,
         std::vector<grid_place_and_mode> row_column_transposed_tuples) const;

     /**
      * @return True if a block emplacement leads to the merging of at least two colors.
      *
      * The method marks all the CellInfos of the color that has to be split and which do belong
      * to the list of grid_place s (that part which is supposed to receive the replacement)
      * with ReColoringFlag::MERGE - which is then used by the actual color merging method.
      *
      * TODO Implement!
      */
     bool does_replacement_require_color_merge(std::vector<grid_place_and_mode> row_column_tuples);

     /*!
      * Find the first place in the matrix (oredr from left to right, top to bottom)
      * where a certain color appears.
      *
      * @param[in] color_to_find the color to find
      *
      * @param[out] block_row the row where the first such block was found
      * @param[out] block_column the column where the first such block was found
      */
     void find_first_appearance_of_color(size_t color_to_find,
                                         size_t& block_row , size_t& block_column) const;

     /**
      * Get the next place in the grid after block_row, block_column.
      * Throws an std::logic_error exception if the last index pair
      * is given as input.
      *
      * @param[in, out] block_row the current block row
      * @param[in, out] block_column the current block column
      *
      * @throws std::logic_error exception if the last index pair
      * is given as input
      */
     void get_next_cell_grid_index(size_t& block_row, size_t& block_column) const;

     /*
      * Check if a given index pair is the last one in the cell_info_grid_.
      *
      * @param[in] block_row the current block row
      * @param[in] block_column the current block column
      *
      * @return true if this is the index pair of the lower right corner
      */
     bool is_last_index_pair(size_t block_row, size_t block_column) const;

     /*
      * Mark every cell of color_to_mark whose position appears in row_column_tuples
      * with the recoloring flag "SPLIT". After that, call split_color to perform the
      * actual splitting of the color.
      *
      * @param[in] color_to_mark Cells of color color_to_mark are marked when their
      * position belongs to row_column_tuples.
      * @param[in] row_column_transposed_tuples Cells of color color_to_mark are marked when their
      * position belongs to row_column_tuples.
      */
     void mark_for_color_split(size_t color_to_mark,
                               std::vector<grid_place_and_mode> row_column_transposed_tuples);

     /*
      * Merges all those colors marked with ReColoringFlag::MERGE into one.
      * Assumes conditions 1, 2 and 3 to hold and maintains them.
      *
      * TODO Implement!
      */
     void merge_colors();

     /*
      * Splits a color in two, assuming that conditions 1,2 and 3 hold and maintaining them.
      * Also relies on one part of the cells of the color to be split are marked with
      * ReColoringFlag::SPLIT. (Make sure to call mark_for_color_split in advance).
      *
      * @param[in] color_to_split The color which has to be split in two
      */
     void split_color(size_t color_to_split);


     // not yet adapted members of BlockMatrix

  public:

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
     void add_scaled(const ColoredBlockMatrix &A, double factor = 1.0);



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

     /** @brief return this ColoredBlockMatrix as one TMatrix
      *
      * This returns a merged version of this martix. That means this matrix then
      * exists twice, as blocks and as a combined matrix. Changes to one of them
      * does not affect the other.
      *
      * Usually this is used to pass this matrix to a solver.
      */
     std::shared_ptr<TMatrix> get_combined_matrix();


     /// @brief total number of rows (> n_block_rows)
     unsigned int n_total_rows() const;
     /// @brief total number of columns(> n_block_columns)
     unsigned int n_total_cols() const;
     /// @brief total number of entries
     unsigned int n_total_entries() const;

     /** @brief return the TMatrix located in the r-th block row and c-th block
      *         column
      */
     const TMatrix& block(const unsigned int r, const unsigned int c) const;

     /** @brief print some information on this ColoredBlockMatrix */
     void info(size_t verbose) const;
};



#endif /* USER_PROJECTS_COLOREDBLOCKMATRIX_H_ */
