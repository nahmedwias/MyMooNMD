/** ************************************************************************
*
* @class      ColoredBlockFEMatrix
* @brief      extends ColoredBlockMatrix by handling active degrees of freedom
*
*             A BlockMatrix of subtype BlockFEMatrix stores FEMatrices
*             instead of simple algebraic TMatrices. Thus it can access
*             information on finite element spaces and active degrees of
*             freedoms and exploits this algoritmically. The majority of
*             BlockMatrices in ParMooN are in fact BlockFEMatrices.
*
*             Each cell row is associated with a certain FE test space
*             and each cell column with a certain FE ansatz space.
*
*             @note The possibility to store matrices as transposed
*             leads to a confusion of Test- and Ansatzspace between the
*             BlockFEMatrix and the stored FEMatrices. If an FEMatrix is
*             stored as transposed somewhere and one requests the testspace
*             of the matrix stored there, it gives a different space than
*             when requiring the testspacc of the cell where it is stored.
*             Because only the BlockFEMatrix knows, whether a block is stored
*             as tranposed or not, its "answer" is the correct one.
*             For the moment the only advice we can give is: do not request
*             test- or ansatzspace of an FEMatrix stored as a block in a
*             BlockFEMatrix, but instead request the test- or ansatzspace
*             of the cell where it is stored.
*
*             The current implementation of BlockFE Matrix has restrictions
*             due to the way ParMooN handles non-active degrees of freedom.
*             It is restricted to:
*               - a symmetrical test- and ansatzraum structure - i.e. the
*                 rowwise sequence of testspaces is THE SAME as the columnwise
*                 sequence of ansatzspaces
*               - whose non-active dofs stem only from Dirichlet boundaries,
*                 i.e. no hanging nodes! (which is the second reason for non-
*                 active degrees of freedom in ParMooN)
*             The reason for these restriction is in the specific handling of,
*             non-active dofs, i.e. dofs from Dirichlet boundary conditions and
*             hanging nodes. The issue is, that the corresponding rows are set in
*             each FEMatrix seperately, depending on its testspace. The FEMatrices
*             do not know about the space dimension they belong to, which leads
*             to a discrepancy between the BlockFEMatrix' global non-active row
*             and the FEMatrices' local non-active rows. We do not know how to entangle
*             that relation at the moment, only in the case given above it is quite
*             easily done: of all blocks in one row only the non-active entries of the
*             block on the diagonal are of relevance. This circumstance is paid
*             respect in the methods
*               - apply (only indirectly)
*               - apply_scaled_add
*               - check_vector_fits_pre_image (actives check)
*               - check_vector_fits_image (actives check)
*               - get_combined_matrix.
*             Should you intend to extend the class to the handling of hanging
*             nodes and/or non-symmetric test- and ansatzspaces, these are the
*             methods to changes (and of course implementing a new constructor
*             and/or remove the "no-hanging-nodes" invariant.)
*
*             TODO Implement move constructor/assignment.
*
*             TODO Put up the 3D Version (basically means wrapping all statements
*             which contain a TFESpace2D in ifdefs and adding a TFESpace3D alternative).
*
*             TODO Named constructors for the typically used types:
*             Convection-diffusion (1x1...)
*             5 NSTypes in 2D
*             5 NSTypes in 3D

* @author     Naveed Ahmed, Clemens Bartsch, Ulrich Wilbrandt
* @date       2015/12/08
*
* @ruleof0
*
****************************************************************************/

#ifndef USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_
#define USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_

#include <ColoredBlockMatrix.h>
#include <FEMatrix.h>
#include <FESpace.h>

class ColoredBlockFEMatrix : public ColoredBlockMatrix
{
  public:
    /**
     * @brief Creates a ColoredBlockFEMatrix which is filled with fitting zero blocks.
     *
     * The number of block rows and the number of block columns
     *  is the length of spaces.
     *
     * @param spaces FE spaces, which represent the test spaces rowwise
     * as well as the ansatz spaces columnwise. Each space
     * applies to all blocks of a particular row as testspace and
     * to all blocks of a particualr column as ansatz space.
     */
    ColoredBlockFEMatrix(std::vector< const TFESpace2D*  > spaces);


    /**
     * Add the actives of a certain FEMatrix to several blocks at once.
     * This just figures out whether the adding will work, whether the
     * coloring scheme must be adapted and then delegates to FEMatrix
     * to perform the actual adding of actives.
     */
    void add_scaled_actives(
        const FEMatrix& summand, double factor,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states);

    /** @brief compute y = Ax, paying respect to non-active degrees of freedom.
     *
     * @param[in] x the BlockVector which is multiplied by this matrix
     * @param[out] y result of matrix-vector-multiplication
     */
    virtual void apply(const BlockVector& x, BlockVector& y) const override;

    /** @brief Compute y = y + a * Ax
     *
     * Add the matrix-vector product "Ax", scaled by "a", to y. "A" is this
     * matrix. Non-Active degrees of freedom are taken into account, i.e. the
     * matrix is treated just the way it would be treated as if the non-
     * active rows of the blocks would have been assembled globally correct.
     *
     *
     *
     * @param x the BlockVector which is multiplied by this matrix
     * @param y Gets added the result of the scaled matrix-vector-multiplication "aAx"
     * @param a optional factor, defaults to 1.0
     */
    virtual void apply_scaled_add(const BlockVector & x, BlockVector & y,
                          double a = 1.0) const override;

    /**
     * Used as a developmental tool to discover slicing,
     * should not be used anymore when the class is finished.
     * Checks if all TMatrix smart pointers stored in the
     * base class can be painlessly casted into smart pointers
     * to FEMatrix. Complains if not so, but does not throw.
     */
    void check_pointer_types();

    /*! Check whether a BlockVector b is fit to be the
     * rhs b of the equation Ax=b. (Including actives check.)
     * @param[in] b Rhs b of the equation Ax=b.
     */
    virtual void check_vector_fits_image(const BlockVector& b) const override;

    /*! Check whether a BlockVector x is fit to be the
     * factor x of the equation Ax=b. (Including actives check.)
     * @param[in] x Factor x in the equation Ax=b.
     */
    virtual void check_vector_fits_pre_image(const BlockVector& x) const override;

    // Getters.
    /**
     * Get the ansatz space of a certain grid cell. Since the ansatzspace is
     * identical for all cells in a column, the input cell_row is not actually needed.
     * It is only requested so that the user does not have to think about
     * whether it's rowwise or columnwise constant every time she calls this method...
     *
     * @param[in] cell_row The cell row of the grid cell.
     * @param[in] cell_column The cell column of the grid cell.
     * @return The ansatzspace, which is the same for the entire column.
     */
    const TFESpace2D& get_ansatz_space(size_t cell_row, size_t cell_column) const;

    /** @brief return this ColoredBlockMatrix as one TMatrix
     *
     * TODO Implement!
     *
     * This returns a merged version of this matix. Note that the merged
     * matrix does not get stored internally, for it cannot easily be kept
     * up to date, but recreated on every call.
     *
     * Usually this is used to pass this matrix to a solver.
     */
    virtual std::shared_ptr<TMatrix> get_combined_matrix() const override;

    /**
     * Get the test space of a certain grid cell. Since the testspace is
     * identical for all cells in a row, the input cell_column is not needed.
     * It is only requested so that the user does not have to think about
     * whether it's rowwise or columnwise constant every time she calls this method...
     *
     * @param[in] cell_row The cell row of the grid cell.
     * @param[in] cell_column The cell column of the grid cell.
     * @return The testspace, which is the same for the entire row.
     */
    const TFESpace2D& get_test_space(size_t cell_row, size_t cell_column) const;


    const TFESpace2D& get_column_space(size_t cell_column) const
    {
      return *ansatz_spaces_columnwise_.at(cell_column);
    }

    const TFESpace2D& get_row_space(size_t cell_row) const
    {
      return *test_spaces_rowwise_.at(cell_row);
    }

    /**
     * Overrides the method from the base class
     * and does nothing but print an error, when called. This ensures, that
     * no TMatrix is put into a BlockFEMatrix.
     */
    virtual void replace_blocks(
        const TMatrix& new_block,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states) override;

    /**
     * Replace the blocks in a set of cells with an FEMatrix.
     * The new block's test- and ansatzspace must fit those of
     * the grid positions where it should be placed (resp. transposed
     * mode).
     * Also, if the testspace of one of the cells has non-active
     * degrees of freedom, the matrix may not be stored as transposed
     * there.
     * If either of these conditions is violated, the method will quit
     * the program with an error.
     *
     * @param[in] new_block The new block, a copy will be made.
     * @param[in] cell_positions The positions where the new block is supposed to go.
     * @param[in] transposed_states The transposed mode of the block in the cell positions.
     */
    virtual void replace_blocks(
        const FEMatrix& new_block,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states);

    /*!
     *  Scale active degrees of freedom in the blocks at the cells
     *  whose positions are given by cell_positions. The method just figures
     *  out which blocks are to be treated and whether the coloring
     *  scheme has to be adapted, then delegates the active-scaling to the
     *  FEMatrix class.
     *
     *  @param[in] factor The scaling factor
     *  @param[in] cell_positions Those cells whose blocks are to scale.
     */
    void scale_blocks_actives(
        double factor,
        const std::vector<std::vector<size_t>>& cell_positions );


    // Special member functions.

    /*! Copy constructor.
     *  @note Is a bit expensive because it copies each stored
     *  matrix twice - once in the base class copy constructor
     *  and once in its body to get the cast to FEMatrix right.
     */
    ColoredBlockFEMatrix(const ColoredBlockFEMatrix&);

    ///! Deleted move constructor. (Should be implemented in time)
    ColoredBlockFEMatrix(ColoredBlockFEMatrix&&) = delete;

    /*! Copy assignment operator. Uses copy-and-swap.
     */
    ColoredBlockFEMatrix& operator=(ColoredBlockFEMatrix);

    /** Swap function used for copy-and swap in copy assignment.
     * @param[in,out] first The object to be swapped with second.
     * @param[in,out] second The object to be swapped with first.
     */
    friend void swap(ColoredBlockFEMatrix& first, ColoredBlockFEMatrix& second);

    //! Deleted move assignment operator. (Should be implemented in time)
    ColoredBlockFEMatrix& operator=(ColoredBlockFEMatrix&&) = delete;

    /// @brief Destructor. Tidies up nice and clean.
    virtual ~ColoredBlockFEMatrix() = default;

  protected:
    /// Store pointers to the testspaces rowwise. (TODO could be changed to weak_ptr)
    std::vector<const TFESpace2D* > test_spaces_rowwise_;
    /// Store pointers to the ansatzspaces columnwise.
    std::vector<const TFESpace2D* > ansatz_spaces_columnwise_;

  private:
    /**
     * Try if a given TMatrix can be cast to an FEMatrix (this meaning the object actually is
     * an FEMatrix) and if so, make an FEMatrix copy of it and return a smart pointer to that,
     * so that it can be stored as a shared pointer to TMatrix in the baseclass.
     *
     * This quirky method is needed due to the baseclass storing TMatrices only,
     * but this class dealing with FEMatrices.
     */
    virtual std::shared_ptr<TMatrix> create_block_shared_pointer(const TMatrix& block) override;

    /**
     * Actual implementation of add scaled actives method, whose interface is given
     * in the public part.
     */
    void add_scaled_actives(
        const FEMatrix& summand, double scaling_factor,
        std::vector<grid_place_and_mode> row_column_transpose_tuples);

    /**
     * Actual implementation of the scale actives method, whose interface is given
     * in the public part.
     */
    void scale_blocks_actives( double scaling_factor,
                       std::vector<grid_place_and_mode> row_column_transpose_tuples);


};



#endif /* USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_ */
