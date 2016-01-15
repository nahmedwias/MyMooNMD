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
*             The ColoredBlockFEMatrix assumes that with each cell row
*             and cell column exactly one FESpace is associated.
*             Everything else might lead to problems in the current set up.
*             TODO Enforce this as invariant when emplacing blocks?!
*
* @author     Naveed Ahmed, Clemens Bartsch, Ulrich Wilbrandt
* @date       2015/12/08
*
*
*
****************************************************************************/

#ifndef USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_
#define USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_

#include <ColoredBlockMatrix.h>
#include <FEMatrix.h>

class ColoredBlockFEMatrix : public ColoredBlockMatrix
{
  public:
    /**
     * @brief Creates a ColoredBlockFEMatrix which is filled with fitting zero blocks.
     *
     * The number of block rows is the length of cell_row_numbers,
     * the number of block columns the length of cell_column_numbers.
     *
     * @note Since there is no default constructor for FEMatrix
     * an object created with this constructor starts its life with
     * storing only TMatrix 0 blocks - this means, that it is not usable
     * as BlockFEMatrix until all cells are filled with actual FEMatrices!
     *
     * @param cell_row_numbers holds the number of rows for all matrices in one block row
     * @param cell_column_numbers holds the number of columns for all matrices in one block column
     */
    ColoredBlockFEMatrix(std::vector<size_t> cell_row_numbers,
                         std::vector<size_t> cell_column_numbers);


    /*
     * Add the actives of a certain FEMatrix to several blocks at once.
     */
    void add_scaled_actives(
        const FEMatrix& summand, double factor,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states);

    /*
     * Used as a developmental tool to discover slicing,
     * should not be used anymore when the class is finished.
     * Checks if all TMatrix smart pointers stored in the
     * base class can be painlessly casted into smart pointers
     * to FEMatrix. Complains if not so, but does not throw.
     */
    void check_pointer_types();

    //! Return the number of active dofs in a certain cell row
    size_t get_n_row_actives( size_t cell_row ) const;

    //! Return the number of active dofs in a certain cell column
    size_t get_n_column_actives( size_t cell_column ) const;

    /*
     * Overrides the method from the base class
     * and does nothing but print an error, when called. This ensures, that
     * no TMatrix is put into a BlockFEMatrix (except from those that are there
     * when it is created).
     */
    virtual void replace_blocks(
        const TMatrix& new_block,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states) override;

    /*
     * Replace the blocks in a set of cells with an FEMatrix.
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
     *  whose positions are given by cell_positions.
     *
     *  @param[in] factor The scaling factor
     *  @param[in] cell_positions those cells which are to scale.
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

    /*! Copy assignment operator. Uses base class copy-and-swap
     * but derived class copy constructor by default.
     */
    ColoredBlockFEMatrix& operator=(const ColoredBlockFEMatrix&) = default;

    //! Deleted move assignment operator. (Should be implemented in time)
    ColoredBlockFEMatrix& operator=(ColoredBlockFEMatrix&&) = delete;

    /// @brief Destructor. Tidies up nice and clean.
    virtual ~ColoredBlockFEMatrix() = default;

  private:
    virtual std::shared_ptr<TMatrix> create_block_shared_pointer(const TMatrix& block) override;

    void add_scaled_actives(
        const FEMatrix& summand, double scaling_factor,
        std::vector<grid_place_and_mode> row_column_transpose_tuples);

    void scale_blocks_actives( double scaling_factor,
                       std::vector<grid_place_and_mode> row_column_transpose_tuples);


};



#endif /* USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_ */
