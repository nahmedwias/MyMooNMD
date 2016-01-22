/** ************************************************************************
*
* @class      ColoredBlockFEMatrix
* @brief      extends ColoredBlockMatrix by handling active degrees of freedom
*
*             A BlockMatrix of subtype BlockFEMatrix stores FEMatrices
*             instead of simple algebraic TMatrices. Thus it can access
*             information on finite element spaces and active degrees of
*             freedoms and exploits this algorithmically. The majority of
*             BlockMatrices in ParMooN will in fact be BlockFEMatrices.
*
*             Each cell row is associated with a certain FE test space
*             and each cell column with a certain FE ansatz space.
*
*             The possibility to store matrices as transposed
*             leads to a confusion of Test- and Ansatzspace between the
*             BlockFEMatrix and the stored FEMatrices. If an FEMatrix is
*             stored as transposed somewhere and one requests the testspace
*             of the matrix stored there, it gives a different space than
*             when requiring the testspace of the cell where it is stored.
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
*            TODO Write and test move constructor and assignment. Clang and gcc treat
*            implicit generation and linking differently, and then move constructor
*            and assignment differently, too. Somehow both compilers create a default move ctor,
*            which does not make use of copying.
*            But the move assignment is replaced by copies in gcc and gives a compiler error
*            in clang (declared default, implicitely deleted, not replaced by copies)...
*
*            TODO Put up the 3D Version (basically means wrapping all statements
*            which contain a TFESpace2D in ifdefs and adding a TFESpace3D alternative).
*
* @author     Clemens Bartsch
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
    //constructors
    /**
     * @brief Creates a ColoredBlockFEMatrix which is filled with fitting zero blocks.
     * The TStructure of each of these zero blocks is the empty zero-map TStructure.
     *
     * The number of block rows and the number of block columns
     * is the length of spaces.
     *
     * @param[in] spaces FE spaces, which represent the test spaces rowwise
     * as well as the ansatz spaces columnwise. Each space
     * applies to all blocks of a particular row as testspace and
     * to all blocks of a particular column as ansatz space.
     */
    ColoredBlockFEMatrix(std::vector< const TFESpace2D*  > spaces);

    /**
     * Default constructor. Constructs emtpy matrix which maps zero space to zero space.
     */
    ColoredBlockFEMatrix();

    // named constructors for block fe matrices often used in ParMooN
    //TODO All named constructors should further reduce the number of TMatrix-Copies made
    // in their body!

    /**
     * @brief Named constructor for a block Matrix used in 2D convection-
     * diffusion problems.
     *
     * Constructs a 1x1 block matrix which holds one single block with
     * fitting TStructure.
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param space Ansatz- equals testspace.
     * @return A newly constructed BlockFEMatrix for CD2D problems.
     */
    static ColoredBlockFEMatrix CD2D( const TFESpace2D& space );

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 1.
     * The matrix takes the block structure
     *
     * ( A  0  B1T )
     * ( 0  A  B2T )
     * ( B1 B2 0   )
     *
     * where B1T and B2T are explicitly stored (and marked non-transposed).
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 1.
     */
    static ColoredBlockFEMatrix NSE2D_Type1( const TFESpace2D& velocity, const TFESpace2D& pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 2.
     * The matrix takes the block structure
     *
     * ( A  0  B1T )
     * ( 0  A  B2T )
     * ( B1 B2 0   )
     *
     * where B1T and B2T are explicitly stored (and marked non-transposed).
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 2.
     */
    static ColoredBlockFEMatrix NSE2D_Type2( const TFESpace2D& velocity, const TFESpace2D& pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 3.
     * The matrix takes the block structure
     *
     * ( A11 A12 B1^T )
     * ( A21 A22 B2^T )
     * ( B1  B2  0    )
     *
     * where B1^T and B2^T are are not explicitly stored.
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 3.
     */
    static ColoredBlockFEMatrix NSE2D_Type3( const TFESpace2D& velocity, const TFESpace2D& pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 4.
     * The matrix takes the block structure
     *
     * ( A11 A12 B1T )
     * ( A21 A22 B2T )
     * ( B1  B2  0   )
     *
     * where B1^T and B2^T are explicitly stored (and marked non-transposed)..
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 4.
     */
    static ColoredBlockFEMatrix NSE2D_Type4( const TFESpace2D& velocity, const TFESpace2D& pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 14.
     * The matrix takes the block structure
     *
     * ( A11 A12 B1^T )
     * ( A21 A22 B2^T )
     * ( B1  B2  C    ),
     *
     * where B1^T and B2^T are explicitly stored (and marked non-transposed).
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 14.
     */
    static ColoredBlockFEMatrix NSE2D_Type14( const TFESpace2D& velocity, const TFESpace2D& pressure);

    /**
     * Named constructor for a matrix for Darcy type problems in 2D.
     * Creates a 2x2 block matrix of the form
     *  ( A  B1' )
     *  ( B2 C   ),
     * where A is velo-velo coupling, B1' velo-pressure,
     * B2 pressure-velo and C pressure-pressure coupling.
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for Darcy problems in 2D.
     */
    static ColoredBlockFEMatrix Darcy2D( const TFESpace2D& velocity, const TFESpace2D& pressure);

    //public methods

    /**
     * Add the actives of a certain FEMatrix to several blocks at once.
     * This just figures out whether the adding will work, whether the
     * coloring scheme must be adapted and then delegates to FEMatrix
     * to perform the actual adding of actives.
     *
     * @param[in] summand The FEMatrix to be added.
     * @param[in] factor The factor by which to scale it.
     * @param[in] cell_positions Where to add the matrix.
     * @param[in] transposed_states In which transposed state to do it.
     */
    void add_matrix_actives(
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
     * @param x the BlockVector which is multiplied by this matrix
     * @param y Gets added the result of the scaled matrix-vector-multiplication "aAx"
     * @param a optional factor, defaults to 1.0
     */
    virtual void apply_scaled_add(const BlockVector & x, BlockVector & y,
                          double a = 1.0) const override;

    /**
     * Used as a developmental tool to discover slicing,
     * there should be no reason to use it anymore when the class is finished.
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

    /**
     * This method is the main interface to ParMooN solving procedures which
     * depend on the block structure.
     * It traverses the whole matrix from left to right, top to bottom, and fills
     * the return vector with pointers to all blocks, regardless of the storage
     * structure.
     *
     * The classes which make use of this method must decide how to proceed
     * with the gained blocks. They must make use of a priori knowledge
     * about the block structure of this matrix, or of methods which check
     * the block structure.
     *
     * @note The method returns shared_pointers which is not optimal,
     * for the solvers should not share ownership. But since most solvers
     * in ParMooN still expect raw pointer and one can't get those from weak
     * pointers, we chose sharedp pointers as return values.
     *
     * @return A list of pointers to the blocks. Left to right, top to bottom,
     * all blocks appear as often as they appear in cells.
     */
    std::vector<std::shared_ptr<const FEMatrix>> get_blocks() const;

    /**
     * This is an awful method which returns non-const shared pointers to all blocks,
     * ordered from left to right, top to bottom.
     * Of course calling this method means you can break your entire
     * ColoredBlockFEMatrix, because whenever you changes something in one
     * block, you will have no idea which other cells are affected due to
     * storing the same block.
     *
     * At the moment we make use of it, because many of the solvers request
     * non-const pointers to TMatrices. This must be changed by and by.
     *
     * FIXME Fix the solvers (esp. Multigrid!), then delete this method immediately.
     */
    std::vector<std::shared_ptr<FEMatrix>> get_blocks_TERRIBLY_UNSAFE();

    /**
     * This method is the main interface to ParMooN assembling procedures.
     * It traverses the whole matrix from left to right, top to bottom, and fills
     * the return vector with each new block it finds.
     *
     * This results in an output vector which contains each stored matrix
     * exactly once. That is the only thing the method cares about: that
     * no stored block is handed out twice.
     *
     * The classes which make use of this method must decide how to proceed
     * with the gained blocks. They must make use of a priori knowledge
     * about the block structure of this matrix, or of methods which check
     * the block structure.
     *
     * By default the method skips all those blocks which have a zero-map
     * Structure. To include them, give "true" as an input parameter.
     *
     * @param[in] include_zeroes Whether to include blocks with zero-map
     * Structure or skip them. Defaults to "false" (skipping zero-map blocks).
     * @return A list of pointers to the blocks. Each (relevant) blocks appears once.
     */
    std::vector<std::shared_ptr<FEMatrix>> get_blocks_uniquely(bool include_zeroes = false);

    /**
     * Acts just as get_blocks_uniquely(bool), but on only those
     * matrix cells whose coordinates are given in input vector "cells".
     *
     * @param[in] cells The positions of the cells whose blocks we are interested in.
     * should be a vector of double pairs, obviously, and not out of range.
     * @param[in] include_zeroes Whether to include blocks with zero-map
     * Structure or skip them. Defaults to "false" (skipping zero-map blocks).
     * @return A list of pointers to the matrix blocks. Each (relevant) blocks appears once.
     */
    std::vector<std::shared_ptr<FEMatrix>> get_blocks_uniquely(std::vector<std::vector<size_t>> cells,
                                               bool include_zeroes = false);

    /// @return The column (means: ansatz-)space of a certain cell column.
    const TFESpace2D& get_column_space(size_t cell_column) const
    {
      return *ansatz_spaces_columnwise_.at(cell_column);
    }

    /** @brief return this ColoredBlockMatrix as one TMatrix
     *
     * This returns a merged version of this matix. Note that the merged
     * matrix does not get stored internally, for it cannot easily be kept
     * up to date, but recreated on every call.
     *
     * Treats Dirichlet rows correctly and globally, regardless of
     * what the particular blocks hold in their Dirichlet rows.
     *
     * Usually this is used to pass this matrix to a solver.
     */
    virtual std::shared_ptr<TMatrix> get_combined_matrix() const override;

    /**
     * This method returns the number of actives of a certain cell column's
     * ansatz-space. It is needed only for the templated constructor of
     * BlockVector, and wil be removed as soon as ParMooN has a new
     * handling of non-actives.
     *
     * @param cell_column The cell column.
     * @return The number of actives of the column's ansatzspace.
     */
    size_t get_n_column_actives(size_t cell_column) const;

    /**
     * This method returns the number of actives of a certain cell row's
     * test-space. It is needed only for the templated constructor of
     * BlockVector, and wil be removed as soon as ParMooN has a new
     * handling of non-actives.
     *
     * @param cell_row The cell row.
     * @return The number of actives of the row's testspace.
     */
    size_t get_n_row_actives(size_t cell_row) const;


    /// @return The row (means: test-)space of a certain cell row.
    const TFESpace2D& get_row_space(size_t cell_row) const
    {
      return *test_spaces_rowwise_.at(cell_row);
    }

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

    ///! Leave moving to the compiler. TODO Investigate, whether this only performs correct
    /// because it is replaced by copies - if so, implement it by hand!
    ColoredBlockFEMatrix(ColoredBlockFEMatrix&&) = default;

    /// Copy assignment operator. Uses copy-and-swap.
    ColoredBlockFEMatrix& operator=(ColoredBlockFEMatrix);

    /** Swap function used for copy-and swap in copy assignment.
     * @param[in,out] first The object to be swapped with second.
     * @param[in,out] second The object to be swapped with first.
     */
    friend void swap(ColoredBlockFEMatrix& first, ColoredBlockFEMatrix& second);

    /// Leave moving to the compiler. TODO Investigate, whether this only performs correct
    /// because it is replaced by copies - if so, implement it by hand!
    /// With clang, this will not work because it is implicitely deleted,
    /// clang will not replace it by copying, while gcc will!
    ColoredBlockFEMatrix& operator=(ColoredBlockFEMatrix&&) = default;

    /// @brief Destructor. Tidies up nice and clean.
    virtual ~ColoredBlockFEMatrix() = default;

  protected:
    /// Store pointers to the testspaces rowwise. (TODO could be changed to weak_ptr)
    std::vector<const TFESpace2D* > test_spaces_rowwise_;
    /// Store pointers to the ansatzspaces columnwise. (TODO could be changed to weak_ptr)
    std::vector<const TFESpace2D* > ansatz_spaces_columnwise_;

  private:
    /**
     * Actual implementation of add scaled actives method, whose interface is given
     * in the public part.
     */
    void add_scaled_actives(
        const FEMatrix& summand, double scaling_factor,
        std::vector<grid_place_and_mode> row_column_transpose_tuples);

    /**
     * Try if a given TMatrix can be cast to an FEMatrix (this meaning the object actually is
     * an FEMatrix) and if so, make an FEMatrix copy of it and return a smart pointer to that,
     * so that it can be stored as a shared pointer to TMatrix in the baseclass.
     *
     * This quirky method is needed due to the base class storing TMatrices only,
     * but this class dealing with FEMatrices.
     */
    virtual std::shared_ptr<TMatrix> create_block_shared_pointer(const TMatrix& block) override;

    /**
     * Actual implementation of the scale actives method, whose interface is given
     * in the public part.
     */
    void scale_blocks_actives( double scaling_factor,
                       std::vector<grid_place_and_mode> row_column_transpose_tuples);


};



#endif /* USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_ */
