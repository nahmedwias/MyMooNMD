/** ************************************************************************
*
* @class      BlockFEMatrix
* @brief      extends BlockMatrix by handling active degrees of freedom
*
*             A BlockMatrix of subtype BlockFEMatrix stores FEMatrices
*             instead of simple algebraic TMatrices. Thus it can access
*             information on finite element spaces and active degrees of
*             freedoms and exploits this algoritmically. The majority of
*             BlockMatrices in ParMooN are in fact BlockFEMatrices.
*
*             \todo This class is still in its infancy.
*
*
* @author     Naveed Ahmed, Clemens Bartsch, Ulrich Wilbrandt
* @date       2015/12/08
*
*
*
****************************************************************************/

#ifndef INCLUDE_MATRIX_BLOCKFEMATRIX_H_
#define INCLUDE_MATRIX_BLOCKFEMATRIX_H_

#include <BlockMatrix.h>

class BlockFEMatrix : public BlockMatrix
{
    protected:
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

    private:
    const TFESpace2D *TestSpace;
    
    const TFESpace2D *AnsatzSpace;
    
    const TFESpace2D *PressureSpace;
    
    public:
    /** @brief construct a BlockFEMatrix suitable for the given problem type
     *
     * All blocks contain zeros only.
     *
     * This constructor essentially creates a BlockPattern using the given
     * arguments.
     */
    BlockFEMatrix(const Problem_type, unsigned int space_dimension,
                bool mass_matrix = false);

    /** @brief a constructor
     * construct a matrix used for the pressure robust methods
     */
    BlockFEMatrix(const TFESpace2D* testSpace, const TFESpace2D* ansatzSpace, 
                  const TFESpace2D* pressureSpace);
    
    /**
     * @brief adding a scaled matrix to this matrix, but only the active entries
     *
     * This does exactly the same as add_scaled, except that nonactive entries
     * are not changed.
     */
    void add_scaled_active(const BlockMatrix &A, double factor = 1.0);

    /**
     * @brief scale all active entries of this matrix
     *
     * That means for each submatrix all active entries are scaled.
     *
     * Possibly existing special matrices are not changed.
     */
    void scale_active(double factor);

    /**
     * @brief Return the number of active dofs in a specific block row.
     * @param[in] index The block row index.
     */
    size_t get_n_row_actives( size_t index ) const;

    /**
     * @brief Return the number of active dofs in a specific block column.
     * @param[in] index The block column index.
     */
    size_t get_n_column_actives( size_t index ) const;
};



#endif /* INCLUDE_MATRIX_BLOCKFEMATRIX_H_ */
