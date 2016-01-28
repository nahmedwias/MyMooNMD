/** ************************************************************************ 
*
* @class    BlockMatrixDarcy2D
* @brief    stores the information of a 2D Darcy system matrix 
*           
*           During the contructor all objects are created. However the matrices
*           are only stored in the base class BlockMatrix as TMatrix. They are
*           created as TSquareMatrix2D and TMatrix2D, so that a cast is 
*           possible, see e.g. the method get_A_block() of this class. This way
*           it is possible to access the finite element spaces as well. 
*           
*           This class does not store pointers to the matrices and the finite
*           element spaces.
* 
* @author   Ulrich Wilbrandt
* @date     15.03.15
 ************************************************************************  */


#ifndef __SYSTEMMATDARCY2D__
#define __SYSTEMMATDARCY2D__

#include <LocalAssembling2D.h>
#include <SquareMatrix2D.h>
#include <Matrix2D.h>
#include <BlockVector.h>
#include <BlockFEMatrix.h>
#include <array>

/**class for 2D scalar system matrix */
class BlockMatrixDarcy2D : public BlockFEMatrix
{
  protected:
    
    /** @brief Boundary value */ 
    std::array<BoundValueFunct2D * const, 2> boundary_values;
    
  public:
    /** constructor */
     BlockMatrixDarcy2D(const TFESpace2D& velocity, const TFESpace2D& pressure,
                        BoundValueFunct2D * const * const BoundValue);
    
    /** destrcutor */
    ~BlockMatrixDarcy2D();
    
    /** @brief copy constructor */
    BlockMatrixDarcy2D(BlockMatrixDarcy2D&) = delete;
    
    /** @brief move constructor */
    BlockMatrixDarcy2D(BlockMatrixDarcy2D&&) = delete;
    
    /** assemble the system matrix */
    void Assemble(LocalAssembling2D& la, BlockVector& rhs);
    
    /** @brief compute y = factor* A*x 
     *
     * write the matrix-vector product "Ax" scaled by a factor to y. 
     * "A" is this matrix. Both "A" and "x" remain unchanged. If the factor is
     * 0.0, the vector y is only reset without performing the actual 
     * multiplication.
     *
     * @param x the vector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication and scaling
     * @param factor optional scaling factor, default to 1.0
     */
    void apply(const double *x, double *y, double factor = 1.0) const;

    /** @brief compute y = y + a * Ax 
     *
     * add the matrix-vector product "Ax", scaled by "a", to y.
     * "A" is this matrix.
     * 
     * This function can be used to compute the residual r = b - Ax.
     *
     * @param x the vector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication and scaling
     * @param factor optional scaling   factor, default to 1.0
     */
    void apply_scaled_add(const double *x, double *y, double factor = 1.0) 
      const;

    /** @brief return the velocity-velocity block
     * 
     * This is created as a TSquareMatrix2D in the constructor of this class. 
     * Therefore the following cast works.
     */
    TSquareMatrix2D * get_A_block()
    { return (TSquareMatrix2D*)this->BlockMatrix::blocks[0].get(); }
    
    /** @brief return the velocity-velocity block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    const TSquareMatrix2D * get_A_block() const
    { return (TSquareMatrix2D*)this->BlockMatrix::blocks[0].get(); }
    
    /** @brief return the pressure-perssure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    TSquareMatrix2D * get_C_block()
    { return (TSquareMatrix2D*)this->BlockMatrix::blocks[3].get(); }
    
    /** @brief return the pressure-pressure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    const TSquareMatrix2D * get_C_block() const
    { return (TSquareMatrix2D*)this->BlockMatrix::blocks[3].get(); }
    
    /** @brief return the velocity-perssure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    TMatrix2D * get_BT_block()
    { return (TMatrix2D*)this->BlockMatrix::blocks[1].get(); }
    
    /** @brief return the velocity-pressure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    const TMatrix2D * get_BT_block() const
    { return (TMatrix2D*)this->BlockMatrix::blocks[1].get(); }
    
    /** @brief return the pressure-perssure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    TMatrix2D * get_B_block()
    { return (TMatrix2D*)this->BlockMatrix::blocks[2].get(); }
    
    /** @brief return the pressure-pressure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    const TMatrix2D * get_B_block() const
    { return (TMatrix2D*)this->BlockMatrix::blocks[2].get(); }
    
    /** @brief return the finite element space for the velocity */
    const TFESpace2D * get_velocity_space() const
    { return this->get_A_block()->GetFESpace2D(); }
    
    /** @brief return the finite element space for the pressure */
    const TFESpace2D * get_pressure_space() const
    { return this->get_C_block()->GetFESpace2D(); }
};

#endif // __SYSTEMMATDARCY2D__
