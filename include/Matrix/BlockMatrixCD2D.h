/** ************************************************************************ 
*
* @class    BlockMatrixCD2D
* @brief    stores the information of a 2D scalar system matrix 
*           
*           During the contructor all objects are created. However the matrix
*           is only stored in the base class BlockMatrix as TMatrix. It is
*           created as TSquareMatrix2D, so that a cast is possible, see e.g. the
*           method get_matrix() of this class. This way it is possible to access
*           the finite element space as well. 
*           
*           This class does not store pointers to the matrix and the finite
*           element space.
*           
* @author   Ulrich Wilbrand, 
* @date     11.09.2015
 ************************************************************************  */


#ifndef __SYSTEMMATSCALAR2D__
#define __SYSTEMMATSCALAR2D__

#include <BlockVector.h>
#include <SquareMatrix2D.h>
#include <LocalAssembling2D.h>

/**class for 2D scalar system matrix */
class BlockMatrixCD2D : public BlockMatrix
{
  protected:
    
    /** @brief Boundary value */
    BoundValueFunct2D * const boundary_values;
    
  public:
    /** constructor */
     BlockMatrixCD2D(const TFESpace2D &fespace, 
                     BoundValueFunct2D * const BoundValue,
                     bool mass_matrix = false);

    /** destrcutor */
    ~BlockMatrixCD2D();
    
    /** assemble the system matrix */
    void Assemble(LocalAssembling2D& la, BlockVector& sol, BlockVector& rhs);

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
    
    /** @brief return the matrix
     * 
     * This is created as a TSquareMatrix2D in the constructor of this class. 
     * Therefore the following cast works.
     */
    TSquareMatrix2D * get_matrix()
    { return (TSquareMatrix2D*)this->BlockMatrix::blocks[0].get(); }
    
    /** @brief return the matrix
     * 
     * See also the description of TSquareMatrix2D * get_matrix().
     */
    const TSquareMatrix2D * get_matrix() const
    { return (TSquareMatrix2D*)this->BlockMatrix::blocks[0].get(); }
    
    /** @brief return the finite element space */
    const TFESpace2D * get_fe_space() const
    { return this->get_matrix()->GetFESpace(); }
    
    
    /** @brief return the fe_space
     * 
     * This is implemented such that 
     * BlockVector::BlockVector<BlockMatrixCD2D>(const BlockMatrixCD2D& , bool)
     * works.
     */
    const TFESpace2D * get_space_of_block(unsigned int b, bool test) const
    { return this->get_fe_space(); }
};

#endif // __SYSTEMMATSCALAR2D__
