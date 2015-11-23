/** ************************************************************************ 
*
* @class    BlockMatrixCD3D
* @brief    Stores the information of a 3D scalar system matrix
*
*           During the constructor all objects are created. However the matrix
*           is only stored in the base class BlockMatrix as TMatrix. It is
*           created as TSquareMatrix2D, so that a cast is possible, see e.g. the
*           method get_matrix() of this class. This way it is possible to access
*           the finite element space as well.
*
*           This class does not store pointers to the matrix and the finite
*           element space.
*
*			The class is initially a copy of BlockMatrixCD2D.C, with "2" replaced.
*
* @author   Clemens Bartsch
* @date     2015/10/26
 ************************************************************************  */


#ifndef __SYSTEMMATSCALAR3D__
#define __SYSTEMMATSCALAR3D__

#include <BlockVector.h>
#include <SquareMatrix3D.h>
#include <LocalAssembling3D.h>

/**@brief Class for 3D scalar system matrix. */
class BlockMatrixCD3D : public BlockMatrix
{
  protected:

    /** @brief Boundary value */
    BoundValueFunct3D * const boundaryValues_;
    
  public:
    /** @brief Construct a matrix from an FE space and a boundary value function.
     * @param[in] fespace The FESpace from which this matrix is constructed (should be const but can't...).
     * @param[in] BoundValue A function pointer to the used boundary values.
     * @param[in] mass_matrix Whether this matrix will contain a mass matrix part
     */
     BlockMatrixCD3D(const TFESpace3D &feSpace,
                     BoundValueFunct3D * const BoundValue,
                     bool containsMassMatrix = false);
    
    /** assemble the system matrix TODO CB Sashi wants this to be not a class member!*/
    void assemble(const LocalAssembling3D& la, BlockVector& sol, BlockVector& rhs);


    //TODO Both the apply methods are terribly unsafe and probably deprecated.

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
     * This is created as a TSquareMatrix3D in the constructor of this class.
     * Therefore the following cast works.
     */
    TSquareMatrix3D * get_matrix()
    { return (TSquareMatrix3D*)this->BlockMatrix::blocks[0].get(); }

    /** @brief return the matrix
     *
     * See also the description of TSquareMatrix3D * get_matrix().
     */
    const TSquareMatrix3D * get_matrix() const
    { return (TSquareMatrix3D*)this->BlockMatrix::blocks[0].get(); }

    /** @brief return the finite element space */
    const TFESpace3D * get_fe_space() const
    { return this->get_matrix()->GetFESpace3D(); }

    /** @brief return the finite element space, but non-const */
    const TFESpace3D * get_fe_space()
    { return this->get_matrix()->GetFESpace3D(); }


    /** @brief return the fe_space
     *
     * This is implemented such that
     * BlockVector::BlockVector<BlockMatrixCD3D>(const BlockMatrixCD3D& , bool)
     * works.
     */
    const TFESpace3D * get_space_of_block(unsigned int b, bool test) const
    { return this->get_fe_space(); }

};

#endif // __SYSTEMMATSCALAR3D__
