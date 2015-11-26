/** ************************************************************************ 
*
* @class     BlockMatrixNSE2D
* @brief     stores the information of a 2D NSE system matrix 
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History   Added methods (Sashi, 23.08.14)
*            made this class a derived class fo BlockMatrix2D (Ulrich, 19.03.2015)
*            further simplifications (Ulrich, 25.03.2015)
* ************************************************************************  */


#ifndef __SYSTEMMATNSE2D__
#define __SYSTEMMATNSE2D__

#include <SquareMatrix2D.h>
#include <Matrix2D.h>
#include <BlockVector.h>
#include <LocalAssembling2D.h>
#include <array>

/**class for 2D  NSE system matrix */
class BlockMatrixNSE2D : public BlockMatrix
{
  protected:
    
    /** @brief Boundary value */
    std::array<BoundValueFunct2D * const, 3> boundary_values;
    
  public:
    /** constructor */
     BlockMatrixNSE2D(const TFESpace2D& velocity, const TFESpace2D& pressure,
                      BoundValueFunct2D * const * const BoundValue);
     
    /** destrcutor */
    ~BlockMatrixNSE2D();
    
    /** assemble the system matrix */
    void Assemble(LocalAssembling2D& la, BlockVector& rhs);
    
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear(LocalAssembling2D& la);
    
    /** solve the system matrix */
    void Solve(double *sol, double *rhs);
    
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
    void apply_scaled_add(const double *x, double *y, double factor = 1.) const;
    
    
    /** @brief return the test or ansatz space for a given block
     * 
     * @param b the index of the block whose test/ansatz space is returned
     * @param test true to return the test space, false to return the ansatz 
     *             space
     */
    const TFESpace2D * get_space_of_block(unsigned int b, bool test) const;
    
    /** @brief return the velocity-velocity block
     * 
     * There are four A-blocks: A11 (i=0), A12 (i=1), A21 (i=2), A22 (i=3).
     * 
     * They are created as a TSquareMatrix2D in the constructor of this class. 
     * Therefore the following cast works. However in case of NSTYPE 1 or 2, 
     * only the diagonal blocks are of type TSquareMatrix2D, the off diagonal
     * ones are not. Therefore if `i` is 0 or 3 this method will always return
     * what is expected. If `i` is 2 or 3 and the NSTYPE is 1 or 2 an exception
     * will be thrown.
     */
    TSquareMatrix2D * get_A_block(unsigned int i = 0);
    
    /** @brief return the velocity-velocity block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    const TSquareMatrix2D * get_A_block(unsigned int i = 0) const;
    
    /** @brief return the pressure-perssure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     * 
     * Note that the C block is only defined for NSTYPE 14, otherwise an 
     * exception is thrown.
     */
    TSquareMatrix2D * get_C_block();
    
    /** @brief return the pressure-pressure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     * 
     * Note that the C block is only defined for NSTYPE 14, otherwise an 
     * exception is thrown.
     */
    const TSquareMatrix2D * get_C_block() const;
    
    /** @brief return the velocity-perssure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     * 
     * Note that in case of NSTYPE 1 or 3 these blocks are not explicitly stored
     * and this method will throw an exception in that case.
     */
    TMatrix2D * get_BT_block(unsigned int i = 0);
    
    /** @brief return the velocity-pressure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     * 
     * Note that in case of NSTYPE 1 or 3 these blocks are not explicitly stored
     * and this method will throw an exception in that case.
     */
    const TMatrix2D * get_BT_block(unsigned int i = 0) const;
    
    /** @brief return the pressure-perssure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    TMatrix2D * get_B_block(unsigned int i = 0);
    
    /** @brief return the pressure-pressure block
     * 
     * See also the description of TSquareMatrix2D * get_A_block().
     */
    const TMatrix2D * get_B_block(unsigned int i = 0) const;
    
    /** @brief return the finite element space for the velocity */
    const TFESpace2D * get_velocity_space() const
    { return this->get_A_block()->GetFESpace2D(); }
    
    /** @brief return the finite element space for the pressure */
    const TFESpace2D * get_pressure_space() const
    { return this->get_B_block()->GetTestSpace2D(); }
};

#endif // __SYSTEMMATNSE2D__
