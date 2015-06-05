/** ************************************************************************ 
*
* @class     TSystemMatScalar2D
* @brief     stores the information of a 2D scalar system matrix 
* @author    Sashikumaar Ganesan, 
* @date      08.08.14
* @History   Added methods (Sashi, 22.08.14)
 ************************************************************************  */


#ifndef __SYSTEMMATSCALAR2D__
#define __SYSTEMMATSCALAR2D__

#include <SystemMat2D.h>
#include <LocalAssembling2D.h>

/**class for 2D scalar system matrix */
class TSystemMatScalar2D : public SystemMat2D
{
  protected:
    
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[1];

     /** Boundary value */   
    BoundValueFunct2D *BoundaryValues[1];
    
  public:
    /** constructor */
     TSystemMatScalar2D(TFESpace2D *fespace);

    /** destrcutor */
    ~TSystemMatScalar2D();
    
    /** Initilize the discrete forms and the matrices */
    void Init(BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue);
 
 
    /** assemble the system matrix */
    void Assemble(LocalAssembling2D& la, double *sol, double *rhs);

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
    void apply_scaled_add(const double *x, double *y, double factor =
        1.0) const;
};

#endif // __SYSTEMMATSCALAR2D__
