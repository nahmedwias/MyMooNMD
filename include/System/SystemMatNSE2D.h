/** ************************************************************************ 
*
* @class     TSystemMatNSE2D
* @brief     stores the information of a 2D NSE system matrix 
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History   Added methods (Sashi, 23.08.14)
*            made this class a derived class fo SystemMat2D (Ulrich, 19.03.2015)
*            further simplifications (Ulrich, 25.03.2015)
* ************************************************************************  */


#ifndef __SYSTEMMATNSE2D__
#define __SYSTEMMATNSE2D__

#include <SystemMat2D.h>
#include <LocalAssembling2D.h>

/**class for 2D  NSE system matrix */
class TSystemMatNSE2D : public SystemMat2D
{
  protected:
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[2];
    
    /** Boundary values */
    BoundValueFunct2D *BoundaryValues[2];
    
  public:
    /** constructor */
     TSystemMatNSE2D(TFEVectFunct2D *Velocity, TFEFunction2D *p);
     
    /** destrcutor */
    ~TSystemMatNSE2D();
    
    /** methods */
    /** Initilize the discrete forms and the matrices */
    void Init(BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue);
    
    /** assemble the system matrix */
    void Assemble(LocalAssembling2D& la, double *sol, double *rhs);
    
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear(LocalAssembling2D& la, double *sol, double *rhs);
    
    /** get the resudual of the NSE system */
    void GetResidual(double *sol, double *rhs, double *res);
    
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
    void apply_scaled_add(const double *x, double *y,
                                  double factor = 1.0) const;
    
    unsigned int n_rows() const; // number of block rows
    unsigned int n_cols() const; // number of block columns
    unsigned int n_total_rows() const; // total number of rows (> nRows)
    unsigned int n_total_cols() const; // total number of columns (> nCols)
};

#endif // __SYSTEMMATNSE2D__
