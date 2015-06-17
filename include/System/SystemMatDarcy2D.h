/** ************************************************************************ 
*
* @class     SystemMatDarcy2D
* @brief     stores the information of a 2D Darcy system matrix 
* @author    Ulrich Wilbrandt
* @date      15.03.15
 ************************************************************************  */


#ifndef __SYSTEMMATDARCY2D__
#define __SYSTEMMATDARCY2D__

#include <LocalAssembling2D.h>
#include <SystemMat2D.h>

/**class for 2D scalar system matrix */
class SystemMatDarcy2D : public SystemMat2D
{
  protected:
    
    /** Boundary conditon (one for u.n and one for pressure) */
    BoundCondFunct2D *BoundaryConditions[2];
    
    /** Boundary value */ 
    BoundValueFunct2D *BoundaryValues[2];
    
  public:
    /** constructor */
     SystemMatDarcy2D(TFESpace2D *velocity, TFESpace2D* pressure,
                       BoundValueFunct2D **BoundValue);
    
    /** destrcutor */
    ~SystemMatDarcy2D();
    
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
    void apply_scaled_add(const double *x, double *y,
                                  double factor = 1.0) const;
    
    unsigned int n_rows() const; // number of block rows
    unsigned int n_cols() const; // number of block columns
    unsigned int n_total_rows() const; // total number of rows (> nRows)
    unsigned int n_total_cols() const; // total number of columns (> nCols)
    unsigned int n_total_entries() const; // total number of entries
};

#endif // __SYSTEMMATDARCY2D__
