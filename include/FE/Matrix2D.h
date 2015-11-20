// =======================================================================
// @(#)Matrix2D.h        1.2 11/20/98
// 
// Class:       TMatrix2D
//
// Purpose:     store a  matrix (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#ifndef __MATRIX2D__
#define __MATRIX2D__

#include <Structure.h>
#include <Matrix.h>

class TMatrix2D : public TMatrix
{
  public:
    /** generate the matrix */
    TMatrix2D(const TFESpace2D * testspace, const TFESpace2D * ansatzspace);
    
    TMatrix2D(const TMatrix2D & m) = default;
    
    /** destructor: free Entries array */
    ~TMatrix2D() = default;

    /** @brief scale all active rows */
    TMatrix2D& operator*=(double alpha);

    /** @brief set all Dirichlet rows to zero.
     * 
     *  That means all rows where the test space has nonactive degrees of 
     *  freedom are set to zero.
     */
    void reset_non_active();
    
    /** @brief set all non-Dirichlet rows to zero.
      * 
      *  That means all rows where the test space has active degrees of 
      *  freedom are set to zero.
      */
    void reset_active();
    
    /** @brief scale this matrix by a factor
     * 
     * Only rows corresponding to active d.o.f are scaled. Other rows remain
     * unscaled.
     */
    void scale_active(double factor = 1.0);
    
    /** @brief adding a scaled matrix to this matrix
     * 
     * This is only done for those rows which correspond to active degrees of 
     * freedom.
     * 
     * The summation is index-wise, i.e. A(i,j) += factor*m(i.j), where A is 
     * this matrix. 
     * 
     * Note that this only works if the sparsity structure is the same for this
     * matrix and m.
     */
    void add_active(const TMatrix2D& m, double factor = 1.0);

    /** @brief add two matrices A and B
     * 
     * note: only active DOF are added
     * note: only works for matrices with the same sparsity pattern 
     */
    friend TMatrix2D& operator+(const TMatrix2D & A, const TMatrix2D & B);
    /** @brief substract matrices A and B, i.e. C = A - B
     * 
     * note: only active DOF are substracted
     * note: only works for matrices with the same sparsity pattern
     */
    friend TMatrix2D& operator-(const TMatrix2D & A, const TMatrix2D & B);
    
    /** @brief C = A*alpha 
     * 
     * C will consist of the active entries of A scaled by alpha. C will have
     * zero entries in nonactive rows.
     */
    friend TMatrix2D& operator*(const TMatrix2D & A, const double alpha);
    /** @brief C = alpha*A 
     * 
     * Same as TMatrix2D& operator*(const TMatrix2D & A, const double alpha);
     */
    friend TMatrix2D& operator*(const double alpha, const TMatrix2D & A);

    /** @brief y = A*x
     * 
     * note: only active DOF are multiplied, others are just copied from x
     * note: the user has to delete y
     */
    friend double* operator*(const TMatrix2D & A, const double* x);
};


#endif // __MATRIX2D__
