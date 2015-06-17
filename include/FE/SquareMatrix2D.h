// =======================================================================
// @(#)SquareMatrix2D.h        1.3 11/20/98
// 
// Class:       TSquareMatrix2D
//
// Purpose:     store a square matrix (ansatz = test space) in 2d
//
// Author:      Gunar Matthies
//
// History:     10.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUAREMATRIX2D__
#define __SQUAREMATRIX2D__

#include <SquareMatrix.h>
#include <SquareStructure2D.h>

class TSquareMatrix2D : public TSquareMatrix
{
  protected:
    /** matrix strcuture */
    TSquareStructure2D *structure;

  public:
    /** generate the matrix */
    TSquareMatrix2D(TSquareStructure2D *squarestructure);

    /** generate an empty nxn matrix */
    explicit TSquareMatrix2D(int n);
    
    /** fill empty matrix, 
     you can either call the constructor TSquareMatrix2D(TSquareStructure2D*);
     or use TSquareMatrix2D(); and then  void SetStructure(TSquareStructure2D*);
    */
    void SetStructure(TSquareStructure2D *squarestructure);
    
    /** destructor: free Entries array */
    ~TSquareMatrix2D();

    /** return FESpace */
    TFESpace2D *GetFESpace() const
    { return structure->GetFESpace(); }

    /** return used matrix structure */
    TSquareStructure2D *GetMatrixStructure() const
    { return structure; }
    TSquareStructure2D *GetStructure() const
    { return structure; }
    
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
    void add_active(const TSquareMatrix2D& m, double factor = 1.0);
    
    
    /** @brief scale matrix by scalar (only active entries) */
    TSquareMatrix2D& operator*=(double alpha);
    /** @brief add another matrix to this one (only active entries) */
    TSquareMatrix2D& operator+=(const TSquareMatrix2D& rhsMat);
    /** @brief add another matrix to this one (only active entries) */
    TSquareMatrix2D& operator+=(const TSquareMatrix2D* rhsMat)
    { *this += *rhsMat; return *this; }
    
    /** @brief copy matrix 'rhs' to this */
    TSquareMatrix2D& operator=(const TSquareMatrix2D& rhs);

    /** @brief add to matrices A and B
     * 
     * note: only active DOF are added
     * note: only works for matrices with the same sparsity pattern
    */
    friend TSquareMatrix2D& operator+(const TSquareMatrix2D & A, 
                                      const TSquareMatrix2D & B);
    
    /**  @brief C= A*alpha 
     * 
     * note: only active DOF are multiplied, others are just copied
     * note: the user is responsible to delete the newly created matrix C
     */
    friend TSquareMatrix2D& operator*(const TSquareMatrix2D & A,
                                      const double alpha);
    /**  @brief C= A*alpha 
     * 
     * Same as TSquareMatrix2D& operator+(const TSquareMatrix2D & A, 
     *                                    const TSquareMatrix2D & B);
     * 
     * note: only active DOF are multiplied, others are just copied
     * note: the user is responsible to delete the newly created matrix C
     */
    friend TSquareMatrix2D& operator*(const double alpha,
                                      const TSquareMatrix2D & A);

    /** @brief y= A*x
     * 
     * note: only active DOF are multiplied, others are just copied from x
     */
    friend double* operator*(const TSquareMatrix2D & A, const double* x);
};

#endif
