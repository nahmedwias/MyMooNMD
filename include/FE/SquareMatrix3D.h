// =======================================================================
// @(#)SquareMatrix3D.h        1.3 11/20/98
// 
// Class:       TSquareMatrix3D
//
// Purpose:     store a square matrix (ansatz = test space) in 3d
//
// Author:      Gunar Matthies
//
// History:     10.08.1998 start implementation
//
// =======================================================================

#ifndef __SQUAREMATRIX3D__
#define __SQUAREMATRIX3D__

#include <SquareMatrix.h>

class TSquareMatrix3D : public TSquareMatrix
{
  public:
    /** generate the matrix */
    TSquareMatrix3D(TStructure *squarestructure);

//     /** generate an empty nxn matrix */
//     explicit TSquareMatrix3D(int n);
    
    /** destructor: free Entries array */
    ~TSquareMatrix3D();

    /** return FESpace */
    const TFESpace3D *GetFESpace() const
    { return structure->GetFESpace3D(); }

    
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
    void add_active(const TSquareMatrix3D& m, double factor = 1.0);
};

#endif
