// =======================================================================
// @(#)Matrix3D.h        1.2 11/20/98
// 
// Class:       TMatrix3D
//
// Purpose:     store a  matrix (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#ifndef __MATRIX3D__
#define __MATRIX3D__

#include <Structure.h>
#include <Matrix.h>

class TMatrix3D : public TMatrix
{
  public:
    /** generate the matrix */
    TMatrix3D(TStructure *structure);

    /** destructor: free Entries array */
    ~TMatrix3D();

    /** @brief set all Dirichlet rows to zero. That means all rows where the 
     * test space has nonactive degrees of freedom. 
     */
    void resetNonActive();
    
    
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
    void add_active(const TMatrix3D& m, double factor = 1.0);
};

#endif
