 // =======================================================================
// @(#)SquareMatrix2D.h        
// 
// Class:       TSquareMatrix2D
//
// Purpose:     store a square matrix (ansatz = test space) in 1d
//
// Author:      Sashikumaar Ganesan
//
// History:     17.05.2007 start implementation
//
// =======================================================================

#ifndef __SQUAREMATRIX1D__
#define __SQUAREMATRIX1D__

#include <SquareMatrix.h>

class TSquareMatrix1D : public TSquareMatrix
{
  protected:
    /** matrix strcuture */
    TStructure *structure;

  public:
    /** generate the matrix */
    TSquareMatrix1D(TStructure *squarestructure);

    /** destructor: free Entries array */
    ~TSquareMatrix1D();

    /** return FESpace */
    const TFESpace1D *GetFESpace() const
    { return structure->GetFESpace1D(); }

    /** return used matrix structure */
    TStructure *GetMatrixStructure() const
    { return structure; }

};

#endif

