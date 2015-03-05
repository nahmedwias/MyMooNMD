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
    /** FESpace */
    TFESpace2D *FESpace;

    /** matrix strcuture */
    TSquareStructure2D *Structure;

  public:
    /** generate the matrix */
    TSquareMatrix2D(TSquareStructure2D *squarestructure);

    /** destructor: free Entries array */
    ~TSquareMatrix2D();

    /** return FESpace */
    TFESpace2D *GetFESpace()
    { return FESpace; }

    /** return used matrix structure */
    TSquareStructure2D *GetMatrixStructure()
    { return Structure; }

};

#endif
