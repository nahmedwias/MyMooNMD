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
#include <SquareStructure3D.h>

class TSquareMatrix3D : public TSquareMatrix
{
  protected:
    /** FESpace */
    TFESpace3D *FESpace;

    /** matrix strcuture */
    TSquareStructure3D *Structure;

  public:
    /** generate the matrix */
    TSquareMatrix3D(TSquareStructure3D *squarestructure);

    /** destructor: free Entries array */
    ~TSquareMatrix3D();

    /** return FESpace */
    TFESpace3D *GetFESpace()
    { return FESpace; }

    /** return used matrix structure */
    TSquareStructure3D *GetMatrixStructure()
    { return Structure; }

};

#endif
