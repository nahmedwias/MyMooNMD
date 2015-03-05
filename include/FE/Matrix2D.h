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

#include <Structure2D.h>
#include <Matrix.h>

class TMatrix2D : public TMatrix
{
  protected:
    /** matrix structure */
    TStructure2D *Structure;

  public:
    /** generate the matrix */
    TMatrix2D(TStructure2D *structure);

    /** destructor: free Entries array */
    ~TMatrix2D();

    TStructure2D *GetStructure()
    { return Structure; }

};

void AllocateMatricesNSE_2D(int mg_level,
			    TFESpace2D *velocity_space, 
			    TFESpace2D *pressure_space,
			    TSquareStructure2D *&sqstructureA, 
			    TSquareStructure2D *&sqstructureC, 
			    TStructure2D *&structureB, 
			    TStructure2D *&structureBT,
			    TSquareMatrix2D *&sqmatrixA,
			    TSquareMatrix2D *&sqmatrixA11,
			    TSquareMatrix2D *&sqmatrixA12,
			    TSquareMatrix2D *&sqmatrixA21,
			    TSquareMatrix2D *&sqmatrixA22,
			    TSquareMatrix2D *&sqmatrixC,
			    TMatrix2D *&matrixB1,
			    TMatrix2D *&matrixB2,
			    TMatrix2D *&matrixB1T,
			    TMatrix2D *&matrixB2T,
			    TSquareMatrix2D **MatricesA,
			    TSquareMatrix2D **MatricesA11,
			    TSquareMatrix2D **MatricesA12,
			    TSquareMatrix2D **MatricesA21,
			    TSquareMatrix2D **MatricesA22,
			    TSquareMatrix2D **MatricesC,
			    TMatrix2D **MatricesB1,			    
			    TMatrix2D **MatricesB2,			    
			    TMatrix2D **MatricesB1T,			    
			    TMatrix2D **MatricesB2T);


#endif
