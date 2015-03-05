// =======================================================================
// @(#)Matrix2D.C        1.2 11/20/98
// 
// Class:       TMatrix2D
//
// Purpose:     store a  matrix2D (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#include <Database.h>
#include <Matrix2D.h>
#include <SquareMatrix2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <string.h>

TMatrix2D::TMatrix2D(TStructure2D *structure)
 : TMatrix(structure)
{
  Structure = structure;
  N_Rows = structure->GetN_Rows();
  N_Columns = structure->GetN_Columns();
  N_Entries = structure->GetN_Entries();
  KCol = structure->GetKCol();
  RowPtr = structure->GetRowPtr();

  Entries = new double[N_Entries];
  memset(Entries, 0, N_Entries*SizeOfDouble);

  HangingN_Entries = structure->GetHangingN_Entries();
  HangingRowPtr = structure->GetHangingRowPtr();
  HangingKCol = structure->GetHangingKCol();

}

TMatrix2D::~TMatrix2D()
{
}

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
			    TMatrix2D **MatricesB2T)
{
    // matrix structures
    structureB = new TStructure2D(pressure_space, velocity_space);
    structureBT = new TStructure2D(velocity_space, pressure_space);
    sqstructureA = new TSquareStructure2D(velocity_space);
    sqstructureA->Sort();
    
    // allocate matrices
    switch(TDatabase::ParamDB->NSTYPE)
    {
        case 1:
	    matrixB1 = new TMatrix2D(structureB);
	    matrixB2 = new TMatrix2D(structureB);
	    
	    MatricesB1[mg_level] = matrixB1;
	    MatricesB2[mg_level] = matrixB2;
	    
	    sqmatrixA = new TSquareMatrix2D(sqstructureA);
	    
          MatricesA[mg_level] = sqmatrixA;
          break;

        case 2:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          matrixB1T = new TMatrix2D(structureBT);
          matrixB2T = new TMatrix2D(structureBT);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;

          sqmatrixA = new TSquareMatrix2D(sqstructureA);

          MatricesA[mg_level] = sqmatrixA;
          break;

        case 3:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;

          sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix2D(sqstructureA);

          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
          break;

        case 4:
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          matrixB1T = new TMatrix2D(structureBT);
          matrixB2T = new TMatrix2D(structureBT);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;

          sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix2D(sqstructureA);
	 
          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
	 
          break;
        case 14:
	    sqstructureC = new TSquareStructure2D(pressure_space);
	    sqstructureC->Sort();
          matrixB1 = new TMatrix2D(structureB);
          matrixB2 = new TMatrix2D(structureB);
          matrixB1T = new TMatrix2D(structureBT);
          matrixB2T = new TMatrix2D(structureBT);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;

          sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix2D(sqstructureA);
          sqmatrixC = new TSquareMatrix2D(sqstructureC);

          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
          MatricesC[mg_level] = sqmatrixC;
          break;
      }
    return;
}
