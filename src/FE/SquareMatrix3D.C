// =======================================================================
// @(#)SquareMatrix3D.C        1.2 11/20/98
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

#include <SquareMatrix3D.h>
#include <string.h>

TSquareMatrix3D::TSquareMatrix3D(TSquareStructure3D *squarestructure)
  : TSquareMatrix(squarestructure)
{
  Structure = squarestructure;

  N_Rows = squarestructure->GetN_Rows();
  N_Columns = squarestructure->GetN_Columns();
  N_Entries = squarestructure->GetN_Entries();
  KCol = squarestructure->GetKCol();
  RowPtr = squarestructure->GetRowPtr();

  Entries = new double[N_Entries];
  // memset(Entries, 0, N_Entries*SizeOfDouble);
  memset(Entries, 0, ((long)N_Entries)*SizeOfDouble);

  HangingN_Entries = squarestructure->GetHangingN_Entries();
  HangingRowPtr = squarestructure->GetHangingRowPtr();
  HangingKCol = squarestructure->GetHangingKCol();

  FESpace = squarestructure->GetFESpace();
  if (FESpace!=NULL)
  {
    ActiveBound = FESpace->GetActiveBound();
    HangingBound = FESpace->GetHangingBound();
  }
  else
  {
    ActiveBound = N_Rows;
    HangingBound = N_Rows;
  }

  ColOrder = squarestructure->GetColOrder();
}

TSquareMatrix3D::~TSquareMatrix3D()
{
  delete [] Entries;
  Entries = NULL;
}
