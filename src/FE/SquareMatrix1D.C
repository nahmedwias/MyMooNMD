// =======================================================================
// @(#)SquareMatrix2D.C        1.2 11/20/98
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

#include <SquareMatrix1D.h>
#include <string.h>

TSquareMatrix1D::TSquareMatrix1D(TSquareStructure1D *squarestructure)
  : TSquareMatrix(squarestructure)
{
  Structure = squarestructure;

  N_Rows = squarestructure->GetN_Rows();
  N_Columns = squarestructure->GetN_Columns();
  N_Entries = squarestructure->GetN_Entries();
  KCol = squarestructure->GetKCol();
  RowPtr = squarestructure->GetRowPtr();

  Entries = new double[N_Entries];
  memset(Entries, 0, N_Entries*SizeOfDouble);

  HangingN_Entries = squarestructure->GetHangingN_Entries();
//   HangingRowPtr = squarestructure->GetHangingRowPtr();
//   HangingKCol = squarestructure->GetHangingKCol();

  FESpace = squarestructure->GetFESpace();

//   ActiveBound = FESpace->GetActiveBound();
  ActiveBound = N_Rows;

  ColOrder = squarestructure->GetColOrder();

//   cout << " ColOrder " << ColOrder <<endl;
}

TSquareMatrix1D::~TSquareMatrix1D()
{
//  if(Entries) delete [] Entries;
}
 
