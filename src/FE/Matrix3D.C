// =======================================================================
// @(#)Matrix3D.C        1.2 11/20/98
// 
// Class:       TMatrix3D
//
// Purpose:     store a  matrix3D (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#include <Matrix3D.h>
#include <string.h>

TMatrix3D::TMatrix3D(TStructure3D *structure)
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
}

TMatrix3D::~TMatrix3D()
{
}
