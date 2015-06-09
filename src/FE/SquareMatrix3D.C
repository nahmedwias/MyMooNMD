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
#include <LinAlg.h>

TSquareMatrix3D::TSquareMatrix3D(TSquareStructure3D *squarestructure)
  : TSquareMatrix(squarestructure), structure(squarestructure)
{
}

TSquareMatrix3D::TSquareMatrix3D(int n) 
 : structure(new TSquareStructure3D(n)), TSquareMatrix(structure)
{
}

TSquareMatrix3D::~TSquareMatrix3D()
{
}

void TSquareMatrix3D::reset_non_active()
{
  int n_active_rows = this->structure->GetFESpace()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  int index_nonactive = rowPtr[n_active_rows];
  int n_nonactive_entries = rowPtr[structure->GetN_Rows()] - index_nonactive;
  memset(Entries + index_nonactive, 0.0, n_nonactive_entries * SizeOfDouble);
}

void TSquareMatrix3D::reset_active()
{
  int n_active_rows = this->structure->GetFESpace()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  // numer of entries in active rows
  int n_active = rowPtr[n_active_rows];
  memset(this->Entries, 0.0, n_active * SizeOfDouble);
}

void TSquareMatrix3D::scale_active(double factor)
{
  if(factor == 1.0)
    return; // no scaling
  if(factor == 0.0)
    this->reset_active();
  
  // number of active rows
  int n_active_rows = this->structure->GetFESpace()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  // numer of entries in active rows
  int n_active = rowPtr[n_active_rows];
  Dscal(n_active, factor, this->Entries);
}

void TSquareMatrix3D::add_active(const TSquareMatrix3D& m, double factor)
{
  if(this->structure != m.GetStructure() // compare pointers
     && (*(this->structure)) != (*(m.GetStructure()))) // compare objects
  {
    ErrMsg("TMatrix::add : the two matrices do not match.");
    throw("TMatrix::add : the two matrices do not match.");
  }
  
  // number of active rows
  int n_active_rows = this->structure->GetFESpace()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  // numer of entries in active rows
  int n_active = rowPtr[n_active_rows];
  Daxpy(n_active, factor, m.GetEntries(), this->Entries);
}
