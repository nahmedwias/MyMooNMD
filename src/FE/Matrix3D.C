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
#include <LinAlg.h>

TMatrix3D::TMatrix3D(std::shared_ptr<TStructure> structure)
 : TMatrix(structure)
{
}

TMatrix3D::~TMatrix3D()
{
}


void TMatrix3D::resetNonActive()
{
  int n_active = this->structure->GetTestSpace()->GetN_DegreesOfFreedom()
                -this->structure->GetTestSpace()->GetN_Dirichlet();
  int * rowPtr = this->structure->GetRowPtr();
  int index_nonactive = rowPtr[n_active];
  int n_nonactive_entries = rowPtr[this->structure->GetN_Rows()]
                           - index_nonactive;
  memset(this->GetEntries() + index_nonactive, 0.0,
         n_nonactive_entries * SizeOfDouble);
}


void TMatrix3D::reset_non_active()
{
  int n_active_rows = this->structure->GetTestSpace()->GetActiveBound();
  int * rowPtr = this->structure->GetRowPtr();
  int index_nonactive = rowPtr[n_active_rows];
  int n_nonactive_entries = rowPtr[structure->GetN_Rows()] - index_nonactive;
  memset(this->GetEntries() + index_nonactive, 0.0,
         n_nonactive_entries * SizeOfDouble);
}

void TMatrix3D::reset_active()
{
  int n_active_rows = this->structure->GetTestSpace()->GetActiveBound();
  int * rowPtr = this->structure->GetRowPtr();
  // number of entries in active rows
  int n_active = rowPtr[n_active_rows];
  memset(this->GetEntries(), 0.0, n_active * SizeOfDouble);
}

void TMatrix3D::scale_active(double factor)
{
  if(factor == 1.0)
    return; // no scaling
  if(factor == 0.0)
    this->reset_active();
  
  // number of active rows
  int n_active_rows = this->structure->GetTestSpace()->GetActiveBound();
  int * rowPtr = this->structure->GetRowPtr();
  // number of entries in active rows
  int n_active = rowPtr[n_active_rows];
  Dscal(n_active, factor, this->GetEntries());
}

void TMatrix3D::add_active(const TMatrix3D& m, double factor)
{
  if(this->GetStructure() != m.GetStructure()) // compare objects
  {
    ErrMsg("TMatrix::add : the two matrices do not match.");
    throw("TMatrix::add : the two matrices do not match.");
  }
  
  // number of active rows
  int n_active_rows = this->structure->GetTestSpace()->GetActiveBound();
  int * rowPtr = this->structure->GetRowPtr();
  // number of entries in active rows
  int n_active = rowPtr[n_active_rows];
  Daxpy(n_active, factor, m.GetEntries(), this->GetEntries());
}
