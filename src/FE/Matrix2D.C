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
#include <string.h>
#include <LinAlg.h>

TMatrix2D::TMatrix2D(TStructure *structure)
 : TMatrix(structure)
{
}

TMatrix2D::~TMatrix2D()
{
}


void TMatrix2D::reset_non_active()
{
  int n_active_rows = this->structure->GetTestSpace2D()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  int index_nonactive = rowPtr[n_active_rows];
  int n_nonactive_entries = rowPtr[structure->GetN_Rows()] - index_nonactive;
  memset(Entries + index_nonactive, 0.0, n_nonactive_entries * SizeOfDouble);
}

void TMatrix2D::reset_active()
{
  int n_active_rows = this->structure->GetTestSpace2D()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  // number of entries in active rows
  int n_active = rowPtr[n_active_rows];
  memset(this->Entries, 0.0, n_active * SizeOfDouble);
}

void TMatrix2D::scale_active(double factor)
{
  if(factor == 1.0)
    return; // no scaling
  if(factor == 0.0)
    this->reset_active();
  
  // number of active rows
  int n_active_rows = this->structure->GetTestSpace2D()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  // number of entries in active rows
  int n_active = rowPtr[n_active_rows];
  Dscal(n_active, factor, this->Entries);
}

void TMatrix2D::add_active(const TMatrix2D& m, double factor)
{
  if(this->structure != m.GetStructure() // compare pointers
     && (*(this->structure)) != (*(m.GetStructure()))) // compare objects
  {
    ErrMsg("TMatrix::add : the two matrices do not match.");
    throw("TMatrix::add : the two matrices do not match.");
  }
  
  // number of active rows
  int n_active_rows = this->structure->GetTestSpace2D()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  // number of entries in active rows
  int n_active = rowPtr[n_active_rows];
  Daxpy(n_active, factor, m.GetEntries(), this->Entries);
}

TMatrix2D& TMatrix2D::operator*=(double alpha)
{
  this->scale_active(alpha);
  return *this;
}


// add to matrices A and B
// note: only active DOF are added
// note: only works for matrices with the same sparsity pattern
TMatrix2D& operator+(const TMatrix2D & A, const TMatrix2D & B)
{
  double *AEntries, *BEntries, *CEntries;
  if (A.GetStructure() == B.GetStructure()) 
  {
    TMatrix2D *C = new TMatrix2D(A.GetStructure());
    AEntries = A.GetEntries();
    BEntries = B.GetEntries();
    CEntries = C->GetEntries();
    for (int i=0; i<A.GetN_Entries(); i++) 
    {
      CEntries[i] = AEntries[i] + BEntries[i];
    }
    return *C;
  } else 
  {
    cout << " Matrix2D: ERROR: can't add Matrices "
         << " with different sparse structures " << endl;
    exit(1);
  }
}

TMatrix2D& operator*(const TMatrix2D & A, const double alpha)
{
  double *AEntries, *CEntries;
  TMatrix2D *C = new TMatrix2D(A.GetStructure());
  AEntries = A.GetEntries();
  CEntries = C->GetEntries();

  //TFESpace2D *fespace = A.GetStructure()->GetAnsatzSpace2D();
  const TFESpace2D *fespace = A.GetStructure()->GetTestSpace2D();
  int nDOFActive = fespace->GetN_ActiveDegrees();
  //int nDOF = fespace->GetN_DegreesOfFreedom();

  // multiply each active entry by alpha and write it into matrix C
  for (int i=0; i<nDOFActive; i++) 
    {
      CEntries[i] = alpha*AEntries[i];
    }
  return *C;
}

TMatrix2D& operator*(const double alpha, const TMatrix2D & A)
{// just to allow alpha*A as well (in addition to A*alpha)
  return A*alpha;
}

double* operator*(const TMatrix2D & A,const double* x)
{
  double *AEntries;
  int *ARowPtr,*AColIndex;
  AEntries = A.GetEntries();
  ARowPtr = A.GetRowPtr();
  AColIndex = A.GetKCol();

  const TFESpace2D *fespace = A.GetStructure()->GetTestSpace2D();
  int nDOFActive = fespace->GetN_ActiveDegrees();
  int nDOF = fespace->GetN_DegreesOfFreedom();

  double *y=new double[nDOF];
  double value;
  int index;

  for(int i=0;i<nDOFActive;i++)
  {
    value = 0;
    for (int j=ARowPtr[i]; j<ARowPtr[i+1]; j++)
    {
      index = AColIndex[j];
      value += AEntries[j] * x[index];
    }
    y[i] = value;
  }
  
  for (int i=nDOFActive; i<nDOF; i++)
  {
    y[i]=x[i];
  }
  return y;
}
