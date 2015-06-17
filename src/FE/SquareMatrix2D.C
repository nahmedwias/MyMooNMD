// =======================================================================
// @(#)SquareMatrix2D.C        1.2 11/20/98
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

#include <SquareMatrix2D.h>
#include <string.h>
#include <LinAlg.h>
#include <stdlib.h>

TSquareMatrix2D::TSquareMatrix2D(TSquareStructure2D *squarestructure)
  : TSquareMatrix(squarestructure), structure(squarestructure)
{
}


// set all class variables
void TSquareMatrix2D::SetStructure(TSquareStructure2D *squarestructure)
{
  structure = squarestructure;
  this->TSquareMatrix::SetStructure((TSquareStructure*) squarestructure);
}


TSquareMatrix2D::TSquareMatrix2D(int n) 
 : structure(new TSquareStructure2D(n)), TSquareMatrix(structure)
{
}

TSquareMatrix2D::~TSquareMatrix2D()
{
}


void TSquareMatrix2D::reset_non_active()
{
  int n_active_rows = this->structure->GetFESpace()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  int index_nonactive = rowPtr[n_active_rows];
  int n_nonactive_entries = rowPtr[structure->GetN_Rows()] - index_nonactive;
  memset(Entries + index_nonactive, 0.0, n_nonactive_entries * SizeOfDouble);
}

void TSquareMatrix2D::reset_active()
{
  int n_active_rows = this->structure->GetFESpace()->GetN_ActiveDegrees();
  int * rowPtr = this->structure->GetRowPtr();
  // numer of entries in active rows
  int n_active = rowPtr[n_active_rows];
  memset(this->Entries, 0.0, n_active * SizeOfDouble);
}

void TSquareMatrix2D::scale_active(double factor)
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

void TSquareMatrix2D::add_active(const TSquareMatrix2D& m, double factor)
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

TSquareMatrix2D& TSquareMatrix2D::operator=(const TSquareMatrix2D& rhs)
{
  // compare structures (first pointers, then structures themselves)
  if(structure != rhs.structure && *structure != *rhs.structure)
  {
    OutPut("WARNING: TSquareMatrix2D& operator= Matrices have different "
           << "fe space or structure. Are you sure this is what you want?" << endl);
    structure = rhs.structure;
  }
  // no further tests, make sure you know these two matrices have the same 
  // structure. We cannot compare the two TSquareStructure2D, because during 
  // intitialization every matrix might get its own TSquareStructure2D.
  
  // copy matrix entries.
  for(int i=0; i<structure->GetN_Entries(); i++)
    Entries[i] = rhs.Entries[i];
  
  return *this;
}


TSquareMatrix2D& TSquareMatrix2D::operator*=(double alpha)
{
  this->scale_active(alpha);
  return *this;
}


TSquareMatrix2D& TSquareMatrix2D::operator+=(const TSquareMatrix2D & rhsMat)
{

  this->add_active(rhsMat, 1.0);
  return *this;
}

// overloaded operators 
// add to matrices A and B
// note: only active DOF are added
// note: only works for matrices with the same sparsity pattern
TSquareMatrix2D& operator+(const TSquareMatrix2D & A, const TSquareMatrix2D & B)
{
  double *AEntries, *BEntries, *CEntries;
  if (A.GetMatrixStructure() == B.GetMatrixStructure()) 
  {
    // create bew TSquareMatrix2D on heap (otherwise return did not work)
    TSquareMatrix2D *C = new TSquareMatrix2D(A.GetMatrixStructure());
    AEntries = A.GetEntries();
    BEntries = B.GetEntries();
    CEntries = C->GetEntries();
    for (int i=0; i<A.GetActiveBound(); i++) 
    {
      CEntries[i] = AEntries[i] + BEntries[i];
    }
    return *C;
  } else 
  {
    cout << " SquareMatrix2D: ERROR: can't add Matrices "
       << " with different sparse structures " << endl;
    exit(1);
  }
}

// C= A*alpha 
// note: only active DOF are multiplied, others are just copied
TSquareMatrix2D& operator*(const TSquareMatrix2D & A,const double alpha)
{
  double *AEntries, *CEntries;
  TSquareMatrix2D *C = new TSquareMatrix2D(A.GetMatrixStructure());
  AEntries = A.GetEntries();
  CEntries = C->GetEntries();
  // multiply each active entry by alpha and write it into matrix C
  for (int i=0; i<A.GetActiveBound(); i++) 
  {
    CEntries[i] = alpha*AEntries[i];
  }
  // non active entries are just copied
  for (int i=A.GetActiveBound(); i<A.GetN_Entries(); i++) 
  {
    CEntries[i] = AEntries[i];
  }
  return *C;
}
TSquareMatrix2D& operator*(const double alpha,const TSquareMatrix2D & A)
{// just to allow alpha*A as well (in addition to A*alpha)
  return A*alpha;
}



double* operator*(const TSquareMatrix2D & A,const double* x)
{
  double *AEntries;
  int *ARowPtr,*AColIndex;
  AEntries = A.GetEntries();
  ARowPtr = A.GetRowPtr();
  AColIndex = A.GetKCol();

  int nDOFActive = A.GetActiveBound();
  int nDOF = A.GetFESpace()->GetN_DegreesOfFreedom();
  double *y=new double[nDOF];
  double value;
  int index;

  // multiply each active entry by alpha and write it into y
  for(int i=0;i<nDOFActive;i++) {
    value = 0;
    for (int j=ARowPtr[i]; j<ARowPtr[i+1]; j++) {
      index = AColIndex[j];
      value += AEntries[j] * x[index];
    }
    y[i] = value;
  }
  
  for (int i=nDOFActive; i<nDOF; i++) {
    y[i]=x[i];
  }
  
  return y;
}


