#include <FEMatrix.h>
#include <MooNMD_Io.h>
#include <algorithm>


FEMatrix::FEMatrix(const TFESpace1D * space)
 : TMatrix(std::make_shared<TStructure>(space)), 
   AnsatzSpace1D(space), AnsatzSpace2D(nullptr), AnsatzSpace3D(nullptr),
   TestSpace1D(space), TestSpace2D(nullptr), TestSpace3D(nullptr)
{
  
}

FEMatrix::FEMatrix(const TFESpace2D * space)
: TMatrix(std::make_shared<TStructure>(space)),
  AnsatzSpace1D(nullptr), AnsatzSpace2D(space), AnsatzSpace3D(nullptr),
  TestSpace1D(nullptr), TestSpace2D(space), TestSpace3D(nullptr)
{
  
}

#ifdef __3D__
FEMatrix::FEMatrix(const TFESpace3D * space)
: TMatrix(std::make_shared<TStructure>(space)),
  AnsatzSpace1D(nullptr), AnsatzSpace2D(nullptr), AnsatzSpace3D(space),
  TestSpace1D(nullptr), TestSpace2D(nullptr), TestSpace3D(space)
{
  
}
#endif // 3D

FEMatrix::FEMatrix(const TFESpace2D * testspace, const TFESpace2D * ansatzspace, bool is_empty)
 : TMatrix(std::make_shared<TStructure>(testspace, ansatzspace, is_empty)),
   AnsatzSpace1D(nullptr), AnsatzSpace2D(ansatzspace), AnsatzSpace3D(nullptr),
   TestSpace1D(nullptr), TestSpace2D(testspace), TestSpace3D(nullptr)
{
  
}

#ifdef __3D__
FEMatrix::FEMatrix(const TFESpace3D * testspace, const TFESpace3D * ansatzspace)
: TMatrix(std::make_shared<TStructure>(testspace, ansatzspace)),
  AnsatzSpace1D(nullptr), AnsatzSpace2D(nullptr), AnsatzSpace3D(ansatzspace),
  TestSpace1D(nullptr), TestSpace2D(nullptr), TestSpace3D(testspace)
{
 
}
#endif // 3D

void FEMatrix::resetActive()
{
  // numer of entries in active rows
  int nActive = this->structure->getNActiveEntries();
  std::fill(this->entries.begin(), this->entries.begin()+nActive, 0.0);
}

void FEMatrix::resetNonActive()
{
  // numer of entries in active rows
  int nActive = this->structure->getNActiveEntries();
  std::fill(this->entries.begin() + nActive, this->entries.end(), 0.0);
}

void FEMatrix::scaleActive(double factor)
{
  if(factor == 1.0)
    return; // no scaling
  if(factor == 0.0)
    this->resetActive();
  
  // numer of entries in active rows
  int nActive = this->structure->getNActiveEntries();
  std::for_each(this->entries.begin(), this->entries.begin() + nActive,
                 [factor](double & a){ a = a*factor; } );
}

void FEMatrix::addActive(const FEMatrix& m, double factor)
{
  if(this->GetStructure() != m.GetStructure()) // compare objects
  {
    ErrThrow("FEMatrix::add : the two matrices do not match.");
  }
  
  // numer of entries in active rows
  int nActive = this->structure->getNActiveEntries();
  std::transform(this->entries.begin(), this->entries.begin() + nActive,
                 m.entries.begin(), this->entries.begin(), 
                 [factor](const double & a, const double & b)
                 { return a + factor * b; } );
}

void FEMatrix::multiplyActive(const double* x, double* y, double factor) const
{
  int nActive= this->GetActiveBound();
  const int * rowPtr = this->GetRowPtr();
  const int * colIndex = this->GetKCol();
  
  for(int i=0; i<nActive; ++i)
  {
    double val=0.;
    for(int j=rowPtr[i]; j<rowPtr[i+1]; ++j)
      val += this->entries[j]*x[colIndex[j]];
    y[i] += factor*val;
  }  
}


int FEMatrix::GetActiveBound() const
{
  return structure->GetActiveBound();
}

const TFESpace1D *FEMatrix::GetTestSpace1D() const
{
  return TestSpace1D;
}

const TFESpace1D *FEMatrix::GetAnsatzSpace1D() const
{
  return AnsatzSpace1D;
}

const TFESpace2D *FEMatrix::GetTestSpace2D() const
{
  return TestSpace2D;
}

const TFESpace2D *FEMatrix::GetAnsatzSpace2D() const
{
  return AnsatzSpace2D;
}

#ifdef __3D__
const TFESpace3D *FEMatrix::GetTestSpace3D() const
{
  return TestSpace3D;
}

const TFESpace3D *FEMatrix::GetAnsatzSpace3D() const
{
  return AnsatzSpace3D;
}
#endif // 3D

const TFESpace *FEMatrix::GetTestSpace() const
{
  if(TestSpace1D)
  {
    return TestSpace1D;
  }
  else if(TestSpace2D)
  {
    return TestSpace2D;
  }
  else
  {
    return TestSpace3D;
  }
}

const TFESpace *FEMatrix::GetAnsatzSpace() const
{
  if(AnsatzSpace1D)
  {
    return AnsatzSpace1D;
  }
  else if(AnsatzSpace2D)
  {
    return AnsatzSpace2D;
  }
  else
  {
    return AnsatzSpace3D;
  }
}

const TFESpace1D *FEMatrix::GetFESpace1D() const
{
  if(this->structure->isSquare())
    return TestSpace1D;
  else
    ErrThrow("accessing FESpace for non-square matrix, but which one?");
}

const TFESpace2D *FEMatrix::GetFESpace2D() const
{
  if(this->structure->isSquare())
    return TestSpace2D;
  else
    ErrThrow("accessing FESpace for non-square matrix, but which one?");
}

#ifdef __3D__
const TFESpace3D *FEMatrix::GetFESpace3D() const
{
  if(this->structure->isSquare())
    return TestSpace3D;
  else
    ErrThrow("accessing FESpace for non-square matrix, but which one?");
}
#endif // 3D

