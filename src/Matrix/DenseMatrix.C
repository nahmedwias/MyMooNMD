/*!
 * DenseMatrix.C
 *
 * Implements class DenseMatrix declared in DenseMatrix.h.
 *
 *  @date: 2015/06/03
 *  @author: Clemens Bartsch
 */

#include <MooNMD_Io.h>
#include <Constants.h> //needed for SizeOfDouble macro

#include <cassert> //not the most sophisticated solution, but fit four our case
#include <cstring>
#include <stdexcept>

#include <DenseMatrix.h>

// Extern declaration of Lapack methods.
extern "C" {
// LU factorize...
void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info);
// ... and solve!
void dgetrs_(char* transa, int* n, int* nrhs, double *a, int* lda,
            int *ipiv, double *b, int* ldb, int *info);
}

DenseMatrix::DenseMatrix(const size_t nRows, const size_t nColumns) :
	nRows_(nRows), nColumns_(nColumns){

	// The way we intend to used the matrix, leading dimension is number of rows.
	leadingDimension_ = nRows_;

	// Allocate the entries_ array and fill it with zeroes. "()".
	entries_ = new double[leadingDimension_*nColumns_]();

	// The memory needed for LU factorization is not allocated yet.
	pivotsLU_= nullptr;
	entriesLU_ = nullptr;
}


DenseMatrix::DenseMatrix(
    size_t nRows, size_t nColumns, size_t leadingDimension) :
    nRows_(nRows), nColumns_(nColumns), leadingDimension_(leadingDimension)
{
  if(leadingDimension < nRows){
    // Leading dimension may not be lesser than number of rows.
    ErrThrow("Error in DenseMatrix leading dimension must be greater"
        " or equal number of rows.");
  }
  // Allocate the entries_ array and fill it with zeroes. The latter is achieved due to "()".
  entries_ = new double[leadingDimension_*nColumns_]();

  // The memory needed for LU factorization is not allocated yet.
  pivotsLU_= nullptr;
  entriesLU_ = nullptr;
}

DenseMatrix::~DenseMatrix(){
	delete[] entries_;
	delete[] pivotsLU_;
	delete[] entriesLU_;
}

// TODO
DenseMatrix::DenseMatrix(const DenseMatrix &obj) :
		nRows_(obj.getNRows()), nColumns_(obj.getNColumns()),
		leadingDimension_(obj.getLeadingDimension())
{
	// Allocate the entries_ array.
	entries_ = new double[leadingDimension_*nColumns_];
	// And fill it with the entries from obj.
	memcpy(entries_,obj.getEntries(),leadingDimension_*nColumns_*SizeOfDouble);

	pivotsLU_= nullptr;
	entriesLU_ = nullptr;

	if(obj.pivotsLU_)
	{
	  pivotsLU_ = new int[nColumns_];
  	memcpy(pivotsLU_,obj.getPivotsLU(),nColumns_*SizeOfInt);
	}

	if(obj.entriesLU_)
	{
    entriesLU_ = new double[leadingDimension_*nColumns_];
    memcpy(entriesLU_,obj.getEntries(),leadingDimension_*nColumns_*SizeOfDouble);
	}
}

DenseMatrix& DenseMatrix::operator=(DenseMatrix other)
{
  //do a swap with the copy constructed object "other"
  swap(*this, other);

  return *this;
}

void DenseMatrix::setEntry(size_t r, size_t c, double value){

  // Throw if wrong line or column number is entered.
  assert(r < nRows_ && c < nColumns_);

	//Write the value to the right place.
	entries_[c*leadingDimension_ + r] = value;
}


double DenseMatrix::getEntry( size_t r, size_t c) const {
	// Throw if wrong line or column number is entered.
	assert (r < nRows_ && c < nColumns_);

	//Get the entry from the right place.
	return entries_[c*leadingDimension_ + r];
}

double DenseMatrix::getLUEntry(const size_t r, size_t c) const {

  if(!entriesLU_)
    ErrThrow("You have to factorize before you can view the factorization!");

  // Throw if wrong line or column number is entered.
  assert (r < nRows_ && c < nColumns_);

  //Get the entry from the right place.
  return entriesLU_[c*leadingDimension_ + r];
}

/*!
 * @brief Solve a linear system with a FortranStyle solver routine.
 * @param[in,out] rhsToSolution An array which contains the rhs on input and the solution calculated by FortranStyle on output.
 *
 * The memory space for the in/out array and the correct array length must be taken care of outside the class.
 */
void DenseMatrix::solve(double* rhsToSolution) const{

	// Throw if the matrix is not quadratic.
	if(!(nRows_ == nColumns_ && nColumns_ == leadingDimension_))
		throw std::runtime_error("DenseMatrix::solve: Solver is only operating when "
				"nRows_ == nColumns_ == leadingDimension_ so far.");

	//Solve, assuming the LU decomposition was performed before.
	if(!entriesLU_ || !pivotsLU_)
	{
	  ErrThrow("DenseMatrix::decomposeLU must be called"
	      " before DenseMatrix::solve");
	}

	// Call the solve routine.
	int info;
	char control[1] = {'n'};
	int nC = nColumns_;
	int LDA = leadingDimension_;
	int nrhs = 1;
	dgetrs_(control , &nC  , &nrhs, entriesLU_,&LDA,
	       pivotsLU_, rhsToSolution, &nC,&info);

  if(info != 0)
    ErrThrow("LAPACK dgetrs (solve) failed with info = ", info);
}

void DenseMatrix::decomposeLU(){

  // get memory for LU fact if necessary and copy entries into it
  if(!entriesLU_)
    entriesLU_ = new double[leadingDimension_*nColumns_]();
  if(!pivotsLU_)
    pivotsLU_= new int[nColumns_]();

  memcpy(entriesLU_, entries_, leadingDimension_*nColumns_*SizeOfDouble);

	//Do the LU decomp and store in entries_ and pivotsLU_
  int info = 0;
  int nR = nRows_;
  int nC = nColumns_;
  int LDA = leadingDimension_;
	dgetrf_(&nR, &nC, entriesLU_, &LDA, pivotsLU_, &info);

	if(info != 0)
	  ErrThrow("LAPACK dgetrf (factorize) failed with info = ", info);
}

void DenseMatrix::print(std::string name) const
{
	Output::info("DenseMatrix: ", name);
	for (size_t i = 0; i < nRows_ ;++i){
		OutPut("Row " << i << " ");
		for(size_t j =0; j < nColumns_; ++j){
			OutPut(getEntry(i,j) << " ");
		}
		OutPut(endl);
	}
}

void DenseMatrix::printLU(std::string name) const
{
  assert(entriesLU_ && pivotsLU_);

  Output::info("DenseMatrix LU decomp: ", name);
  for (size_t i = 0; i < nRows_ ;++i)
  {
    OutPut("Row " << i << " ");
    for(size_t j =0; j < nColumns_; ++j)
    {
      OutPut(getLUEntry(i,j) << " ");
    }
    OutPut(endl);
  }

  Output::info("DenseMatrix LU pivots: ", name);
  for(size_t j =0; j < nColumns_; ++j)
    OutPut(pivotsLU_[j] << " ");
  OutPut(endl);
}

double DenseMatrix::norm() const
{
  double sum = 0;
  for(size_t k = 0; k < leadingDimension_*nColumns_ ;++k)
  {
    sum += entries_[k]*entries_[k];
  }
  return sqrt(sum);
}

void swap(DenseMatrix& first, DenseMatrix& second)
{
  //std::swap all members
  std::swap(first.nRows_, second.nRows_);
  std::swap(first.nColumns_, second.nColumns_);
  std::swap(first.leadingDimension_, second.leadingDimension_);

  std::swap(first.entries_, second.entries_);
  std::swap(first.entriesLU_, second.entriesLU_);
  std::swap(first.pivotsLU_, second.pivotsLU_);

}
