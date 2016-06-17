/*!
 * FortranStyleMatrix.C
 *
 * Implements class FortranStyleMatrix declared in FortranStyleMatrix.h.
 *
 *  @date: 2015/06/03
 *  @author: Clemens Bartsch
 */

#include <MooNMD_Io.h>
#include <Constants.h>
#include <LinAlg.h>

#include <cstring>
#include <stdexcept>

#include <FortranStyleMatrix.h>

extern "C" {

void dgetrf_(int m, int n, double *a, int lda, int *ipiv, int *info);

void dgetrs_(char transa, int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb, int *info);
}

/*!
 * @brief Construct a nRows x nColumns matrix filled with 0.
 * @param[in] nRows The number of rows.
 * @param[in] nColumns The number of columns.
 */
FortranStyleMatrix::FortranStyleMatrix(const size_t nRows, const size_t nColumns) :
	nRows_(nRows), nColumns_(nColumns){
	// The way we intend to used the matrix, leading dimension equals number of rows.
	leadingDimension_ = nRows_;
	// Allocate the entries_ array and fill it with zeroes. The latter is achieved by "()".
	entries_ = new double[leadingDimension_*nColumns_]();
	// Allocate the LU pivot array and fill it with zeroes.
	pivotsLU_= new int[nColumns_]();
}

/*!
 * @brief Construct a nRows x nColumns matrix filled with 0 and manually sets leading dimension.
 * @param[in] nRows The number of rows.
 * @param[in] nColumns The number of columns.
 */
FortranStyleMatrix::FortranStyleMatrix(const size_t nRows, const size_t nColumns, const size_t leadingDimension) :
	nRows_(nRows), nColumns_(nColumns), leadingDimension_(leadingDimension){
	if(leadingDimension < nRows){
		// Leading dimension may not be lesser than number of rows.
		leadingDimension_ = nRows_;
		OutPut("Error in FortranStyle Matrix: leading dimension must be greater or equal number of rows, I did set leadingDimension_ = nRows_.");
	}
	// Allocate the entries_ array and fill it with zeroes. The latter is achieved due to "()".
	entries_ = new double[leadingDimension_*nColumns_]();
	// Allocate the pivotsLU_ array and fill it with zeroes. The latter is achieved due to "()".
	pivotsLU_= new int[nColumns_]();
}

//! @brief Standard destructor.
FortranStyleMatrix::~FortranStyleMatrix(){
	//Delete the entries_ array.
	delete[] entries_; entries_=nullptr;
	delete[] pivotsLU_; pivotsLU_=nullptr;
}

/*!
 * @brief Copy constructor. Performs deep copy.
 * @param[in] obj The object to be copied.
 */
FortranStyleMatrix::FortranStyleMatrix(const FortranStyleMatrix &obj) :
		nRows_(obj.getNRows()), nColumns_(obj.getNColumns()), leadingDimension_(obj.getLeadingDimension())
{
	// Allocate the entries_ array.
	entries_ = new double[leadingDimension_*nColumns_];
	// And fill it with the entries from obj.
	memcpy(entries_,obj.getEntries(),leadingDimension_*nColumns_*SizeOfDouble);

	// Allocate the LU pivots array.
	pivotsLU_= new int[nColumns_];
	//And fill it.
	memcpy(pivotsLU_,obj.getPivotsLU(),nColumns_*SizeOfInt);
}

/*! @brief Copy assignment. Performs deep copy.
 *  @param[in] obj The object to be copied.
 */
FortranStyleMatrix& FortranStyleMatrix::operator=( const FortranStyleMatrix& obj ){
    nRows_ = obj.getNRows();
    nColumns_ = obj.getNColumns();
    leadingDimension_= obj.getLeadingDimension();

    // Allocate the entries_ array.
    entries_ = new double[leadingDimension_*nColumns_];
    // And fill it with the entries from obj.
    memcpy(entries_,obj.getEntries(),leadingDimension_*nColumns_*SizeOfDouble);

	// Allocate the LU pivots array.
	pivotsLU_= new int[nColumns_];
	//And fill it.
	memcpy(pivotsLU_,obj.getPivotsLU(),nColumns_*SizeOfInt);

    return *this;
}

/*!
 * @brief Write a value into FortranStyle style matrix.
 * @param[in] rowNumber The row Index of the new entry. Begin counting with index 0.
 * @param[in] columnNumber The column Index of the new entry. Begin counting with index 0.
 *
 * @throws Runtime error if row or column number are out of bound.
 *
 * TODO This is probably slow.
 */
void FortranStyleMatrix::setEntry(const size_t rowNumber, const size_t columnNumber, const double valueToWrite){
	// Throw if wrong line or column number is entered.
	if (rowNumber >= nRows_ || columnNumber >= nColumns_)   {
		throw std::runtime_error("FortranStyleMatrix::setEntry: Out of matrix dimension error.");
	}
	//Write the value to the right place.
	entries_[columnNumber*leadingDimension_ + rowNumber] = valueToWrite;
}

/*!
 * @brief Get a value from a FortranStyle style matrix.
 * @param[in] rowNumber The row Index of the entry to get. Begin counting with index 0.
 * @param[in] columnNumber The column Index of the entry to get. Begin counting with index 0.
 * @return The entry at the specified place.
 *
 * @throws Runtime error if row or column number are out of bound.
 */
double FortranStyleMatrix::getEntry(const size_t rowNumber, const size_t columnNumber) const {
	// Throw if wrong line or column number is entered.
	if (rowNumber >= nRows_ || columnNumber >= nColumns_)   {
		throw std::runtime_error("FortranStyleMatrix::getEntry: Out of matrix dimension error.");
	}
	//Get the entry from the right place.
	return entries_[columnNumber*leadingDimension_ + rowNumber];
}

/*!
 * @brief Solve a linear system with a FortranStyle solver routine.
 * @param[in,out] rhsToSolution An array which contains the rhs on input and the solution calculated by FortranStyle on output.
 *
 * The memory space for the in/out array and the correct array length must be taken care of outside the class.
 */
void FortranStyleMatrix::solve(double* rhsToSolution) const{
	//OutPut("Solve FStyleMatrix." << endl);

	// Throw if the matrix is not quadratic.
	if(!(nRows_ == nColumns_ && nColumns_ == leadingDimension_))
		throw std::runtime_error("FortranStyleMatrix::solve: Solver is only operating when "
				"nRows_ == nColumns_ == leadingDimension_ so far.");

	//Solve, assuming the LU decomposition was performed before.
	int info;
	   dgetrs_('n' , nColumns_ , 1, entries_,leadingDimension_, pivotsLU_, rhsToSolution, nColumns_,&info);
}

/*!
 * Do and store the LU decomposition directly in entries_ to save memory.
 */
void FortranStyleMatrix::decomposeLU(){
	//little aux variable needed by the fortran method, 0 if succesful
	int info;
	//Do the LU decomp and store in entries_ and pivotsLU_
	   dgetrf_(nRows_, nColumns_, entries_, leadingDimension_, pivotsLU_, &info);
}

void FortranStyleMatrix::print() const{
	OutPut("FortranStyleMatrix = " << endl);
	for (unsigned int i = 0; i < nRows_ ;++i){
		OutPut("Row " << i << " ");
		for(unsigned int j =0; j < nColumns_; ++j){
			OutPut(getEntry(i,j) << " ");
		}
		OutPut(endl);
	}
}
