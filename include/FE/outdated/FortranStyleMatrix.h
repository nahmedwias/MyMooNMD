/*!
 * FortranStyleMatrix.h
 *
 *  @date: 2015/06/03
 *  @author: Clemens Bartsch
 */

#ifndef FORTRANSTYLEMATRIX_H_
#define FORTRANSTYLEMATRIX_H_

#include <cstddef>

/*! @brief Wraps up a matrix in a fashion the FortranStyle solver needs it.
 *
 * This class is intended to make the handling of FortranStyle style matrices more convenient by wrapping up a FORTRAN style matrix.
 *
 */
class FortranStyleMatrix{

public:
	//! @brief Construct a nRows x nColumns matrix filled with 0.
	FortranStyleMatrix(const size_t nRows, const size_t nColumns);

	//! @brief Construct a nRows x nColumns matrix filled with 0 and manually set leading dimension.
	FortranStyleMatrix(const size_t nRows, const size_t nColumns, const size_t leadingDimension);

	//! @brief Copy constructor. Performs deep copy.
	FortranStyleMatrix(const FortranStyleMatrix &obj);

	//! @brief Copy assignment. Performs deep copy.
	FortranStyleMatrix& operator=( const FortranStyleMatrix& obj );

	//! @brief Standard destructor.
	~FortranStyleMatrix();

	//! @brief Write a value into FortranStyle style matrix.
	void setEntry(const size_t rowNumber, const size_t columnNumber, const double valueToWrite);

	//! @brief Get a value from a FortranStyle style matrix.
	double getEntry(const size_t rowNumber, const size_t columnNumber) const;

	//! Do and store the LU decomposition.
	void decomposeLU();

	//! @brief Solve a linear system with a FortranStyle solver routine.
	void solve(double* rhsToSolution) const;

	void print() const;

	//! Getters
	size_t getNRows() const{
		return nRows_;
	}

	size_t getNColumns() const{
		return nColumns_;
	}

	size_t getLeadingDimension() const{
		return leadingDimension_;
	}

	//! @brief returns a pointer to the entries_ array. Use cautiously.
	double* getEntries() const{
		return entries_;
	}
private:
	//! @brief returns a pointer to the pivotsLU_ array. Use cautiously.
	int* getPivotsLU() const{
		return pivotsLU_;
	}

	/*! An array which contains the entries of the matrix columnwise. Using FORTRAN style, each column consists of
	* leadingDimension_ subsequent locations in the array, but only the first nRows_ locations of these are used.
	* That is why we require leadingDimension >= nRows_ (usually =).
	* After call of decomposeLU() this stores the LU factorization.
	*/
	double* entries_;

	//! The LAPACK LU decomposition calculates and stores pivot entries, which have to be passed to the solver.
	int* pivotsLU_;

	//! The number of rows.
	size_t nRows_;

	//! The number of columns.
	size_t nColumns_;

	/*! The leading dimension, for a detailed explanation see e.g.:
	* http://www-01.ibm.com/support/knowledgecenter/SSFHY8_5.3.0/com.ibm.cluster.essl.v5r3.essl100.doc/am5gr_leaddi.htm
	*/
	size_t leadingDimension_;
};



#endif /* FORTRANSTYLEMATRIX_H_ */
