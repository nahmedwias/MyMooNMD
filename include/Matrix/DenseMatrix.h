/*!
 * DenseMatrix.h
 * Provides a handy wrapper for general dense band matrices stored
 * as packed matrix.
 *
 * Note that row and column indexing starts at 0 (C++ style), as in the rest of
 * ParMooN.
 *
 * It is not suited to do anything but being set, decomposed and solved.
 * Should you require more functionality, you'll have to take the trouble and
 * implement it.
 *
 * @date: 2015/06/03, import to ParMooN: 2016/06/21
 * @author: Clemens Bartsch
 */

#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include <cstddef>
#include <string>

class DenseMatrix
{

public:
    /*!
     * @brief Construct a nRows x nColumns matrix filled with 0.
     * @param[in] nRows The number of rows.
     * @param[in] nColumns The number of columns.
     */
  DenseMatrix(size_t nRows, size_t nColumns);

  /*!
   * @brief Construct a nRows x nColumns matrix filled with 0
   * and manually set leading dimension.
   * @param[in] nRows The number of rows.
   * @param[in] nColumns The number of columns.
   */
	DenseMatrix(size_t nRows, size_t nColumns, size_t leadingDimension);

	//! @brief Write a value into the matrix.
	void setEntry(size_t rowNumber, size_t columnNumber, double valueToWrite);

  //! Do and store the LU decomposition.
  void decomposeLU();

  /// Solve a linear system with LAPACK's dgetrs routine.
  /// Must make sure that decomposeLU has been called before - the method
  /// will complain if not so.
  void solve(double* rhsToSolution) const;

	/// Get a value at a certain position-
  /// this is merely used for debugging and testing.
	double getEntry(size_t rowNumber, size_t columnNumber) const;

	/// Get a value at a certain position from the LU matrix -
	/// this is merely used for debugging and testing.
	double getLUEntry(const size_t r, size_t c) const;

	/// Print the matrix the console.
	/// The method is not perfect and is intended for testing and debugging only.
	void print(std::string name) const;

	void printLU(std::string name) const;

	/// Calculates the matrix' Frobenius norm. We implement only Frobenius,
	/// because this one is independent on d.o.f. numbering and therefore
	/// handy for debugging in MPI.
	double norm() const;

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


	//Special member functions.

  //! Copy constructor. Performs deep copy.
  DenseMatrix(const DenseMatrix &obj);

  //! Move constructor.
  DenseMatrix(DenseMatrix&&);

  //! Assignment operator, using copy-and-swap idiom for deep copy.
  DenseMatrix& operator=(DenseMatrix other);

  //! Standard destructor.
  ~DenseMatrix();

private:

  //! @brief returns a pointer to the entries_ array. Use cautiously.
  double* getEntries() const{
    return entries_;
  }

	//! @brief returns a pointer to the pivotsLU_ array. Use cautiously.
	int* getPivotsLU() const{
		return pivotsLU_;
	}

	friend void swap(DenseMatrix& one, DenseMatrix& other);

	//! The number of rows.
	size_t nRows_;

	//! The number of columns.
	size_t nColumns_;

	/*! The leading dimension, for a detailed explanation see e.g.:
	* http://www-01.ibm.com/support/knowledgecenter/SSFHY8_5.3.0/com.ibm.cluster.essl.v5r3.essl100.doc/am5gr_leaddi.htm
	*/
	size_t leadingDimension_;

  /*!
   * An array which contains the entries of the matrix columnwise (FORTRAN!).
   * Each column consists of leadingDimension_ subsequent locations in the
   * array, but only the first nRows_ locations of these are used.
   * That is why we require leadingDimension >= nRows_ (usually =).
   */
  double* entries_;

  /// This should store the entries of the LU factorization, if computed.
  double* entriesLU_;

  //! The LAPACK LU decomposition calculates and stores pivot entries,
  //! which have to be passed to the solver.
  int* pivotsLU_;

};




#endif /* DENSEMATRIX_H_ */
