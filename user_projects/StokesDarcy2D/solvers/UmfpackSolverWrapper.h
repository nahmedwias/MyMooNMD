/** ***********************************************************************
*
* @name       UmfpackSolver
* @brief      Wrapper for the Umfpack Solver
*
*             Provides a Wrapper to easily access the Umfpack Solver
*             found at http://www.cise.ufl.edu/research/sparse/umfpack/
*
* @author     Johannes Neumann
* @date       10.01.2012
*
***************************************************************************/

#ifndef UMFPACKSOLVERWRAPPER_H_
#define UMFPACKSOLVERWRAPPER_H_

#include <iostream>
#include <Matrix.h>

using namespace std;

/*
 * @brief	Wrapper class for the umfpack solver
 */

class UmfpackSolverWrapper {
public:

	/**
	 * @brief Creates an empty solver. Call setMatrix to initialize!
	 *
	 * @param	mtype	type of the matrix (ref. @b PARDISO_MTYPE_*)
	 */
	UmfpackSolverWrapper() {};

	/**
	 *  @brief Initializes a new Solver with the given matrix in CSR. The
	 *  indices are shifted to conform with Fortran.
	 *
	 *  @param 	ia 		Row pointers (length of @a n)
	 *  	   	ja 		column positions (length of @a values
	 *  		a		non zero values of the matrix
	 *  		n		dimension of the matrix (i.e. size)
	 *  		mtype 	type of the matrix (ref. @b PARDISO_MTYPE_*)
	 *
	 *  @returns		Pointer to the new initialized solver class
	 */
	UmfpackSolverWrapper(
			int* ia,
			int* ja,
			double* a,
			int n);

	/**
	 *  @brief Initializes a new Solver with the given matrix in CSR. The
	 *  indices are shifted to conform with Fortran.
	 *
	 *  @param 	matrix	the matrix A where Ax=b
	 *  		mtype 	type of the matrix (ref. @b PARDISO_MTYPE_*)
	 *
	 *  @returns		Pointer to the new initialized solver class
	 */
	UmfpackSolverWrapper(TMatrix* matrix);

	/**
	 * @brief	Releases all Matrices and shifts back indices.
	 */
	~UmfpackSolverWrapper();

	/**
	 *  @brief (Re-)Sets the Matrix A for an (possibly empty) solver object.
	 *
	 *  @param 	matrix	the matrix A where Ax=b
	 */
	void setMatrix(TMatrix* matrix);

	/**
	 *  @brief Substitutes the values of the matrix A with new entries.
	 *  The structure (i.e. the sparsity pattern) of the matrix cannot be
	 *  changed.
	 *
	 *  @param 	entries		the new entries for the matrix A
	 */
	void setEntries(double* entries);

	/**
	 *  @brief Solves the equation A*(solution)=rhs for solution.
	 *  The computed solution is stored in the provided array solution.
	 *
	 *  @param 	rhs			the right-hand side of the problem Ax=b
	 *  		solution	vector to store the solution into
	 */
	void solve(double* rhs, double* solution);


	/**
	 * @brief prints some information about this instance.
	 */
	void printDebug();

	/*
	 * @brief simply solves for matrix*sol = rhs
	 */
	static void solve(TMatrix *matrix, double *rhs, double *sol);


	/**
	 * @brief Error class that is thrown, if Umfpack encounters some problem.
	 */
	class UmfpackSolverWrapperError : public std::runtime_error {
	public:

		/**
		 * @brief Creates a new Error with a message.
		 *
		 * @param msg	Detailed description of the error.
		 */
		UmfpackSolverWrapperError(const std::string& msg = "")
			: runtime_error("ERROR (Umfpack): " + msg) {
			cerr << endl << this->what() << endl << endl;
		}
	};

private:
	/**
	 * @brief Performs a symbolic facorization.
	 */
	void symFactorize();

	/**
	 * @brief Performs a numerical factorization.
	 */
	void numFactorize();

	/**
	 * @brief Determines the number of processors.
	 *
	 * @returns The number of processors of the current system.
	 */
	static int get_nproc(void);

	/**
	 * @brief Handles the error codes and throws errors accordingly with
	 * detailed messages.
	 *
	 * @param ierror	error code (ref. @b PARDISO_ERROR_*)
	 * @param line		The line number of the occurrence.
	 */
	void handle_error(int ierror, int line);

	/**
	 *  @brief Initializes a new Solver with the given matrix in CSR. The
	 *  indices are shifted to conform with Fortran.
	 *
	 *  @param 	ia 		Row pointers (length of @a n)
	 *  	   	ja 		column positions (length of @a values
	 *  		a		non zero values of the matrix
	 *  		n		dimension of the matrix (i.e. size)
	 *  		mtype 	type of the matrix (ref. @b PARDISO_MTYPE_*)
	 */
	void initSolver(int* ia, int* ja, double* a, int n);

	void reorderMatrix();

	int* ia; 		/// row indices
	int* ja;		/// column pointer
	double* a;		/// non zero values
	int n_eq;		/// size of matrix
	int nja;		/// number of non zero elements

	void* symbolic;
	void* numeric;


	//dummy
	double dzero;
	int izero;
};

#endif /* UMFPACKSOLVERWRAPPER_H_ */
