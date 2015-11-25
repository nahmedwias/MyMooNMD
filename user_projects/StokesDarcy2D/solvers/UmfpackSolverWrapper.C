/*
 * PardisoSolver.cpp
 *
 *  Created on: Nov 22, 2012
 *      Author: neumann
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "UmfpackSolverWrapper.h"
#include <omp.h>
#include <Matrix.h>
#include <Database.h>
#include <umfpack.h>

void UmfpackSolverWrapper::reorderMatrix() {
//    OutPut("\tPARDISO: reordering of the columns will be performed"<<endl);
//    OutPut("\tPARDISO:   no back ordering implemented !!!"<<endl);
	int begin, end, l;
	double value;

    for(int i=0;i<n_eq;i++) {
        begin=ia[i];
        end=ia[i+1];
        for(int j=begin;j<end;j++) {
            for(int k=j+1;k<end;k++) {
                if(ja[j] > ja[k]) {
                    l = ja[j];     	value = a[j];
                    ja[j] = ja[k]; 	a[j] = a[k];
                    ja[k] = l;       	a[k] = value;
                }
            }
        }
    }
}

void UmfpackSolverWrapper::initSolver(int* ia, int* ja, double* a, int n) {

	this->a = a;
	this->ia = ia;
	this->ja = ja;
	this->n_eq = n;
	this->nja = ia[n];

	//reorderMatrix();

	symFactorize();
	numFactorize();
}

UmfpackSolverWrapper::UmfpackSolverWrapper(
		int* ia,
		int* ja,
		double* a,
		int n) {

	initSolver(ia,ja,a,n);
}

UmfpackSolverWrapper::UmfpackSolverWrapper(TMatrix* matrix) {
	initSolver(matrix->GetRowPtr(),matrix->GetKCol(),matrix->GetEntries(),
			matrix->GetN_Rows());
}

UmfpackSolverWrapper::~UmfpackSolverWrapper() {
	umfpack_di_free_symbolic(&symbolic);
	umfpack_di_free_numeric(&numeric);
}


void UmfpackSolverWrapper::symFactorize()
{
	int error = umfpack_di_symbolic(n_eq,n_eq,ia,ja,a,&symbolic,NULL,NULL);
	handle_error(error,__LINE__);
}

void UmfpackSolverWrapper::numFactorize() {
  double Info[UMFPACK_INFO];
  double Control[UMFPACK_CONTROL];
  umfpack_di_defaults(Control);
  
	int error = umfpack_di_numeric(ia,ja,a,symbolic,&numeric,Control,Info);
	handle_error(error,__LINE__);
	if(TDatabase::ParamDB->SC_VERBOSE>1)
	  Output::print<1>("  estimated condition number ", 1/Info[UMFPACK_RCOND]);
}

void UmfpackSolverWrapper::solve(double* rhs, double* solution) {
  //todo: why transposed (UMFPACK_At) instead of normal (UMFPACK_A) ?!?
	int error = umfpack_di_solve(UMFPACK_At,ia,ja,a,solution,rhs,numeric,NULL,NULL);
  handle_error(error,__LINE__);
}

/*
 * 		Setter/Getter
 */

void UmfpackSolverWrapper::setMatrix(TMatrix* matrix) {

	//TODO CB 2015/30/06 This method causes SIGSEGV trouble when called upon object created with default constructor!
	umfpack_di_free_symbolic(&symbolic);
	umfpack_di_free_numeric(&numeric);

	initSolver(matrix->GetRowPtr(),matrix->GetKCol(),matrix->GetEntries(),
			matrix->GetN_Rows());
}

void UmfpackSolverWrapper::setEntries(double* entries) {
	umfpack_di_free_numeric(&numeric);

	this->a = entries;
	numFactorize();
}

/*
 *		HELPER
 */


void UmfpackSolverWrapper::printDebug() {
	cout << "UMFPACK (INFO):" << endl;
	cout << "\t" << UMFPACK_VERSION << endl;
}

void UmfpackSolverWrapper::solve(TMatrix *matrix, double *rhs, double *sol) {
  try {
	  UmfpackSolverWrapper *solver = new UmfpackSolverWrapper(matrix);
	  solver->solve(sol,rhs);
	  delete solver;
  } catch (UmfpackSolverWrapper::UmfpackSolverWrapperError* e) {
	printf("\n\n### ERROR_UMFPACK: %s ###\n\n",e->what());
	throw;
  }
}


int UmfpackSolverWrapper::get_nproc(void)
{
  return 1;
}

void UmfpackSolverWrapper::handle_error(int ierror, int line)
{
	if (ierror==UMFPACK_OK) return;

	switch(ierror)
	{
	//WARNINGS
	case UMFPACK_WARNING_singular_matrix:
		cout << "Warning (UMFPACK): Matrix is singular!\n";
    throw ;
		break;
	case UMFPACK_WARNING_determinant_underflow:
		cout << "Warning (UMFPACK): Determinant smaller than eps\n";
		break;
	case UMFPACK_WARNING_determinant_overflow:
		cout << "Warning (UMFPACK): Determinant is larger than IEEE Inf\n";
		break;
	//ERRORS
	case UMFPACK_ERROR_out_of_memory:
		throw new UmfpackSolverWrapperError("Out of Memory");
		break;
	case UMFPACK_ERROR_invalid_Numeric_object:
		throw new UmfpackSolverWrapperError("Invalid numeric factorization object");
		break;
	case UMFPACK_ERROR_invalid_Symbolic_object:
		throw new UmfpackSolverWrapperError("Invalid symbolic factorization object");
		break;
	case UMFPACK_ERROR_argument_missing:
		throw new UmfpackSolverWrapperError("Argument Missing.");
		break;
	case UMFPACK_ERROR_n_nonpositive:
		throw new UmfpackSolverWrapperError("Matrix dimensions not positive.");
		break;
	case UMFPACK_ERROR_invalid_matrix:
		throw new UmfpackSolverWrapperError("Invalid Matrix Structure.");
		break;
	case UMFPACK_ERROR_different_pattern:
		throw new UmfpackSolverWrapperError("Different sparse pattern.");
		break;
	case UMFPACK_ERROR_invalid_system:
		throw new UmfpackSolverWrapperError("Invalid system provided with sys.");
		break;
	case UMFPACK_ERROR_invalid_permutation:
		throw new UmfpackSolverWrapperError("Invalid permutation vector.");
		break;
	case UMFPACK_ERROR_file_IO:
		throw new UmfpackSolverWrapperError("Fille IO error.");
		break;
	case UMFPACK_ERROR_internal_error:
		throw new UmfpackSolverWrapperError("Internal error.");
		break;

	default:
		throw new UmfpackSolverWrapperError("Error number " + ierror); break;
	break;
	}
}
