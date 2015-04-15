// =======================================================================
// @(#)ParDirectSolver.h
//
// Class:      TParDirectSolver
// Purpose:    Solve equation system by OpenMP or MPI based direct solvers
//
// Author:     Sashikumaar Ganesan (30.06.09)
//
// History:    Start of implementation 30.06.09 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __PARDIRECTSOLVER__
#define __PARDIRECTSOLVER__

#include <SquareMatrix.h>
#include <Matrix.h>
#include <SquareMatrix2D.h>
#include <Matrix2D.h>
#include <SquareMatrix3D.h>
#include <Matrix3D.h>


/** general class for all parallel direct solvers */

class TParDirectSolver
{
  protected:
// =======================================================================
//  internal adress pointers of pardiso solver
// =======================================================================
    /** Internal solver memory pointer pt  **/
    void **pt;

    /** controll parameter  **/
    int *iparm;

    /** number of threads  used in pardiso **/
    int num_threads;

    /** default matrix type of pardiso is real unsymmetric*/
    int mtype;

    /** more information (default values) for solver **/
    int maxfct;         /* Maximum number of numerical factorizations.  */
    int mnum;         /* Which factorization to use. */
    int msglvl;         /* Print statistical information  */
    int nrhs;        /* number of rhs */
    int  idum;             /* integer dummy  */

    int N_Eqn; 
    int *RowPtr;
    int *KCol;
    int N_Entries;

    double ddum;             /* double dummy  */

  private:

    /** copying given parameters into inner storage places */
    int Init_SMP_Default(int InN_Eqn, int *InRowPtr, int *InKCol,
                         int InN_Entries);



  public:
    /** constructor */
     TParDirectSolver(int InN_Eqn, int *InRowPtr, int *InKCol, int InN_Entries);

     ~TParDirectSolver();

    /** methods **/

     void ParDirectSolver_SMP_Terminate(double *Values);

     void SymbolicFactorize(double *Entries);

     void FactorizeAndSolve(TSquareMatrix *matrix, double *rhs, double *sol);

     void Solve(TSquareMatrix *matrix, double *rhs, double *sol);

#ifdef __3D__

   void FactorizeAndSolve(TSquareMatrix3D *sqmatrixA11,TSquareMatrix3D *sqmatrixA12,
                            TSquareMatrix3D *sqmatrixA13, TSquareMatrix3D *sqmatrixA21,
                            TSquareMatrix3D *sqmatrixA22,TSquareMatrix3D *sqmatrixA23,
                            TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
                            TSquareMatrix3D *sqmatrixA33,
                            double *sol, double *rhs);

   void Solve(TSquareMatrix3D *sqmatrixA11,TSquareMatrix3D *sqmatrixA12,
                            TSquareMatrix3D *sqmatrixA13, TSquareMatrix3D *sqmatrixA21,
                            TSquareMatrix3D *sqmatrixA22,TSquareMatrix3D *sqmatrixA23,
                            TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
                            TSquareMatrix3D *sqmatrixA33,
                            double *sol, double *rhs);

#endif

};
#endif
