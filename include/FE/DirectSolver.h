// =======================================================================
// @(#)DirectSolver.h
// 
// Purpose:     solve equation system by direct solver
//
// Author:      Gunar Matthies (06.09.05)
//
// History:     start of implementation 06.09.05 (Gunar Matthies)
//
// =======================================================================

#ifndef __DIRECTSOLVER__
#define __DIRECTSOLVER__

#include <SquareMatrix2D.h>
#include <Matrix2D.h>

#ifdef __3D__
#include <SquareMatrix3D.h>
#include <Matrix3D.h>
#endif

class TDirectSolver
{
  public:
    TDirectSolver() {};
    
    virtual ~TDirectSolver() {};
    
#ifdef __3D__
    /** NSTYPE 2 */
    virtual void SetMatrix(TSquareMatrix3D *sqmatrixA,
		   TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
		   TMatrix3D *matrixB1, TMatrix3D *matrixB2, TMatrix3D *matrixB3) = 0;
		   
    /** NSTYPE 4 */
    virtual void SetMatrix(TSquareMatrix3D *sqmatrixA11, TSquareMatrix3D *sqmatrixA12,
		      TSquareMatrix3D *sqmatrixA13, TSquareMatrix3D *sqmatrixA21,
		      TSquareMatrix3D *sqmatrixA22, TSquareMatrix3D *sqmatrixA23,
		      TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
		      TSquareMatrix3D *sqmatrixA33,
		      TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
		      TMatrix3D *matrixB1, TMatrix3D *matrixB2, TMatrix3D *matrixB3) = 0;
#endif
    
    virtual void Analyse() = 0;
    virtual void Factorize() = 0;
    virtual void Solve(double *sol, double *rhs) = 0;
    virtual void FactorizeSolve(double *sol, double *rhs) = 0;
    
};

/** solve equation system */

void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol);
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag);

void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp);
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag);

void DirectSolverLong(TSquareMatrix *matrix, double *rhs, double *sol);

void DirectSolver(TSquareMatrix2D *sqmatrixA,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol);

void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
                  double *rhs1, double *rhs2, double *sol1, double *sol2, int rb_flag=3);

/*void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag=3);*/

void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag);

//****************************************************************************/
// for NSTYPE == 1
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol);

void DirectSolver(TSquareMatrix2D *sqmatrixA, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag);

//****************************************************************************/
// for NSTYPE == 2
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  TMatrix2D *matrixC,
                  double *rhs, double *sol);

//****************************************************************************/
// for NSTYPE == 4
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol);

//****************************************************************************/
// for NSTYPE == 14
//****************************************************************************/
void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
		  TSquareMatrix2D *sqmatrixC,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol);

#ifdef __3D__
//****************************************************************************/
// for NSTYPE == 2
//****************************************************************************/
void DirectSolver(TSquareMatrix3D *sqmatrixA,
                  TMatrix3D *matrixB1T, TMatrix3D *matrixB2T,
                  TMatrix3D *matrixB3T,
                  TMatrix3D *matrixB1,  TMatrix3D *matrixB2,
                  TMatrix3D *matrixB3,
                  double *rhs, double *sol);
//****************************************************************************/
// for NSTYPE == 4
//****************************************************************************/
void DirectSolver(TSquareMatrix3D *sqmatrixA11, TSquareMatrix3D *sqmatrixA12,
		  TSquareMatrix3D *sqmatrixA13,
                  TSquareMatrix3D *sqmatrixA21, TSquareMatrix3D *sqmatrixA22,
		  TSquareMatrix3D *sqmatrixA23,
                  TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
		  TSquareMatrix3D *sqmatrixA33,
                  TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T, 
                  TMatrix3D *matrixB1,  TMatrix3D *matrixB2, TMatrix3D *matrixB3,
                  double *rhs, double *sol, int flag);

void DirectSolver(TSquareMatrix3D **sqmatrices, int n_row, int n_column,
		  double *sol, double *rhs);

#endif
#endif
