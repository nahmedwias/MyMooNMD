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

/**
 * 
 */

#include <Matrix.h>
#include <BlockVector.h>

// forward declaration
class ColoredBlockMatrix;

class DirectSolver
{
  public:
    enum class DirectSolverTypes {umfpack, pardiso};
    /**
     * @brief compute the factorization of a matrix, ready to call solve
     * 
     * This calls the other (private) constructor with 
     * matrix.get_combined_matrix().
     *
     * @param  matrix the matrix A where Ax=b
     */
    DirectSolver(const ColoredBlockMatrix& matrix, DirectSolverTypes type);
    
    /** @brief release all memory, the matrix is not touched */
    ~DirectSolver();

    /**
     * @brief Solves the equation A*(solution)=rhs for solution.
     * 
     * This calls the other (private) method solve.
     *
     * @param rhs the right-hand side of the problem Ax=b
     * @param solution vector to store the solution into
     */
    void solve(BlockVector& rhs, BlockVector& solution);
    
  private:
    /** @brief type of direct solver used */
    DirectSolverTypes type;
    
    /** @brief the matrix of the linear equation A*x=b */
    std::shared_ptr<TMatrix> matrix;
    
    /** @brief storage for umfpack direct solver */
    //@{
    void* symbolic;
    void* numeric;
    //@}
    
    /**
     * @brief compute the factorization of a matrix, ready to call solve
     * 
     * The indices are shifted to conform with Fortran in case of the Pardiso 
     * solver. This is why this methods has to change the matrix and therefore
     * does not take a const TMatrix.
     *
     * @param  matrix  the matrix A where Ax=b
     */
    DirectSolver(std::shared_ptr<TMatrix> matrix, DirectSolverTypes type);
    
    /**
     * @brief Solves the equation A*(solution)=rhs for solution.
     * 
     * The computed solution is stored in the provided array solution.
     *
     * @param rhs the right-hand side of the problem Ax=b
     * @param solution vector to store the solution into
     */
    void solve(double* rhs, double* solution);
    
    /** @brief true if indices start with 1 (Fortran style, for pardiso) */
    bool isFortranShifted;
    
    void symetric_factorize();
    void numeric_factorize();
    
    /**
     * @brief shifts back and forth the index to comply with fortran.
     * 
     * @Note this changes the matrix. It is not usable as usual until this is
     * called a second time.
     */
    void fortranShift();
};


#include <SquareMatrix2D.h>
#include <Matrix2D.h>

#ifdef __3D__
#include <SquareMatrix3D.h>
#include <Matrix3D.h>
#endif

/** solve equation system */

void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol);
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag);

void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp);
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag);

void DirectSolverLong(TSquareMatrix *matrix, double *rhs, double *sol);

void DirectSolver(const TSquareMatrix2D *sqmatrixA,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  double *rhs, double *sol);

void DirectSolver(const TSquareMatrix2D *sqmatrixA11, const TSquareMatrix2D *sqmatrixA12,
                  const TSquareMatrix2D *sqmatrixA21, const TSquareMatrix2D *sqmatrixA22,
                  double *rhs1, double *rhs2, double *sol1, double *sol2, int rb_flag=3);

/*void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
                  TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
                  TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                  TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag=3);*/

void DirectSolver(const TSquareMatrix2D *sqmatrixA11, const TSquareMatrix2D *sqmatrixA12,
                  const TSquareMatrix2D *sqmatrixA21, const TSquareMatrix2D *sqmatrixA22,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag);

//****************************************************************************/
// for NSTYPE == 1
//****************************************************************************/
void DirectSolver(const TSquareMatrix2D *sqmatrixA,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  double *rhs, double *sol);

void DirectSolver(const TSquareMatrix2D *sqmatrixA,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  double *rhs, double *sol, int rb_flag);

//****************************************************************************/
// for NSTYPE == 2
//****************************************************************************/
void DirectSolver(const TSquareMatrix2D *sqmatrixA,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  const TMatrix2D *matrixC,
                  double *rhs, double *sol);

//****************************************************************************/
// for NSTYPE == 4
//****************************************************************************/
void DirectSolver(const TSquareMatrix2D *sqmatrixA11, const TSquareMatrix2D *sqmatrixA12,
                  const TSquareMatrix2D *sqmatrixA21, const TSquareMatrix2D *sqmatrixA22,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  double *rhs, double *sol);

//****************************************************************************/
// for NSTYPE == 14
//****************************************************************************/
void DirectSolver(const TSquareMatrix2D *sqmatrixA11, const TSquareMatrix2D *sqmatrixA12,
                  const TSquareMatrix2D *sqmatrixA21, const TSquareMatrix2D *sqmatrixA22,
                  const TSquareMatrix2D *sqmatrixC,
                  const TMatrix2D *matrixB1T, const TMatrix2D *matrixB2T,
                  const TMatrix2D *matrixB1,  const TMatrix2D *matrixB2,
                  double *rhs, double *sol);

//****************************************************************************/
// for Darcy
//****************************************************************************/
void DirectSolver(const TSquareMatrix2D *sqmatrixA, const TSquareMatrix2D *sqmatrixC,
                  const TMatrix2D *matrixBT, const TMatrix2D *matrixB,
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
