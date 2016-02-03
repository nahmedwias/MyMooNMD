/** ************************************************************************ 
*
* @class DirectSolver
* @brief solving a linear system using a direct solver
*
* Given a BlockMatrix the constructor of this class computes a factorization of
* that matrix. Then the method DirectSolver::solve enables the actual solving
* step which can be called multiple time, reusing the computed factorization.
* 
* Note that there is no way to update the matrix, i.e. to recompute a 
* factorization. In case your matrix changes, you should also create a new
* object of this class DirectSolver.
* 
* @ruleof0
*
****************************************************************************/

#ifndef __DIRECTSOLVER__
#define __DIRECTSOLVER__

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
    
    /** @brief This class is not copy constructible */
    DirectSolver(const DirectSolver&) = delete;

    /** @brief move constructor */
    DirectSolver(DirectSolver&&);

    /** @brief This class is not copy assignable */
    DirectSolver& operator=(const DirectSolver&) = delete;

    /** @brief move assignment operator */
    DirectSolver& operator=(DirectSolver&&);
    
    /** @brief release all memory */
    ~DirectSolver();

    /**
     * @brief Solves the equation A*(solution)=rhs for solution.
     * 
     * This calls the other (private) method solve.
     *
     * @param rhs the right-hand side of the problem Ax=b
     * @param solution vector to store the solution into
     */
    void solve(const BlockVector& rhs, BlockVector& solution);
    
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
    
    /** @brief true if indices start with 1 (Fortran style, for pardiso) */
    bool isFortranShifted;
    
    
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
    void solve(const double* rhs, double* solution);
    
    /** @brief compute symbolic factorization */
    void symetric_factorize();
    /** @brief compute numeric factorization (requires symbolic factorization) 
     */
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

[[deprecated("use the class DirectSolver instead")]]
void DirectSolver_old(TSquareMatrix *matrix, double *rhs, double *sol);
[[deprecated("use the class DirectSolver instead")]]
void DirectSolver_old(TSquareMatrix2D *sqmatrixA11,
                      TSquareMatrix2D *sqmatrixA12,
                      TSquareMatrix2D *sqmatrixA21,
                      TSquareMatrix2D *sqmatrixA22,
                      TMatrix2D *matrixB1T, TMatrix2D *matrixB2T, 
                      TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                      double *rhs, double *sol, int rb_flag);
[[deprecated("use the class DirectSolver instead")]]
void DirectSolver_old (TSquareMatrix2D *sqmatrixA, 
                       TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
                       double *rhs, double *sol, int rb_flag);
#ifdef __3D__
[[deprecated("use the class DirectSolver instead")]]
void DirectSolver_old(TSquareMatrix3D *sqmatrixA11,
                      TSquareMatrix3D *sqmatrixA12,
                      TSquareMatrix3D *sqmatrixA13,
                      TSquareMatrix3D *sqmatrixA21,
                      TSquareMatrix3D *sqmatrixA22,
                      TSquareMatrix3D *sqmatrixA23,
                      TSquareMatrix3D *sqmatrixA31,
                      TSquareMatrix3D *sqmatrixA32,
                      TSquareMatrix3D *sqmatrixA33,
                      TMatrix3D *matrixB1T, TMatrix3D *matrixB2T,
                      TMatrix3D *matrixB3T, TMatrix3D *matrixB1,
                      TMatrix3D *matrixB2, TMatrix3D *matrixB3,
                      double *rhs, double *sol, int flag);
#endif // __3D__
#endif // __DIRECTSOLVER__
