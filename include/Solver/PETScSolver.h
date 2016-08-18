/** ************************************************************************ 
*
* @class PETScSolver
* @brief solving a linear system using PETSc
****************************************************************************/

#ifndef __PETSCSOLVER__
#define __PETSCSOLVER__

#include <petscksp.h> // defines type 'Mat' below and much more
#include <ParameterDatabase.h>

// forward declaration
class BlockMatrix;
class BlockFEMatrix;
class BlockVector;

class PETScSolver
{
  public:
	PETScSolver(const BlockFEMatrix& matrix, const ParameterDatabase& db);

    /** @brief constructor: an internal copy is created and stored. Ready to 
     * solve.
     */
    PETScSolver(const BlockMatrix& matrix, const ParameterDatabase& db);
    
    /** @brief This class is not copy constructible */
    PETScSolver(const PETScSolver&) = delete;

    /** @brief move constructor */
    PETScSolver(PETScSolver&&);

    /** @brief This class is not copy assignable */
    PETScSolver& operator=(const PETScSolver&) = delete;

    /** @brief move assignment operator */
    PETScSolver& operator=(PETScSolver&&);
    
    /** @brief release all memory */
    ~PETScSolver();

    /**
     * @brief Solves the equation A*(solution)=rhs for solution.
     * 
     * @param rhs the right-hand side of the problem Ax=b
     * @param solution vector to store the solution into
     */
    void solve(const BlockVector& rhs, BlockVector& solution);
    
  private:
    /** @brief the matrix of the linear equation A*x=b, PETSc format */
    Mat petsc_mat;
};
#endif // __PETSCSOLVER__

