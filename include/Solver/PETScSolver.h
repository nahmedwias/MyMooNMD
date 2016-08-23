/** ************************************************************************ 
*
* @class PETScSolver
* @brief solving a linear system using PETSc
*
*
* \page page_using_petsc Using PETSc
* \brief Some examples how to use PETSc
*
* PETSc solves your problems with preconditioned Krylov subspace methods.
* You can have an overview of all Krylov methods and preconditioners at
* http://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html
*
* The default Krylov method is GMRES with ILU as preconditioner.
*
* PETSc is configured by setting two things in your input
* *.dat file:
*   - set 'solver_type: petsc'
*   - set 'petsc_arguments:' to the usual PETSc command line
*     arguments as found in the PETSc documentation
*     https://www.mcs.anl.gov/petsc/
*
*
* \section solving_scalar Solving scalar problems
* set 'petsc_arguments:' as
* \code{.sh}
* -ksp_type richardson -ksp_monitor -pc_type lu -pc_factor_mat_solver_package umfpack
* \endcode
*
*
* \section solving_saddle_point Solving saddle point problems
* set 'petsc_arguments:' as
* \code{.sh}
* -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_type schur -fieldsplit_0_ksp_atol 1.0e-13 -fieldsplit_0_ksp_rtol 0. -fieldsplit_1_ksp_atol 1.0e-13 -fieldsplit_1_ksp_rtol 0.
* \endcode
*
*
*
****************************************************************************/

#ifndef __PETSCSOLVER__
#define __PETSCSOLVER__

#include <petscksp.h> // defines type 'Mat' below and much more
#include <ParameterDatabase.h>
#include <vector>
#include <memory>

// forward declaration
class FEMatrix;
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
    std::vector<Mat> sub_petsc_mats;

    /**
     * Convert a FEBlock of a BlockFEMatrix into a sub matrix from PETSc.
     *
     * @param feblock[in] shared pointer of FEMatrix
     * @param sub_mat[in,out] PETSc matrix with the values of the feblock
     */
    void FEBlock2PETScBlock(std::shared_ptr<const FEMatrix> feblock,
    						            bool isOnDiag,
    						            Mat &sub_mat);
};
#endif // __PETSCSOLVER__

