/** ************************************************************************ 
*
* @class     BlockMatrixCD3D
* @brief     stores the information of a 3D scalar system matrix 
* @author    Sashikumaar Ganesan 
* @date      23.01.15
* @History    
 ************************************************************************  */


#ifndef __SYSTEMMATSCALAR3D__
#define __SYSTEMMATSCALAR3D__

//#include <AssembleMat3D.h>

#ifdef _MPI
//#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>

#include <ParDirectSolver.h>
#endif
#include <BlockMatrix3D.h>
#include <ItMethod.h>

#ifdef _OMPONLY
#include <ParDirectSolver.h>
#endif

/** class for 3D scalar system matrix */
class BlockMatrixCD3D : public BlockMatrix3D
{
  protected:
    /** Boundary condition and Boundary Value */
    BoundCondFunct3D *BoundaryConditions[1];
    
    /** Boundary condition and Boundary Value */
    BoundValueFunct3D *BoundaryValues[1];
    
    /** instance of the Assemble class */
    //TAssembleMat3D **AMatRhsAssemble;
    
  public:
    /** Constructor*/
    BlockMatrixCD3D(TFESpace3D *fespace);
    
    /** destrcutor */
    ~BlockMatrixCD3D();
    
    /** Initilize the discrete forms and the matrices */
    void Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond,
              BoundValueFunct3D *BoundValue);
    
    /** assemble the system matrix */
    void Assemble(CoeffFct3D *BilinearCoeffs, double *sol, double *rhs);
    
    /** @brief compute y = factor* A*x 
     *
     * write the matrix-vector product "Ax" scaled by a factor to y. 
     * "A" is this matrix. Both "A" and "x" remain unchanged. If the factor is
     * 0.0, the vector y is only reset without performing the actual 
     * multiplication.
     *
     * @param x the vector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication and scaling
     * @param factor optional scaling factor, default to 1.0
     */
    void apply(const double *x, double *y, double factor = 1.0) const;

    /** @brief compute y = y + a * Ax 
     *
     * add the matrix-vector product "Ax", scaled by "a", to y.
     * "A" is this matrix.
     * 
     * This function can be used to compute the residual r = b - Ax.
     *
     * @param x the vector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication and scaling
     * @param factor optional scaling   factor, default to 1.0
     */
    void apply_scaled_add(const double *x, double *y, double factor = 1.0)
      const;
};

#endif
