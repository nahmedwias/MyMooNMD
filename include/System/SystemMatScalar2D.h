/** ************************************************************************ 
*
* @class     TSystemMatScalar2D
* @brief     stores the information of a 2D scalar system matrix 
* @author    Sashikumaar Ganesan, 
* @date      08.08.14
* @History   Added methods (Sashi, 22.08.14)
 ************************************************************************  */


#ifndef __SYSTEMMATSCALAR2D__
#define __SYSTEMMATSCALAR2D__

#include <SystemMat2D.h>
#include <LocalAssembling2D.h>

/**class for 2D scalar system matrix */
class TSystemMatScalar2D : public SystemMat2D
{
  protected:
    
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[1];

     /** Boundary value */   
    BoundValueFunct2D *BoundaryValues[1];
    
  public:
    /** constructor */
     TSystemMatScalar2D(TFESpace2D *fespace);

    /** destrcutor */
    ~TSystemMatScalar2D();
    
    /** Initilize the discrete forms and the matrices */
    void Init(BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue);
 
    /** assemble the system matrix */
    //void Assemble(TAuxParam2D *aux, double *sol, double *rhs);
    void Assemble(LocalAssembling2D& la, double *sol, double *rhs);

    /** solve the system matrix */
    void Solve(double *sol, double *rhs);
};

#endif // __SYSTEMMATSCALAR2D__
