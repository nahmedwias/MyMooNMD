/** ************************************************************************ 
*
* @class     TSystemMatNSE2D
* @brief     stores the information of a 2D NSE system matrix 
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History   Added methods (Sashi, 23.08.14)
*            made this class a derived class fo SystemMat2D (Ulrich, 19.03.2015)
*            further simplifications (Ulrich, 25.03.2015)
* ************************************************************************  */


#ifndef __SYSTEMMATNSE2D__
#define __SYSTEMMATNSE2D__

#include <SystemMat2D.h>
#include <LocalAssembling2D.h>

/**class for 2D  NSE system matrix */
class TSystemMatNSE2D : public SystemMat2D
{
  protected:
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[2];
    
    /** Boundary values */
    BoundValueFunct2D *BoundaryValues[2];
    
  public:
    /** constructor */
     TSystemMatNSE2D(TFEVectFunct2D *Velocity, TFEFunction2D *p);
     
    /** destrcutor */
    ~TSystemMatNSE2D();
    
    /** methods */
    /** Initilize the discrete forms and the matrices */
    void Init(BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue);
    
    /** assemble the system matrix */
    void Assemble(LocalAssembling2D& la, double *sol, double *rhs);
    
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear(LocalAssembling2D& la, double *sol, double *rhs);
    
    /** get the resudual of the NSE system */
    void GetResidual(double *sol, double *rhs, double *res);
    
    /** solve the system matrix */
    void Solve(double *sol, double *rhs);
};

#endif // __SYSTEMMATNSE2D__
