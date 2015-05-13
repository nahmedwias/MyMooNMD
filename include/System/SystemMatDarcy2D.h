/** ************************************************************************ 
*
* @class     TSystemMatDarcy2D
* @brief     stores the information of a 2D Darcy system matrix 
* @author    Ulrich Wilbrandt
* @date      15.03.15
 ************************************************************************  */


#ifndef __SYSTEMMATDARCY2D__
#define __SYSTEMMATDARCY2D__

#include <LocalAssembling2D.h>
#include <SystemMat2D.h>

/**class for 2D scalar system matrix */
class TSystemMatDarcy2D : public SystemMat2D
{
  protected:
    
    /** Boundary conditon (one for u.n and one for pressure) */
    BoundCondFunct2D *BoundaryConditions[2];
    
    /** Boundary value */ 
    BoundValueFunct2D *BoundaryValues[2];
    
  public:
    /** constructor */
     TSystemMatDarcy2D(TFESpace2D *velocity, TFESpace2D* pressure,
                       BoundValueFunct2D **BoundValue);
    
    /** destrcutor */
    ~TSystemMatDarcy2D();
    
    /** assemble the system matrix */
    void Assemble(LocalAssembling2D& la, double *sol, double *rhs);
    
    /** solve the system matrix */
    void Solve(double *sol, double *rhs);
};

#endif // __SYSTEMMATDARCY2D__
