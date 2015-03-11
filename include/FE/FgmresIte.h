// =======================================================================
// @(#)ItMethod.h        1.6 10/18/99
// 
// Class:       TFgmresIte
// Purpose:     defines the fixed point iteration
//
// Author:      Volker John
//
// History:     start of implementation 24.10.2000
//
// =======================================================================
#ifndef __FGMRESITE__
#define __FGMRESITE__

#include <ItMethod.h>
// #include <ParFECommunicator2D.h>
#include <ParFECommunicator3D.h>

/** iteration method */
class TFgmresIte : public TItMethod
{
  protected:

  /** arrays for flexible gmres depending on restart */
  double *s;
  double *cosi;
  double *sn;

  /** arrays for flexible gmres depending on number of dof */
  double *d;
  double *z;
  
  /** matrices for flexible gmres depending on restart */
  double **H;
  double **v;
  double **zv;

#ifdef _MPI   
   #ifdef __3D__
      TParFECommunicator3D *ParComm;
   #else
      TParFECommunicator2D *ParComm;
   #endif
#endif  
      
  public:
    /** constructor */
    TFgmresIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
               int n_aux, int N_Unknowns, int scalar);
 
#ifdef _MPI
    TFgmresIte(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
               int n_aux, int N_Unknowns, int scalar,
  #ifdef  __3D__			       
			       TParFECommunicator3D *parcomm
  #else			       
                               TParFECommunicator2D *parcomm
  #endif
	      );
#endif
    /** destructor */
    virtual ~TFgmresIte();
    
    /** iterate routine */
    int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
                double *rhs);    
};
#endif
