// =======================================================================
// @(#)ItMethod.h        1.6 10/18/99
// 
// Class:       TItMethod
// Purpose:     defines an iteration method
//
// Author:      Volker John
//
// History:     start of implementation 24.10.2000
//
// =======================================================================

#ifndef __ITMETHOD__
#define __ITMETHOD__

#include <Constants.h>

/** iteration method */
class TItMethod
{
  protected:
    /** routine for matrix vector product */
    MatVecProc *matvec;

    /** routine for defect */
    DefectProc *matvecdefect;

    /** preconditioner */
    TItMethod *prec;

    /** number of degrees of freedom  */
    int N_DOF;

    /** id of the system type  */
    int system_id;

    /** maximal number of preconditioner iterations */
    int prec_maxit;

    /** absolute tolerance for stopping */
    double res_norm_min;

    /** relative tolerance for stopping */
    double red_factor;

    /** limit divergence */
    double div_factor;

    /** maximal number of iterations */
    int maxit;

    /** minimal number of iterations */
    int minit;

    /** restart */
    int restart;

    /** array for defect */
    double *defect;

    /** auxiliary arrays */
    double **AuxArray; 

    /** number of auxiliary arrays */
    int N_Aux;

  public:
    /** constructor */
    TItMethod(MatVecProc *MatVec, DefectProc *Defect, TItMethod *Prec,
              int n_aux, int N_Unknowns);

    /** destructor */
    virtual ~TItMethod();

    /** iterate */
    /** input **A, **B */
    /** output (solution) *sol */
    /** input (rhs) *rhs */
    virtual int Iterate(TSquareMatrix **A, TMatrix **B, double *sol, 
                        double *rhs) = 0;

    /** return system id */
    int GetSystemId()
      {return system_id;};

    /** return maximal number of preconditioner iterations */
    int GetPrecMaxit()
      {return prec_maxit;};

    /** return absolute tolerance for stopping */
    double GetResNormMin()
      {return res_norm_min;};

    /** return relative tolerance for stopping */
    double GetRedFactor()
      {return red_factor;};

    /** return tolerance for divergence */
    double GetDivFactor()
      {return div_factor;};

    /** return maximal number of iterations */
    int GetMaxit()
      {return maxit;};

    /** return minimal number of iterations */
    int GetMinit()
      {return minit;};

    /** return restart */
    int GetRestart()
      {return restart;};

};

#endif
