// ======================================================================
// @(#)Darcy2DMixed.h        15.03.2015
//
// common declaration for all Darcy problems
// ======================================================================

#ifndef __DARCY_MIXED_2D__ 
#define __DARCY_MIXED_2D__

/** the local assembling routines. Each of them corresponds to one 
 * LocalAssembling2D_type */
int GetSignOfThisDOF(int N_DOF, int DOF);
// ======================================================================
// Standard Galerkin with Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM)
// elements
void BilinearAssembleDarcyGalerkin(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues,
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs);

#endif // __DARCY_MIXED_2D__
