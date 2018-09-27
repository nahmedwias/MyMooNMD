#ifndef __DARCY_MIXED_2D__ 
#define __DARCY_MIXED_2D__

/** the local assembling routines. Each of them corresponds to one 
 * LocalAssembling2D_type */

// ======================================================================
// Standard Galerkin with Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM)
// elements
template <int d>
void BilinearAssembleDarcyGalerkin(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues,
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs);

#endif // __DARCY_MIXED_2D__
