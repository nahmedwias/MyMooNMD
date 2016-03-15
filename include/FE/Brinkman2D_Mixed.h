// ======================================================================
// @(#)Brinkman2D_Mixed.h        1.2 04/13/00
//
// common declaration for all Navier-Stokes problems
// ======================================================================

#ifndef __Brinkman2DMixed__
#define __Brinkman2DMixed__

#include <Enumerations.h>
// ======================================================================
// Type 4, Standard Galerkin for Brinkman
// ======================================================================
void BrinkmanType4Galerkin(double Mult, double *coeff,
                     double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);
// ======================================================================
#endif
