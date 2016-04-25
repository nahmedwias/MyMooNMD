// ======================================================================
// @(#)Brinkman2D_Mixed.h        1.2 04/13/00
//
// common declaration for all Brinkman problems
// ======================================================================

#ifndef __Brinkman2DMixed__
#define __Brinkman2DMixed__

#include <Enumerations.h>

// ======================================================================
// Type 1, Standard Galerkin for Brinkman in [p div u] formulation
// ======================================================================
void BrinkmanType1Galerkin(double Mult, double *coeff,
                     double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);
// ======================================================================

// ======================================================================
// Type 1.2, Standard Galerkin for Brinkman in [p div u] formulation with stabilization (PSPG) for P1/P1
// ======================================================================
void BrinkmanType1GalerkinStab(double Mult, double *coeff,
                           double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);
// ======================================================================

// ======================================================================
// Type 1.3, Standard Galerkin for Brinkman in [p div u] formulation with stabilization (PSPG) for P2/P2
// ======================================================================
void BrinkmanType1GalerkinStab2(double Mult, double *coeff,
                               double *param, double hK,
                               double **OrigValues, int *N_BaseFuncts,
                               double ***LocMatrices, double **LocRhs);
// ======================================================================

// ======================================================================
// Type 2, Standard Galerkin for Brinkman in [u gradp] formulation
// ======================================================================
void BrinkmanType2Galerkin(double Mult, double *coeff,
                           double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);
// ======================================================================
#endif
