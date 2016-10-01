/** ************************************************************************
 *
 * @brief     common declaration for all Brinkman problems
 *
 *
 * @author    Alfonso Caiazzo & Laura Blank
 * @date      18.05.16
 ************************************************************************  */

// TO DO:
// 1) PSPGSTAB mit Nitsche gekoppelt definieren (boundary integrals aus Stabilisierung nicht vergessen)
// 2) Formulation using the deformation tensor

#ifndef __Brinkman2DMixed__
#define __Brinkman2DMixed__

#include <Enumerations.h>

// ======================================================================
// Type 1, Standard Galerkin for Brinkman in [p div v] formulation
// ======================================================================
/** @brief brinkman weak formulation
 effective_viscosity (gradu,gradv)-(p,div v)+\frac{viscosity}{permeability}(u,v)=(f,v)
 (q, div u)=(g,q)
 */
void BrinkmanType1Galerkin(double Mult, double *coeff,
                           double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);



// ======================================================================
// Type 1b, Standard Galerkin for Scaled Brinkman in [p div v] formulation
// ======================================================================
/** @brief brinkman weak formulation
 t^2(gradu,gradv)-(p,div v)+(u,v)=(f,v)
 (div u, q)=(g,q)
 */
void BrinkmanType1bGalerkin(double Mult, double *coeff,
                            double *param, double hK,
                            double **OrigValues, int *N_BaseFuncts,
                            double ***LocMatrices, double **LocRhs);



// ======================================================================
// // Type 2, Standard Galerkin for Brinkman in [u gradp] formulation (macht nicht wirklich Sinn!!!!!!)
// ======================================================================
/** @brief brinkman weak formulation
 t^2(gradu,gradv)-(grad p,v)+(u,v)=(f,v)
 (u, grad q)=(g,q)
 */
void BrinkmanType2Galerkin(double Mult, double *coeff,
                           double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);



// ======================================================================
// Type 1.2, Standard Galerkin for Scaled Brinkman in [p div v] formulation with stabilization (PSPG) for P1/P1
// ======================================================================
/** @brief brinkman weak formulation
 t^2 (gradu,gradv)-(p,div v)+(u,v)=(f,v)
 (q, div u)=(g,q)
 PSPGStab \frac{h_k^2}{t^2+h_k^2}(-t^2 \Delta u + \nabla p + u,-t^2 \Delta v + \nabla q + v ) = PSPGStab \frac{h_k^2}{t^2+h_k^2} (f,-t^2 \Delta v + \nabla q + v )
 */
// for P1/P1 Delta u=0
void BrinkmanType1GalerkinResidualStab(double Mult, double *coeff,
                                       double *param, double hK,
                                       double **OrigValues, int *N_BaseFuncts,
                                       double ***LocMatrices, double **LocRhs);



// ======================================================================
// Type 1.3, Standard Galerkin for Scaled Brinkman in [p div v] formulation with stabilization (PSPG) for P2/P2
// ======================================================================
/** @brief brinkman weak formulation
 t^2(gradu,gradv)-(p,div v)+(u,v)=(f,v)
 (q, div u)=(g,q)
 PSPGStab \frac{h_k^2}{t^2+h_k^2}(-t^2 \Delta u + \nabla p + u,-t^2 \Delta v + \nabla q + v ) = PSPGStab \frac{h_k^2}{t^2+h_k^2} (f,-t^2 \Delta v + \nabla q + v )
 */
void BrinkmanType1GalerkinResidualStab2(double Mult, double *coeff,
                                        double *param, double hK,
                                        double **OrigValues, int *N_BaseFuncts,
                                        double ***LocMatrices, double **LocRhs);




#endif
