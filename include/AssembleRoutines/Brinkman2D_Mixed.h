/** ************************************************************************
 * @brief     Modules of Brinkman problems in 2D
 * @author    Alfonso Caiazzo & Laura Blank
 * @date      18.05.16
 ************************************************************************  */

// TO DO: Formulation using the deformation tensor

#ifndef __Brinkman2DMixed__
#define __Brinkman2DMixed__

#include <Enumerations.h>

// ======================================================================
// Standard Galerkin for Brinkman in (p, div v) formulation
// ======================================================================
/** @brief Weak formulation of the Brinkman problem:
  mueff * (grad u, grad v) - (p, div v) + sigma * (u, v) = (f, v)
  (q, div u) = (g, q)
 */
void BrinkmanType1Galerkin(double Mult, double *coeff,
    double *param, double hK,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs);

// ======================================================================
// Standard Galerkin for Brinkman in (grad p, v) formulation (Check if this function is meaningfull at all)
// ======================================================================
/** @brief Weak formulation of the Brinkman problem:
  mueff * (grad u, grad v) - (grad p, v) + sigma * (u, v) = (f, v)
  (u, grad q) = (g,q)
 */
void BrinkmanType2Galerkin(double Mult, double *coeff,
    double *param, double hK,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs);

// ======================================================================
// GLS-stabilization terms restricted to P1/P1 for the Brinkman in [p div v] formulation
// ======================================================================
/** @brief GLS method terms for P1/P1 for the Brinkman Problem:
  alpha * \frac{h_K^2}{mueff + sigma l_K^2} * (-mu_eff \Delta u + \nabla p + sigma u, -mueff \Delta v + \nabla q + sigma v ) = alpha * \frac{h_K^2}{mueff + sigma l_K^2} (f, -mueff \Delta v + \nabla q + sigma v )
 */
// for P1/P1 it is Delta u = 0
void BrinkmanType1GalerkinResidualStabP1(double Mult, double *coeff,
    double *param, double hK,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs);

// ======================================================================
// GLS-stabilization terms restricted to P2/P2 without P1/P1 terms for the Brinkman in [p div v] formulation↲
// ======================================================================
/** @brief GLS method terms for P2/P2 without P1/P1 terms for the Brinkman Problem:
  alpha * \frac{h_K^2}{mueff + sigma l_K^2} * [ (-mu_eff \Delta u, -mueff \Delta v + \nabla q + sigma v ) + ( \nabla p + sigma u, -mueff \Delta v) ] = alpha * \frac{h_K^2}{mueff + sigma l_K^2} (f, -mueff \Delta v)
 */
// for Pk/Pk, k>1, it is Delta u != 0
void BrinkmanType1GalerkinResidualStabP2(double Mult, double *coeff,
    double *param, double hK,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs);

// ======================================================================
// Grad-Div-Stabilization
// ======================================================================
/** @brief Grad-Div Stabilization terms:
  \delta * \frac{h_K^2}{mueff + sigma l_K^2} * (div u, div v) = \delta * \frac{h_K^2}{mueff + sigma l_K^2} * (g, div v)
 */
void BrinkmanGradDivStab(double Mult, double *coeff,
    double *param, double hK,
    double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs);


#endif
