/** ************************************************************************
 *
 * @brief     common declaration for all Brinkman problems
 *
 *
 * @author    Alfonso Caiazzo & Laura Blank
 * @date      06.10.2016
 ************************************************************************  */

// TO DO:
// 1) PSPGSTAB mit Nitsche gekoppelt definieren (boundary integrals aus Stabilisierung nicht vergessen)
// 2) Formulation using the deformation tensor

#ifndef __Brinkman3DMixed__
#define __Brinkman3DMixed__

#include <Enumerations.h>

// ======================================================================
// Type 1, Standard Galerkin for Brinkman in [p div v] formulation
// ======================================================================
/** @brief brinkman weak formulation for NSTYPE 4
 effective_viscosity (gradu,gradv)-(p,div v)+\frac{viscosity}{permeability}(u,v)=(f,v)
 (q, div u)=(g,q)
 */
void Brinkman3DType1Galerkin(double Mult,
                             double *coeff,
                             double *param,
                             double hK,
                             double **OrigValues,
                             int *N_BaseFuncts,
                             double ***LocMatrices,
                             double **LocRhs);

// ======================================================================
// Type 2, Standard Galerkin for Brinkman in [p div v] formulation
// ======================================================================
/** @brief brinkman weak formulation for NSTYPE 4
 effective_viscosity (gradu,gradv)-(p,div v)+\frac{viscosity}{permeability}(u,v)=(f,v)
 (q, div u)=(g,q)
 */
void Brinkman3DType2Galerkin(double Mult,
                             double *coeff,
                             double *param,
                             double hK,
                             double **OrigValues,
                             int *N_BaseFuncts,
                             double ***LocMatrices,
                             double **LocRhs);


// ======================================================================
// Type 1.2, Standard Galerkin for Scaled Brinkman in [p div v] formulation with stabilization (PSPG) for P1/P1
// ======================================================================
/** @brief brinkman weak formulation
 t^2 (gradu,gradv)-(p,div v)+(u,v)=(f,v)
 (q, div u)=(g,q)
 PSPGStab \frac{h_k^2}{t^2+h_k^2}(-t^2 \Delta u + \nabla p + u,-t^2 \Delta v + \nabla q + v ) = PSPGStab \frac{h_k^2}{t^2+h_k^2} (f,-t^2 \Delta v + \nabla q + v )
 */
// for P1/P1 Delta u=0
void Brinkman3DType1GalerkinResidualStabP1(double Mult, double *coeff,
                                           double *param, double hK,
                                           double **OrigValues, int *N_BaseFuncts,
                                           double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1.2, Standard Galerkin for Scaled Brinkman in [p div v] formulation with stabilization (PSPG) for P2/P2
// ======================================================================
/** @brief brinkman weak formulation
 t^2 (gradu,gradv)-(p,div v)+(u,v)=(f,v)
 (q, div u)=(g,q)
 PSPGStab \frac{h_k^2}{t^2+h_k^2}(-t^2 \Delta u + \nabla p + u,-t^2 \Delta v + \nabla q + v ) = PSPGStab \frac{h_k^2}{t^2+h_k^2} (f,-t^2 \Delta v + \nabla q + v )
 */
void Brinkman3DType1GalerkinResidualStabP2(double Mult, double *coeff,
                                           double *param, double hK,
                                           double **OrigValues, int *N_BaseFuncts,
                                           double ***LocMatrices, double **LocRhs);








void ResidualStabPkPk_for_Brinkman3DType1Galerkin(double Mult, double *coeff,
                                                double *param, double hK,
                                                double **OrigValues, int *N_BaseFuncts,
                                                double ***LocMatrices, double **LocRhs);



void GradDivStab_for_Brinkman3DType1Galerkin(double Mult, double *coeff,
                                             double *param, double hK,
                                             double **OrigValues, int *N_BaseFuncts,
                                             double ***LocMatrices, double **LocRhs);


//// ======================================================================
//// // Type 2, Standard Galerkin for Brinkman in [u gradp] formulation (macht nicht wirklich Sinn!!!!!!)
//// ======================================================================
///** @brief brinkman weak formulation
// t^2(gradu,gradv)-(grad p,v)+(u,v)=(f,v)
// (u, grad q)=(g,q)
// */
//void Brinkman3DType2Galerkin(double Mult, double *coeff,
//                           double *param, double hK,
//                           double **OrigValues, int *N_BaseFuncts,
//                           double ***LocMatrices, double **LocRhs);
//
//
//

//// ======================================================================
//// Type 1.3, Standard Galerkin for Scaled Brinkman in [p div v] formulation with stabilization (PSPG) for P2/P2
//// ======================================================================
///** @brief brinkman weak formulation
// t^2(gradu,gradv)-(p,div v)+(u,v)=(f,v)
// (q, div u)=(g,q)
// PSPGStab \frac{h_k^2}{t^2+h_k^2}(-t^2 \Delta u + \nabla p + u,-t^2 \Delta v + \nabla q + v ) = PSPGStab \frac{h_k^2}{t^2+h_k^2} (f,-t^2 \Delta v + \nabla q + v )
// */
//void Brinkman3DType1GalerkinResidualStabP2(double Mult, double *coeff,
//                                        double *param, double hK,
//                                        double **OrigValues, int *N_BaseFuncts,
//                                        double ***LocMatrices, double **LocRhs);
//
//


#endif
