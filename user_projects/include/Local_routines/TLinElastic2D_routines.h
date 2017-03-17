/*
 * TLinElastic2D_routines.h
 *
 * common declaration for routines for TLinElastic2D
 *
 *  Created on: Mar 17, 2017
 *      Author: alia
 */

#ifndef __TLINELASTIC2D_ROUTINES__
#define __TLINELASTIC2D_ROUTINES__



// ======================================================================
// Stiffness, mass matrix and rhs
// ======================================================================
void TimeLinearElasticityWholeSystem(double Mult, double *coeff,
                                     double *param, double hK,
                                     double **OrigValues, int *N_BaseFuncts,
                                     double ***LocMatrices, double **LocRhs);


// ======================================================================
// Stiffness matrix
// ======================================================================
void TimeLinearElasticityStiffness(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs);

// ======================================================================
// Mass matrix
// ======================================================================
void TimeLinearElasticityMass(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs);

// ======================================================================
// Rhs
// ======================================================================
void TimeLinearElasticityRhs(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs);


// ======================================================================
// Parameter Fct
// ======================================================================
void ParameterFunction(double *in, double *out);


#endif /* __TLINELASTIC2D_ROUTINES__ */
