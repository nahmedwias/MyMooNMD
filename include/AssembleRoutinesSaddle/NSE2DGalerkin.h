#ifndef NSE2DGALERKIN_H
#define NSE2DGALERKIN_H

#include <Enumerations.h>

//===========================================================================================
// Local assembling routines for the Navier-Stokes equations:
// All nstypes
//===========================================================================================
void NSType1Galerkin(double Mult, double *coeff, double *param, double hK, 
                     double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs);

void NSType2Galerkin(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs);

void NSType1_2NLGalerkin(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs);
//===========================================================================================
void NSType3Galerkin(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs);

void NSType4Galerkin(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                       double **LocRhs);


void NSType3_4NLGalerkin(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs);
//===========================================================================================
void NSType3GalerkinDD(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                       double **LocRhs);

void NSType4GalerkinDD(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs);

void NSType3_4NLGalerkinDD(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs);

#endif // NSE2DGALERKIN_H
