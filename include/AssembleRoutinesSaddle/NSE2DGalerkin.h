#ifndef NSE2DGALERKIN_H
#define NSE2DGALERKIN_H

#include <Enumerations.h>

//===========================================================================================
// Local assembling routines for the Navier-Stokes equations:
// All nstypes
//===========================================================================================
void NSLaplaceGradGradSingle(double Mult, double *coeff, double *param,
                             double hK, double**OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs);
void NSLaplaceGradGrad(double Mult, double *coeff, double *param,
                       double hK, double**OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs);
void NSLaplaceDeformation(double Mult, double *coeff, double *param, double hK,
                          double**OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices,  double **LocRhs);
void NSDivergenceBlocks(double Mult, double *coeff, double *param, double hK,
                        double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                        double **LocRhs);
void NSGradientBlocks(double Mult, double *coeff, double *param, double hK,
                      double**OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs);
void NSRightHandSide(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);

void NSNonlinearTermSingle(double Mult, double *coeff, double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);
void NSNonlinearTerm(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);

void NSParamsVelo(double *in, double *out);

#endif // NSE2DGALERKIN_H
