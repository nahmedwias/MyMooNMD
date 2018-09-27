#ifndef NSE2DSUPG_H
#define NSE2DSUPG_H

#include <Enumerations.h>

//===========================================================================================
// Local assembling routines for the Navier-Stokes equations:
// All nstypes
//===========================================================================================
void NSType2SUPG(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs);

void NSType2NLSUPG(double Mult, double *coeff, double *param, double hK,
                   double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                   double **LocRhs);
//===========================================================================================
void NSType4SUPG(double Mult, double *coeff, double *param, double hK, 
                 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                 double **LocRhs);
void NSType4SUPGDD(double Mult, double *coeff, double *param, double hK, 
                   double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                   double **LocRhs);

void NSType4NLSUPG(double Mult, double *coeff, double *param, double hK,
                   double **OrigValues, int *N_BaseFuncts, 
                   double ***LocMatrices, double **LocRhs);

void NSType4NLSUPGDD(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);



double SUPG_Parameter(double hK, double eps, double b1, double b2, double c);

#endif // NSE2DSUPG_H
