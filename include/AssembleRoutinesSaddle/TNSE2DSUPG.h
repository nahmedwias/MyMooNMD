#ifndef TNSE2DSUPG_H
#define TNSE2DSUPG_H

//==============================================================================
//==============================================================================
void TimeNSType4SUPG(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType4NLSUPG(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType4RHSSUPG(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);
//===============================================================================
//===============================================================================
void TimeNSType14SUPG(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType14NLSUPG(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType14RHSSUPG(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);

//===============================================================================
#endif // TNSE2DSUPG_H
