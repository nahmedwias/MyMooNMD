#ifndef TNSE3DSUPG_H
#define TNSE3DSUPG_H

//==============================================================================
//==============================================================================
void TimeNSType4SUPGDD3D(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType4NLSUPGDD3D(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType4RHSSUPG3D(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);
//===============================================================================
void TimeNSType4Params_SUPG(double *in, double *out);
//===============================================================================
void TimeNSType14SUPGDD3D(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType14NLSUPGDD3D(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType14RHSSUPGDD3D(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);

//===============================================================================
void TimeNSType14Params_SUPG(double *in, double *out);
//===============================================================================
#endif // TNSE3DSUPG_H
