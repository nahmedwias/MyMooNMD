#ifndef TNSE3DRESBASEDVMS
#define TNSE3DRESBASEDVMS
//===============================================================================
// everything for the full residual-based vms method: 
// Either only in the test function previous time solution will be used
// or also in the ansatz (IMEX) both variants are implemented here
//===============================================================================
void TimeNSType4Residual_VMSDD3D(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType4NLResidual_VMSDD3D(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType4RHS_Residual_VMS(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);

//===============================================================================
void TimeNSParams_Type4Residual_VMS3D(double* in, double* out);
//===============================================================================

void TimeNSType14L_NL_Residual_VMS3D(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType14RHS_Residual_VMS3D(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSParams_Type14Residual_VMS3D(double* in, double* out);
#endif // TNSE3DRESBASEDVMS
