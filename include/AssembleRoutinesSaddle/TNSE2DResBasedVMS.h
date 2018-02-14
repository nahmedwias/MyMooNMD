#ifndef TNSE2DRESBASEDVMS
#define TNSE2DRESBASEDVMS
//===============================================================================
// everything for the full residual-based vms method: 
// Either only in the test function previous time solution will be used
// or also in the ansatz (IMEX) both variants are implemented here
//===============================================================================
void TimeNSType4Residual_VMS(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType4NLResidual_VMS(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType4RHS_Residual_VMS(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);

//===============================================================================
void TimeNSParams_Residual_VMS(double* in, double* out);
//===============================================================================
void TimeNSType14Residual_VMS(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType14NLResidual_VMS(double Mult, double *coeff, double *param, double hK, 
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType14RHS_Residual_VMS(double Mult, double *coeff, 
  double *param, double hK, double **OrigValues, int *N_BaseFuncts, 
  double ***LocMatrices, double **LocRhs);


#endif // TNSE2DRESBASEDVMS
