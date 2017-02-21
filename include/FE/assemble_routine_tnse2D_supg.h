#ifndef ASSEMBLE_ROUTINE_TNSE2D_SUPG_H
#define ASSEMBLE_ROUTINE_TNSE2D_SUPG_H

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

void TimeNSType14ParamsSUPG(double *in, double *out);
//===============================================================================
//===============================================================================
void TimeNSType4SUPGExtr(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);

void TimeNSType4NLSUPGExtr(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);

void TimeNSType4RHSSUPGExtr(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);

void TimeNSType4SUPGExtrParam(double *in, double *out);
//===============================================================================

//===============================================================================
void TimeNSType14SUPGExtr(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);

void TimeNSType14NLSUPGExtr(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);

void TimeNSType14RHSSUPGExtr(double Mult, double* coeff, double* param, double hK, 
 double** OrigValues, int* N_BaseFuncts, double*** LocMatrices, double** LocRhs);
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
#endif // ASSEMBLE_ROUTINE_TNSE2D_SUPG_H
