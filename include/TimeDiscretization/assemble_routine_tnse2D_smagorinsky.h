#ifndef ASSEMBLE_ROUTINE_TNSE2D_SMAGORINSKY_H
#define ASSEMBLE_ROUTINE_TNSE2D_SMAGORINSKY_H

void TimeNSType4SmagorinskyDD(double Mult, double *coeff, double *param, double hK,
double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType3_4NLSmagorinskyDD(double Mult, double *coeff, double *param, double hK,
double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSRHSSmagorinskyExplicit(double Mult, double *coeff, double *param, double hK,
double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

#endif // ASSEMBLE_ROUTINE_TNSE2D_SMAGORINSKY_H
