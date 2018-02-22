#ifndef TNSE3DPROJBASEDNVMS_H
#define TNSE3DPROJBASEDNVMS_H


void TimeNSType4VMS_ProjectionDD3D(double Mult, double *coeff, double *param, double hK, 
     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

void TimeNSType3_4NLVMS_ProjectionDD3D(double Mult, double *coeff, double *param, double hK, 
     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

// ========================================================================
// parameters: u, grad u, G^H
// ========================================================================
void TimeNSParamsVelo_GradVelo_LargeScale3D(double *in, double *out);

#endif // TNSE3DPROJBASEDNVMS_H
