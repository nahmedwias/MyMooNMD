#ifndef ASSEMBLE_ROUTINE_TNSE2D_PBVMS_H_
#define ASSEMBLE_ROUTINE_TNSE2D_PBVMS_H_

//===============================================================================
// everything for the full Projection-Based VMS method
//===============================================================================
void TimeNSParamsVelo_GradVelo_LargeScale2D(double *in, double *out);

void TimeNSType4VMS_ProjectionDD2D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

void TimeNSType3_4NLVMS_ProjectionDD2D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);


#endif /* ASSEMBLE_ROUTINE_TNSE2D_PBVMS_H_ */
