#ifndef INCLUDE_ASSEMBLEROUTINESSADDLE_TIME_NSE_LOCAL_ASSEMBLING_ROUTINES_H
#define INCLUDE_ASSEMBLEROUTINESSADDLE_TIME_NSE_LOCAL_ASSEMBLING_ROUTINES_H

/** Store only the local assembling routines that are additional from the 
 * time stepping schemes, e.g., the mass matrices, additional matrices in the 
 * case of stabilization schemes. Note that: we can still use the local assembly
 * of velocity-velocity coupling from the steady case but due 
 * to the ordering of matrices in the Assemble function, rectangular blocks 
 * will have different number's because of the additional mass matrices. Keeping 
 * this in mind, velocity-pressure coupling blocks are copied here with correct numbering. 
 */

template<int d>
void NSMassMatrixSingle(double Mult, double *coeff, double *, double,
  double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);

template<int d>
void NSMassMatrix(double Mult, double *, double *, double, 
  double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);

template<int d> 
void NSLaplaceGradGradSingleSmagorinsky(double Mult, double *coeff, double *param, double hK,
  double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);

template<int d> 
void NSLaplaceGradGradSmagorinsky(double Mult, double *coeff, double *param, double hK,
  double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);

template<int d> 
void NSLaplaceDeformationSmagorinsky(double Mult, double *coeff, double *param, double hK,
  double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);

template<int d> 
void NSLaplaceDeformationVariationalMS(double Mult, double *, double *param, double hK,
  double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);

template<int d>
void NSLumpMassMatrix(double Mult, double *, double *, double ,
  double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);

template<int d>
void NSVariationlMS_GMatrices(double Mult, double *, double *, double ,
  double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);

template<int d>
void NSVariationlMS_GTildeMatrices(double Mult, double *, double *param, double hK,
  double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **);

template<int d>
void NSParamVelGradSmagorinsky(const double *in, double *out);

template<int d>
void NSParamsVariationalMSLargeScale(const double *in, double *out);



//===============================================================================
// nonlinear terms
//===============================================================================

#endif
