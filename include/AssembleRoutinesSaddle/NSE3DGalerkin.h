#ifndef NSE3DGALERKIN_H
#define NSE3DGALERKIN_H

//===========================================================================================
// Local assembling routines for the Navier-Stokes equations:
// All nstypes
//===========================================================================================
void NSLaplaceGradGradSingle(double Mult, double *coeff, double *param,
                             double hK, double**OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs);
void NSLaplaceGradGrad(double Mult, double *coeff, double *param,
                       double hK, double**OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs);
void NSLaplaceDeformation(double Mult, double *coeff, double *param, double hK,
                          double**OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices,  double **LocRhs);
void NSDivergenceBlocks(double Mult, double *coeff, double *param, double hK,
                        double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                        double **LocRhs);
void NSGradientBlocks(double Mult, double *coeff, double *param, double hK,
                      double**OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs);
void NSRightHandSide(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);

void NSNonlinearTermSingle(double Mult, double *coeff, double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);
void NSNonlinearTerm(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);

// delete all down except Ulrich


void NSType1Galerkin(double Mult, double *coeff, double *param, double hK, 
                     double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs);

void NSType2Galerkin(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs);

void NSType1_2NLGalerkin(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs);
//===========================================================================================
void NSType3Galerkin(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs);

void NSType4Galerkin(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                       double **LocRhs);


void NSType3_4NLGalerkin(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs);
//===========================================================================================
void NSType3GalerkinDD(double Mult, double *coeff, double *param, double hK, 
                       double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                       double **LocRhs);

void NSType4GalerkinDD(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                     double **LocRhs);

void NSType3_4NLGalerkinDD(double Mult, double *coeff, double *param, double hK,
                         double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs);

// ========================================================================
void NSParamsVelo3D(double *in, double *out);


// check this with Ulrich
void mat_coriolis(double Mult, double *coeff, double *param, double hK, 
                  double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                  double **LocRhs);

#endif // NSE3DGALERKIN_H
