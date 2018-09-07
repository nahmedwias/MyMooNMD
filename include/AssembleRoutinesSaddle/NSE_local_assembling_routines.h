#ifndef INCLUDE_ASSEMBLEROUTINESSADDLE_NSE_LOCAL_ASSEMBLING_ROUTINES_H
#define INCLUDE_ASSEMBLEROUTINESSADDLE_NSE_LOCAL_ASSEMBLING_ROUTINES_H

///////////////////////////////////////////////////////////////////////////////
// standard terms (needed for Galerkin)

template <int d>
void NSLaplaceGradGradSingle(double Mult, double *coeff, double *param,
                             double hK, double**OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs);
template <int d>
void NSLaplaceGradGrad(double Mult, double *coeff, double *param,
                       double hK, double**OrigValues, int *N_BaseFuncts,
                       double ***LocMatrices, double **LocRhs);
template <int d>
void NSLaplaceDeformation(double Mult, double *coeff, double *param, double hK,
                          double**OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices,  double **LocRhs);
template <int d>
void NSDivergenceBlocks(double Mult, double *coeff, double *param, double hK,
                        double**OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                        double **LocRhs);
template <int d>
void NSGradientBlocks(double Mult, double *coeff, double *param, double hK,
                      double**OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs);
template <int d>
void NSRightHandSide(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);

template <int d>
void NSNonlinearTermSingle(double Mult, double *coeff, double *param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);
template <int d>
void NSNonlinearTerm(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs);

///////////////////////////////////////////////////////////////////////////////
// more special terms

template <int d>
void NSCoriolis(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                double **LocRhs);

///////////////////////////////////////////////////////////////////////////////
// stabilization terms

// PSPG (Pressure stabilization Petrov-Galerkin)
double compute_PSPG_delta(double hK, double nu);
template <int d>
void NSPSPG(double Mult, double *coeff, double *param, double hK,
            double **OrigValues, int *N_BaseFuncts,
            double ***LocMatrices, double **LocRhs);
template <int d>
void NSPSPG_RightHandSide(double Mult, double *coeff, double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs);

// symmetric GLS (Galerkin least-squares) method
double compute_GLS_delta(double hK, double nu);
template <int d>
void NSsymmGLS(double Mult, double *coeff, double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);
template <int d>
void NSsymmGLS_RightHandSide(double Mult, double *coeff, double *param,
                             double hK, double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs);

// non-symmetric GLS (Galerkin least-squares) method
template <int d>
void NSnonsymmGLS(double Mult, double *coeff, double *param, double hK,
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs);
template <int d>
void NSnonsymmGLS_RightHandSide(double Mult, double *coeff, double *param,
                                double hK, double **OrigValues,
                                int *N_BaseFuncts, double ***LocMatrices,
                                double **LocRhs);

///////////////////////////////////////////////////////////////////////////////
// routines to pass parameters (values of fe functions) to local assembling 
// routines

template <int d>
void NSParamsVelocity(double *in, double *out);


#endif // INCLUDE_ASSEMBLEROUTINESSADDLE_NSE_LOCAL_ASSEMBLING_ROUTINES_H
