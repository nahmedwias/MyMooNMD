#ifndef INCLUDE_ASSEMBLEROUTINESSADDLE_NSE_LOCAL_ASSEMBLING_ROUTINES_H
#define INCLUDE_ASSEMBLEROUTINESSADDLE_NSE_LOCAL_ASSEMBLING_ROUTINES_H

///////////////////////////////////////////////////////////////////////////////
// standard terms, linear (needed for Galerkin)

template<int d>
void NSResistanceMassMatrixSingle(double Mult, double *coeff,
				  double *param,
				  double hK, double **OrigValues,
				  int *N_BaseFuncts,
				  double ***LocMatrices, double **LocRhs);

template<int d>
void NSResistanceMassMatrix(double Mult, double *coeff,
			    double *param,
			    double hK, double **OrigValues,
			    int *N_BaseFuncts,
			    double ***LocMatrices, double **LocRhs);

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
                        double **LocRhs, int sign=1);
template <int d>
void NSGradientBlocks(double Mult, double *coeff, double *param, double hK,
                      double**OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs);
template <int d>
void NSRightHandSide(double Mult, double *coeff, double *param, double hK,
                     double **OrigValues, int *N_BaseFuncts,
                     double ***LocMatrices, double **LocRhs, int sign = 1);

///////////////////////////////////////////////////////////////////////////////
// local assemblings for the nonlinear term
template <int d>
void NSNonlinearTerm_convective_Single(
  double Mult, double *coeff,double *param, double hK,double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template <int d>
void NSNonlinearTerm_convective(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template <int d>
void NSNonlinearTerm_skew_symmetric_Single(
  double Mult, double *coeff,double *param, double hK,double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template <int d>
void NSNonlinearTerm_skew_symmetric(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
template <int d>
void NSNonlinearTerm_rotational(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
// Energy, Momentum, and Angular momentum Conserving (EMAC)
template <int d>
void NSNonlinearTerm_emac(
  double Mult, double *coeff, double *param, double hK, double **OrigValues,
  int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);

///////////////////////////////////////////////////////////////////////////////
// more special terms

template <int d>
void NSCoriolis(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices,
                double **LocRhs);

///////////////////////////////////////////////////////////////////////////////
// stabilization terms

// Grad-Div
template <int d>
void NSGradDiv(double Mult, double *coeff, double *param, double hK,
            double **OrigValues, int *N_BaseFuncts,
            double ***LocMatrices, double **LocRhs, double delta0);
template <int d>
void NSGradDiv_RightHandSide(double Mult, double *coeff, double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs, double delta0);

// PSPG (Pressure stabilization Petrov-Galerkin)
double compute_PSPG_delta(double delta0, double hK, double nu);
template <int d>
void NSPSPG(double Mult, double *coeff, double *param, double hK,
            double **OrigValues, int *N_BaseFuncts,
            double ***LocMatrices, double **LocRhs, double delta0);
template <int d>
void NSPSPG_RightHandSide(double Mult, double *coeff, double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs, double delta0);

// symmetric GLS (Galerkin least-squares) method
double compute_GLS_delta(double delta0, double hK, double nu, double sigma=0.);

template <int d>
void NSsymmGLS(double Mult, double *coeff, double *param, double hK,
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs, double delta0);
template <int d>
void NSsymmGLS_RightHandSide(double Mult, double *coeff, double *param,
                             double hK, double **OrigValues, int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs,
                             double delta0);

// non-symmetric GLS (Galerkin least-squares) method
template <int d>
void NSnonsymmGLS(double Mult, double *coeff, double *param, double hK,
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs, double delta0);
template <int d>
void NSnonsymmGLS_RightHandSide(double Mult, double *coeff, double *param,
                                double hK, double **OrigValues,
                                int *N_BaseFuncts, double ***LocMatrices,
                                double **LocRhs, double delta0);

template <int d>
void NS_BrezziPitkaeranta(double Mult, double *coeff, double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs,
                          double delta0);

///////////////////////////////////////////////////////////////////////////////
// routines to pass parameters (values of fe functions) to local assembling 
// routines

template <int d>
void NSParamsVelocity(double *in, double *out);


#endif // INCLUDE_ASSEMBLEROUTINESSADDLE_NSE_LOCAL_ASSEMBLING_ROUTINES_H
