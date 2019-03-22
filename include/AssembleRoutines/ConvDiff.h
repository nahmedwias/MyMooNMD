// ======================================================================
// @(#)ConvDiff.h        12/06/26
//
// common declaration for all convection diffusion problems
// ======================================================================

#ifndef __CONVDIFF__
#define __CONVDIFF__

#include <Collection.h>
#include "Constants.h"

template <int d>
double Mesh_size_in_convection_direction_without_storing(
  double hK, std::array<double, d> b);

template <int d>
double Mesh_size_in_convection_direction(double hK, std::array<double, d> b);

template <int d>
double Compute_SDFEM_delta(double hK, double eps, std::array<double, d> b,
                           double react, double linfb);

template <int d>
double Compute_SOLD_sigma(double hK, double eps, std::array<double, d> b,
                          double c, double f, double linfb, double deltaK,
                          double *param, double residual, int residual_computed,
                          int time_dependent_problem);

/** coercivity constant of problem 
 * used eg for residual based estimator of Verf"uhrt 2005 */
double EstimateCoercivityConstant(TCollection *Coll, 
#ifdef __2D__
                                  const CoeffFct2D& Coeff
#else // 3D
                                  const CoeffFct3D& Coeff
#endif
                                 );

/** it should be i>=0 and i<= 31, otherwise error*/
void SetSoldParameters(int i);


#ifdef __2D__
void EdgeStabilization(TFESpace2D *fespace,  TFEFunction2D *u, 
                       const CoeffFct2D& Coeffs, double *rhs,
                       int time_dependent,
                       double *time_step, TFEFunction2D *old_u);
#endif

double ComputeAlpha(double hK);

//=============================================================================
// local assembling routines:

template<int d>
void BilinearAssembleGalerkin(double Mult, double *coeff, double* param, 
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);
// SUPG/SDFEM
// SDFEM - Streamline Diffusion Finite Element Method,
// SUPG - Streamline Upwind Petrov Galerkin
template<int d>
void BilinearAssemble_SD(double Mult, double *coeff, double* param,
                         double hK, double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs);
// Galerking least squares
template<int d>
void BilinearAssemble_GLS(double Mult, double *coeff, double* param,
                         double hK, double **OrigValues, int *N_BaseFuncts,
                         double ***LocMatrices, double **LocRhs);

//=============================================================================
// local error computing routines:
template<int d>
void conv_diff_l2_h1_linf_error(int N_Points, std::array<double*, d> xyz,
                                double *AbsDetjk, const double *Weights,
                                double hK, double **Der, double **Exact,
                                double **coeffs, double *LocError);

#endif // __CONVDIFF__
