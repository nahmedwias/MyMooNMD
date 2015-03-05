// ======================================================================
// @(#)ConvDiff3D_Routines.h        2007/05/09
//
// common routines for all convection diffusion problems
// ======================================================================

#ifndef __CONVDIFF3D_ROUTINES__
#define __CONVDIFF3D_ROUTINES__

double Compute_SDFEM_delta(double hK, double eps, double b1, double b2, double b3,
double react, double linfb);

double Compute_SOLD_sigma(double hK, double eps, double b1,
			  double b2, double b3, double c, double f,
			  double linfb, double deltaK, double *param,
			  double residual, int residual_computed, 
			  int time_dependent_problem);

#endif
