// ======================================================================
// @(#)ConvDiff2D_Routines.h        2006/02/22
//
// common routines for all convection diffusion problems
// ======================================================================

#ifndef __CONVDIFF2D_ROUTINES__
#define __CONVDIFF2D_ROUTINES__

double Mesh_size_in_convection_direction(double hK, double b1, double b2);

double Compute_SDFEM_delta(double hK, double eps, double b1, double b2,
double react, double linfb);

double Compute_SOLD_sigma(double hK, double eps, double b1,
			  double b2, double c, double f, 
			  double linfb, double deltaK, double *param,
			  double residual, int residual_computed,
			  int time_dependent_problem);

void EdgeStabilization(TFESpace2D *fespace, 
		       TFEFunction2D *u, 
		       CoeffFct2D *Coeffs,
		       double *rhs,
		       int time_dependent,
		       double *time_step,
		       TFEFunction2D *old_u);

/** coercivity constant of problem */
double EstimateCoercivityConstant(TCollection *Coll,
				  CoeffFct2D *Coeff);

void JumpTermsForAdjointProblemP1(TFESpace2D *fespace,
				TFEFunction2D *u,
				CoeffFct2D *Coeffs,
				BoundCondFunct2D *BoundaryConditions,
				double *rhs);

void JumpTermsForAdjointProblem(TFESpace2D *fespace,
				TFEFunction2D *u,
				CoeffFct2D *Coeffs,
				BoundCondFunct2D *BoundaryConditions,
				double *rhs);


#endif
