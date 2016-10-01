#ifndef _TNSE3D_Routines_
#define _TNSE3D_Routines_

// ======================================================================
// @(#)TNSE3D_Routines.h        1.5 11/24/99
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

// ======================================================================
// compute turbulent viscosity for LES 
// ======================================================================
double TurbulentViscosity3D(double hK, double* gradU, double* u, 
			    double* uConv, double* x, double* y, double* z,
                            double proj_space);

// ======================================================================
// compute stabilization for div--div term
// ======================================================================
double DivDivStab3D(double u1, double u2, double u3, double hK, double eps);  

// ======================================================================
/// compute the Frobenious norm tensor
// ======================================================================
double frobeniusNormTensor(double* u, double* gradu,double *uConv,  int proj_space = 0);
// ======================================================================
/// compute turbulence viscosity using Smagorinsky model
// ======================================================================
double turbulentViscosity3D(double hK, double* u, double* gradu, double* uConv, 
                            double* x, double* y, double* z);
#endif