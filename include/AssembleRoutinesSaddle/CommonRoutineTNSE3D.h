#ifndef COMMONROUTINETNSE3D_H
#define COMMONROUTINETNSE3D_H

// ======================================================================
// compute turbulent viscosity for LES 
// ======================================================================
double TurbulentViscosity3D(double hK, double* gradU, double* u, 
			    double* uConv, double* x, double* y, double* z,
                            double proj_space);
// ========================================================================
void stabilization_parameters_equal_order(double Mult, double* u, double* coeff, 
                                          double* params);

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
                            double* x, double* y, double* z, double proj_space);

#endif // COMMONROUTINETNSE3D_H
