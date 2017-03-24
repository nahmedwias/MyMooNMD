#ifndef COMMONROUTINETNSE2D_H
#define COMMONROUTINETNSE2D_H

// ========================================================================
void stabilization_parameters_equal_order(double Mult, double* u, double* coeff, 
                                          double* params);

double frobeniusNormTensor(double *gradu);

double turbulentviscosity(double hK, double* gradU, double* u, double* uConv);

#endif // COMMONROUTINETNSE2D_H
