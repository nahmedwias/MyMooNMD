#ifndef COMMONROUTINETNSE2D_H
#define COMMONROUTINETNSE2D_H

// ========================================================================
void stabilization_parameters_equal_order(double Mult, double* u, double* coeff, 
                                          double* params);

void stab_params(double hk, double *u, double nu, double K, double *param);
#endif // COMMONROUTINETNSE2D_H
