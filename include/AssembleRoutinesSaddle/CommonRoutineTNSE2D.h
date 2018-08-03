#ifndef COMMONROUTINETNSE2D_H
#define COMMONROUTINETNSE2D_H

double turbulentviscosity(double hK, double* gradU, double* u, double* uConv);


// stabilization parameaters for equal order finite elements. 
void stabilization_parameters_equal_order(double Mult, double* u, double* coeff,
                                          double* params);
// ========================================================================
// parameter routines
// ========================================================================


// ========================================================================
// parameters: u1old, u2old
// ========================================================================

void TimeNSParams2(double *in, double *out);

// ========================================================================
// parameters: u1old, u2old, combination of previous time solutions
// ========================================================================
void TimeNSType14ParamsSUPG(double *in, double *out);

// ========================================================================
// coletti, without g_\delta \ast u
// ========================================================================
void TimeNSParamsVelo_GradVelo(double *in, double *out);



#endif // COMMONROUTINETNSE2D_H
