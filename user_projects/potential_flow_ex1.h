#ifndef POTENTIAL_FLOW_EX3_H
#define POTENTIAL_FLOW_EX3_H

#include <math.h>

double DIMENSIONLESS_VISCOSITY;
void ExampleFile()
{
  Output::print("Example: potential_flow_ex1.h");  
}
// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    values[0] = 0; values[1] = 0; values[2] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;values[1] = 0; values[2] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = (-157*sqrt(0.00017886802095311101))/12. +
              2*sqrt(0.04024530471444998)*pow(x,3)*y +
              2*sqrt(0.04024530471444998)*pow(x,2)*pow(y,2) +
              2*sqrt(0.04024530471444998)*pow(y,4);

  values[1] = 6*sqrt(0.04024530471444998)*pow(x,2)*y +
              4*sqrt(0.04024530471444998)*x*pow(y,2);

  values[2] = 2*sqrt(0.04024530471444998)*pow(x,3) +
              4*sqrt(0.04024530471444998)*pow(x,2)*y +
              8*sqrt(0.04024530471444998)*pow(y,3);
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}


void U1BoundValue(int BdComp, double Param, double &value)
{
  value=0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
 value=0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = DIMENSIONLESS_VISCOSITY;
  double valp[4];

  for(int i=0;i<n_points;i++)
  {
    ExactP(X[i], Y[i], valp);
    coeffs[i][0] = nu;
    coeffs[i][1] = valp[1];
    coeffs[i][2] = valp[2];
  }
}

#endif // POTENTIAL_FLOW_EX3_H
