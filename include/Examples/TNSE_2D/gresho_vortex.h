#ifndef GRESHO_VORTEX_H
#define GRESHO_VORTEX_H

#include <math.h>

double DIMENSIONLESS_VISCOSITY;


void ExampleFile()
{
  Output::print("Example: TNSE_2D/gresho_vortex.h");
}

// exact solution 
void ExactU1(double x, double y, double* values)
{
  const double r = sqrt(x*x + y*y);
  
  if((fabs(r) >= 0) &&  (fabs(r) <= 0.2))
    values[0] = -5.*y;
  else if((fabs(r) >= 0.2) &&  (fabs(r) <= 0.4))
    values[0] = -2.*y/r + 5.*y;
  else
    values[0] = 0.;
}
void ExactU2(double x, double y, double* values)
{
  const double r = sqrt(x*x + y*y);
  
  if((fabs(r) >= 0) &&  (fabs(r) <= 0.2))
    values[0] = 5.*x;
 else if((fabs(r) >= 0.2) &&  (fabs(r) <= 0.4))
    values[0] = 2.*x/r - 5.*x;
  else
    values[0] = 0.;
}
void ExactP(double x, double y, double* values)
{
  values[0]=0.;
  values[1]=0.;
  values[2]=0.;
  values[3]=0.;
}
void InitialU1(double x, double y, double* values)
{
  const double r = sqrt(x*x + y*y);
  
  if(r <= 0.2)
    values[0] = -5.*y;
  else if(r <= 0.4)
    values[0] = -2.*y/r + 5.*y;
  else
    values[0] = 0.;
}

void InitialU2(double x, double y, double* values)
{
  const double r = sqrt(x*x + y*y);
  
  if(r <= 0.2)
    values[0] = 5.*x;
 else if(r <= 0.4)
    values[0] = 2.*x/r - 5.*x;
  else
    values[0] = 0.;
}

void InitialP(double x, double y, double* values)
{
  values[0] = 0.;
}

void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  value=0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
 value=0;
}

void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double init_u1[1], init_u2[1];
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = DIMENSIONLESS_VISCOSITY;
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
    
    InitialU1(X[i], Y[i], init_u1);
    InitialU2(X[i], Y[i], init_u2);
    coeffs[i][3] = init_u1[0];
    coeffs[i][4] = init_u2[0];
  }
}
#endif // GRESHO_VORTEX_H