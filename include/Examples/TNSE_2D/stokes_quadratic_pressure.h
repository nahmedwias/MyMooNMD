// transient Stokes: potential flow example
// U = grad phi, phi = x*x - y*y

void ExampleFile()
{
  Output::print<1>("Example: stokes_quadratic_pressure.h") ;
}

#include<iostream>
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 2.*t*x;
}

void InitialU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -2.*t*y;
}

void InitialP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = (x*x - y*y);
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = 2.*t*x;
  values[1] = 2.*t;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -2.*t*y;
  values[1] =  0;
  values[2] = -2*t;
  values[3] =  0;
}

void ExactP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] =  (x*x - y*y);
  values[1] =  2.*x;
  values[2] = -2.*y;
  values[3] =  0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double x, y;
  switch(BdComp)
  {
    case 0:
      x = Param; y = 0.;
      break;
    case 1: 
      x = 1.; y = Param;
      break;
    case 2:
      x = 1. - Param; y = 1.;
      break;
    case 3:
      x = 0; y = 1. - Param;
      break;
    default: 
      ErrThrow("wrong boundary part number ", BdComp);
      break;
  }
  double v[4];
  ExactU1(x,y,v);
  value = v[0];
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double x, y;
  switch(BdComp)
  {
    case 0:
      x = Param; y = 0.;
      break;
    case 1: 
      x = 1.; y = Param;
      break;
    case 2:
      x = 1. - Param; y = 1.;
      break;
    case 3:
      x = 0; y = 1. - Param;
      break;
    default: 
      ErrThrow("wrong boundary part number ", BdComp);
      break;
  }
  double v[4];
  ExactU2(x,y,v);
  value = v[0];
  
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = 1/TDatabase::ParamDB->RE_NR;
  double u1[4], u2[4], p[4];
  
  for(int i=0;i<n_points;i++)
  {
    ExactU1(X[i], Y[i], u1);
    ExactU2(X[i], Y[i], u2);
    ExactP(X[i], Y[i], p);
    
    coeffs[i][0] = nu;    
    // Stokes
    coeffs[i][1] =  2*X[i] - nu*u1[3] + p[1];
    coeffs[i][2] = -2*Y[i] - nu*u2[3] + p[2];
  }
}


