#ifndef COSINE_SINE_H
#define COSINE_SINE_H

#include <math.h>
// u1 = cos(t) * ( sin(\pi * x -0.7) *sin( \pi * y +0.2 ) );
// u2 = cos(t) * ( cos(\pi * x -0.7) *cos( \pi * y +0.2 ) );
// p =  cos(t) * ( sin(x) * cos(y) * (cos(1) * sin(1) - sin(1)) );
//
double DIMENSIONLESS_VISCOSITY;
void ExampleFile()
{
  Output::print("Example: TNSE_2D/cosine_sin.h") ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = cos(Pi*t)*( sin(Pi * x -0.7) *sin( Pi * y +0.2 ) );
}

void InitialU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0]= cos(Pi*t)*( cos(Pi * x -0.7) *cos( Pi * y +0.2 ) );
}

void InitialP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = cos(Pi*t)*( sin(Pi*x) * cos(Pi*y));
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double u  = cos(Pi*t)*( sin(Pi * x -0.7) *sin( Pi * y +0.2 ) );
  double ux = cos(Pi*t)*( Pi * cos(Pi * x -0.7) *sin( Pi * y +0.2 ) );
  double uy = cos(Pi*t)*( sin(Pi * x -0.7) * Pi* cos( Pi * y +0.2 ) );  
  double ut = -Pi*sin(Pi*t)*( sin(Pi * x -0.7) *sin( Pi * y +0.2 ) );
  
  values[0] = u;
  values[1] = ux;
  values[2] = uy;
  values[3] = -2.*Pi*Pi*u;
  values[4] = ut;
}

void ExactU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double u = cos(Pi*t)*( cos(Pi * x -0.7) *cos( Pi * y +0.2 ) );
  double ux= - cos(Pi*t)*( Pi * sin(Pi * x -0.7) *cos( Pi * y +0.2 ) );
  double uy= - cos(Pi*t)*( cos(Pi * x -0.7) *Pi *sin( Pi * y +0.2 ) ); 
  double ut = -Pi*sin(Pi*t)*( cos(Pi * x -0.7) *cos( Pi * y +0.2 ) );
  
  values[0] = u;
  values[1] = ux;
  values[2] = uy;
  values[3] = -2*Pi*Pi*u;
  values[4] = ut;
}

void ExactP(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  values[0] =  cos(Pi*t)*( sin(Pi*x) * cos(Pi*y));
  values[1] =  cos(Pi*t)*Pi*( cos(Pi*x) * cos(Pi*y) );
  values[2] = - cos(Pi*t)*Pi*( sin(Pi*x) * sin(Pi*y) );
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
  double x=0, y=0;
  switch(BdComp)
  {
    case 0: x = Param; y=0; break;
    case 1: x = 1; y = Param; break;
    case 2: x = 1-Param; y = 1; break;
    case 3: x = 0; y = 1-Param; break;
    default: cout << "wrong boundary part number" << endl;
      break;
  }
  double u[4];
  ExactU1(x,y,u);
  value= u[0];
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double x=0, y=0;

  switch(BdComp)
  {
    case 0: x = Param; y=0; break;
    case 1: x = 1; y = Param; break;
    case 2: x = 1-Param; y = 1; break;
    case 3: x = 0; y = 1-Param; break;
    default: cout << "wrong boundary part number" << endl;
      break;
  }
  double u[4];
  ExactU2(x,y,u);
  value= u[0];
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = DIMENSIONLESS_VISCOSITY;
  double u1[5], u2[5], p[5];
  
  for(int i=0;i<n_points;i++)
  {  
    ExactU1(X[i], Y[i], u1);
    ExactU2(X[i], Y[i], u2);
    ExactP(X[i], Y[i], p);
    
    coeffs[i][0] = nu;

    coeffs[i][1] = u1[4] - nu * u1[3] + p[1];
    coeffs[i][2] = u2[4] - nu * u2[3] + p[2];    
  }
}

#endif // COSINE_SINE_H
