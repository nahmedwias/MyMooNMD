#ifndef CORIOLIS_BETAPLANE
#define CORIOLIS_BETAPLANE

#include <math.h>
// u1 = cos(t) * ( sin(\Pi * x -0.7) *sin( \Pi * y +0.2 ) );
// u2 = cos(t) * ( cos(\Pi * x -0.7) *cos( \Pi * y +0.2 ) );
// p =  cos(t) * ( sin(x) * cos(y) * (cos(1) * sin(1) - sin(1)) );
//
double DIMENSIONLESS_VISCOSITY;
void ExampleFile()
{
  Output::print("Example: TNSE_2D/cosine_sine_poly.h") ;
}



// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{ 
  values[0] = 1.0;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 0.;  
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 0.;
}

void ExactP(double x, double y, double *values)
{
  values[0] =  y*y - 1./3.;
  values[1] =  0.;
  values[2] =  2.*y;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double val[4];
  ExactU1(x,y,val);
  values[0] = val[0];
}

void InitialU2(double x, double y, double *values)
{
  double val[4];
  ExactU2(x,y,val);
  
  values[0] = val[0];
}

void InitialP(double x, double y, double *values)
{
  double val[4];
  ExactP(x,y,val);
  
  values[0] = val[0];
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
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  for(int i=0;i<n_points;i++)
  {  
    double x = X[i];
    double y = Y[i];
    
    coeffs[i][0] = nu;
    coeffs[i][1] = 0.;
    coeffs[i][2] = 0.;
    coeffs[i][3] = 2.*y; // angular velocity (third component)
  }
}

#endif // COSINE_SINE_H
