// h = x*y^3 - x^3 y
// u = t*grad h

void ExampleFile()
{
  Output::print<1>("Example: Navier_Stokes_Test_code1.h") ;
}

#include<iostream>
#include <math.h>
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(y*y*y - 3.*x*x*y);
}

void InitialU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(3.*x*x*y- x*x*x);
}

void InitialP(double x, double y, double *values)
{
  values[0] = pow(x, 6)+ 3*pow(x,4)*y*y-pow(x,3)*y +3*x*x* pow(y,4)
              +x *pow(y,3)+pow(y,6);
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  // {y*y*y-3*x*x*y);
  values[0] = t*(y*y*y - 3.*x*x*y);
  values[1] = t*(-6.*x*y);
  values[2] = t*3.*(y*y-x*x);
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  // 3 x y^2-x^3
  values[0] = t*(3.*x*x*y- x*x*x); // 
  values[1] = t*(6.*x*y -3.*x*x);
  values[2] = t*(3*x*x);
  values[3] = 6.*t*(y-x);
}

void ExactP(double x, double y, double *values)
{
  // values[0] = x*y*y*y - x*x*x*y;
  // x (x (x (x (12. y^2-0.5 x^2)-y)-10.5 y^4)+y^3)+y^6
  values[0] = x*(x*(x*(x*(12.*y*y-0.5*x*x)-y)-10.5*pow(y,4))+pow(y,3))+pow(y,6);
  // y (y (y (1.-21. x y)+48. x^3)-3. x^2)-3. x^5
  values[1] = y*(y*(y*(1.-21.*x*y)+48.*pow(x,3))-3.*pow(x,2))-3.*pow(x,5);
  // x (x (x (24. x y-1.)-42. y^3)+3. y^2)+6. y^5
  values[2] = x*(x*(x*(24.*x*y-1.)-42.*pow(y,3))+3.*y*y)+6.*pow(y,5);
  // not needed
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
  double t = TDatabase::TimeDB->CURRENTTIME;
  double u1[4], u2[4], p[4];
  
  for(int i=0;i<n_points;i++)
  {
    ExactU1(X[i], Y[i], u1);
    ExactU2(X[i], Y[i], u2);
    ExactP(X[i], Y[i], p);
    
    coeffs[i][0] = nu;
    // Navier-Stokes
    coeffs[i][1] = 1./t*u1[0] - nu*u1[3] + p[1];
    coeffs[i][2] = 1./t*u2[0] - nu*u2[3] + p[2];
    
//    if(TDatabase::ParamDB->PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
//    {
//      coeffs[i][1] += u1[0]*u1[1] + u2[0]*u1[2]; // f1
//      coeffs[i][2] += u1[0]*u2[1] + u2[0]*u2[2]; // f2
//    }
  }
}


