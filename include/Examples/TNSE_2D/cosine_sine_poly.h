#ifndef COSINE_SINE_POLY
#define COSINE_SINE_POLY

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
  double t = TDatabase::TimeDB->CURRENTTIME;
  double xp2=x*x;
  double yp2=y*y;
  double xm1=(1.0-x);
  double ym1=(1.0-y);
  double xm1p2=xm1*xm1;
  double ym1p2=ym1*ym1;
  double tempsx = sin(2*Pi*x);
  double tempcy = cos(2*Pi*y);
  double tempcx = cos(2*Pi*x);
  double tempsy = sin(2*Pi*y);
  
  double u1 = 2.0*xm1p2*xp2*ym1p2*y-2.0*xm1p2*xp2*ym1*yp2;
  double u2 = -4*Pi*tempsx*tempsx*tempcy*tempsy;
  double u = cos(Pi*t)*128*u1 + sin(Pi*t)*u2;
  double du1 = 4. * x * (1. - 3. * x + 2. * xp2) * y * (1. - 3. * y + 2. * yp2);
  double du2 = 2. * (-1. + x)*(-1. + x) * xp2 * (1. - 6. * y + 6. * yp2);
  double du2_1 = -16*Pi*Pi*tempcx*tempcy*tempsx*tempsy;
  double du2_2 = -8*Pi*Pi*tempsx*tempsx*(tempcy*tempcy - tempsy*tempsy);
  
  values[0] = u;
  values[1] = cos(Pi*t)*128*du1+sin(Pi*t)*du2_1;
  values[2] = cos(Pi*t)*128*du2+sin(Pi*t)*du2_2;
  values[3] = 0.;
  values[4] = 0.;  
}

void ExactU2(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double xp2=x*x;
  double yp2=y*y;
  double xm1=(1.0-x);
  double ym1=(1.0-y);
  double xm1p2=xm1*xm1;
  double ym1p2=ym1*ym1;
  double tempsx = sin(2*Pi*x);
  double tempcy = cos(2*Pi*y);
  double tempcx = cos(2*Pi*x);
  double tempsy = sin(2*Pi*y);
  
  double v1 = -2.0*xm1p2*x*ym1p2*yp2+2.0*xm1*xp2*ym1p2*yp2;
  double v2 = 4*Pi*tempcx*tempsx*tempsy*tempsy;
  double u = cos(Pi*t)*128*v1 + sin(Pi*t)*v2; 
  double dv1 = -2. * (1. - 6. * x + 6. * xp2) * (-1. + y)*(-1. + y) * yp2;
  double dv2 = -4. * x * (1. - 3. * x + 2. * xp2) * y * (1. - 3. * y + 2. * yp2);
  double dv2_1 = 8*Pi*Pi*tempsy*tempsy*(tempcx*tempcx - tempsx*tempsx);
  double dv2_2 = 16*Pi*Pi*tempcx*tempcy*tempsx*tempsy ;
  
  values[0] = u;
  values[1] = cos(Pi*t)*128*dv1+sin(Pi*t)*dv2_1;
  values[2] = cos(Pi*t)*128*dv2+sin(Pi*t)*dv2_2;
  values[3] = 0.;
  values[4] = 0.;
}

void ExactP(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  values[0] =  128.*(x*x*x + y*y*y - 0.5);
  values[1] =  128.*(3.*x*x);
  values[2] =  128.*(3.*y*y);
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
    coeffs[i][0] = nu;
    
    double x = X[i];
    double y = Y[i];
    
    double xp2=x*x;
    double yp2=y*y;
    double yp3=yp2*y;
    double yp4=yp2*yp2;
    double xp3=xp2*x;
    double xp4=xp2*xp2;
    double xm1=(1.-x);
    double ym1=(1.-y);
    double xm1p2=xm1*xm1;
    double ym1p2=ym1*ym1;
    
    double f1 = 128.*cos(Pi*t)*nu*(-4.*(2.*y - 1.)*(3.*xp4 - 6.*xp3 + 6.*xp2*yp2 
                - 6.*xp2*y + 3.*xp2 - 6.*x*yp2 + 6*x*y + yp2 - y));
    double f2 = 128.*cos(Pi*t)*nu*4*(2*x - 1)*(3.*yp4 - 6.*yp3 + 6.*yp2*xp2
                - 6.*yp2*x + 3.*yp2 - 6*y*xp2 + 6.*y*x + xp2 - x);
    f1 = f1 + sin(Pi*t)*nu*16.*Pi*Pi*Pi*(-1. + 2.*cos(4.*Pi*x))*sin(4.*Pi*y);
    f2 = f2 - sin(Pi*t)*nu*16.*Pi*Pi*Pi*(-1. + 2.*cos(4.*Pi*y))*sin(4.*Pi*x);
    f1 = f1 + 128.*3.*x*x;
    f2 = f2 + 128.*3.*y*y;
    double tempsx = sin(2.*Pi*x);
    double tempcy = cos(2.*Pi*y);
    double tempcx = cos(2.*Pi*x);
    double tempsy = sin(2.*Pi*y);
    double u1 = 2.0*xm1p2*xp2*ym1p2*y-2.0*xm1p2*xp2*ym1*yp2;
    double v1 = -2.0*xm1p2*x*ym1p2*yp2+2.0*xm1*xp2*ym1p2*yp2;
    double u2 = -4*Pi*tempsx*tempsx*tempcy*tempsy;
    double v2 = 4.*Pi*tempcx*tempsx*tempsy*tempsy;
    f1 = f1 - Pi*sin(Pi*t)*128.*u1 + Pi*cos(Pi*t)*u2;
    f2 = f2 - Pi*sin(Pi*t)*128.*v1 + Pi*cos(Pi*t)*v2;
    coeffs[i][1] = f1;
    coeffs[i][2] = f2;
  }
}

#endif // COSINE_SINE_H
