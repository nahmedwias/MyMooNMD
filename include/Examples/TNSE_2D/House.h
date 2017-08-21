#ifndef HOUSE_H
#define HOUSE_H

#include <math.h>

double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::print("Example: CowHouse.h");
}

void InitialU1(double x, double y, double *values)
{
  values[0] = 0.;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0.;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0.;
}

void ExactU1(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}

void BoundCondition(int i, double t, BoundCond &cond)
{
  
  if(i == 17)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
  if(i>27)
  {
    Output::print("cannot assign a boundary condition to component ",i);
    exit(-4711);
  }
  
}

void U1BoundValue(int BdComp, double Param, double &value)
{

  // logarithmic profile [provided by D.Janke]
  // u(z) = u_ref*log(z/z_0)/log(z_r/z_0)
  double U_REF,Z_R,Z_0,HEIGHT;
  U_REF = 6.46; //m/s
  Z_0 = 0.07; //m
  Z_R = 10; //m
  HEIGHT = 100; //m
  
  // top open boundary
  if(BdComp == 18)
  {
    value = U_REF/log(Z_R/Z_0)*log(HEIGHT/Z_0);
  }
  else if(BdComp == 19) // inflow boundary
  {
    double y = HEIGHT*(1-Param);
    if (y<=Z_0) {
      value = 0;
    } else {
      value = U_REF/log(Z_R/Z_0)*log(y/Z_0);
    }
    //cout << y << "," << value << ";" << endl;
  }
  else // no-slip
  {
    value = 0.;
  }
  
  if(BdComp>27)
  {
    Output::print("wrong boundary", BdComp);
    exit(-4711);
  }
  

}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0.;
}

void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  for(int i=0; i<n_points; i++)
  {
    double *coeff = coeffs[i];

    coeff[0] = 1.81*1e-5/1.225; //m^2/s = [kg/(m.s)] / [kg/m^3]
    coeff[1] = 0.;
    coeff[2] = 0.;
  }
}


#endif // HOUSE_H
