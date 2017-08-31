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
    // interpolate experimental data (50m downstream)
    // value = ...TODO...;
  }
  else if(BdComp == 19) // inflow boundary
  {   
    double y = HEIGHT*(1-Param);
    if (y<=Z_0) {
      value = 0;
    } else {
      value = U_REF/log(Z_R/Z_0)*log(y/Z_0);
    }
    /*
    // interpolate experimental data (50m downstream)
    std::vector<double> yval(15),uval(15);
    yval[0] = 0.; uval[0] = 0.;
    yval[1] = 1.; uval[1] = 3.17;
    yval[2] = 2.; uval[2] = 3.41;
    yval[3] = 3.; uval[3] = 3.59;
    yval[4] = 4.; uval[4] = 3.78;
    yval[5] = 5.; uval[5] = 3.91;
    yval[6] = 7.5; uval[6] = 4.18;
    yval[7] = 10.; uval[7] = 4.35;
    yval[8] = 12.5; uval[8] = 4.54;
    yval[9] = 15.; uval[9] = 4.61;
    yval[10] = 20.; uval[10] = 4.8;
    yval[11] = 30.; uval[11] = 5.18;
    yval[12] = 40.; uval[12] = 5.41;
    yval[13] = 50.; uval[13] = 5.75;
    yval[14] = 60.; uval[14] = 6.04;
    int i_range = -1;
    for (unsigned int i=0; i<14; i++)
    {
      if ( (y>=yval[i]) && (y<yval[i+1]) )
	i_range = i;
    }
    if (i_range>=0)
    {
      value = uval[i_range] +
	(uval[i_range+1]-uval[i_range])*
	(y - yval[i_range])/(yval[i_range+1] - yval[i_range]);
    } else {
      Output::print(" ** ERROR: point at the inflow now found ** ");
      exit(1);
    }
    */
    
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
