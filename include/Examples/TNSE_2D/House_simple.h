#ifndef HOUSE_SIMPLE_H
#define HOUSE_SIMPLE_H

#include <math.h>

double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::print("Example: House_simple.h");
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
  
  if(i == 1)
    cond = NEUMANN;
  else
  {
    cond = DIRICHLET;
    for (int k=0;k<TDatabase::ParamDB->n_slip_boundary;k++)
    {
      int slip_boundary = TDatabase::ParamDB->slip_boundary_id[k];
	if (i == slip_boundary)
	{
	  cond = NEUMANN;
	  //	  Output::print(" component ",i," -> SLIP BC");
	}	    
    }	  
  }
  
  if(i>3)
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

  // extrapolation constants
  double A1,B1,C1;
  A1 = 1.1;
  B1 = 2.0;
  C1 = 0.66;
    
  // top open boundary
  if(BdComp == 2)
  {
    //value = U_REF/log(Z_R/Z_0)*log(HEIGHT/Z_0);
    // interpolate experimental data (50m downstream)
    value = A1*log(B1*HEIGHT) + C1;
    
  }
  else if(BdComp == 3) // inflow boundary
  {   
    double y = HEIGHT*(1-Param);
    /*if (y<=Z_0) {
      value = 0;
    } else {
      value = U_REF/log(Z_R/Z_0)*log(y/Z_0);
      }*/
    
    // interpolate experimental data (50m downstream)
    std::vector<double> yval(15),uval(15);
    yval[0] = 0.; uval[0] = 3.2;     
    yval[1] = 1.; uval[1] = 3.52;
    yval[2] = 2.; uval[2] = 3.77;
    yval[3] = 3.; uval[3] = 3.94;
    yval[4] = 4.; uval[4] = 4.01;
    yval[5] = 5.; uval[5] = 4.10;
    yval[6] = 7.5; uval[6] = 4.24;
    yval[7] = 10.; uval[7] = 4.41;
    yval[8] = 12.5; uval[8] = 4.55;
    yval[9] = 15.; uval[9] = 4.67;
    yval[10] = 20.; uval[10] = 4.82;
    yval[11] = 30.; uval[11] = 5.18;
    yval[12] = 40.; uval[12] = 5.46;
    yval[13] = 50.; uval[13] = 5.71;
    yval[14] = 60.; uval[14] = 5.96;
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
    } else { // point between 60 and 100
      value = A1*log(B1*y) + C1;
      
    }
    
    
  }
  else  if(BdComp == 0)
  {
    value = 3.2;
  }
  else // no-slip
  {
    value = 0.;
  }
  
  if(BdComp>3)
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


#endif // HOUSE_SIMPLE_H
