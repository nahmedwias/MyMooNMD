// This is a test to to calculate the pressure of a drop (circle)
// submitted only to surface tension, lying in a fluid at rest.
// The circle has density 1000kg/m3, and radius 0.02m
// the fluid in background has 500kg/m3, and
// the surface tension tau is 0.02364 N/m
// In this test, velocity is 0, only DeltaP should be non zero
// It should be equal to tau*kappa = tau/R, according to Laplace result.
// Initial solution is known. The problem is stationary.
// No gravity. Geometry should be UnitSquare.
// BENCHMARK OF EQUILIBRIUM ROD, FROM BRACKBILL ET AL. 1996



void ExampleFile()
{
  Output::info<3>("Example: 43_DropPressureCSF_CD_NSE.h");
}

double get_nu()
{
  return 1;
}

constexpr bool rhs_depends_on_time = false;
constexpr bool coefficients_depend_on_time = true;

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double Param, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int i, double Param, double &value)
{
  value = 1; //background liquid everywhere, drop in the middle
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  double pi = 3.14159265358979;
  double x0 = TDatabase::ParamDB->P4;
  double y0 = TDatabase::ParamDB->P5; // center of unit square=center of drop
  double R = TDatabase::ParamDB->P6; // radius of drop, 2cm
//  if ( (x-x0)*(x-x0)+(y-y0)*(y-y0) <= R*R )
//    values[0] = 0;
//  else
//    values[0] = 1;

  double phi = R*R - (x-x0)*(x-x0) - (y-y0)*(y-y0); // level set of circle
//  double a=0.2; double b=0.1;
//  double phi = 1-((x-x0)*(x-x0)/(a*a))-((y-y0)*(y-y0)/(b*b)); // level set of ellipse
  double eps = 0.0;
  if (phi < -eps)
    values[0] = 1;
  else if ( phi > eps)
    values[0] = 0;
  else
    values[0] =0.5*(1 - (phi/eps) - (1/pi)*sin(phi*pi/eps));

//  /* BELOW, SQUARE FOR NON-EQUILIBRIUM ROD */
//  if (x <= 0.045 && x >= 0.015)
//  {
//    if (y <= 0.045 && y >= 0.015)
//      values[0] = 0;
//    else
//      values[0] = 1;
//  }
//  else
//    values[0] = 1;

}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
//  double eps=1/TDatabase::ParamDB->RE_NR;
//  double a=1, b=2, c=1;
  double u_x,u_y;
  int i;
  double *coeff;
//  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    x = X[i];
//    y = Y[i];

    u_x = parameters[i][0];
    u_y = parameters[i][1];

    coeff[0] = 0;
    coeff[1] = u_x;//0;
    coeff[2] = u_y;//0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
  double R = 0.02; // radius of drop, 2cm
  double x0 = 0.5; double y0 = 0.5; // center of unit square=center of drop
  if ( (x*x-x0*x0) + (y*y-y0*y0) <= R*R )
    values[0] = 1;
  else
    values[0] = 0;
}

