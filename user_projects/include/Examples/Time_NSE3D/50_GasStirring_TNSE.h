// Time-dependent 3D Navier-Stokes problem
// Corresponds to the gas stirring in a cylindrical ladle

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;
// coordinates of the plugs p1 and p2
//  double p1x = -0.566,  p1y = +0.5;
//  double p2x = -0.468,  p2y = -0.27;
double p1x = 0,  p1y = 0;
double p2x = 0,  p2y = 0;
double p_radius = 0.05;
double height = 0.4;
double free_zone = 0.2; // check that DRIFT_Z=height+free_zones
double inflow_velocity = 0.1;

void ExampleFile()
{
  Output::print<1>("Example: 50_GasStirring_TNSE.h");
}

// ========================================================================
// initial conditions
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  double in_plug1 = p_radius*p_radius-(x-p1x)*(x-p1x)-(y-p1y)*(y-p1y);
//  double in_plug2 = p_radius*p_radius-(x-p2x)*(x-p2x)-(y-p2y)*(y-p2y);

  if ( (z == 0) && (in_plug1 >= 0 ))// || in_plug2 >= 0) )
    values[0] = 0; //gas inflow
  else
    values[0] = 0;
}

void InitialP(double x, double y, double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// kind of boundary condition (for FE space needed) and values
// ========================================================================
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  double in_plug1 = p_radius*p_radius-(x-p1x)*(x-p1x)-(y-p1y)*(y-p1y);
//  double in_plug2 = p_radius*p_radius-(x-p2x)*(x-p2x)-(y-p2y)*(y-p2y);
  double total_height = height + free_zone; // = geometry height

  if ( (z == 0) && (in_plug1 >= 0) )// || in_plug2 > 0) )
    cond = DIRICHLET; //gas inflow
  else if (z == total_height)
    cond = DIRICHLET;  // outflow for the gas, in 2d it is Dirichlet ... weird
  else
    cond = DIRICHLET;

//  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  double in_plug1 = p_radius*p_radius-(x-p1x)*(x-p1x)-(y-p1y)*(y-p1y);
//  double in_plug2 = p_radius*p_radius-(x-p2x)*(x-p2x)-(y-p2y)*(y-p2y);

//  double t = TDatabase::TimeDB->CURRENTTIME;
//  bool start = 0;
//  if (t > 0)
//    start = 1;

  if ( (z == 0) && (in_plug1 >= 0))// || in_plug2 > 0) )
    value = 0; //gas inflow
  else
    value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  const double nu = DIMENSIONLESS_VISCOSITY;
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = nu;
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = -10; // f3 = gravity in z-direction
  }
}
