// Time-dependent 3D Navier-Stokes problem, Poiseuille Problem
// 
// u(x,y,z) = (t*y,2*t*x,0)
// p(x,y,z) = 0

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("Example","linear_space_time.h");
}

// ========================================================================
// initial conditions
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*y;
}

void InitialU2(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = 2*t*x;
}

void InitialU3(double x, double y, double z, double *values)
{
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
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*y;
  values[1] = 0;
  values[2] = t;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = 2*t*x;
  values[1] = 2*t;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// ========================================================================
// kind of boundary condition (for FE space needed) and values
// ========================================================================
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  value = t*y;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  value = 2*t*x;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
    value = 0 ;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  const double nu = DIMENSIONLESS_VISCOSITY;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double x, y;
  if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
  {
    for(int i = 0; i < n_points; i++)
    {
      coeffs[i][0] = nu;
      coeffs[i][1] = 0; // f1
      coeffs[i][2] = 0; // f2
      coeffs[i][3] = 0; // f3
    }
  }
  else
  {
    for(int i = 0; i < n_points; i++)
    {
      x = X[i];
      y = Y[i];

      coeffs[i][0] = nu;
      coeffs[i][1] = y+2*t*t*x; // f1
      coeffs[i][2] = 2*(x+t*t*y); // f2
      coeffs[i][3] = 0; // f3
    }
  }
}
