// Time-dependent 3D Navier-Stokes problem
// Corresponds to the dam break in 3D

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::print<1>("Example: 40_DamBreak3D_TNSE.h");
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
