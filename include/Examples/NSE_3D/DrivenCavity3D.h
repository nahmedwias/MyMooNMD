/**
 * @file The driven cavity benchmark example in 3D.
 *
 * The boundary data is adapted to the [0.1]^3 unit cube example, which is
 * availabe as default geometry in ParMooN. It will throw an error if you
 * try running it on any other domain - just to make you aware of that fact.
 */

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("EXAMPLE"," DrivenCavity3D.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double, double, double z, double &value)
{
  double tol = 1e-10;

  if (fabs(1-z)<tol)
  {
    value = 1.0;
  }
  else
    value = 0.0 ;
}

// value of boundary condition
void U2BoundValue(double, double, double, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(double, double, double, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *, double *, double *, double **,
               double **coeffs)
{
  static double eps = DIMENSIONLESS_VISCOSITY;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // f3
    coeff[4] = 0; // g
    coeff[5] = 0; // sigma
  }
}
