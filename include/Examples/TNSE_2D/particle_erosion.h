#ifndef PARTICLE_EROSION_H
#define PARTICLE_EROSION_H

double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::print("Example: particle_erosion.h" );
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void InitialP(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 0;
}

void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

void U1BoundValue_diff(double x, double y, double z, double &value)
{
  value = 0;
}

void U2BoundValue_diff(double x, double y, double z, double &value)
{
  value = 0;
}

void U3BoundValue_diff(double x, double y, double z, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  const double eps = DIMENSIONLESS_VISCOSITY;
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;
    coeffs[i][4] = 0;
    coeffs[i][5] = 0;
    coeffs[i][6] = 0;
  }
}

#endif // PARTICLE_EROSION_H
