// time dependent Navier-Stokes problem 3D, ansatz
//
//u1(t,x,y,z) = 2*t
//u2(t,x,y,z) = 3
//u3(t,x,y,z) = 4
//p(t,x,y,z) = 0

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("Example", "Bsp0.h");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 2*t;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
values[0] = 3;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
values[0] = 4;
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
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 2*t;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void ExactU2(double x, double y, double z, double *values)
{
values[0] = 3;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void ExactU3(double x, double y, double z, double *values)
{
values[0] = 4;
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
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 2*t;
}

void U2BoundValue(double x, double y, double z, double &value)
{
value = 3;
}

void U3BoundValue(double x, double y, double z, double &value)
{
value = 4;
}

void U1BoundValue_diff(double x, double y, double z, double &value)
{
value = 2;
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
    coeffs[i][1] = 2;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;
    coeffs[i][4] = 0;
    coeffs[i][5] = 0;
    coeffs[i][6] = 0;
  }
}
