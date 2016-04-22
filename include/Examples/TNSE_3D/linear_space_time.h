// Time-dependent 3D Navier-Stokes problem, Poiseuille Problem
// 
// u(x,y,z) = (t*y,2*t*y,0)
// p(x,y,z) = 0

void ExampleFile()
{
  Output::print<1>("Example: linear_space_time.h");
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
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
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
  if ( y == 1 )
    cond = NEUMANN;
  else
    cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
    value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  if ( (y == 0) && (x >= 0.4) && (x <= 0.6) && (z >= 0.4) && (z <= 0.6))
    value = 10*t;
//  else if ( y == 1 )
//    value = 10*t;
  else
    value = 0;
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
  const double nu = 1/TDatabase::ParamDB->RE_NR;
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
      coeffs[i][0] = nu;
      coeffs[i][1] = 0; // f1
      coeffs[i][2] = -10; // f2
      coeffs[i][3] = 0; // f3
    }
  }
}
