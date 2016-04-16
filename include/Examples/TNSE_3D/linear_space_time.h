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
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0*t;
}

void InitialU2(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0*t;
}

void InitialU3(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0*t;
}

void InitialP(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0*t;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0*t*x;;
  values[1] = 0*t;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0*t;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0*t;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}
// void ExactNull(double x, double y, double z, double *values)
// {
//   values[0] =0;
//   values[1] =0;
//   values[2] =0;
//   values[3] =0;
//   values[4] =0;
// }

// ========================================================================
// kind of boundary condition (for FE space needed) and values
// ========================================================================
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
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
    value = 0 ;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double nu = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff; // x, y, z;
  double t = TDatabase::TimeDB->CURRENTTIME;

 if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      
      coeff[0] = nu;
      coeff[1] = 0; // f1
      coeff[2] = 0; // f2
      coeff[3] = 0; // f3
    }
  }
  else
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      
//      x = X[i];
//      y = Y[i];
//      z = Z[i];
      coeff[0] = nu;
      coeff[1] = 0*t; // f1
      coeff[2] = -10; // f2
      coeff[3] = 0; // f3
    }
  }
    
}
