// Time-dependent 3D Navier-Stokes problem, Poiseuille Problem
// 
// u(x,y,z) = (t*y,2*t*x,0)
// p(x,y,z) = 0

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::print<1>("Example: linear_space_time.h");
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
  values[5] = y;// time derivative
}

void ExactU2(double x, double y,  double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = 2*t*x;
  values[1] = 2*t;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 2*x;// time derivative
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;// time derivative
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
// initial conditions
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double u[6];
  ExactU1(x,y,z,u);
  values[0] = u[0] ;
}

void InitialU2(double x, double y, double z, double *values)
{
  double u[6];
  ExactU2(x,y,z,u);
  values[0] = u[0] ;
}

void InitialU3(double x, double y, double z, double *values)
{
  double u[6];
  ExactU3(x,y,z,u);
  values[0] = u[0] ;
}

void InitialP(double x, double y, double z, double *values)
{
  values[0] = 0;
}
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
  double u[6];
  ExactU1(x,y,z,u);
  value = u[0] ;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  double u[6];
  ExactU2(x,y,z,u);
  value = u[0] ;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
    double u[6];
    ExactU3(x,y,z,u);
    value = u[0] ;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  const double nu = DIMENSIONLESS_VISCOSITY;
  double x, y, z;
  double u1[6], u2[6], u3[6], p[6];
  
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
      z = Z[i];
      ExactU1(x,y,z,u1);
      ExactU2(x,y,z,u2);
      ExactU3(x,y,z,u3);
      ExactP(x,y, z,p);

      coeffs[i][0] = nu;
      
      coeffs[i][1] = u1[5]-nu*u1[4] + (u1[0]*u1[1] + u2[0]*u1[2] + u3[0]*u1[3] ) + p[1];
      coeffs[i][2] = u2[5]-nu*u2[4] + (u1[0]*u2[1] + u2[0]*u2[2] + u3[0]*u2[3] ) + p[2];
      coeffs[i][3] = u3[5]-nu*u3[4] + (u1[0]*u3[1] + u2[0]*u3[2] + u3[0]*u3[3] ) + p[3];
    }
  }
}
