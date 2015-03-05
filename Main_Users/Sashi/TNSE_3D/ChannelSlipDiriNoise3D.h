// flow in a semi infinite channel

#define NOISE_FACTOR 0.00

void ExampleFile()
{
  OutPut("Example: ChannelSlipDiriNoise3D.h, NOISE "<< NOISE_FACTOR << endl);
}

// ========================================================================
// exact solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  values[0] = y*(1-y)*4;
  values[0] += -NOISE_FACTOR*4*Pi*sin(4*Pi*y);
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = -NOISE_FACTOR*3*Pi*sin(3*Pi*z);
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = NOISE_FACTOR*3*Pi*cos(3*Pi*x);
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 4*y*(1-y);
  values[1] = 0;
  values[2] = 4-8*y;
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
  double eps = 1/TDatabase::ParamDB->RE_NR;

  values[0] = -8*eps*(x-10);
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  // inflow and lower and upper wall 
  if ((fabs(x)<1e-6)||(fabs(y)<1e-6) ||(fabs(y-1)<1e-6))
  //if ((fabs(x)<1e-6))
    cond  = DIRICHLET;
  else
  {
    // outflow
    if (fabs(x-10)<1e-6)
    {
      cond = NEUMANN;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
    }
    else      
      // boundary conditions on the outer walls
    {
      cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
      TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
      //cond  = DIRICHLET;
    }
  }
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  //if( (fabs(x)<1e-6) || (fabs(x-2.5)<1e-6) )
  if((fabs(x)<1e-6))
  {
    value = 4*y*(1-y);
  }
  else
  {
    value = 0;
  }

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
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
      
    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
  }
}

