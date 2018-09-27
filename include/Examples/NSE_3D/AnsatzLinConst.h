// Navier-Stokes problem, solution in ansatz space
// velocity pw linear, pressure constant
// 
// u(x,y) = (y+z,5x-3z,-x-2y)^T
// p(x,y) = 0

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("EXAMPLE","AnsatzLinConst.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = y+z;
  values[1] = 0;
  values[2] = 1;
  values[3] = 1;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 5*x-3*z;
  values[1] = 5;
  values[2] = 0;
  values[3] = -3;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = -x-2*y;
  values[1] = -1;
  values[2] = -2;
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
// void ExactNull(double x, double y, double z, double *values)
// {
//   values[0] =0;
//   values[1] =0;
//   values[2] =0;
//   values[3] =0;
//   values[4] =0;
// }

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  value = y+z;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 5*x-3*z;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
    value = -x-2*y ;     
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = DIMENSIONLESS_VISCOSITY;
  int i;
  double *coeff;

 if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      
      coeff[0] = eps;
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
      coeff[0] = eps;
      coeff[1] = 0; // f1
      coeff[2] = 0; // f2
      coeff[3] = 0; // f3
        /*coeff[1] = 4*x-2*y-3*z; // f1
      coeff[2] = 3*x+11*y+5*z; // f2
      coeff[3] = -10*x-y+5*z; // f3*/
    }
  }
    
}
