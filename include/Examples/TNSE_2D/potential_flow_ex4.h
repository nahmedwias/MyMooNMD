
double DIMENSIONLESS_VISCOSITY;
void ExampleFile()
{
  OutPut("Example: potential_flow_ex4.h" << endl) ;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 97;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0]= 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
}

void ExactP(double x, double y, double *values)
{
  values[0] = x*x*x + y*y*y -1./2;
  values[1] = 3.*x*x;
  values[2] = 3.*y*y;
}
// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  value=0.;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
 value=0.;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = DIMENSIONLESS_VISCOSITY;
  for(int i=0;i<n_points;i++)
  {    
    coeffs[i][0] = nu;    
    coeffs[i][1] = 3.*X[i]*X[i];
    coeffs[i][2] = 3.*Y[i]*Y[i];
  }
}


