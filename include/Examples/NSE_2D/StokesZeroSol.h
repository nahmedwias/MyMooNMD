// Navier-Stokes problem with sine and cosine functions
// 

void ExampleFile()
{
  Output::print<1>("Example: StokesZeroSol.h with INTERNAL_PROBLEM_IDENTITY ", 
                   TDatabase::ParamDB->PROBLEM_TYPE);
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  // psi = x^2*(1-x)^2*y^2(1-y)^2
    values[0] = 0.;
    values[1] = 0.;
    values[2] = 0.;
    values[3] = 0.;
}

void ExactU2(double x, double y, double *values)
{
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;
    values[3] = 0.0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 2*x*x*(1-x)*y*(1-y) - 1./36;
  values[1] = (4*x*(1-x) - 2*x*x)*y*(1-y);
  values[2] = 2*x*x*(1-x)*(1-2*y); 
  values[3] = 0.;
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
  value = 0.;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0.;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  const double nu=1./TDatabase::ParamDB->RE_NR;  
  double u1[4], u2[4], p[4];
  
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = nu;
    
    ExactU1(X[i], Y[i], u1);
    ExactU2(X[i], Y[i], u2);
    ExactP(X[i], Y[i], p);
    
    coeffs[i][1] = -u1[3] + p[1]; // f1
    coeffs[i][2] = -u2[3] + p[2]; // f2    
    
    if(TDatabase::ParamDB->PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
    {
      ErrThrow("flow problem type: " , TDatabase::ParamDB->PROBLEM_TYPE ," not tested yet");
    }
    coeffs[i][3] = u1[1] + u2[2]; // g (divergence)
  }  
}
