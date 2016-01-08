// Navier-Stokes problem with sine and cosine functions
// 

void ExampleFile()
{
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = STOKES;
  Output::print<1>("Example: quadratic_pressure.h with INTERNAL_PROBLEM_IDENTITY ", 
                   TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY);
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU2(double x, double y, double *values)
{
  values[0] = -(2*x*y*y*(x - 1)*(x - 1)*(y-1)*(y-1) + x*x*y*y*(2*x - 2)*(y-1)*(y-1) );
  
  values[1] = -(2*x*x*y*y*(y-1)*(y-1) + 2*y*y*(x - 1)*(x - 1)*(y-1)*(y-1) 
               + 4*x*y*y*(2*x - 2)*(y-1)*(y-1) );
  
  values[2] = -(2*x*y*y*(2*y - 2)*(x - 1)*(x - 1) + 2*x*x*y*(2*x - 2)*(y-1)*(y-1) 
            + x*x*y*y*(2*x - 2)*(2*y - 2) + 4*x*y*(x - 1)*(x - 1)*(y-1)*(y-1) );
  
  values[3] = -(2*x*x*y*y*(2*x - 2) + 4*x*(x - 1)*(x - 1)*(y-1)*(y-1) 
            + 2*x*x*(2*x - 2)*(y-1)*(y-1) + 6*y*y*(2*x - 2)*(y-1)*(y-1) 
            + 4*x*y*y*(x - 1)*(x - 1) + 12*x*y*y*(y-1)*(y-1) 
            + 4*x*x*y*(2*x - 2)*(2*y - 2) + 8*x*y*(2*y - 2)*(x - 1)*(x - 1) );
}

void ExactU1(double x, double y, double *values)
{
  values[0]= 2*x*x*y*(x - 1)*(x - 1)*(y - 1)*(y-1) + x*x*y*y*(2*y - 2)*(x - 1)*(x - 1);
  
  values[1]= 2*x*y*y*(2*y - 2)*(x - 1)*(x - 1) + 2*x*x*y*(2*x - 2)*(y - 1)*(y-1) 
           + x*x*y*y*(2*x - 2)*(2*y - 2) + 4*x*y*(x - 1)*(x - 1)*(y - 1)*(y-1) ;
  
  values[2]= 2*x*x*y*y*(x - 1)*(x - 1) + 2*x*x*(x - 1)*(x - 1)*(y - 1)*(y-1) 
           + 4*x*x*y*(2*y - 2)*(x - 1)*(x - 1);
           
  values[3]= 2*x*x*y*y*(2*y - 2) + 4*y*(x - 1)*(x - 1)*(y - 1)*(y-1) 
           + 6*x*x*(2*y - 2)*(x - 1)*(x - 1) + 2*y*y*(2*y - 2)*(x - 1)*(x - 1) 
           + 12*x*x*y*(x - 1)*(x - 1) + 4*x*x*y*(y - 1)*(y-1) 
           + 4*x*y*y*(2*x - 2)*(2*y - 2) + 8*x*y*(2*x - 2)*(y - 1)*(y-1);  
}

void ExactP(double x, double y, double *values)
{
  values[0]= 2*x*x*(1-x)*y*(1-y) - 1./36;
  
  values[1] = (4*x*(1-x) - 2*x*x)*y*(1-y);
  
  values[2] = 2*x*x*(1-x)*(1-2*y); 
  
  values[3] = 0.;
}

void InitialU1(double x, double y, double *values)
{
  values[0]=0;
}

void InitialU2(double x, double y, double *values)
{
  values[0]=0;
}

void InitialP(double x, double y, double *values)
{
  values[0]=0;
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
  value = 0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  const double nu=1./TDatabase::ParamDB->RE_NR;
  double val1[4];
  double val2[4];
  double val3[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = nu;
    
    ExactU1(X[i], Y[i], val1);
    ExactU2(X[i], Y[i], val2);
    ExactP(X[i], Y[i], val3);
    
    coeffs[i][1] = -nu*val1[3] + val3[1]; // f1
    coeffs[i][2] = -nu*val2[3] + val3[2]; // f2
    
    if(TDatabase::ParamDB->PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
    {
      coeffs[i][1] += val1[0]*val1[1] + val2[0]*val1[2]; // f1
      coeffs[i][2] += val1[0]*val2[1] + val2[0]*val2[2]; // f2
    }
    coeffs[i][3] = val1[1] + val2[2]; // g (divergence)
  }
  
}
