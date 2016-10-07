// Navier-Stokes problem with sine and cosine functions
// h = x*y^3 - x^3 * y
double DIMENSIONLESS_VISCOSITY;
void ExampleFile()
{
  TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 5;
  Output::print<1>("Example: Bsp_PR1.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = y*y*y-3*x*x*y;
  
  values[1] = -6.*x*y;
  
  values[2] = 3*y*y - 3.*x*x;
  
  values[3] = 0.;
}

void ExactU2(double x, double y, double *values)
{
  values[0]= 3.*x*y*y - x*x*x;
  
  values[1]= 3.*y*y - 3.*x*x;
  
  values[2]= 6.*x*y;
           
  values[3]= 0.;  
}

void ExactP(double x, double y, double *values)
{
  values[0]= -0.5*((y*y*y - 3.*x*x*y)*(y*y*y - 3.*x*x*y) 
             +(3*x*y*y - x*x*x)*(3*x*y*y - x*x*x)) + 12./35.;
  
  values[1] = -0.5*(6*pow(x,5)+y*y*(12*x*x*x+6*x*y*y));
  
  values[2] = -0.5*(x*x*(6*x*x*y+12.*y*y*y)+6*pow(y,5)); 
  
  values[3] = 0.0;
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
  double x, y;
  switch(BdComp)
  {
    case 0:
      x = Param; y = 0.;
      break;
    case 1: 
      x = 1.; y = Param;
      break;
    case 2:
      x = 1. - Param; y = 1.;
      break;
    case 3:
      x = 0; y = 1. - Param;
      break;
    default: 
      ErrThrow("wrong boundary part number ", BdComp);
      break;
  }
  double v[4];
  ExactU1(x,y,v);
  value = v[0];
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double x, y;
  switch(BdComp)
  {
    case 0:
      x = Param; y = 0.;
      break;
    case 1: 
      x = 1.; y = Param;
      break;
    case 2:
      x = 1. - Param; y = 1.;
      break;
    case 3:
      x = 0; y = 1. - Param;
      break;
    default: 
      ErrThrow("wrong boundary part number ", BdComp);
      break;
  }
  double v[4];
  ExactU2(x,y,v);
  value = v[0];
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  const double nu=DIMENSIONLESS_VISCOSITY;
  double val1[4];
  double val2[4];
  double val3[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = nu;
    
    ExactU1(X[i], Y[i], val1);
    ExactU2(X[i], Y[i], val2);
    ExactP(X[i], Y[i], val3);
    
    coeffs[i][1] = 0.; // f1
    coeffs[i][2] = 0.; // f2
        
    coeffs[i][3] = 0;// val1[1] + val2[2]; // g (divergence)
  }
  
}
