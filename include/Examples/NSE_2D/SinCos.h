// Navier-Stokes problem with sine and cosine functions
// 

void ExampleFile()
{
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = OSEEN_PROBLEM;
  OutPut("Example: SinCos.h with INTERNAL_PROBLEM_IDENTITY " << 
	 TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY <<  endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    values[0] = sin(Pi*x);
    values[1] = Pi*cos(Pi*x);
    values[2] = 0;
    values[3] = -Pi*Pi*sin(Pi*x);
}

void ExactU2(double x, double y, double *values)
{
    values[0] = -Pi*y*cos(Pi*x);
    values[1] = Pi*Pi*y*sin(Pi*x);
    values[2] = -Pi*cos(Pi*x);
    values[3] = Pi*Pi*Pi*y*cos(Pi*x);
}

void ExactP(double x, double y, double *values)
{
  values[0] = sin(Pi*x)*cos(Pi*y);
  values[1] = Pi*cos(Pi*x)*cos(Pi*y);
  values[2] = -Pi*sin(Pi*x)*sin(Pi*y);
  values[3] = -Pi*Pi*sin(Pi*x)*cos(Pi*y)-Pi*Pi*sin(Pi*x)*cos(Pi*y);
}

void InitialU1(double x, double y, double *values)
{
  values[0] = sin(Pi*x);
}

void InitialU2(double x, double y, double *values)
{
  values[0] = -Pi*y*cos(Pi*x);
}

void InitialP(double x, double y, double *values)
{
  values[0] = sin(Pi*x)*cos(Pi*y);
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=sin(Pi*Param);
            break;
    case 1: value=sin(Pi);
            break;
    case 2: value=sin(Pi*(1-Param));
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=Pi*Param;
            break;
    case 2: value=-Pi*cos(Pi*(1-Param));
            break;
    case 3: value=-Pi*(1-Param);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
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
