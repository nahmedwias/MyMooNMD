// Navier-Stokes problem with sine and cosine functions
// 

void ExampleFile()
{
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = OSEEN_PROBLEM;
  Output::print<1>("Example: SinCos.h with INTERNAL_PROBLEM_IDENTITY ", 
                   TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY);
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
  //cond = DIRICHLET;
  //TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
  cond = (i == 3) ? NEUMANN : DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = 1./TDatabase::ParamDB->RE_NR;
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  switch(BdComp)
  {
    case 0:
      if(cond == DIRICHLET)
        value = sin(Pi*Param); // Dirichlet
      else
        value = lt == 0 ? 0. : -nu*Pi*Pi*sin(Pi*Param); // Neumann
      break;
    case 1:
      if(cond == DIRICHLET)
        value = 0.; // Dirichlet
      else
        value = lt == 0 ? -nu*Pi : -2*nu*Pi; // Neumann
      break;
    case 2:
      if(cond == DIRICHLET)
        value = sin(Pi*(1-Param)); // Dirichlet
      else
        value = lt == 0 ? 0. : nu*Pi*Pi*sin(Pi*(1-Param)); // Neumann
      break;
    case 3:
      if(cond == DIRICHLET)
        value = 0.; // Dirichlet
      else
        value = lt == 0 ? -nu*Pi : -2*nu*Pi; // Neumann
      break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
  return;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = 1./TDatabase::ParamDB->RE_NR;
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  switch(BdComp)
  {
    case 0:
      if(cond == DIRICHLET)
        value = 0; // Dirichlet
      else
      {
        value = nu*Pi*cos(Pi*Param) * (lt == 0 ? 1. : 2.); // Neumann
        value += sin(Pi*Param);// Neumann (from pressure)
      }
      break;
    case 1:
      if(cond == DIRICHLET)
        value = Pi*Param; // Dirichlet
      else
        value = 0.; // Neumann
      break;
    case 2:
      if(cond == DIRICHLET)
        value = -Pi*cos(Pi*(1-Param)); // Dirichlet
      else
      {
        value = -nu*Pi*cos(Pi*(1-Param)) * (lt == 0 ? 1. : 2.); // Neumann
        value += sin(Pi*(1-Param));// Neumann (from pressure)
      }
      break;
    case 3:
      if(cond == DIRICHLET)
        value = -Pi*(1-Param); // Dirichlet
      else
        value = 0.; // Neumann
      break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
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
