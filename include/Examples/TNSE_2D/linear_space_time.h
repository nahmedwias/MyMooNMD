// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (t*y, 2*t*y)
// p(x,y) = 0

void ExampleFile()
{
  OutPut("Example: linear_space_time.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*y;
}

void InitialU2(double x, double, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 2*t*x;
}

void InitialP(double, double, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*y;
  values[1] = 0;
  values[2] = t;
  values[3] = 0;
}

void ExactU2(double x, double, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 2*t*x;
  values[1] = 2*t;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = 0;
    break;
    case 1:
      value = t*Param;
    break;
    case 2:
      value = t;
    break;
    case 3:
      value = t*(1-Param);
    break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = 2*t*Param;
    break;
    case 1:
      value = 2*t;
    break;
    case 2:
      value = 2*t*(1-Param);
    break;
    case 3:
      value = 0;
    break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double **, double **coeffs)
{
  static double nu = 1/TDatabase::ParamDB->RE_NR;
  double t = TDatabase::TimeDB->CURRENTTIME;
  int i;
  double *coeff, x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;

    coeff[1] = y+2*t*t*x;
    coeff[2] = 2*(x+t*t*y);
  }
}
