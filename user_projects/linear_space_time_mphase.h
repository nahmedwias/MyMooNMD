
void ExampleFile()
{
  OutPut("Example: lsp_tnse_tcd2D.h" << endl) ; //linear_space_time
}

double get_nu()
{
  return 1;
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = false;



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

void InitialT(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x+3*t*y;
  values[1] = 2;
  values[2] = 3*t;
  values[3] = 0;
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
  void ExactT(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x+3*t*y;
  values[1] = 2;
  values[2] = 3*t;
  values[3] = 0;
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

void TBoundValue(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = 1+2*Param;
    break;
    case 1:
      value = 3+3*Param*t;
    break;
    case 2:
      value = 1+3*t+2*(1-Param);
    break;
    case 3:
      value = 1+3*t*(1-Param);
    break;
  } // endswitch
}
// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double **param, double **coeffs)
{
  static double nu = 1/TDatabase::ParamDB->RE_NR;
  double eps=1/TDatabase::ParamDB->RE_NR;
  double t = TDatabase::TimeDB->CURRENTTIME;
  int i;
  double x, y;
  double a=1, b=2, c=1;

  for(i=0;i<n_points;i++)
  {
    x = X[i];
    y = Y[i];

    coeffs[i][0] = nu;
    coeffs[i][1] = y+2*t*t*x;      //f1 ns; 
    coeffs[i][2] = 2*(x+t*t*y);   //f2 ns;
    
    coeffs[i][3] =eps;
    coeffs[i][4] =param[i][0];
    coeffs[i][5] =param[i][1];     
  }
}
