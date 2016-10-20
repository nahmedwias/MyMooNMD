// Navier-Stokes problem, SinCos (Bsp1)
// 
// u(x,y) = (sin(t)*sin(Pi*x)*sin(Pi*y), sin(t)*cos(Pi*x)*cos(Pi*y))
// p(x,y) = sin(t)*(sin(Pi*x)+cos(Pi*y)-2/Pi)

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;



void ExampleFile()
{
  OutPut("Example: 2_TestExampleForTCDTestExample1.h , Velocity Field (1,-1)" << endl) ;
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 1;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = -1;
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
  values[0] = 1;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = -1;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
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
  value = 1;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = -1;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = REYNOLDS_number;
  int i;
  double *coeff, x, y; 

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];
//    double rho = parameters[i][2];
//    double mu  = parameters[i][3];

    coeff[0] = nu;

/*
    // Stokes
    coeff[1] =0;
    coeff[2] =0;
*/

    // Navier-Stokes
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
  }
}
