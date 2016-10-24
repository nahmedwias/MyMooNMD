// Navier-Stokes problem for the ROTATING SEMI-CIRCLE TEST
// THIS IS USED TO TEST THE BENCHMARK MULTIPHASE EXAMPLE
// WHERE A SEMI_CIRCLE (REPRESENTED BY CONVECTION EQUATION)
// IS CONVECTED BY A CONSTANT ROTATING VELOCITY FIELD
// The following lines should generate this constant rotating
// velocity field (=omega.r centered in the unit square box)
// omega is 2 rd/sec
// GRAVITY IS 0.

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;



void ExampleFile()
{
  Output::info<3>("Example: 6_TestTNSE2D.h ") ;
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 2.0*(y-0.5);
  values[1] = 0;
  values[2] = 2.0;
  values[3] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = -2.0*(x-0.5);
  values[1] = -2.0;
  values[2] = 0;
  values[3] = 0;
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
  values[0] = 2.0*(y-0.5);
  values[1] = 0;
  values[2] = 2.0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = -2.0*(x-0.5);
  values[1] = -2.0;
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
  if (BdComp == 0) value = 2.0*(0-0.5);
  if (BdComp == 1) value = 2.0*(Param-0.5);
  if (BdComp == 2) value = 2.0*(1-0.5);
  if (BdComp == 3) value = 2.0*(0.5-Param);
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  if (BdComp == 0) value = -2.0*(Param-0.5);
  if (BdComp == 1) value = -2.0*(1-0.5);
  if (BdComp == 2) value = -2.0*(0.5-Param);
  if (BdComp == 3) value = -2.0*(0-0.5);
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = REYNOLDS_number;
  int i;
  double *coeff;
//  double  x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    x = X[i];
//    y = Y[i];
//    double rho = parameters[i][2];
//    double mu  = parameters[i][3];

    coeff[0] = nu;    // coefficient in front of viscosity term

/*
    // Stokes
    coeff[1] =0;
    coeff[2] =0;
*/

    // Navier-Stokes
    coeff[1] = 0;     // f1
    coeff[2] = 0;     // f2
    coeff[3] = 0;
  }
}
