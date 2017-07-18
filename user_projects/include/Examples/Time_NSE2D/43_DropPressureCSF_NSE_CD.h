// Navier-Stokes problem for the drop pressure 43_TCD2D
// Should calculate the drop pressure in a fluid at rest.
// The drop is in middle and subject to surface tension only.
// initial velocity =0 everyhwere. No gravity
// See paper Brackbille et al. 1996

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;



void ExampleFile()
{
  Output::info<3>("Example: 43_DropPressureCSF_NSE_CD.h ") ;
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
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
  values[0] = 0;
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
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
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
  cond = NEUMANN;
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

    coeff[0] = nu;

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
