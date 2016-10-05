// Navier-Stokes problem, Driven cavity
//
// u(x,y) = unknown
// p(x,y) = unknown

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;
double USER_parameter3;

void ExampleFile()
{
  Output::info<1>("Example", "DrivenCavity.h");
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
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
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
    break;
    case 1: value = 0;
    break;
    case 2: if(Param<0.00001 || Param>0.99999)
      value = 0;
    else
      value = 1;
    break;
    case 3: value = 0;
    break;
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
  if(BdComp>3) cout << "wrong boundary part number" << endl;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  /* Note: coeff is used in the local assembling routine in the
   * following way: coeff0 is Re-1, it is a factor in the diffusive
   * term and it is used only when parameter "dimensional_nse" is false
   * Coeff1 and Coeff2 are the right hand side. Coeff3 and 4 are
   * the density and the dynamic viscosity, rho and mu, used as factors
   * in the corresponding terms of the equations. They are active only
   * when "dimensional_nse" is true (in which case, coeff0 is not taken
   * into account anymore in the assemble routine).
   * Param are values from a FE Function used in the assemble routine:
   * param0 and 1 are the x- and y- velocity component at the point i
   * while param2 and 3 are rho and mu at this point.
   */

  double nu = REYNOLDS_number; // this is actually Re-1
  int i;
  double val1[4];   // U1-Exact function and its derivatives
  double val2[4];   // U2-Exact function and its derivatives
  double val3[4];   // P-Exact function and its derivatives
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    double rho = parameters[i][2];
    double mu  = parameters[i][3];

    ExactU1(x[i], y[i], val1);
    ExactU2(x[i], y[i], val2);
    ExactP (x[i], y[i], val3);

    // ATTENTION: IT IS NOT CONSISTENT! THE RIGHT HAND SIDE HAS TO BE
    // COMPUTED WITH RHO AND MU IN CASE OF DIMENSIONAL NSE!
    coeff[0] = nu;
    coeff[1] = -nu*val1[3] + val3[1]; // f1
    coeff[2] = -nu*val2[3] + val3[2]; // f2

    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
    {
      coeff[1] += val1[0]*val1[1] + val2[0]*val1[2]; // f1
      coeff[2] += val1[0]*val2[1] + val2[0]*val2[2]; // f2
    }

    coeff[3] = rho;     // rho
    coeff[4] = mu;      // mu

  }
}
