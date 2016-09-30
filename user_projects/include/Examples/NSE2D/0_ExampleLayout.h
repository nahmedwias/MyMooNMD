// THIS IS AN EXAMPLE LAY-OUT USED FOR DIMENSIONAL NSE
// COPY-PASTE IT AND CHANGE IT TO MAKE YOUR OWN EXAMPLE

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;
double USER_parameter3;

void ExampleFile()
{
  Output::info<1>("Example", "Example_Layout.h");
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
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
    break;
    case 1: value=0;
    break;
    case 2: value=0;
    break;
    case 3: value=0;
    break;
    default: cout << "wrong boundary part number" << endl;
    break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
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

  double eps = REYNOLDS_number; // this is actually Re-1
  int i;
  double *coeff;
  //  double rho = USER_parameter1;
  //  double mu  = USER_parameter2;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0;                    // f1
    coeff[2] = 0;                    // f2
    coeff[3] = parameters[i][2];     // rho
    coeff[4] = parameters[i][3];     // mu
    // cout << coeff[3] << " " << coeff[4] << endl;
  }
}
