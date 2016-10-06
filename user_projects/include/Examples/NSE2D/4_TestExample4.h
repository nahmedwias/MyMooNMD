// Navier-Stokes problem, CREATED EXAMPLE TO GENERATE A (-y,x)
// VELOCITY FIELD
// u(x,y) = (-y,x)
// p(x,y) = 0

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;
double USER_parameter3;

void ExampleFile()
{
  Output::info<1>("Example", "Velocity for two interior layers");
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = -y;
  values[1] = 0;
  values[2] = -1;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = x;
  values[1] = 1;
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
    case 1: value = -Param;
    break;
    case 2: value = -1;
    break;
    case 3: value = -1+Param;
    break;
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = Param;
    break;
    case 1: value = 1;
    break;
    case 2: value = 1-Param;
    break;
    case 3: value = 0;
    break;
    default: cout << "wrong boundary part number" << endl;
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

    /* Attention, when dimensional_nse is false,
    // the reynolds number should be input by user. Then the next line
    // should be commented. When dimensional_nse is true, it means
    // that only rho and mu are of importance, and then
    // the next line should be used. In practice, we do it the following way:
    // the user ensures that rho, mu and reynolds number are
    // correctly input in the input file, so that nu is always
    // mu/rho. That's why the next line is always commented and
    // the program relies on the value given in "reynolds_number". */
    // nu = mu/rho;
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
