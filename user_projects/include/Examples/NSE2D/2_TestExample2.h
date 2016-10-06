// Navier-Stokes problem with sine and cosine functions
//

// some variables from user input
double REYNOLDS_number;
double USER_parameter1;
double USER_parameter2;
double USER_parameter3;

void ExampleFile()
{
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = OSEEN_PROBLEM;
    Output::info<1>("Example", "SinCos.h with INTERNAL_PROBLEM_IDENTITY ",
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

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void transform(const int BdComp, const double Param, double& x, double& y,
               double& nx, double& ny)
{
  switch(BdComp)
  {
    case 0:
      x = Param;
      y = 0.;
      nx = 0.;
      ny = -1.;
      break;
    case 1:
      x = 1.;
      y = Param;
      nx = 1.;
      ny = 0.;
      break;
    case 2:
      x = 1. - Param;
      y = 1.;
      nx = 0.;
      ny = 1.;
      break;
    case 3:
      x = 0.;
      y = 1. - Param;
      nx = -1.;
      ny = 0.;
      break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = 1./TDatabase::ParamDB->RE_NR;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double u1[4];
  double u2[4];
  double p[4];
  ExactU1(x, y, u1);
  ExactU2(x, y, u2);
  ExactP(x, y, p);
  if(cond == DIRICHLET)
    value = u1[0];
  else
  {
    // NEUMANN
    value = nu * (nx * u1[1] + ny * u1[2]) - p[0] * nx;
    if(lt == 1)
    {
      value *= 0.5;
      value += 0.5 * nu * (nx * u1[1] + ny * u2[1]);
    }
  }
  return;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
   // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
   const int lt = TDatabase::ParamDB->LAPLACETYPE;
   const double nu = 1./TDatabase::ParamDB->RE_NR;
   // find out boundary condition at the evaluation point on the boundary
   BoundCond cond;
   BoundCondition(BdComp, Param, cond);
   // find coordinates and normal of evaluation point on the boundary
   double x, y, nx, ny;
   transform(BdComp, Param, x, y, nx, ny);
   // evaluate the exact solution at the given point
   double u1[4];
   double u2[4];
   double p[4];
   ExactU1(x, y, u1);
   ExactU2(x, y, u2);
   ExactP(x, y, p);
   if(cond == DIRICHLET)
     value = u2[0];
   else
   {
     // NEUMANN
     value = nu * (nx * u2[1] + ny * u2[2]) - p[0] * ny;
     if(lt == 1)
     {
       value *= 0.5;
       value += 0.5 * nu * (nx * u1[2] + ny * u2[2]);
     }
   }
   return;
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
    coeff[1] = -nu*val1[3] + val3[1]; // f1 non dimensional case
    coeff[2] = -nu*val2[3] + val3[2]; // f2 non dimensional case

    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
    {
      coeff[1] += val1[0]*val1[1] + val2[0]*val1[2]; // f1 nondimensional
      coeff[2] += val1[0]*val2[1] + val2[0]*val2[2]; // f2 nondimensional
    }

    coeff[3] = rho;     // rho
    coeff[4] = mu;      // mu

  }
}
