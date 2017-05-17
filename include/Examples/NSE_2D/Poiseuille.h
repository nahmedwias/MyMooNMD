/** Navier-Stokes problem, Poiseuille-Problem
 *
 * u(x,y) = u_max*(1 - t^2), t is boundary parameter
 * p(x,y) = ??
 *
 * Inflow boundary is identified as 3.
 *
 * TODO This example is adapted to go along with ASA_crystallizer example.
 *  - diffusion coefficient?
 *  - exact solution, pressure solution?
 *  - do we need a right hand side?
 */


int VELOCITY_CODE=0; //can be controlled byu input parameter "velocity_code"
double u_max[4] = {0.102958, 0.1825, 0.241916, 0.26738}; // maximum inflow velocity (m/s)
double diffusion = 0.001; // diffusion coefficient - unit of your choice

void ExampleFile()
{
  Output::print<1>("Example: Poiseuille.h (adapted to ASA crystallizer example)");
  Output::dash("with velocity_code ", VELOCITY_CODE, " and u_max ", u_max[VELOCITY_CODE]);
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
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
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  switch(BdComp)
  {
    case 0: cond = DIRICHLET;
            break;
    case 1: cond = NEUMANN;
            break;
    case 2: cond = DIRICHLET;
            break;
    case 3: cond = DIRICHLET;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double u = u_max[VELOCITY_CODE];

  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 0;
            break;
    case 2: value = 0;
            break;
    case 3: value = u * (4 * Param*(1-Param));;
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

void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  for(int i=0;i<n_points;i++)
  {
    double * coeff = coeffs[i];

    coeff[0] = diffusion; //nu
    coeff[1] = 0; //f1
    coeff[2] = 0; //f2
  }
}
