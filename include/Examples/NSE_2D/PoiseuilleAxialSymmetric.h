/** Navier-Stokes problem, Poiseuille-Problem in axisymmetric formulation.
 * Note that the geometry must be as such, that the boundary parts are
 * numbered as:
 *
 * Wall boundary is 0.
 * Outflow boundary is 1.
 * Spurious boundary is 2. (Symmetry axis)
 * Inflow boundary is 3.
 *
 * AND that x_cart=z_cyl, y_cart=r_cyl (lying cylinder)
 *
 */


int VELOCITY_CODE=0; //can be controlled by input parameter "velocity_code"
double u_max[4] = {0.102958, 0.1825, 0.241916, 0.26738}; // maximum inflow velocity (m/s)
double diffusion = 0.0001; // diffusion coefficient - unit of your choice

void ExampleFile()
{
  Output::print<1>("Example: PoiseuilleAxialSymmetric.h (adapted to ASA crystallizer example)");
  Output::dash("with velocity_code ", VELOCITY_CODE, " and u_max ", u_max[VELOCITY_CODE]);
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactUZ(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactUR(double x, double y, double *values)
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
void boundary_conditions_axial(int BdComp, double t, BoundCond &cond)
{
  switch(BdComp)
  {
    case 0: cond = DIRICHLET;
            break;
    case 1: cond = NEUMANN;
            break;
    case 2: cond = NEUMANN;
            break;
    case 3: cond = DIRICHLET;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void boundary_conditions_radial(int BdComp, double t, BoundCond &cond)
{
  switch(BdComp)
  {
    case 0: cond = DIRICHLET; //y=R_max (wall boundary)
            break;
    case 1: cond = NEUMANN;   //outflow
            break;
    case 2: cond = DIRICHLET; //y=0 (spurios boundary)
            break;
    case 3: cond = DIRICHLET; //inflow
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void UZBoundValue(int BdComp, double Param, double &value)
{

	double u = u_max[VELOCITY_CODE];

  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value = u * (4 * 0.5 * (1 - Param) *(1- 0.5 * (1 - Param)));
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void URBoundValue(int BdComp, double Param, double &value)
{

  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 0;
            break;
    case 2: value = 0;
            break;
    case 3: value = 0;
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
