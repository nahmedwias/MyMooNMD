/**
 * @file Batch crystallizer example - flow only.
 * Experiment was reported by Viktoria Wiedmeyer (VW) in 2017.
 */


void ExampleFile()
{
  Output::info<1>("EXAMPLE"," WiedmeyerBatchCrystallizer.h");
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;
}

namespace FluidProperties
{
double eta = 0.0014;     // ( kg /(m*s) ) the dynamic viscosity, (of a Kalialaun solution, VW)
double rho = 1050; // ( kg / m^3  ) the density (of a Kalialaun solution, VW)

double u_infty = 1;    // (m/s) the characteristic velocity of the fluid
double l_infty = 1;    // (m) the characteristic length scale of the tube

double mass_flux = 56.0/3600; // (kg/s) the mass flux at in- and outflow (of the experiment, as reported by VW)

double u_avg_in = 0.047157;   //m/s
double u_max_in = 2*u_avg_in; //m/s, assuming HP-inflow
double u_avg_out = 0.00083834702534; //m/s

// note: in the coefficients function the de-dimensionalized diffusion
// coefficient will be calculated as:
//      eps = (eta/rho) / (u_infty*l_infty);
}

//Boundary parts of the Geometry, numbered "bottom up"
enum class BoundaryPart {BOTTOM, WALL, TOP};

//Find out which boundary part a boundary point lies on. Measures are in m.
BoundaryPart determine_boundary_part(double x, double y, double z)
{
	double tol = 1e-8;
	//catch first those values that lie on the wall parts of top and bottom
	if(fabs(sqrt(x*x + y*y) - 0.01) < tol || fabs(sqrt(x*x + y*y) - 0.075) < tol)
	{
		return BoundaryPart::WALL;
	}
	if (fabs(z-0) < tol )
	{
		return BoundaryPart::BOTTOM;
	}
	else if (fabs(z-0.5) < tol )
	{
		return BoundaryPart::TOP;
	}
	else //everything "complicated" is on the wall boundary
		return BoundaryPart::WALL;
}

// kind of boundary condition
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
	if (determine_boundary_part(x,y,z) == BoundaryPart::BOTTOM)
		cond = DIRICHLET;
	else if (determine_boundary_part(x,y,z) == BoundaryPart::WALL)
		cond = DIRICHLET;
	else //TOP
		cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
	value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
	using namespace FluidProperties;
	if (determine_boundary_part(x,y,z) == BoundaryPart::BOTTOM) //inflow
	{//HP inflow profile
		double R = 0.01;
		double r = sqrt(x*x + y*y);
		value = u_max_in * (1 - r*r/R*R);
	}
	else if (determine_boundary_part(x,y,z) == BoundaryPart::WALL)	//no-slip
		value = 0;
	else //TOP - outflow
		value = u_avg_out; //eigene Berechnung - konstante Velo ueber Querschnitt
}

// ========================================================================
// coefficients for Stokes form
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  using namespace FluidProperties;

  static double eps = (eta/rho) / (u_infty*l_infty);
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // f3		TODO Put some gracity forcing here?
    coeff[4] = 0; // g
  }
}

// ========================================================================
// exact solution is unknown, everything set to 0.
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}
