/**
 * @file Batch crystallizer example - flow only.
 * Experiment was reported by Viktoria Wiedmeyer (VW) in 2017.
 */

bool TIME_DEPENDENT = false;

namespace FluidProperties
{
double eta = 0.0014;     // ( kg /(m*s) ) the dynamic viscosity, (of a Kalialaun solution, VW)
double rho = 1050; // ( kg / m^3  ) the density (of a Kalialaun solution, VW)

double u_infty = 1;    // (m/s) the characteristic velocity of the fluid
double l_infty = 1;    // (m) the characteristic length scale of the tube

double r_in = 0.01;   //cm the radius of the inlet
double r_out = 0.075; //cm the radius of the outlet

double mass_flow_rate = 0; // (kg/s) the mass flow rate at in- and outflow
double u_avg_in  = 0; //0.047157;   //m/s
double u_max_in  = 0; //m/s, assuming HP-inflow
double u_avg_out  = 0; //m/s

bool gravity;

void set_gravity(bool val){gravity = val;}
void set_mass_flow_rate(double mfr)
{
	mass_flow_rate = mfr/3600; // (kg/s) the mass flow rate at in- and outflow
	u_avg_in = mass_flow_rate / (rho * M_PI * r_in * r_in ); //0.047157;   //m/s
	u_max_in = 2*u_avg_in; //m/s, assuming HP-inflow
	u_avg_out = mass_flow_rate / (rho * M_PI * r_out*r_out ); //m/s
}

// note: in the coefficients function the de-dimensionalized diffusion
// coefficient will be calculated as:
//      eps = (eta/rho) / (u_infty*l_infty);
}

void ExampleFile()
{
  Output::info<1>("EXAMPLE"," WiedmeyerBatchCrystallizer.h");
  std::string gravity_string = FluidProperties::gravity ? "gravity enabled" : "gravity disabled";
  Output::info<1>("EXAMPLE"," With ", gravity_string,
		  " and mass flow rate ", FluidProperties::mass_flow_rate * 3600, " kg/h.");
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;
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
		value = u_max_in * (1 - (r*r)/(R*R));
	      if( TIME_DEPENDENT )
	      {
	        double t = TDatabase::TimeDB->CURRENTTIME;
	        if(t < 1) //within first second of the simulated time
	        	//multiply inflow with t ("anstroemen")
	        	value  *=t;
	      }
	}
	else if (determine_boundary_part(x,y,z) == BoundaryPart::WALL)	//no-slip
		value = 0;
	else //TOP - outflow
		value = u_avg_out; // constant velocity
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
    coeff[3] = FluidProperties::gravity ? -9.81 : 0; // f3 - gravity forcing, if enabled
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

// In case this should be instationary.
void InitialU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}
