/**
 * @file Implementation of the twisted pipe flow example,
 * as declared in twisted_pipe_flow.h.
 */


#include <Database.h>
#include <MooNMD_Io.h>
#include <CoiledPipe.h>

#ifdef _MPI
#include <mpi.h>
#endif

namespace twisted_pipe_flow{ //must be included within the same namespace as in Exa
  #include <twisted_pipe_flow.h>
};

namespace FluidProperties
{
double eta = 1.11e-5;     // ( kg /(cm*s) ) the dynamic viscosity, (of a Kalialaun solution, as reported by V.Wiedmeyer)
double rho = 1100 * 1e-6; // ( kg / cm^3  ) the density (of a Kalialaun solution, as reported by V.Wiedmeyer)

double u_infty = 1;    // (cm/s) the characteristic velocity of the fluid
double l_infty = 1;    // (cm) the characteristic length scale of the tube

double vol_flux = 7.2; // (cm^3/s) the volume flux at inflow (of the experiment, as reported by V.Wiedmeyer)

// note: in the coefficients function the de-dimensionalized diffusion
// coefficient will be calculated as:
//      eps = (eta/rho) / (u_infty*l_infty);
}

void twisted_pipe_flow::ExampleFile()
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &my_rank);
#else
  int my_rank = 0;
#endif

  // set global parameters which are too deeply rooted in the
  // code to be easily removed at the moment
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1356; //FIXME Remove this from the entire code!
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;

  if(my_rank == 0)
    Output::print(" > Example: twisted_pipe_flow.h.");

  using namespace FluidProperties;

  double eps = (eta/rho) / (u_infty*l_infty);
  if(my_rank == 0)
    Output::print(" > with diffusion coefficient eps = ", eps);

  double inflow_velo[5];
  ExactU1(- CoiledPipe::GeoConsts::l_inflow, 0, 0, inflow_velo);
  if(my_rank == 0)
    Output::print(" > and center inflow velocity u1: ", inflow_velo[0]);

}


// ========================================================================
// Exact solutions (unknown, set to zero in most parts of the domain)
// ========================================================================
void twisted_pipe_flow::ExactU1(double x, double y,  double z, double *values)
{

  CoiledPipe::PIPE_PIECE piece = CoiledPipe::in_which_piece_is_this_point(x,y,z);
  switch(piece)
  {
    case CoiledPipe::PIPE_PIECE::INFLOW_FACE:
    case CoiledPipe::PIPE_PIECE::INFLOW_PIECE:
    {
      using namespace FluidProperties;
      using namespace CoiledPipe::GeoConsts;
      //in straight inflow part set Hagen-Poiseuille Profile as exact solution

      double cross_section = Pi * tube_radius * tube_radius;
      double u_avr = vol_flux / cross_section; // reminder: u_max = 2*u_average for Hagen-Poiseuille flow

      double r2 = y * y + z * z ; //inflow piece aligned along x-axis!
      double R2 = tube_radius *tube_radius;

      // dimensionalized flow in x direction:
      double u_1 = 2 * u_avr * ( 1 -  r2 / R2);
      double u_1_y = -2 * 2 * u_avr * y / R2; //y derivative
      double u_1_z = -2 * 2 * u_avr * z / R2;; //z derivative

      if(true) //TODO Only for time-dependent case - find a way to find out, if its time-dependent or not!
      {
        double t = TDatabase::TimeDB->CURRENTTIME;
        if(t < 1.0) //within first second of the simulated time
        {//multiply inflow with t ("anstroemen")
          u_1   *=t;
          u_1_y *=t;
          u_1_z *=t;
        }
      }

      // ...dedimensionalize and assign the results
      values[0] = u_1 / u_infty;
      values[1] = 0;
      values[2] = u_1_y / u_infty;
      values[3] = u_1_z / u_infty;
      values[4] = 0;
      break;
    }
    default:
      values[0]= 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      break;
  }
}


void twisted_pipe_flow::ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void twisted_pipe_flow::ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void twisted_pipe_flow::ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


// ========================================================================
// Initial solutions
// ========================================================================
/// Initial solution in U1 direction - Hagen-Poiseuille in inflow piece.
void twisted_pipe_flow::InitialU1(double x, double y, double z, double *values)
{
  ExactU1(x, y, z, values);
}

/// Initial solution in U2 direction. No flow, nowhere.
void twisted_pipe_flow::InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

/// Initial solution in U3 direction. No flow, nowhere.
void twisted_pipe_flow::InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
}

/// Initial pressure. No pressure, nowhere.
void twisted_pipe_flow::InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}


// ========================================================================
// Boundary conditions and values
// ========================================================================
void twisted_pipe_flow::BoundCondition(double x, double y, double z, BoundCond &cond)
{
  CoiledPipe::PIPE_PIECE piece = CoiledPipe::in_which_piece_is_this_point(x,y,z);
  switch(piece)
  {
    case CoiledPipe::PIPE_PIECE::OUTFLOW_FACE:
      cond = NEUMANN;
      break;
    default: //everywhere else:
      cond = DIRICHLET;
      break;
  }
}


void twisted_pipe_flow::U1BoundValue(double x, double y, double z, double &value)
{

  CoiledPipe::PIPE_PIECE piece = CoiledPipe::in_which_piece_is_this_point(x,y,z);
  switch(piece)
  {
    case CoiledPipe::PIPE_PIECE::INFLOW_FACE:
      double values[5];
      ExactU1(x,y,z,values);
      value = values[0];
      break;

    default:
      value =0;
      break;
  }
}

void twisted_pipe_flow::U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

void twisted_pipe_flow::U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// ========================================================================
// Coefficients for weak form: eps, f1, f2, f3, g
// ========================================================================
void twisted_pipe_flow::LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  using namespace FluidProperties;

  double eps =  (eta / rho) / (l_infty * u_infty);

  for(int i=0; i < n_points; i++)
  {
    double* coeff = coeffs[i];

    coeff[0] = eps;       // momentum diffusion
    coeff[1] = 0;         // f1
    coeff[2] = 0;         // f2
    coeff[3] = 0;         // f3
    coeff[4] = 0;         // g
  }
}
