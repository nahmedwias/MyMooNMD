/**
 * @brief A test program for the solution of Time dependent Navier-Stokes
 * 
 * This test program is intended to check whether different solvers for TNSE2D
 * @author 
 * @history 14.05.2016
 * 
 */
#include <Time_NSE2D.h>
#include <Database.h>
#include <Domain.h>
#include <Multigrid.h> // newly implemented by clemens

void compare(const Time_NSE2D &tnse2d, std::array<double, int(4)> errors, double tol)
{
  std::array<double, int(6)> computed_errors;
  computed_errors = tnse2d.get_errors();
  
  // check the L2-error of the velcoity
  if( fabs(computed_errors[0]-errors[0]) > tol )
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the H1-error of the velcoity
  if( fabs(computed_errors[1] - errors[1]) > tol )
  {
    ErrThrow("H1 norm of velocity: ", computed_errors[1], "  ", errors[1]);
  }    
  // check the L2-error of the pressure
  if( fabs(computed_errors[2] - errors[2]) > tol)
  {
    ErrThrow("L2 norm of pressure: ", computed_errors[2], "  ", errors[2]);
  }
  // check the H1-error of the pressure
  if(fabs(computed_errors[3] - errors[3]) > tol )
  {
    ErrThrow("H1 norm of pressure: ", computed_errors[3], "  ", errors[3]);
  }
}

void check(TDomain& domain, int example, int velocity_order, int pressure_order, 
           int nstype, int laplace_type, int nonlinear_form, int time_disc, 
           std::array<std::array<double, int(4)>,4> errors)
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  
  db["problem_type"] = 5;
  Output::print("******* Check ",example, " ",velocity_order, " ",
		pressure_order, " ",nstype," *******"); // and more stuff

  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  TDatabase::ParamDB->NSTYPE = nstype;
  
}

void set_solver_globals(std::string solver_name, ParameterDatabase& db)
{
  db["solver_type"] = std::string("iterative");
  db["direct_solver_type"] = std::string("umfpack");
  db["iterative_solver_type"] = std::string("fgmres");
  db["preconditioner"] = std::string("no_preconditioner");
  db["residual_tolerance"] = 1.0e-13;
  
  if (solver_name.compare("lsc") == 0)
  {
    db["preconditioner"] = "least_squares_commutator";
    db["nonlinloop_epsilon"] = 1e-12;
    // just to not distract 'NSE3D::check_parameters'
    TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE = 20;
  }
  else if (solver_name.compare("multigrid") == 0)
  {
    db.merge(Multigrid::default_multigrid_database());
    db["preconditioner"] = "multigrid";
    db["refinement_n_initial_steps"] = 1;
    //control nonlinear loop
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_maxit"] = 5;
    // New multigrid parameters
    db["multigrid_n_levels"] = 2;
    db["multigrid_cycle_type"] = "V";
    db["multigrid_smoother"] = "nodal_vanka";
    db["multigrid_smoother_coarse"] = "nodal_vanka";
    db["multigrid_correction_damp_factor"] = 0.8;
    db["multigrid_n_pre_smooth"] = 2;
    db["multigrid_n_post_smooth"] = 2;
    db["multigrid_coarse_residual"] = 1.0e-1;
    db["multigrid_coarse_max_n_iterations"] = 5;
    db["multigrid_vanka_damp_factor"]=0.7;

  }
  else
  {
    throw std::runtime_error("Unknown solver for Time_NSE2D problem!");
  }
}

double get_tolerance(std::string solver_name)
{
  if(solver_name.compare("umfpack") == 0)
    return 1e-9;
  if(solver_name.compare("lsc") == 0)
    return 1e-9;
  if(solver_name.compare("multigrid") == 0)
    return 1e-9;
}

int main(int argc, char* argv[])
{
  
}
