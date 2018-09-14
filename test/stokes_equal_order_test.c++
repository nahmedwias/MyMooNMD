/**
 * @brief A test program for the solving of NSE2D problems.
 *
 * This serves as a test for the solving of NSE2D problems. It is intended to
 * perform NSE2D calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
 * So far only one such test is implemented.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes wrong)
 * the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then this tests
 * shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 *
 * @author Naveed, Ulrich, Clemens
 *
 */

#include <Domain.h>
#include <Database.h>
#include <ParameterDatabase.h>
#include <FEDatabase2D.h>
#include <NSE2D.h>
#include <Example_NSE2D.h>
#include <Multigrid.h>
#include <Chrono.h>
#include <algorithm>
#include "LocalAssembling.h"

double accuracy = 1e-6;

void compare(const NSE2D& stokes, std::array<double, int(5)> errors)
{
  std::array<double, int(5)> computed_errors = stokes.get_errors();
  
  // check the L2-error of the velcoity
  if( fabs(computed_errors[0]-errors[0]) > accuracy )
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the L2-error of the divergence of the velocity
  if( fabs(computed_errors[1] - errors[1]) > accuracy )
  {
    ErrThrow("L2 norm of divergence: ", std::setprecision(14), computed_errors[1], "  ", errors[1]);
  }
  // check the H1-error of the velcoity
  if( fabs(computed_errors[2] - errors[2]) > accuracy )
  {
    ErrThrow("H1 norm of velocity: ", computed_errors[2], "  ", errors[2]);
  }
  // check the L2-error of the pressure
  if( fabs(computed_errors[3] - errors[3]) > accuracy)
  {
    ErrThrow("L2 norm of pressure: ", computed_errors[3], "  ", errors[3]);
  }
  // check the H1-error of the pressure
  if(fabs(computed_errors[4] - errors[4]) > accuracy )
  {
    ErrThrow("H1 norm of pressure: ", computed_errors[4], "  ", errors[4]);
  }
}

void compute(TDomain &domain, ParameterDatabase& db,
             std::array<double, int(5)> errors)
{
  NSE2D stokes(domain, db);
  stokes.assemble();
  stokes.solve();
  stokes.output();
  // compare now the errors
  compare(stokes, errors);
}

void check(TDomain &domain, ParameterDatabase db,
           int velocity_order, int laplace_type,
           int nonlinear_form,
           std::array<double, int(5)> errors)
{
  Output::print("\n\nCalling check with velocity_order=", velocity_order,
                ", laplace_type=", laplace_type, ", and nonlinear_form=",
                nonlinear_form);
  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(Multigrid::default_multigrid_database());
  db["problem_type"] = 5;
  db["solver_type"] = "direct";
  db["iterative_solver_type"] = "fgmres";
  db["residual_tolerance"] = 1.e-12;
  db["preconditioner"] = "sor";
  db["max_n_iterations"] = 2000;
  
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = velocity_order;
  TDatabase::ParamDB->LAPLACETYPE = laplace_type;
  TDatabase::ParamDB->NSE_NONLINEAR_FORM = nonlinear_form;
  
  Chrono timer;
  compute(domain, db, errors);
  timer.restart_and_print("nse2d direct solver,                    velocity "
                          + std::to_string(velocity_order));
  
  // we have to reset the space codes because they are changed in nse2d
  TDatabase::ParamDB->PRESSURE_SPACE = velocity_order;
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  
  db["solver_type"] = "iterative";
  compute(domain, db, errors);
  timer.restart_and_print("nse2d fgmres(lsc preconditioner),       velocity "
                          + std::to_string(velocity_order));
  
  // we have to reset the space codes because they are changed in nse2d
  TDatabase::ParamDB->PRESSURE_SPACE = velocity_order;
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  

  db["preconditioner"] = "multigrid";
  db["multigrid_n_levels"] = db["refinement_n_initial_steps"].get<size_t>();
  //choose smoother on fine grid according to element
  std::vector<int> disc_p = {12,13,14,15,22,23,24};
  if(std::find(disc_p.begin(), disc_p.end(), velocity_order) != disc_p.end())
    db["multigrid_smoother"] = "cell_vanka_store";
  else
    db["multigrid_smoother"] = "patch_vanka_store";

  db["multigrid_type"] = "standard";
  db["multigrid_smoother_coarse"] = "direct_solve";
  db["multigrid_n_pre_smooth"] = 0;
  db["multigrid_n_post_smooth"] = 1;
  db["multigrid_correction_damp_factor"] = 1.0;
  db["multigrid_vanka_damp_factor"] = 1.0;
  compute(domain, db, errors);
  timer.restart_and_print("nse2d fgmres(multigrid preconditioner), velocity "
                          + std::to_string(velocity_order));
}

template <int laplace_type>
void check_one_element(TDomain& domain, ParameterDatabase db,
                       int velocity_order,
                       std::array<double, int(5)> errors)
{
  int nonlinear_form = 0;
  check(domain, db, velocity_order, laplace_type, nonlinear_form, errors);
}

void pspg_on_triangles()
{
  TDatabase Database;
#ifdef __2D__
  TFEDatabase2D FEDatabase;
#else
  TFEDatabase3D FEDatabase;
#endif
  
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(NSE2D::default_NSE_database());
  db.merge(TDomain::default_domain_parameters());

  db["problem_type"].set<size_t>(3); // Stokes
  db["example"] = 2;
  db["refinement_n_initial_steps"] = 2;
  db["boundary_file"] = "Default_UnitSquare";
  db["geo_file"] = "TwoTriangles";
  
  db["reynolds_number"] = 1;
  TDatabase::ParamDB->FLOW_PROBLEM_TYPE=3;
  db["space_discretization_type"] = "pspg";
  db["laplace_type_deformation"] = false;
  TDatabase::ParamDB->NSTYPE = 14;
  TDatabase::ParamDB->LAPLACETYPE = 0;
  
  // possibly parameters in the database
  check_parameters_consistency_NSE(db);
  
  // default construct a domain object
  TDomain domain(db);
  // refine grid
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(unsigned int i=0; i < n_ref; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, int(5)> errors;
  
  //=========================================================================
  Output::print<1>("\nTesting the P1/P1 elements");
  errors = {{ 0.14425380968049, 0.89062163325473, 1.5272900521314,
              1.3353534209894, 8.0717803059975 }};
  check_one_element<0>(domain, db, 1, errors);
  
  //=========================================================================
  Output::print<1>("\nTesting the P2/P2 elements");
  errors = {{ 0.0094075018094777, 0.092498252374659, 0.16396168888847,
              0.10814807138916, 1.9898182882992 }};
  check_one_element<0>(domain, db, 2, errors);

  //=========================================================================
  Output::print<1>("\nTesting the P3/P3 elements");
  errors = {{ 0.0013655488005913, 0.012260506129698, 0.017792417110209,
              0.022236536026839, 0.3674535115527 }};
  check_one_element<0>(domain, db, 3, errors);

  //=========================================================================
  Output::print<1>("\nTesting the P4/P4 elements");
  errors = {{ 5.7237636133624e-05, 0.0015288367247389, 0.0016535091698868,
              0.0017097099552422, 0.072004774733194 }};
  check_one_element<0>(domain, db, 4, errors);
  
  //=========================================================================
  
  db["geo_file"] = "UnitSquare";
}

void pspg_on_quads()
{
  TDatabase Database;
#ifdef __2D__
  TFEDatabase2D FEDatabase;
#else
  TFEDatabase3D FEDatabase;
#endif
  
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(NSE2D::default_NSE_database());
  db.merge(TDomain::default_domain_parameters());
  db["problem_type"].set<size_t>(3); // Stokes
  db["example"] = 2;
  db["refinement_n_initial_steps"] = 2;
  db["boundary_file"] = "Default_UnitSquare";
  db["geo_file"] = "UnitSquare";
  db["reynolds_number"] = 1;
  TDatabase::ParamDB->FLOW_PROBLEM_TYPE=3;
  db["space_discretization_type"] = "pspg";
  db["laplace_type_deformation"] = false;
  TDatabase::ParamDB->NSTYPE = 14;
  TDatabase::ParamDB->LAPLACETYPE = 0;
  
  // possibly parameters in the database
  check_parameters_consistency_NSE(db);
  
  // default construct a domain object
  TDomain domain(db);
  // refine grid
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(unsigned int i=0; i < n_ref; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, int(5)> errors;
  
  //=========================================================================
  Output::print<1>("\nTesting the Q1/Q1 elements");
  errors = {{ 0.14770518128781, 0.5612829754172, 1.1794990495543,
              1.3441843485325, 8.0847736329225 }};
  check_one_element<0>(domain, db, 1, errors);
  
  //=========================================================================
  Output::print<1>("\nTesting the Q2/Q2 elements");
  errors = {{ 0.0067110682559855, 0.065269023656333, 0.11630741869309,
              0.092921606881003, 1.5925769236787 }};
  check_one_element<0>(domain, db, 2, errors);

  //=========================================================================
  Output::print<1>("\nTesting the Q3/Q3 elements");
  errors = {{ 0.0012077200634523, 0.011283264757702, 0.015848841146417,
              0.023611185292096, 0.33408390051687 }};
  check_one_element<0>(domain, db, 3, errors);

  //=========================================================================
  Output::print<1>("\nTesting the Q4/Q4 elements");
  errors = {{ 2.9945160214221e-05, 0.00041913940235019, 0.00061350090960876,
              0.00083988300949934, 0.023863686042054 }};
  check_one_element<0>(domain, db, 4, errors);
}

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  pspg_on_triangles();
  pspg_on_quads();
  return 0;
}
