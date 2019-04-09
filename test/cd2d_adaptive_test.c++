/**
 * @brief A test program for the solving of CD2D problems with adaptively
 * refined grids.
 *
 */
#include "Domain.h"
#include "Database.h"
#include "FEDatabase2D.h"
#include "ConvectionDiffusion.h"
#include "CDErrorEstimator.h"
#include "RefinementStrategy.h"
#include "LoopInfo.h"
#include "Chrono.h"
#include "LocalAssembling.h"
#include "Multigrid.h"


bool write_PS_and_VTU_files = true; // turn this on for debugging purposes

// compare the computed errors in the Darcy2D object with the given ones in 
// the array
void compareErrors(const ConvectionDiffusion<2>& cd2d,
                   std::array<double, 4> errors)
{
  const double eps = 1e-11;
  
  // check the errors
  if( std::abs(cd2d.get_L2_error() - errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 error not correct. ",
             cd2d.get_L2_error() - errors[0]);
  }
  if( std::abs(cd2d.get_H1_semi_error() - errors[1]) > eps )
  {
    ErrThrow("Program 1: H1-semi error not correct. ",
             cd2d.get_H1_semi_error() - errors[1]);
  }
  if( std::abs(cd2d.get_SD_error() - errors[2]) > eps )
  {
    ErrThrow("Program 1: sd error not correct.",
             cd2d.get_SD_error() - errors[2]);
  }
  if( std::abs(cd2d.get_L_inf_error() - errors[3]) > eps )
  {
    ErrThrow("Program 1: L_inf error not correct.",
             cd2d.get_L_inf_error() - errors[3]);
  }
}

// Here the actual computations take place
void check_cd2d(ParameterDatabase& db, int element_code,
                std::array<double, 4> errors)
{
  // default construct a domain object
  TDomain domain(db);
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(db);
  
  TDatabase::ParamDB->ANSATZ_ORDER = element_code;
  int estimator_type = 1;
  bool conform_grid = true;
  CDErrorEstimator<2> estimator(estimator_type, conform_grid,
                                Parameter(db["space_discretization_type"]));
  RefinementStrategy<2> refinementStrategy(db);
  LoopInfo loop_info("adaptive", true, true, 1);
  BlockVector values;
  // for this test we do a fixed number of adaptive refinement steps
  size_t n_adaptive_steps = domain.get_database()["refinement_max_n_adaptive_steps"];
  for(size_t curr_level = 0; curr_level < n_adaptive_steps; ++curr_level)
  {
    Output::print("\nadaptive loop ", curr_level);
    std::ostringstream ostr;
    ostr << "parmoon" << std::setw(2) << std::setfill('0') << curr_level;
    db["output_basename"].set(ostr.str(), false);
    
    ConvectionDiffusion<2> cd2d(domain, db);
    cd2d.assemble();
    cd2d.solve();
    estimator.estimate(cd2d.get_example(), cd2d.get_function());
    estimator.info();
    
    ConvectionDiffusion<2>::FEFunction estimated_error(estimator, values);
    cd2d.add_to_output(&estimated_error);
    cd2d.output();
    
    if(curr_level+1 != n_adaptive_steps)
    {
      loop_info.print(curr_level,
                      estimator.get_estimated_global_error()[estimator_type]);
      refinementStrategy.apply_estimator(estimator);
      domain.RefineByRefinementStrategy(refinementStrategy, true);
    }
    else
    {
      loop_info.finish(curr_level,
                       estimator.get_estimated_global_error()[estimator_type]);
      // compare computed with given errors
      compareErrors(cd2d, errors); // throws upon a difference
    }
    if(write_PS_and_VTU_files)
    {
      std::string name_after = std::string("domain_after")
                               + std::to_string(curr_level) + ".ps";
      domain.PS(name_after.c_str(), It_Finest, 0);
    }
  }
}

// indicate where to refine when using TDomain::RefineByIndicator.
// Cells which intersect with the zero set of this indicator function are 
// regularly refined.
void indicator(double x, double y, double* values)
{
  values[0] = 0.5*x*x + y - 0.5;
}

void test_indicator()
{
  Output::print("testing TDomain::RefineByIndicator");
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(TDomain::default_domain_parameters());
  TDomain domain(db);
  size_t n_adaptive_steps = 8;
  for(size_t curr_level = 0; curr_level < n_adaptive_steps; ++curr_level)
  {
    domain.RefineByIndicator(indicator, true);
    if(write_PS_and_VTU_files)
    {
      std::string name = std::string("domain_indicator_after")
                       + std::to_string(curr_level) + ".ps";
      domain.PS(name.c_str(), It_Finest, 0);
    }
  }
  // finally, solve a simple problem on this domain, only to check
  ConvectionDiffusion<2> cd2d(domain, db);
  cd2d.assemble();
  cd2d.solve();
  cd2d.output();
  
  std::array<double, 4>errors = {{ 0.0012798228687561, 0.040891950610384,
                                   0.040891950610384, 0.0073846377470421 }};
  compareErrors(cd2d, errors);
}


void tests_on_quads(ParameterDatabase& db)
{
  std::array<double, 4> errors;
    
  Output::print("\nstarting with Q1");
  errors = {{ 0.030136750989447, 0.50136473756928, 0.50136473756928, 
              0.086620472046938 }};
  check_cd2d(db, 1, errors);
  
  Output::print("\nstarting with Q2");
  errors = {{ 0.0019271846609828, 0.050971302154831, 0.050971302154831, 
              0.0033739148851507 }};
  check_cd2d(db, 2, errors);
}

void tests_on_triangles(ParameterDatabase& db)
{ 
  std::array<double, 4> errors;

  Output::print("\nstarting with P1");
  errors = {{ 0.47784752238546, 4.1491718167438, 0.14237859736407,
              1.0217899951186 }};
  check_cd2d(db, 1, errors);
  
  Output::print("\nstarting with P2");
  errors = {{ 0.47843700715078, 4.1949720830524, 0.13830367367314,
              1.1101098542477 }};
  check_cd2d(db, 2, errors);
}


// =======================================================================
// main program
// =======================================================================
int main(int, char**)
{
  Output::setVerbosity(2);
  //  declaration of databases
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  test_indicator();
  
  size_t nRefinements = 2; // initial refinements prior to adaptive refinements
  
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(Solver<>::default_solver_database());
  db.merge(Example2D::default_example_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(RefinementStrategy<2>::default_refinement_strategy_database());
  db.merge(TDomain::default_domain_parameters());
  db["refinement_max_n_adaptive_steps"] = 8;
  db["problem_type"] = 0; // problem type is not needed
  TDatabase::ParamDB->PE_NR = 1000;
  db["example"] = 1; // two_interior_layers
  db["residual_tolerance"] = 1.0e-13;
  
  db["output_compute_errors"] = true;
  
  db["boundary_file"] = "Default_UnitSquare";
  db["geo_file"] = "UnitSquare";
  db["refinement_n_initial_steps"] = nRefinements;
  db.add("space_discretization_type", "galerkin", "");
  if(write_PS_and_VTU_files)
    db["output_write_vtu"] = true;
  
  db["refinement_strategy"] = 0;
  db["current_estimator"] = 1;
  
  Chrono time; // measure the time spent for each solver
  
  Output::print("\n\n ----------- direct solver -----------\n");
  db["solver_type"] = "direct";
  db["geo_file"] = "UnitSquare";
  //tests_on_quads(db);
  db["geo_file"] = "TwoTriangles";
  tests_on_triangles(db);
  time.restart_and_print("all tests, direct solver");
  
  
  Output::print("\n\n --------- fgmres+multigrid solver ---------\n");
  db.merge(Multigrid::default_multigrid_database());
  db["iterative_solver_type"] = "fgmres";
  db["preconditioner"] = "multigrid";
  db["multigrid_n_levels"] = 2;
  db["multigrid_cycle_type"] = "W";
  db["multigrid_smoother"] = "jacobi";
  db["multigrid_smoother_coarse"] = "direct_solve";
  db["multigrid_correction_damp_factor"] = 1.0;
  db["multigrid_n_pre_smooth"] = 2;
  db["multigrid_n_post_smooth"] = 2;
  db["geo_file"] = "UnitSquare";
  //tests_on_quads(db);
  db["geo_file"] = "TwoTriangles";
  tests_on_triangles(db);
  time.restart_and_print("all tests, fgmres solver and multigrid "
                         "preconditioner");
  
  
  return 0;
}
