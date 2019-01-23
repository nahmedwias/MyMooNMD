#include <Domain.h>
#include <Database.h>
#include <FESpace3D.h>
#include <FEDatabase3D.h>
#include <FEFunction3D.h>
#include <LinesEval.h>

void analytic_function(double x, double y, double z, double *values)
{
  values[0] = x + y + z;
}

void all_dirichlet_boundary_condition(double, double, double, BoundCond & bc)
{
  bc = DIRICHLET;
}


bool additional_tests(TDomain & domain)
{
  auto lines_eval_db = LinesEval<3>::default_lineseval_parameters();
  lines_eval_db["line_direction"]  = 1;
  lines_eval_db["line_position"] = std::vector<double>({0.,0.,0., 0.,0.5,0.});
  lines_eval_db["line_refinement"] = 1;

  
  try
  {
    lines_eval_db["line_position"] = std::vector<double>({0.,0.,0., 0.5,0.});
    LinesEval<3> lines(domain, lines_eval_db);
    Output::print("there should be an error with invalid 'line_position'");
    return false;
  }
  catch(...)
  {
    // ok, this should have been an error
    lines_eval_db["line_position"] = std::vector<double>({0.,0.,0., 0.,0.5,0.});
  }
  
  
  try
  {
    lines_eval_db["line_direction"] = 5;
    //LinesEval lines(domain, lines_eval_db);
    Output::print("there should be an error with invalid 'line_direction'");
    return false;
  }
  catch(...)
  {
    // ok, this should have been an error
    lines_eval_db["line_direction"] = 1;
  }
  
  
  
  return true;
}

int main()
{
  // declaration of databases
  TDatabase Database;
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(TDomain::default_domain_parameters());
  
  db["refinement_n_initial_steps"] = 2;
  Output::setVerbosity(3);
  
#ifdef __2D__
  ErrThrow("only implemented in 3D");
  //TFEDatabase2D FEDatabase;
  //db["boundary_file"] = "Default_UnitSquare";
  //db["geo_file"] = "TwoTriangles";
#else // 3D
  TFEDatabase3D FEDatabase;
  db["boundary_file"] = "Default_UnitCube";
  db["geo_file"] = "Default_UnitCube_Tetra";
#endif // 2D
  
  TDomain domain(db);
  
  // refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  
  if(!additional_tests(domain))
  {
    Output::print("additional_tests failed");
    return 1;
  }
  
  auto coll = domain.get_grid_collections().front();
  TFESpace3D fespace(coll, "dummy", "", all_dirichlet_boundary_condition, 2);
  int length = fespace.GetN_DegreesOfFreedom();
  std::vector<double> entries(length, 0);
  TFEFunction3D fefunction(&fespace, "testfunction", "", &entries[0], length);
  fefunction.Interpolate(analytic_function);

  auto lines_eval_db = LinesEval<3>::default_lineseval_parameters();
//  lines_eval_db["position_file"];
  lines_eval_db["line_direction"]  = 1;
  std::vector<double> posi = {0.0,0.0,0.0, 0.5,0.0,0.5};
  lines_eval_db["line_position"] = posi;
  lines_eval_db["line_refinement"] = 1;
  db.add_nested_database(lines_eval_db);

  LinesEval<3> lines(domain, lines_eval_db);
  LinesEval<3> lines_test(domain, db);
  
  // code here
  double z_min = 0.;
  double z_max = 1.;
  double avg_0 = 0.5;
  double avg_1 = 1.5;
  
  double z_min_c = 0.;
  double z_max_c = 1.;
  double avg_0_c = lines.GetLine(0).space_average_value(fefunction);
  double avg_1_c = lines.GetLine(1).space_average_value(fefunction);
  
  if(std::abs(z_min - z_min_c) > 1.e-10)
  {
    ErrThrow("wrong value computed, ", z_min, " != ", z_min_c);
    return 1;
  }
  if(std::abs(z_max - z_max_c) > 1.e-10)
  {
    ErrThrow("wrong value computed, ", z_max, " != ", z_max_c);
    return 1;
  }
    if(std::abs(avg_0 - avg_0_c) > 1.e-10)
  {
    ErrThrow("wrong value computed, ", avg_0, " != ", avg_0_c);
    return 1;
  }
    if(std::abs(avg_1 - avg_1_c) > 1.e-10)
  {
    ErrThrow("wrong value computed, ", avg_1, " != ", avg_1_c);
    return 1;
  }
}
