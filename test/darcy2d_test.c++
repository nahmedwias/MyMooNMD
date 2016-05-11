/**
 * @brief A test program for the solving of Darcy2D problems.
 *
 * This serves as a test for the solving of Darcy2D problems. It is intended to
 * perform Darcy2D calculations with different examples in different setups to 
 * test a wide variety of ParMooN core functionality.
 * So far only one such test is implemented.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes 
 * wrong) the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then 
 * this tests shows you how many other program parts are affected by your 
 * changes. If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Darcy2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_CD2D.h>

#include <MainUtilities.h> //for error measuring

// compare the computed errors in the Darcy2D object with the given ones in 
// the array
void compareErrors(const Darcy2D& darcy2d, std::array<double, 5> errors)
{
  const double eps = 2e-9;
  
  // check the errors
  if( fabs(darcy2d.getL2VelocityError() - errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 velocity error not correct. ",
             darcy2d.getL2VelocityError() - errors[0]);
  }
  if( fabs(darcy2d.getL2DivergenceError() - errors[1]) > eps )
  {
    ErrThrow("Program 1: L2 velocity divergence error not correct. ",
             darcy2d.getL2DivergenceError() - errors[1]);
  }
  if( fabs(darcy2d.getH1SemiVelocityError() - errors[2]) > eps )
  {
    ErrThrow("Program 1: H1-semi velocity error not correct. ",
             darcy2d.getH1SemiVelocityError() - errors[2]);
  }
  if( fabs(darcy2d.getL2PressureError() - errors[3]) > eps )
  {
    ErrThrow("Program 1: L2 pressure error not correct.",
             darcy2d.getL2PressureError() - errors[3]);
  }
  if( fabs(darcy2d.getH1SemiPressureError() - errors[4]) > eps )
  {
    ErrThrow("Program 1: H1-semi pressure error not correct.",
             darcy2d.getH1SemiPressureError() - errors[4]);
  }
}

// Here the actual computations take place
void check_darcy2d(TDomain & domain, ParameterDatabase& db, int velocityCode, 
           std::array<double, 5> errors)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocityCode;
  // automatically choose pressure space to get inf-sup stable pair
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  
  Darcy2D darcy2d(domain, db);
  darcy2d.assemble();
  darcy2d.solve();
  darcy2d.output();
  // compare computed with given errors
  compareErrors(darcy2d, errors); // throws upon a difference
}

void tests_on_quads(unsigned int nRefinements, ParameterDatabase& db)
{
  // default construct a domain object
  TDomain domain(db);
  // the domain is initialised with default description and default
  // initial mesh
  domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");

  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }

  std::array<double, 5> errors;
    
  Output::print("\nstarting with RT0 on quads");
  errors = {{ 2.1136884064519, 23.120239110875, 32.350359820686, 
              0.30277518654981, 4.4428829381584 }};
  check_darcy2d(domain, db, 1000, errors);
  
  Output::print("\nstarting with RT1 on quads");
  errors = {{ 0.40548747487874, 4.9461605478514, 12.695872686293, 
              0.063044878850364, 1.9735667622761 }};
  check_darcy2d(domain, db, 1001, errors);
  
  Output::print("\nstarting with RT2 on quads");
  errors = {{ 0.053370236091622, 0.66178203636242, 2.7529898961154,
              0.0084084949151706, 0.43469938937457 }};
  check_darcy2d(domain, db, 1002, errors);
  
  // this is quite slow!
  if(db["solver_type"].is("direct"))
  {
    Output::print("\nstarting with RT3 on quads");
    errors = {{ 0.0052759466028633, 0.065778892133375, 0.39754725046775,
                0.00083489202575657, 0.063038523036227 }};
    check_darcy2d(domain, db, 1003, errors);
  }
  
  Output::print("\nstarting with BDM1 on quads");
  errors = {{ 1.7301620785317, 23.120239110875, 25.376165640701,
              0.32075488021636, 4.4428829381584 }};
  check_darcy2d(domain, db, 1011, errors);
  
  Output::print("\nstarting with BDM2 on quads");
  errors = {{ 0.36989913525384, 8.7083491818683, 8.7386341355059,
              0.1118419830706, 2.6029572935706 }};
  check_darcy2d(domain, db, 1012, errors);
  
  // this is really slow (for the preconditioner
  // semi_implicit_method_for_pressure_linked_equations)
  if(db["solver_type"].is("direct")
    || db["preconditioner"].is("least_squares_commutator"))
  {
    Output::print("\nstarting with BDM3 on quads");
    errors = {{ 0.073739263678954, 2.2160565937362, 2.5477742459801,
                0.028107167396319, 0.99433882928285 }};
    check_darcy2d(domain, db, 1013, errors);
  }
}

void tests_on_triangles(unsigned int nRefinements, ParameterDatabase& db)
{
  // default construct a domain object
  TDomain domain(db);
  // the domain is initialised with default description and default
  // initial mesh
  domain.Init((char*)"Default_UnitSquare", (char*)"TwoTriangles");

  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }

  std::array<double, 5> errors;

  Output::print("\nstarting with RT0 on triangles");
  errors = {{ 2.0193948379924, 19.174978565868, 31.034055967159,
              0.24829883521652, 4.4428829381584 }};
  check_darcy2d(domain, db, 1000, errors);
  
  Output::print("\nstarting with RT1 on triangles");
  errors = {{ 0.45226149917363, 5.7965843062654, 14.561985165855,
              0.074478813057122, 2.3249493136119 }};
  check_darcy2d(domain, db, 1001, errors);
  
  // this is really slow!
  if(db["solver_type"].is("direct")
     || db["preconditioner"].is("least_squares_commutator"))
  {
    Output::print("\nstarting with RT2 on triangles");
    errors = {{ 0.077083425180006, 1.0445467331298, 4.1457189447051,
                0.0038190922876324, 0.81680935765046 }};
    check_darcy2d(domain, db, 1002, errors);
  }
  
  // this is really slow!
  if(db["solver_type"].is("direct"))
  {
    Output::print("\nstarting with RT3 on triangles");
    errors = {{ 0.011582995440399, 0.24199767560786, 0.79591235625871,
                0.0028837387982513, 0.2072548767011 }};
    check_darcy2d(domain, db, 1003, errors);
  }
  
  Output::print("\nstarting with BDM1 on triangles");
  errors = {{ 1.3471909779398, 19.174978565868, 24.552917977967,
              0.25991524351834, 4.4428829381584 }};
  check_darcy2d(domain, db, 1011, errors);
  
  Output::print("\nstarting with BDM2 on triangles");
  errors = {{ 0.21067940108674, 5.7965843062654, 8.2247391949719,
              0.074736732125372, 2.3558018820457 }};
  check_darcy2d(domain, db, 1012, errors);
  
  // this is really slow!
  if(db["solver_type"].is("direct"))
  {
    Output::print("\nstarting with BDM3 on triangles");
    errors = {{ 0.043596263688953, 1.0445467331299, 1.7271188357497,
                0.0035276997665613, 0.82269214135298 }};
    check_darcy2d(domain, db, 1013, errors);
  }
}



// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of databases
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  // Set Database values (this is what is usually done by the input-file)
  TDatabase::ParamDB->PROBLEM_TYPE = 0; // problem type is not needed
  TDatabase::ParamDB->EXAMPLE = 0; // known sin-cos solution
  // velocity space code for Raviart-Thomas (RT) and 
  // Brezzi-Douglas-Marini(BDM) elements:
  // 1000    RT_0
  // 1001    RT_1
  // 1002    RT_2
  // 1003    RT_3
  // 1011    BDM_1
  // 1012    BDM_2
  // 1013    BDM_3
  TDatabase::ParamDB->VELOCITY_SPACE = 1000;
  // automatically choose pressure space to get inf-sup stable pair
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
    // high order quadrature for computing errors
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  TDatabase::ParamDB->SIGMA_PERM = 1.; // permeability
  TDatabase::ParamDB->SOLVER_TYPE = 2; // use direct solver

  unsigned int nRefinements = 2;
  
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_output_database());
  db["problem_type"] = 0; // problem type is not needed
  db["example"] = 0;
  db["residual_tolerance"] = 1.0e-13;
  
  db["output_compute_errors"] = true;
  
  Output::print("\n\n ----------- direct solver -----------\n");
  db["solver_type"] = "direct";
  tests_on_quads(nRefinements, db);
  tests_on_triangles(nRefinements, db);
  
  
  db["solver_type"] = "iterative";
  db["max_n_iterations"] = 1000;
  
  Output::print("\n\n --------- fgmres+lsc solver ---------\n");
  db["preconditioner"] = "least_squares_commutator";
  
  tests_on_quads(nRefinements, db);
  tests_on_triangles(nRefinements, db);
  
  Output::print("\n\n -------- fgmres+simple solver -------\n");
  db["preconditioner"] = "semi_implicit_method_for_pressure_linked_equations";
  
  tests_on_quads(nRefinements, db);
  tests_on_triangles(nRefinements, db);
  
  return 0;
}
