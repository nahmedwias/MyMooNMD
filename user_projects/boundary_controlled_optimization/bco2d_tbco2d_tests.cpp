#include "Time_BoundaryControlledOptimization.hpp"
#include "BoundaryControlledOptimization.hpp"
#include "parmoon_optimization_routines.hpp"
#include <Chrono.h>
#include <Database.h>
#include <Domain.h>
#include <FEDatabase2D.h>
#include <LoopInfo.h>
#include <ParMooN_repository_info.h> // to get SOURCE_DIRECTORY

// Compare optimal solutions to a given one
template <class Opt>
void compare(const Opt& bco_tbco_object,
             std::array<double, 2> reference_solution, double tol)
{
  double solution = bco_tbco_object.get_J_hat();
  std::vector<double> control = bco_tbco_object.get_control_old();
  double norm_control = std::accumulate(control.begin(), control.end(), 0.0,
                                             [](double a, double b)
                                             {return a + b * b;});
  norm_control = std::sqrt(norm_control);

  // check the cost functional and the norm of the control vector
  if( fabs(solution-reference_solution[0]) > tol ||
      fabs(norm_control-reference_solution[1]) > tol ||
      solution != solution) //check for nan!
  {
    ErrThrow("Solution is not correct! Computed solution = ", setprecision(10),solution,
             " , Reference solution = ", reference_solution[0],
             " Computed norm of control vector = ", setprecision(10), norm_control,
             " , Reference = ", reference_solution[1]);
}
}

void set_ref_sol(int test_case,
                 std::array<double, 2>& reference_solution)
{
  // Unfortunately, the solutions are slightly different from one solver
  // to another. This explains the switch here.
  // Note however that they should be in the same order of magnitude.
  switch(test_case)
  {
    case 1: // test 1 - bco with gradient-free
      // this is the solution when nonlinloop_maxiteration = 5, xtol =1.e-5,
      // algorithm_nlopt = LN_COBYLA (gradient-free) = 25, max opt it = 22
      // velo space = 2
      reference_solution[0] = 17.89696331; // min cost functional
      reference_solution[1] = 0.9835083284; // norm of control vector
      break;
    case 2: // test 2 - bco with gradient-free
      // with algo 40 LD_SLSQP (gradient-based)
      // number of iterations needed: 47, computational time = 16.95
      reference_solution[0] = 16.52532739; // min cost functional 16.52440131
      reference_solution[1] = 1.3404218; // norm of control vector
      // note that we take the solutions from (t)bco object, corresponding to
      // the last iterations. They might slightly differ from the solution
      // printed by the nlopt object.
      break;
    case 3: // test 3 - tbco with gradient-free LN_COBYLA (algo 25)
      // the max number of iterations 12 is reached. At this iteration, the solution should be:
      reference_solution[0] = 383.0113357; // min cost functional
      reference_solution[1] = 1.271266404; // norm of control vector
      // note1: by having nl iterations long enough, the solution is similar
      // to stationary
      // note2: the cost functional is not comparable to stationary
      // because it is the sum in time....but the order of magnitude is
      // suspicious
      // note3: the solution at t=0 corresponds to a fully developed
      // Poiseuille flow
      // note4: control dependent on space, but constant in time
      break;
    case 4: // the max number of iterations 9 is reached. At this iteration, the solution should be:
      reference_solution[0] = 392.3244362; // min cost functional
      reference_solution[1] = 1.255678231; // norm of control vector
      // note1: by having nl iterations long enough,
      // the solution is similar to stationary
      break;
    default:
      ErrThrow("Unknown test case.");
  }
}


// this should be done inside the Domain class, but it is not yet.
void refine_domain(TDomain& domain, ParameterDatabase& db, bool write_ps_file)
{
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
  {
    domain.refine_and_get_hierarchy_of_collections(db);
  }
  
  // write grid into an Postscript file
  if(write_ps_file)
  {
    domain.PS("Domain.ps", It_Finest, 0);
  }
}

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    ErrThrow("Please provide an input file with parameters to run the ParMooN ",
             "program 'boundary controlled optimization'");
  }
  
  bool testall = false;
  if (argv[2])
    testall = (std::string(argv[2]).compare("testall") == 0);

  //  declaration of database, you need this in every program
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(TimeDiscretization::default_TimeDiscretization_database());
  parmoon_db.read(argv[1]);

  /** Hard-code for path to geometries **/
  std::string path = parmoon::source_directory + std::string("/data/mesh/backward_facing_step/");
  parmoon_db.add("boundary_file", path+"backward_facing_step.PRM", " ");
  parmoon_db.add("geo_file", path+"backward_facing_step_tria1.mesh", " ");

  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(parmoon_db);

   // set variables' value in TDatabase using argv[1] (*.dat file)
  TDomain domain(parmoon_db);
  refine_domain(domain, parmoon_db, false);
  // Expected solution - if the results is equal to this given solution
  // then the test is passed
  std::array<double, 2> reference_solution;
  double tol = 1.e-7;
  int test_number;

  {
    //=======================================================================
    //============= TEST 1 : (stationary) BCO with gradient-free solver =====
    //=======================================================================
    test_number = 1;
    parmoon_db["max_n_evaluations"] = 22; // force lower iterations because of bug
    set_ref_sol(test_number, reference_solution);
    BoundaryControlledOptimization<2> bco(domain, parmoon_db);
    parmoon_opt::optimize(bco, parmoon_db);
    compare(bco,reference_solution,tol);
  }

  if(testall)
  {
    //=======================================================================
    //============= TEST 2 : (stationary) BCO with gradient-based solver ====
    //== This is an important test to check the adjoint equations ===========
    //=======================================================================
    test_number = 2;
    parmoon_db["algorithm_nlopt"] = 40;
    parmoon_db["max_n_evaluations"] = 100;
    Output::print("Test ", test_number, ", Changing algoritm nlopt to: ",
                  parmoon_db["algorithm_nlopt"]);
    set_ref_sol(test_number, reference_solution);
    BoundaryControlledOptimization<2> bco2(domain, parmoon_db);
    parmoon_opt::optimize(bco2, parmoon_db);
    compare(bco2,reference_solution,tol);
  }

  if(testall)
  {
    //=======================================================================
    //============= TEST 3 : TimeBCO ========================================
    //=======================================================================
    test_number = 3;
    // reset some parameters specifically to test 3
    parmoon_db["algorithm_nlopt"] = 25;
    parmoon_db["max_n_evaluations"] = 12; // force lower iterations because of bug
    Output::print("Test ", test_number, ", Changing algoritm nlopt to: ",
                  parmoon_db["algorithm_nlopt"],
                  ", and changing max_n_evaluations to: ",
                  parmoon_db["max_n_evaluations"]);
    set_ref_sol(test_number, reference_solution);
    Time_BoundaryControlledOptimization<2> tbco(domain, parmoon_db);
    parmoon_opt::optimize(tbco, parmoon_db);
    compare(tbco,reference_solution,tol);
  }

  /* If you want to add time-dependent tests, do not forget to
   * re-initialize these global parameters    */
  TDatabase::TimeDB->TIMESTEPLENGTH=parmoon_db["time_step_length"];
  parmoon_db["time_end"]=0.03;
  TDatabase::TimeDB->CURRENTTIME=  parmoon_db["time_start"];


  {
    Output::print(
     "============= TEST 4 : TimeBCO with time-dependent control ============\n"
     "==== This tests a control which is independent from space, in one =====\n"
     "==== direction only, and dependent of time (n_control = 4 time steps) =\n"
     "=======================================================================\n"
    );
    test_number = 4;
    // reset some parameters specifically to test 4
    parmoon_db["max_n_evaluations"] = 9; // force lower iterations because of bug
    // time-dependent control, only in x-direction
    parmoon_db["control_depends_on_time"] = true;
    parmoon_db["control_depends_on_space"] = false;
    parmoon_db["control_in_x_direction"] = true;
    parmoon_db["control_in_y_direction"] = false;
    set_ref_sol(test_number, reference_solution);
    Time_BoundaryControlledOptimization<2> tbco2(domain, parmoon_db);
    parmoon_opt::optimize(tbco2, parmoon_db);
    compare(tbco2,reference_solution,tol);
  }

  Output::print("TEST SUCCESSFUL!");
  return 0;
}
