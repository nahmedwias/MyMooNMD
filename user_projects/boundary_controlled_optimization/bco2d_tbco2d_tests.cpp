#include "Time_BoundaryControlledOptimization.hpp"
#include "BoundaryControlledOptimization.hpp"
#include "parmoon_optimization_routines.hpp"
#include <Chrono.h>
#include <Database.h>
#include <Domain.h>
#include <FEDatabase2D.h>
#include <LoopInfo.h>

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
      // algorithm_nlopt = LN_COBYLA (gradient-free) = 25, max opt it = 100
      // velo space = 2
      // number of iterations needed: 95, computational time = 18.6s
      reference_solution[0] = 17.42159598; // min cost functional
      reference_solution[1] = 1.006193328; // norm of control vector
      break;
    case 2: // test 2 - bco with gradient-free
      // with algo 40 LD_SLSQP (gradient-based)
      // number of iterations needed: 79, computational time = 16.95
      reference_solution[0] = 16.53840917; // min cost functional 16.5334
      reference_solution[1] = 1.318516539; // norm of control vector 1.31763
      // note that we take the solutions from (t)bco object, corresponding to
      // the last iterations. They might slightly differ from the solution
      // printed by the nlopt object.
      break;
    case 3: // test 3 - tbco with gradient-free LN_COBYLA (algo 25)
      // the max number of iterations 50 is reached (29s), otherwise the test
      // will be too long. At this iteration, the solution should be:
      reference_solution[0] = 374.0993461; // min cost functional
      reference_solution[1] = 1.053401308; // norm of control vector
      // note1: by having nl iterations long enough, the solution is similar
      // to stationary
      // note2: the cost functional is not comparable to stationary
      // because it is the sum in time....but the order of magnitude is
      // suspicious
      // note3: the solution at t=0 corresponds to a fully developed
      // Poiseuille flow
      // note4: control dependent on space, but constant in time
      break;
    default:
      ErrThrow("Unknown test case.");
  }
}


// this should be done inside the Domain class, but it is not yet.
void refine_domain(TDomain& domain, bool write_ps_file)
{
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
  {
    domain.RegRefineAll();
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
  
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);

  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(parmoon_db);

   // set variables' value in TDatabase using argv[1] (*.dat file)
  TDomain domain(parmoon_db, argv[1]);
  refine_domain(domain, false);
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
    set_ref_sol(test_number, reference_solution);
    BoundaryControlledOptimization bco(domain, parmoon_db);
    parmoon_opt::optimize(bco, parmoon_db);
    compare(bco,reference_solution,tol);
  }

  {
    //=======================================================================
    //============= TEST 2 : (stationary) BCO with gradient-based solver ====
    //== This is an important test to check the adjoint equations ===========
    //=======================================================================
    test_number = 2;
    parmoon_db["algorithm_nlopt"] = 40;
    Output::print("Test ", test_number, ", Changing algoritm nlopt to: ",
                  parmoon_db["algorithm_nlopt"]);
    set_ref_sol(test_number, reference_solution);
    BoundaryControlledOptimization bco2(domain, parmoon_db);
    parmoon_opt::optimize(bco2, parmoon_db);
    compare(bco2,reference_solution,tol);
  }


  {
    //=======================================================================
    //============= TEST 3 : TimeBCO ========================================
    //=======================================================================
    test_number = 3;
    // reset some parameters specifically to test 3
    parmoon_db["algorithm_nlopt"] = 25;
    parmoon_db["max_n_evaluations"] = 50;
    Output::print("Test ", test_number, ", Changing algoritm nlopt to: ",
                  parmoon_db["algorithm_nlopt"],
                  ", and changing max_n_evaluations to: ",
                  parmoon_db["max_n_evaluations"]);
    set_ref_sol(test_number, reference_solution);
    Time_BoundaryControlledOptimization tbco(domain, parmoon_db);
    parmoon_opt::optimize(tbco, parmoon_db);
    compare(tbco,reference_solution,tol);
  }

  /* If you want to add time-dependent tests, do not forget to
   * re-initialize these global parameters
     TDatabase::TimeDB->STARTTIME=0;
     TDatabase::TimeDB->TIMESTEPLENGTH=parmoon_db["time_step_length"];
     TDatabase::TimeDB->ENDTIME=0.03;
     TDatabase::TimeDB->CURRENTTIME=  TDatabase::TimeDB->STARTTIME;
   */

  Output::print("TEST SUCCESSFUL!");
  return 0;
}
