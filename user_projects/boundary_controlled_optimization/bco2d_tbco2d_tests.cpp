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
             std::array<double, int(4)> errors, double tol)
{

}

void set_errors(int test_case,
                std::array<std::array<double, int(4)>,3>& errors)
{

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
  std::array<std::array<double, int(4)>,3> errors;
  double tol = 1.e-9;
  int test_number;

  {
    //=======================================================================
    //============= TEST 1 : (stationary) BCO ===============================
    //=======================================================================
    test_number = 1;
    set_errors(test_number, errors);
    BoundaryControlledOptimization bco(domain, parmoon_db);
    parmoon_opt::optimize(bco, parmoon_db);
    compare(bco,errors[0],tol);

  }

  {
    //=======================================================================
    //============= TEST 2 : TimeBCO ========================================
    //=======================================================================
    test_number = 2;
    set_errors(test_number, errors);
    Time_BoundaryControlledOptimization tbco(domain, parmoon_db);
    parmoon_opt::optimize(tbco, parmoon_db);
    compare(tbco,errors[0],tol);
  }

}
