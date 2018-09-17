#include "BoundaryControlledOptimization.hpp"
#include "parmoon_optimization_routines.hpp"
#include "Chrono.h"
#include "Database.h"
#include "Domain.h"
#include "FEDatabase2D.h"
#include "LoopInfo.h"


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
  // start a stopwatch which measures time spent in program parts
  Chrono timer;
  
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase; 
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);
  
  // open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(parmoon_db);
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  
   // set variables' value in TDatabase using argv[1] (*.dat file)
  TDomain domain(parmoon_db, argv[1]);
  refine_domain(domain, parmoon_db["output_write_ps"]);

  BoundaryControlledOptimization bco(domain, parmoon_db);
  
  timer.restart_and_print("all preparations before optimization loop");
  parmoon_opt::print_nlopt_version();
  parmoon_opt::optimize(bco, parmoon_db);
  timer.stop_and_print("optimization loop");
  timer.print_total_time("entire boundary controlled optimization");
}
