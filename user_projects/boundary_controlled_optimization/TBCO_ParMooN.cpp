#include "Time_BoundaryControlledOptimization.hpp"
#include "parmoon_optimization_routines.hpp"
#include <Chrono.h>
#include <Database.h>
#include <Domain.h>
#ifdef __2D__
  #include <FEDatabase2D.h>
#else
  #include <FEDatabase3D.h>
#endif
#include <LoopInfo.h>


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
#ifdef __2D__
  TFEDatabase2D FEDatabase;
#else
  TFEDatabase3D FEDatabase;
#endif
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
  TDomain domain(parmoon_db);
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  if(parmoon_db["output_write_ps"])
  {
    domain.PS("Domain.ps", It_Finest, 0);
  }

#ifdef __2D__
    Time_BoundaryControlledOptimization<2> tbco(domain, parmoon_db);
#else
    Time_BoundaryControlledOptimization<3> tbco(domain, parmoon_db);
#endif
  
  timer.restart_and_print("all preparations before optimization loop");
  parmoon_opt::print_nlopt_version();
  parmoon_opt::optimize(tbco, parmoon_db);
  timer.stop_and_print("optimization loop");
  timer.print_total_time("entire boundary controlled optimization");
}
