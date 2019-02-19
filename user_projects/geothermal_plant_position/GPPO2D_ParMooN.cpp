#include "GeothermalPlantsPositionOptimization.hpp"
#include "parmoon_optimization_routines.hpp"
#include "Chrono.h"
#include "Database.h"
#include "Domain.h"
#include "FEDatabase2D.h"
#include "LoopInfo.h"


// this should be done inside the Domain class, but it is not yet.
/** ************************************************************************ */
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

/** ************************************************************************ */
int main(int argc, char* argv[])
{
  // start a stopwatch which measures time spent in program parts
  Chrono timer;

  //  declaration of database, you need this in every program
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;

  // create a default database
  ParameterDatabase db = GeothermalPlantsPositionOptimization<2>::default_GPPO_database();

    cout << "transversal_dispersion_factor in GPPO2D_ParMooN: " << db["transversal_dispersion_factor"] << endl;

  if(argc < 2)
  {
    Output::print("Please provide an input file with parameters to run the "
            "ParMooN program 'geothermal plants position optimization'");
    std::string default_name = "default_gppo_input.dat";
    db.merge(parmoon_opt::default_optimization_database());
    // Todo: merge with NSE2D_GPPO database and TCD2D_Temperature database)
    db.write(default_name, true);
    Output::print("A default input file to run this program is written to ",
            default_name);
    Output::print("You can now run the program via\n  ", argv[0], " ",
            default_name);
    return 0;
  }
  // read parameters
  db.read(argv[1]);

  // open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(db["outfile"]);
  Output::setVerbosity(db["verbosity"]);

  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(db);
  // write all Parameters to the OUTFILE (not to console) for later reference
  db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

   // set variables' value in TDatabase using argv[1] (*.dat file)
  TDomain domain(db);  // OLD: TDomain domain(db, argv[1]);
  // refine grid
    domain.refine_and_get_hierarchy_of_collections(db); // OLD:   refine_domain(domain, db["output_write_ps"]);


  timer.restart_and_print("Started all preparations before optimization loop.");
  GeothermalPlantsPositionOptimization<2> gppo(domain, db);
  timer.restart_and_print("Finished all preparations before optimization loop.");

  parmoon_opt::print_nlopt_version();

  timer.stop_and_print("Started optimization loop.");
  parmoon_opt::optimize(gppo, db);
  timer.stop_and_print("Finished optimization loop.");
  timer.print_total_time("Finished entire boundary controlled optimization.");
  Output::close_file();
}
