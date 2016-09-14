// =======================================================================
//
// Purpose:
//
// Author:      NA
//
// History:     Start 13.09.2016

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <NSE2D.h>
#include <Example_NSE2D.h>
#include <Chrono.h>
#include <LoopInfo.h>
#include <ParameterDatabase.h>

// ***** MAIN PROGRAM ***** //
int main(int argc, char* argv[])
{
  Chrono stopwatch; // Start a stopwatch for time measurement during execution

  TDatabase Database;  // Initialize User Input Databases.
  ParameterDatabase navierstokes_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  navierstokes_db.read(fs);
  fs.close();
  TFEDatabase2D FEDatabase; // Initialize FE Database.

  TDomain domain(argv[1], navierstokes_db); // Initialize geometry

  // ***** Write parameters to output file ***** //
  Output::set_outfile(navierstokes_db["outfile"]);
  Output::setVerbosity(navierstokes_db["verbosity"]);
  Database.CheckParameterConsistencyNSE();
  Database.WriteParamDB(argv[0]);

  std::list<TCollection* > gridCollections
  = domain.refine_and_get_hierarchy_of_collections(navierstokes_db);
  if(navierstokes_db["output_write_ps"]) domain.PS("Domain.ps", It_Finest, 0);

  domain.print_info("NSE2D domain");  // Output domain info

  // ***** Declaring objects for NSE2D ***** //
  Example_NSE2D example_nse2d(navierstokes_db); // Construct Example for NSE
  NSE2D nse2d(domain, navierstokes_db, example_nse2d); // Construct NSE system
  nse2d.assemble(); // assemble linear term
  nse2d.stopIt(0); // check initial residual

  // ***** Initialize object for iterations residual output ***** //
  LoopInfo loop_info("nonlinear");
  loop_info.print_time_every_step = true;
  loop_info.verbosity_threshold = 1; // full verbosity
  loop_info.print(0, nse2d.getFullResidual());

  stopwatch.print_time("setting up spaces, matrices, linear assemble");
  stopwatch.reset();

  //================================================================
  // ****** NON-LINEAR LOOP ****** //
  for(unsigned int k = 1;; k++)
  {
    nse2d.solve();
    if(navierstokes_db["problem_type"].is(3)) break; // Stokes
    nse2d.assemble_nonlinear_term();

    if(nse2d.stopIt(k)) // Check residuals
    {
      loop_info.finish(k, nse2d.getFullResidual());
      break;
    }
    else loop_info.print(k, nse2d.getFullResidual());
  } // end for k, non linear loop
  stopwatch.print_time("total solving duration: ");
  //================================================================




  nse2d.output();
  Output::close_file();
  return 0;
}
// end main
