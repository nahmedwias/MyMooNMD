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

  Example_NSE2D example_nse2d(navierstokes_db); // Construct Example for NSE
  NSE2D nse2d(domain, navierstokes_db, example_nse2d); // Construct NSE system


  cout << navierstokes_db["solver_type"] << endl;
  Output::print<2>("This is just a sample ", navierstokes_db["example"]);





  stopwatch.print_time("total program duration: ");
  return 0;
}
// end main
