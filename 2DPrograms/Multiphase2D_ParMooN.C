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
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  TFEDatabase2D FEDatabase; // Initialize FE Database.


//  cout << parmoon_db["solver_type"] << endl;
//  Output::print<5>("This is just a sample ", parmoon_db["residual_tolerance"]);






  stopwatch.print_time("total program duration: ");
  return 0;
}
// end main
