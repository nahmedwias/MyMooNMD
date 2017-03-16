// =======================================================================
// Purpose:     Linear Elasticity
//
// Author:      NA
//
// History:     Start 16.03.2017
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Chrono.h>
#include <LoopInfo.h>
#include <ParameterDatabase.h>
#include <TimeDiscRout.h>

#include <Example_TimeLinElastic2D.h>
#include <Time_LinElastic2D.h>

using namespace std;

// ***** MAIN PROGRAM ***** //
int main(int argc, char* argv[])
{
  Chrono        stopwatch;        // Start a stopwatch for time measurement during execution

  TDatabase     Database;         // Initialize User Input Databases.
  TFEDatabase2D FEDatabase;       // Initialize FE Database.

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]); parmoon_db.read(fs);  fs.close();

  TDomain domain(argv[1], parmoon_db);            // Initialize geometry



  /********************************************************************
  * WRITE PARAMETERS TO OUTFILE
  ********************************************************************/
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  Database.WriteParamDB(argv[0]);

  std::list<TCollection* > gridCollections
    = domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  if(parmoon_db["output_write_ps"]) domain.PS("Domain.ps", It_Finest, 0);
  domain.print_info("LinearElastic2D domain");      // Output domain info


  /********************************************************************
   * Creating TLinElastic object
   ********************************************************************/
  SetTimeDiscParameters(0);    // Initialize parameters for time discretization


  Example_TimeLinElastic2D ElasticityExample(parmoon_db);
  Time_LinElastic2D ElasticStructure;


  cout << "HELLO WOOOORLD!" << endl;
  cout << ElasticStructure.test << endl;


  return 0;
}
// end main
