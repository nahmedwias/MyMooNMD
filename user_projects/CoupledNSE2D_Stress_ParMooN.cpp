#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Example_CoupledNS_Stress.h>
#include <TimeDiscretizations.h>
#include <TimeDiscRout.h>
#include <LoopInfo.h>
#include <CoupledNavierStokesStress.h>

using namespace std;

int main(int argc, char* argv[])
{
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase2D;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);
  
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  TDomain Domain(parmoon_db);

  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

  // refine grid
  Domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);

  // set some parameters for time stepping
  SetTimeDiscParameters(0);

  // create an object of TimeNavierStokes<2> class
  Example_CoupledNS_Stress ex(parmoon_db);
  CoupledNavierStokesStress<2> nse2d(Domain, parmoon_db, ex);
  
  nse2d.assemble_linear_terms();
  
  return 0;
}