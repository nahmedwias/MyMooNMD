// =======================================================================
//
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <Darcy.h>

// =======================================================================
// main program
// =======================================================================
int main(int, char* argv[])
{
  double t_start = GetTime();
  //  declaration of database, you need this in every program
  TDatabase Database(argv[1]);
  TFEDatabase3D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(std::string(argv[1]));
  
  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  TDomain domain(parmoon_db);
 
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  
  //=========================================================================
  Darcy<3> darcy3d(domain, parmoon_db);
  darcy3d.assemble();
  darcy3d.solve();
  darcy3d.output();
  //=========================================================================
  
  Output::print<1>("used time: ", GetTime() - t_start, " seconds");
  Output::close_file();
  return 0;
} // end main
