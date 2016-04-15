// =======================================================================
//
// Purpose:  main program for solving a stationary scalar equation using ParMooN
//
// Author:   Sashikumaar Ganesan
//
// History:  Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <CD2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_CD2D.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(argv[1]);

  //set PROBLEM_TYPE to CD if not yet set
  if(parmoon_db["problem_type"].is(0))
    parmoon_db["problem_type"] = 1;
  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
  
  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */
  domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"]);
  
  // refine grid up to the coarsest level
  size_t n_ref = parmoon_db["uniform_refinement_steps"];
  for(size_t i = 0; i < n_ref; i++)
    domain.RegRefineAll();
  
  // write grid into an Postscript file
  if(parmoon_db["write_ps"])
    domain.PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(parmoon_db["WRITE_VTK"].is(1))
    mkdir(parmoon_db["output_directory"].get<std::string>().c_str(), 0777);
   
  //=========================================================================
  CD2D cd2d(domain, parmoon_db);
  cd2d.assemble();
  cd2d.solve();
  cd2d.output();
  //=========================================================================
  //db.info(false);
  Output::close_file();
  return 0;
} // end main
