// =======================================================================
//
// Purpose:     main program for solving a stationary scalar equation using ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Darcy2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <Example_Darcy2D.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  double t_start = GetTime();
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(argv[1]);

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
    domain.PS("domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(parmoon_db["WRITE_VTK"].is(1))
    mkdir(parmoon_db["output_directory"], 0777);
  
  // choose example according to the value of db["example"]
  Example_Darcy2D example;
  
  //=========================================================================
  Darcy2D darcy2d(domain, parmoon_db, example);
  darcy2d.assemble();
  darcy2d.solve();
  darcy2d.output();
  //=========================================================================
  
  Output::print<1>("used time: ", GetTime() - t_start, " seconds");
  Output::close_file();
  return 0;
} // end main
