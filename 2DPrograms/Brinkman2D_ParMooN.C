// =======================================================================
//
// Purpose:     main program for solving a stationary Brinkman equation in ParMooN
//
// Author:      Laura Blank
//
// History:     Implementation started on  14.03.2016

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <Brinkman2D.h>
#include <Output2D.h>
#include <MainUtilities.h>
#include <LocalAssembling2D.h>
#include <Example_Brinkman2D.h>

#include <ParameterDatabase.h>

#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
    //for(UNIFORM_STEPS=1; UNIFORM_STEPS <= 6;++UNIFORM_STEPS)
    //{
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain Domain(argv[1], parmoon_db);
  
  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(parmoon_db["outfile"]);
  
  //// possibly change parameters in the database, if they are not meaningful now
  ////Database.CheckParameterConsistencyNSE();
    
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
    
  /* include the mesh from a mesh generator, for a standard mesh use the
   * build-in function. The GEOFILE describes the boundary of the domain. */
  Domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"]); // call mesh generator
  
  // refine grid
  size_t n_ref =  Domain.get_n_initial_refinement_steps();
  for(size_t i=0; i< n_ref; i++)
    Domain.RegRefineAll();  
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);
    
  Example_Brinkman2D example(parmoon_db["example"]);
    
  //=========================================================================
  // create an object of the Brinkman class
    
  Brinkman2D brinkman2d(Domain, parmoon_db, example);
  brinkman2d.assemble();
  brinkman2d.solve();
  brinkman2d.output();

  //=========================================================================

  Output::close_file();
    return 0;
         // }
} // end main
