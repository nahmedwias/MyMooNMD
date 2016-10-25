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
    //for(refinement_n_initial_steps=1; refinement_n_initial_steps <= 6;++refinement_n_initial_steps)
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
  
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
 
  // refine grid
  size_t n_ref =  Domain.get_n_initial_refinement_steps();
  for(size_t i=0; i< n_ref; i++)
  {
      Domain.RegRefineAll();
  }
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);
    
  Example_Brinkman2D example(parmoon_db);

  //=========================================================================
  // create an object of the Brinkman class
  Brinkman2D brinkman2d(Domain, parmoon_db, example);
    Output::print<>("Datenbank");
    parmoon_db.info();
  brinkman2d.assemble();
  brinkman2d.solve();
    Output::print<>(TDatabase::ParamDB->VELOCITY_SPACE);
        Output::print<>(TDatabase::ParamDB->PRESSURE_SPACE);
  brinkman2d.output();

  //=========================================================================

  Output::close_file();
    return 0;
         // }
} // end main
