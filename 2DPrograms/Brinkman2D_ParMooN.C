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
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain Domain(argv[1]);  
  
 //// //set PROBLEM_TYPE to NSE if not yet set (3 means Stokes, 5 Naver-Stokes)
 //// if(TDatabase::ParamDB->PROBLEM_TYPE!=3 && TDatabase::ParamDB->PROBLEM_TYPE!=5)
 ////   TDatabase::ParamDB->PROBLEM_TYPE = 5;
  
    //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);
  
  //// possibly change parameters in the database, if they are not meaningful now
  ////Database.CheckParameterConsistencyNSE();
    
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
    
  /* include the mesh from a mesh generator, for a standard mesh use the
   * build-in function. The GEOFILE describes the boundary of the domain. */
  Domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); // call mesh generator
  
  // refine grid up to the coarsest level
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    Domain.RegRefineAll();  
  
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain.PS("Domain.ps",It_Finest,0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
    
  Example_Brinkman2D example;
  

  //=========================================================================
  // create an object of the Brinkman class

  Brinkman2D brinkman2d(Domain, example);
  brinkman2d.assemble();
  brinkman2d.solve();
  brinkman2d.output();

  //=========================================================================

  Output::close_file();
    return 0;
         // }
} // end main
