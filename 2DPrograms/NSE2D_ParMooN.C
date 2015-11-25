// =======================================================================
//
// Purpose:     main program for solving a stationary NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 23.08.2014

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <NSE2D.h>
#include <Output2D.h>
#include <MainUtilities.h>
#include <LocalAssembling2D.h>
#include <Example_NSE2D.h>

#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase; 
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain Domain(argv[1]);  
  
  //set PROBLEM_TYPE to NSE if not yet set (3 means Stokes, 5 Naver-Stokes)
  if(TDatabase::ParamDB->PROBLEM_TYPE!=3 && TDatabase::ParamDB->PROBLEM_TYPE!=5)
    TDatabase::ParamDB->PROBLEM_TYPE = 5;
  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);
  
  // possibly change parameters in the database, if they are not meaningful now
  Database.CheckParameterConsistencyNSE();
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
  
  Example_NSE2D example;
  
  // create an object of the Naviert-Stokes class
  NSE2D ns(Domain, example);
  ns.assemble();
  // if solution was not zero up to here, you should call 
  //ns.assemble_nonlinear_term();
  
  ns.stopIt(0);
  Output::print<1>("Nonlinear iteration step   0\t", ns.getResiduals());
  
  //======================================================================
  // nonlinear loop
  // in function 'stopIt' termination condition is checked
  for(unsigned int k = 1;; k++)
  {
    ns.solve();
    
    //no nonlinear iteration for Stokes problem
    if(TDatabase::ParamDB->PROBLEM_TYPE == 3)
      break;
    
    ns.assemble_nonlinear_term();
    
    Output::print<1>("nonlinear iteration step ", setw(3), k, "\t", 
                     ns.getResiduals());
    if(ns.stopIt(k))
      break;
  } // end for k
  
  ns.output();
  
  Output::close_file();
  return 0;
} // end main
