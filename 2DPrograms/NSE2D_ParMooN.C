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
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(argv[1]);
  
  //set PROBLEM_TYPE to NSE if not yet set (3 means Stokes, 5 Naver-Stokes)
  if(!parmoon_db["problem_type"].is(3) && !parmoon_db["problem_type"].is(5))
    parmoon_db["problem_type"] = 5;
  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // possibly change parameters in the database, if they are not meaningful now
  Database.CheckParameterConsistencyNSE();
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
    mkdir(parmoon_db["output_directory"].get<std::string>().c_str(), 0777);
  
  Example_NSE2D example;
  
  // create an object of the Naviert-Stokes class
  NSE2D ns(domain, parmoon_db, example);
  ns.assemble();
  // if solution was not zero up to here, you should call 
  //ns.assemble_nonlinear_term();
  
  ns.stopIt(0);
  
  //======================================================================
  // nonlinear loop
  // in function 'stopIt' termination condition is checked
  for(unsigned int k = 1;; k++)
  {
    Output::print<1>("nonlinear iteration step ", setw(3), k-1, "\t", 
                     ns.getResiduals());
    ns.solve();
    
    //no nonlinear iteration for Stokes problem
    if(parmoon_db["problem_type"].is(3))
      break;
    
    ns.assemble_nonlinear_term();
    
    if(ns.stopIt(k))
      break;
  } // end for k
  
  ns.output();
  
  Output::close_file();
  return 0;
} // end main
