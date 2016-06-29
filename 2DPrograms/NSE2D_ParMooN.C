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
#include <NSE2D.h>
#include <Example_NSE2D.h>
#include <Chrono.h>
#include <LoopInfo.h>
#include <ParameterDatabase.h>


// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // start a stopwatch which measures time spent in program parts
  Chrono chrono_parts;
  
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase; 
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(argv[1], parmoon_db);
  
  //open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // possibly change parameters in the database, if they are not meaningful now
  Database.CheckParameterConsistencyNSE();
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
  
  // refine grid
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
    domain.RegRefineAll();
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    domain.PS("Domain.ps", It_Finest, 0);
  
  Example_NSE2D example(parmoon_db["example"]);
  
  // create an object of the Navier-Stokes class
  NSE2D ns(domain, parmoon_db, example);
  ns.assemble();
  // if solution was not zero up to here, you should call 
  //ns.assemble_nonlinear_term();
  
  ns.stopIt(0);
  
  LoopInfo loop_info("nonlinear");
  loop_info.print_time_every_step = true;
  loop_info.verbosity_threshold = 1; // full verbosity
  loop_info.print(0, ns.getFullResidual());
  
  chrono_parts.print_time("setting up spaces, matrices, linear assemble");
  chrono_parts.reset();
  
  //======================================================================
  // nonlinear loop
  // in function 'stopIt' termination condition is checked
  for(unsigned int k = 1;; k++)
  {
    Output::print(); // new line for a new nonlinear iteration
    ns.solve();
    
    //no nonlinear iteration for Stokes problem
    if(parmoon_db["problem_type"].is(3))
      break;
    
    ns.assemble_nonlinear_term();
    
    if(ns.stopIt(k))
    {
      loop_info.finish(k, ns.getFullResidual());
      break;
    }
    else
      loop_info.print(k, ns.getFullResidual());
  } // end for k
  
  chrono_parts.print_time("solving procedure");
  
  ns.output();
  
  Output::close_file();
  return 0;
} // end main
