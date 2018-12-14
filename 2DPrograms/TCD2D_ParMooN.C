// =======================================================================
//
// Purpose:     main program for scalar equations with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 08.08.2014

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include "TimeConvectionDiffusion.h"

#include <sys/stat.h>
#include <sys/types.h>

#include <TimeDiscRout.h>
#include <SNAPS.h>



using namespace std;

int main(int argc, char* argv[])
{
  double t_start = GetTime();
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;
  
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
  SNAPS snaps( parmoon_db );
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);
  
  TimeConvectionDiffusion<2> tcd(Domain, parmoon_db);
  
  TimeDiscretization& tss = tcd.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  // ======================================================================
  // assemble matrices and right hand side at start time  
  tcd.assemble_initial_time();
  // ======================================================================
  
  double start_time = parmoon_db["time_start"];
  double end_time   = parmoon_db["time_end"];
  TDatabase::TimeDB->CURRENTTIME = start_time;
  tcd.output();
  // store initial condition as snapshot
  snaps.write_data(tcd.get_solution());
  // ======================================================================
  // time iteration
  // ======================================================================
  while(tss.current_time_ < end_time - 1e-10)
  {
    tss.current_step_++;
    // Output::print("mem before: ", GetMemory());
    TDatabase::TimeDB->INTERNAL_STARTTIME = tss.current_time_;
    tss.set_time_disc_parameters();
    double tau = parmoon_db["time_step_length"];
    
    tss.current_time_ += tau;
    TDatabase::TimeDB->CURRENTTIME = tss.current_time_;
    Output::print("\nCURRENT TIME: ", tss.current_time_);
    SetTimeDiscParameters(1);
    
    tcd.assemble();    
    tcd.solve();
    // write the snap shots
    snaps.write_data(tcd.get_solution(), tss.current_step_);
    cout<<" current step: "<< tss.current_step_<< endl;
    if((tss.current_step_-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      tcd.output();
  }
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();
  return 0;
} // end main
