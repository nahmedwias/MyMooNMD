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
#include <TimeConvectionDiffusionROM.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <TimeDiscRout.h>
#include <SnapshotsCollector.h>

using namespace std;

int main(int, char* argv[])
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

  // initialize snapshot writer
  SnapshotsCollector snaps( parmoon_db );
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);
  
  TimeConvectionDiffusion<2> tcd(Domain, parmoon_db);
  // for testing ROM
  //Example_TimeCD2D example = tcd.get_example();
  //TimeConvectionDiffusionROM<2> tcd_rom(parmoon_db, tcd.get_example());
  // --------------------------
  TimeDiscretization& tss = tcd.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();

  // ======================================================================
  // assemble matrices and right hand side at start time  
  tcd.assemble_initial_time();
  // ======================================================================
  
  double start_time = parmoon_db["time_start"];
  TDatabase::TimeDB->CURRENTTIME = start_time;
  tcd.output();

  // store initial condition as snapshot
  if (parmoon_db["write_snaps"])
    snaps.write_data(tcd.get_solution());
  
  // ======================================================================
  // time iteration
  // ======================================================================
  while(!tss.reached_final_time_step())
  {
    tss.current_step_++;
    // Output::print("mem before: ", GetMemory());
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    tss.set_time_disc_parameters();
    tss.current_time_ += tss.get_step_length();
    double tau = parmoon_db["time_step_length"];
    
    TDatabase::TimeDB->CURRENTTIME += tau;
    Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    SetTimeDiscParameters(1);
    
    tcd.assemble();    
    tcd.solve();

    if((tss.current_step_-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      tcd.output();

    // write the snap shots
    if (parmoon_db["write_snaps"])
      snaps.write_data(tcd.get_solution(), tss.current_step_);
    
  }
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();
  return 0;
} // end main
