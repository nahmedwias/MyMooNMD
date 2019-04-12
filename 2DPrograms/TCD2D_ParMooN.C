// =============================================================================
//
// Purpose:     main program for scalar equations with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 08.08.2014

// =============================================================================

#include <Database.h>
#include <Domain.h>
#include <FEDatabase2D.h>
#include <SnapshotsCollector.h>
#include <TimeConvectionDiffusion.h>
#include <TimeConvectionDiffusionPOD.h>
#include <TimeConvectionDiffusionROM.h>
#include <TimeDiscRout.h>

#include <sys/stat.h>
#include <sys/types.h>


using namespace std;

int main(int, char* argv[])
{
  bool do_rom = false;
  double t_start = GetTime();
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);

  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  // ===========================================================================
  // set the database values and generate mesh
  // ===========================================================================
  TDomain Domain(parmoon_db);

  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

  // refine grid
  Domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
  {
    Domain.PS("Domain.ps", It_Finest, 0);
  }

  // ===========================================================================
  // Solve the PDE using FEM and acquire snapshots
  // ===========================================================================
  if( parmoon_db["write_snaps"] || (   !parmoon_db["write_snaps"]
                                    && !parmoon_db["compute_POD_basis"]
                                    && !parmoon_db["ROM_method"]) )
  {
    TimeConvectionDiffusion<2> tcd(Domain, parmoon_db);

    // initialize snapshot writer (only needed if parmoon_db["write_snaps"])
    // the solutions are stored in a special collector which can be later read
    // by the computation of POD the basis
    SnapshotsCollector snaps(parmoon_db);

    TimeDiscretization& tss = tcd.get_time_stepping_scheme();
    tss.current_step_ = 0;
    tss.set_time_disc_parameters();

    // assemble matrices and right hand side at start time
    tcd.assemble_initial_time();

    double start_time = parmoon_db["time_start"];
    TDatabase::TimeDB->CURRENTTIME = start_time;
    tcd.output();

    // store initial condition as snapshot
    if(parmoon_db["write_snaps"])
    {
      snaps.write_data(tcd.get_solution());
    }

    // time iteration
    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;
      // Output::print("mem before: ", GetMemory());
      tss.set_time_disc_parameters();
      tss.current_time_ += tss.get_step_length();
      double tau = parmoon_db["time_step_length"];

      TDatabase::TimeDB->CURRENTTIME += tau;
      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
      SetTimeDiscParameters(1);

      tcd.assemble();
      tcd.solve();

      if((tss.current_step_-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        tcd.output();
      }

      // write the snapshots
      if(parmoon_db["write_snaps"])
      {
        snaps.write_data(tcd.get_solution(), tss.current_step_);
      }
    }
  } // if FEM + snapshots

  // ===========================================================================
  // Compute the POD basis from snapshots
  // ===========================================================================
  if( parmoon_db["compute_POD_basis"] )
  {
    auto collections            = Domain.get_grid_collections();
    TCollection& cellCollection = *collections.front();

    TimeConvectionDiffusionPOD<2> tcd_pod(cellCollection, parmoon_db); // TODO add TYPE (i.e. TCD2D, TNSE2D,...)

    tcd_pod.compute_pod_basis();
    //pod_2d.read_basis();

    tcd_pod.output();
  } // if POD

  // ===========================================================================
  // Solve the PDE using ROM instead of FEM
  // ===========================================================================
  if( parmoon_db["ROM_method"] )
  {
    Output::print("\nROM: not implemented yet");
  } // if ROM

  // ===========================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ===========================================================================
  Output::close_file();
  return 0;
} // end main
