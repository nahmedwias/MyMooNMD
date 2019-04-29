#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Example_TimeNSE2D.h>
#include "TimeNavierStokes.h"
#include <TimeDiscretizations.h>
#include <TimeDiscRout.h>
#include <LoopInfo.h>

using namespace std;

int main(int, char* argv[])
{
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase2D;

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

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);

  // set some parameters for time stepping
  SetTimeDiscParameters(0);

  // create an object of TimeNavierStokes<2> class
  TimeNavierStokes<2> tnse2d(Domain, parmoon_db);
  
  TimeDiscretization& tss = tnse2d.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  
  // assemble everything at the start time
  // this includes assembling of all A's, B's
  // and M's blocks that are necessary
  tnse2d.assemble_initial_time();

  tnse2d.output();

  
  LoopInfo loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1;
  int linear_iteration=0;
  
  TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();
  while(!tss.reached_final_time_step())
  {
    tss.current_step_++;

    // set the time parameters
    tss.set_time_disc_parameters();
    // tau may change depending on the time discretization (adaptive time)
    double tau = tss.get_step_length();
    tss.current_time_ += tss.get_step_length();
    // this is used at several places, e.g., in the example file etc.
    TDatabase::TimeDB->CURRENTTIME += tau;
    Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    
    tnse2d.assemble_matrices_rhs(0);

    LoopInfo loop_info("nonlinear");
    loop_info.print_time_every_step = true;
    loop_info.verbosity_threshold = 1;
    for(unsigned int i=0;; i++)
    {
      if(tnse2d.stop_it(i))
      {
        loop_info.finish(i,tnse2d.get_full_residual());
        linear_iteration +=i;
        loop_info_time.print(linear_iteration, tnse2d.get_full_residual());
        break;
      }
      else
        loop_info.print(i, tnse2d.get_full_residual());

      tnse2d.solve();

      if(tnse2d.imex_scheme())
        break;

      tnse2d.assemble_matrices_rhs(i+1);
    }
    tnse2d.output();
  }
  loop_info_time.finish(linear_iteration, tnse2d.get_full_residual());
  Output::close_file();
  return 0;
}
