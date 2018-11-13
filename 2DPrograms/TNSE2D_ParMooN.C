#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Example_TimeNSE2D.h>
#include <Time_NSE2D.h>
#include <TimeDiscretizations.h>
#include <TimeDiscRout.h>
#include <LoopInfo.h>

using namespace std;

int main(int argc, char* argv[])
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

  // refine grid up to the coarsest level
  size_t n_ref = Domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
    Domain.RegRefineAll();

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);

  // set some parameters for time stepping
  SetTimeDiscParameters(0);

  Example_TimeNSE2D example( parmoon_db );
  // create an object of Time_NSE2D class
  Time_NSE2D tnse2d(Domain, parmoon_db, example);
  
  TimeDiscretization& tss = tnse2d.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  
  // assemble everything at the start time
  // this includes assembling of all A's, B's
  // and M's blocks that are necessary
  tnse2d.assemble_initial_time();

  tnse2d.output(tss.current_step_);

  
  LoopInfo loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1;
  int linear_iteration=0;
  
  double end_time = tss.get_end_time();
  TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    tss.current_step_++;

    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    // set the time parameters
    tss.set_time_disc_parameters();
    double tau = parmoon_db["time_step_length"];
    TDatabase::TimeDB->CURRENTTIME += tau;
    Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    
    tnse2d.assemble_matrices_rhs(0);

    LoopInfo loop_info("nonlinear");
    loop_info.print_time_every_step = true;
    loop_info.verbosity_threshold = 1;
    for(unsigned int i=0;; i++)
    {
      if(tnse2d.stopIte(i))
      {
        loop_info.finish(i,tnse2d.getFullResidual());
        linear_iteration +=i;
        loop_info_time.print(linear_iteration, tnse2d.getFullResidual());
        break;
      }
      else
        loop_info.print(i, tnse2d.getFullResidual());

      tnse2d.solve();

      if(tnse2d.imex_scheme(1))
        continue;

      tnse2d.assemble_matrices_rhs(i+1);
    }
    tnse2d.output(tss.current_step_);
  }
  loop_info_time.finish(linear_iteration, tnse2d.getFullResidual());
  Output::close_file();
  return 0;
}
