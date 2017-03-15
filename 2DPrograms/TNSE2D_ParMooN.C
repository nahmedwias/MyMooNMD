#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Example_TimeNSE2D.h>
#include <Time_NSE2D_Merged.h>
#include <MixingLayerSlip.h>
#include <TimeDiscRout.h>
#include <LoopInfo.h>

using namespace std;

int main(int argc, char* argv[])
{
  double t_start=GetTime();
  TDatabase Database;
  TFEDatabase2D FEDatabase2D;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  //std::string time_disc = parmoon_db["time_discretization"];
  //parmoon_db.merge(ParameterDatabase::default_time_database(time_disc), true);

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  TDomain Domain(argv[1], parmoon_db);

  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

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
  //TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters(0);

  Example_TimeNSE2D example( parmoon_db );
  // create an object of Time_NSE2D class
  Time_NSE2D_Merged tnse2d(Domain, parmoon_db, example);
  if(parmoon_db["example"].is(3))
  {
    Mixing_layer::fill_arrays(tnse2d);
  }
  
  tnse2d.time_stepping_scheme.current_step_ = 0;
  tnse2d.time_stepping_scheme.set_time_disc_parameters();
  
  // assemble everything at the start time
  // this includes assembling of all A's, B's
  // and M's blocks that are necessary
  tnse2d.assemble_initial_time();

  tnse2d.output(tnse2d.time_stepping_scheme.current_step_);

  double end_time = TDatabase::TimeDB->ENDTIME;
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    tnse2d.time_stepping_scheme.current_step_++;

    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    // set the time parameters
    tnse2d.time_stepping_scheme.set_time_disc_parameters();
    double tau = parmoon_db["time_step_length"];
    TDatabase::TimeDB->CURRENTTIME += tau;
    Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    for(unsigned int i=0;; i++)
    {
      tnse2d.assemble_matrices_rhs(i);
      if(tnse2d.stopIte(i))
        break;
      tnse2d.solve();
    }
    tnse2d.output(tnse2d.time_stepping_scheme.current_step_);
    if(parmoon_db["example"].is(3))
    {
      Mixing_layer::compute_mean_velocity(tnse2d);
    }
  }
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();
  return 0;
}
