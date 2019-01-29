/**
 * @brief Main program for solving a 3D time dependent convection diffusion reaction problem
 * @author Naveed Ahmed
 * @history 14.06.16
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include "TimeConvectionDiffusion.h"
#include <TimeDiscRout.h>
#include <Chrono.h>
#include <TimeDiscretizations.h>

#include <sys/stat.h>

using namespace std;

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

int main(int argc, char *argv[])
{
#ifdef _MPI
  MPI_Init(&argc, &argv);
#endif
  {
  Chrono timer;
  // Construct the ParMooN Databases.
  TDatabase Database(argv[1]);
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(TimeDiscretization::default_TimeDiscretization_database());
  parmoon_db.read(argv[1]);
  
  bool i_am_root = true;
#ifdef _MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_rank, size;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &size);
  i_am_root = (my_rank == 0);
#endif
  if(i_am_root)
    Output::print("<<<<< Running ParMooN: Time_CD3D Main Program >>>>>");
#ifdef _MPI
  Output::info("Time_CD3D", "MPI, using ", size, " processes");
#endif
  
  // open outfile, this is where all output is written (additionally to console)
  if(i_am_root)
    Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  TFEDatabase3D feDatabase;
  TDomain domain(parmoon_db);
  
  if(i_am_root)
  {
    parmoon_db.write(Output::get_outfile());
    Database.WriteParamDB(argv[1]);
  }
  

  // Intial refinement
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  //print information on the mesh partition on the finest grid
  domain.print_info("TCD3D domain");
  
  // set some parameters for time stepping
  SetTimeDiscParameters(0);

  timer.restart_and_print("setup(domain, example, database)");
  // create an object of the class Time_CD3D
  TimeConvectionDiffusion<3> tcd3d(domain, parmoon_db);
  timer.restart_and_print("constructing Time_CD3D object");
  
  TimeDiscretization& tss = tcd3d.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.current_time_ = parmoon_db["time_start"];
  
  // assemble the matrices and right hand side at the start time
  tcd3d.assemble_initial_time();
  timer.restart_and_print("initial assembling");
  Chrono timer_solve;
  timer_solve.stop();
  
  double start_time = parmoon_db["time_start"];
  TDatabase::TimeDB->CURRENTTIME = start_time;
  while(!tss.reached_final_time_step())
  {
    tss.current_step_++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = tss.current_time_;
    tss.set_time_disc_parameters();
    SetTimeDiscParameters(1);
    tss.current_time_ += tss.get_step_length();
    // this is used at several places, e.g., in the example file etc.
    TDatabase::TimeDB->CURRENTTIME += tss.get_step_length();
    
    if(i_am_root)
      Output::print<1>("\nCURRENT TIME: ", tss.current_time_);
    tcd3d.assemble();
    timer_solve.start();
    tcd3d.solve();
    timer_solve.restart_and_print("solving in step t=" 
                               +std::to_string(tss.current_time_));
    timer_solve.stop();
    
    tcd3d.output();
    timer.restart_and_print(
      "time step (t=" + std::to_string(tss.current_time_)+ ")");
  }
  
  timer_solve.print_total_time("accumulated solver time");
  
  if(i_am_root)
    Output::print("<<<<< ParMooN Finished: TCD3D Main Program >>>>>");

  timer.print_total_time("TCD3D_ParMooN program");
  Output::close_file();
  }
#ifdef _MPI
  MPI_Finalize();
#endif
}
