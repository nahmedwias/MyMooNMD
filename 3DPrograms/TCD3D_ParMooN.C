/**
 * @brief Main program for solving a 3D time dependent convection diffusion reaction problem
 * @author Naveed Ahmed
 * @history 14.06.16
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <Time_CD3D.h>
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
  TDatabase Database;
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
  TDomain domain(parmoon_db, argv[1]);
  
  if(i_am_root)
  {
    parmoon_db.write(Output::get_outfile());
    Database.WriteParamDB(argv[1]);
  }
  

  // Intial refinement and grabbing of grids for multigrid.
#ifdef _MPI
  int maxSubDomainPerDof = 0;
#endif
  std::list<TCollection* > gridCollections
  = domain.refine_and_get_hierarchy_of_collections(
      parmoon_db
  #ifdef _MPI
      , maxSubDomainPerDof
  #endif
      );
  //print information on the mesh partition on the finest grid
  domain.print_info("TCD3D domain");
  
  // set some parameters for time stepping
  SetTimeDiscParameters(0);

  // Choose example according to the value of "example" in the given database
  Example_TimeCD3D example(parmoon_db);
  
  timer.restart_and_print("setup(domain, example, database)");
  // create an object of the class Time_CD3D
#ifdef _MPI
  Time_CD3D tcd3d(gridCollections, parmoon_db, example, maxSubDomainPerDof);
#else
  Time_CD3D tcd3d(gridCollections, parmoon_db, example);
#endif
  timer.restart_and_print("constructing Time_CD3D object");
  
  // assemble the matrices and right hand side at the start time
  tcd3d.assemble_initial_time();
  int step = 0;
  tcd3d.output(step);
  timer.restart_and_print("initial assembling");
  Chrono timer_solve;
  timer_solve.stop();
  
  double start_time = parmoon_db["time_start"];
  double end_time = parmoon_db["time_end"];
  TDatabase::TimeDB->CURRENTTIME = start_time;
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    step++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    SetTimeDiscParameters(1);

    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    TDatabase::TimeDB->CURRENTTIME += tau;
    
    if(i_am_root)
      Output::print<1>("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    tcd3d.assemble();
    timer_solve.start();
    tcd3d.solve();
    timer_solve.restart_and_print("solving in step t=" 
                               +std::to_string(TDatabase::TimeDB->CURRENTTIME));
    timer_solve.stop();
    tcd3d.descale_stiffness();
    
    tcd3d.output(step);
    timer.restart_and_print(
      "time step (t=" + std::to_string(TDatabase::TimeDB->CURRENTTIME)+ ")");
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
