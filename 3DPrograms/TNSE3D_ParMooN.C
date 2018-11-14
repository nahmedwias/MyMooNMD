/**
 * @brief Main program for solving a 3D time-dependent Navier Stokes equation using ParMooN.
 * @author Najib Alia
 *
 * Implementation started on 2016/04/15.
 *
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <Time_NSE3D.h>
#include <MeshPartition.h>
#include <Chrono.h>
#include <TetGenMeshLoader.h>
#include <LoopInfo.h>
#include <TimeDiscretizations.h>

#include <sys/stat.h>

#include <TimeDiscRout.h>

using namespace std;

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

// main program
// =======================================================================
int main(int argc, char* argv[])
{
  {
#ifdef _MPI
  //Construct and initialise the default MPI communicator.
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  // Hold mpi rank and size ready, check whether the current processor
  // is responsible for output (usually root, 0).
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(my_rank==0)
  {
    Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
    Output::info("Time_NSE3D", "MPI, using ", size, " processes");
  }
#else
  int my_rank = 0;
  Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
  Output::info("Time_NSE3D", "SEQUENTIAL (or OMP...)");
#endif
  double t_start = GetTime();
  //start a stopwatch which measures time spent in program parts
  Chrono timer;

  // Construct the ParMooN Databases.
  TDatabase Database(argv[1]);
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(ParameterDatabase::default_tetgen_database(), true);
  parmoon_db.read(argv[1]);
  
  //open OUTFILE, this is where all output is written to (additionally to console)
  if(my_rank==0)
  {
    Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  }
  Output::setVerbosity(parmoon_db["verbosity"]);
  
#ifdef _MPI
  TDatabase::ParamDB->Comm = comm;
#endif
  TFEDatabase3D feDatabase;

  // Choose and construct example.
  Example_TimeNSE3D example(parmoon_db);

  // Do the parameter check of the Database.
  check_parameters_consistency_NSE(parmoon_db);
  // =====================================================================
  // set the database values and generate mesh
  // =====================================================================
  TDomain domain(parmoon_db);

  //open OUTFILE, this is where all output is written to (additionally to console)
  if(parmoon_db["problem_type"].is(0))
    parmoon_db["problem_type"] = 6;
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);

  if(my_rank==0) //Only one process should do that.
  {
    parmoon_db.write(Output::get_outfile());
    Database.WriteParamDB(argv[0]);
    Database.WriteTimeDB();
  }
  // Initial refinement and grid collection
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
  domain.print_info("TNSE3D domain");
  // set some parameters for time stepping
  SetTimeDiscParameters(0);
  // Construct an object of the Time_NSE3D-problem type.
#ifdef _MPI
  Time_NSE3D tnse3d(gridCollections, parmoon_db, example, maxSubDomainPerDof);
#else
  Time_NSE3D tnse3d(gridCollections, parmoon_db, example);
#endif
  
  TimeDiscretization& tss = tnse3d.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  
  tnse3d.assemble_initial_time();
  tnse3d.output(tnse3d.current_step_);

  int n_substeps = GetN_SubSteps();

  LoopInfo loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1; // full verbosity
  int linear_iterations = 0; 

  timer.restart_and_print("setting up spaces, matrices and initial assembling");
  double end_time = tss.get_end_time();
  TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();  
  //======================================================================
  // time iteration
  //======================================================================
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    // time measuring during every time iteration
    Chrono timer_timeit;

    //tnse3d.current_step_++;
    tss.current_step_++;

    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    for(int j = 0; j < n_substeps; ++j) // loop over substeps in one time iteration
    {
      // setting the time discretization parameters
      SetTimeDiscParameters(1);
     if( tnse3d.current_step_ == 1 && my_rank==0) // a few output, not very necessary
     {
       Output::print<1>("Theta1: ", TDatabase::TimeDB->THETA1);
       Output::print<1>("Theta2: ", TDatabase::TimeDB->THETA2);
       Output::print<1>("Theta3: ", TDatabase::TimeDB->THETA3);
       Output::print<1>("Theta4: ", TDatabase::TimeDB->THETA4);
     }
      // tau may change depending on the time discretization (adaptive time)
      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;

      if (my_rank==0)
        Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      // prepare the right hand side vector - needed only once per time step
      tnse3d.assemble_rhs();

      // assemble the nonlinear matrices
      tnse3d.assemble_nonlinear_term();

      // prepare the matrices for defect computations and solvers
      tnse3d.assemble_system();
      timer_timeit.restart_and_print("preparation of nonlinear iteration");

      // nonlinear iteration
     LoopInfo loop_info("nonlinear");
     loop_info.print_time_every_step = true;
     loop_info.verbosity_threshold = 1; // full verbosity
     for(unsigned int k=0; ; k++)
      {
        tnse3d.compute_residuals();

        if (my_rank==0) // some outputs
        {
          Output::print<1>("\nNONLINEAR ITERATION :", setw(3), k);
          Output::print<1>("Residuals :", tnse3d.get_residuals());
        }

        // checking residuals and stop conditions
        if(tnse3d.stop_it(k))
        {
          loop_info.finish(k, tnse3d.get_full_residual());
          linear_iterations+=k;
          /// @todo provide all parts of the residual
          /// @todo loop_info restricted to the solver only
          loop_info_time.print(linear_iterations, tnse3d.get_full_residual());
	   break;
         }
         else
           loop_info.print(k, tnse3d.get_full_residual());
 
        tnse3d.solve();

        if(tnse3d.imex_scheme(1))
          continue;

        tnse3d.assemble_nonlinear_term();

        tnse3d.assemble_system();

        timer_timeit.restart_and_print("solving and reassembling in the "
                                        "nonlinear iteration " + 
                                        std::to_string(k));
      }  // end of nonlinear loop

      timer_timeit.restart_and_print(
        "solving the time iteration " +
        std::to_string(TDatabase::TimeDB->CURRENTTIME));

      tnse3d.output(tnse3d.current_step_);
      timer_timeit.print_total_time(
        "time step " + std::to_string(TDatabase::TimeDB->CURRENTTIME));
    } // end of subtime loop
  } // end of time loop
  loop_info_time.finish(linear_iterations, tnse3d.get_full_residual());

  timer.print_total_time("whole solving procedure ");

  for (auto i: gridCollections)
    delete i;
  gridCollections.clear();
  
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  
  Output::close_file();
  
}
#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;
}
