/**
 * @brief Main program for solving a 3D stationary Navier Stokes equation using ParMooN.
 * @author Sashikumaar Ganesan, Clemens Bartsch
 *
 * Implementation started on 2015/01/23. Rework since 2015/12/7.
 *
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include "NavierStokes.h"
#include <Chrono.h>
#include <LoopInfo.h>

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

int main(int argc, char* argv[])
{
#ifdef _MPI
  //Construct and initialise the default MPI communicator.
  MPI_Init(&argc, &argv);
#endif
  {
#ifdef _MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(my_rank==0)
    {
      Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
      Output::info("NSE3D", "MPI, using ", size, " processes");
    }
#else
    int my_rank = 0;
    Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
    Output::info("NSE3D", "SEQUENTIAL (or OMP...)");
#endif

    //start a stopwatch which measures time spent in program parts
    Chrono timer;

    // Construct the ParMooN Databases.
    TDatabase Database(argv[1]);
    ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
    parmoon_db.read(argv[1]);

    //open OUTFILE, this is where all output is written to (addionally to console)
    if(my_rank==0)
    {
      Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
    }
    Output::setVerbosity(parmoon_db["verbosity"]);
    
#ifdef _MPI
    TDatabase::ParamDB->Comm = comm;
#endif

    TFEDatabase3D feDatabase;
    TDomain domain(parmoon_db);

    if(my_rank==0) //Only one process should do that.
    {
      parmoon_db.write(Output::get_outfile());
      Database.WriteParamDB(argv[0]);
    }

    // Do the parameter check of the Database.
    check_parameters_consistency_NSE(parmoon_db);

    // Intial refinement
    domain.refine_and_get_hierarchy_of_collections(parmoon_db);

    //print information on the mesh partition on the finest grid
    domain.print_info("NSE3D domain");

    timer.restart_and_print("setup(domain, example, database)");
    // Construct an object of the NSE3D-problem type.
    NavierStokes<3> nse3d(domain, parmoon_db);
    timer.restart_and_print("constructing NSE3D object");
    
    // assemble all matrices and right hand side
    nse3d.assemble_linear_terms();
    nse3d.stop_it(0);  // check initial residuals

    LoopInfo loop_info("nonlinear");
    loop_info.print_time_every_step = true;
    loop_info.verbosity_threshold = 1; // full verbosity
    if(my_rank==0)
      loop_info.print(0, nse3d.get_full_residual());

    timer.restart_and_print("assembling linear terms");

    //Timer which measures time spent in solve() method solely.
    Chrono timer_sol;
    timer_sol.stop();

    //======================================================================
    for(unsigned int k=1;; k++)
    {
      if(my_rank == 0)
        Output::print<3>(); // new line for a new nonlinear iteration
      // solve the system
      timer_sol.start();
      nse3d.solve();
      timer_sol.stop();

      //no nonlinear iteration for Stokes problem
      if(parmoon_db["problem_type"].is(3) || parmoon_db["problem_type"].is(7))
        break;

      nse3d.assemble_nonlinear_term();

      // checking residuals
      if(nse3d.stop_it(k))
      {
        loop_info.finish(k, nse3d.get_full_residual());
        break;
      }
      else
        loop_info.print(k, nse3d.get_full_residual());

    } // end for k
    timer.restart_and_print("nonlinear loop");
    
    timer_sol.print_total_time("solver only");

    nse3d.output();
    timer.restart_and_print("output");
    timer.print_total_time("NSE3D_ParMooN program");

    if(my_rank==0)
      Output::print("<<<<< ParMooN Finished: NSE3D Main Program >>>>>");

    if(my_rank == 0)
      Output::close_file();
  }

#ifdef _MPI
  MPI_Finalize();
#endif
  
  return 0;
}

