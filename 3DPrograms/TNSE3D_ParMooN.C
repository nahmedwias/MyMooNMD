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
#include <Output3D.h>
#include <ChannelFlowRoutines.h>

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
  double t_start = GetTime();
  Output::print("<<<<< Running ParMooN: TNSE3D Main Program >>>>>");

  //start a stopwatch which measures time spent in program parts
  Chrono chrono_parts;

  // Construct the ParMooN Databases.
  TDatabase Database;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(ParameterDatabase::default_tetgen_database(), true);

  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

  // Choose and construct example.
  Example_TimeNSE3D example(parmoon_db);


  // set parameters for particular examples
  if(parmoon_db["example"].is(7))
  {
    ChannelTau180::setParameters(parmoon_db);
  }
  // Do the old parameter check of the Database.
  Database.CheckParameterConsistencyNSE();

#ifdef _MPI
  //Construct and initialize the default MPI communicator and store it.
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  TDatabase::ParamDB->Comm = comm;

  // Hold mpi rank and size ready, check whether the current processor
  // is responsible for output (usually root, 0).
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  int my_rank = 0;
#endif

  TFEDatabase3D feDatabase;
  // =====================================================================
  // set the database values and generate mesh
  // =====================================================================
  // Construct domain, thereby read in controls from the input file.
  TDomain domain(argv[1],parmoon_db);

  //open OUTFILE, this is where all output is written to (additionally to console)
  if(parmoon_db["problem_type"].is(0))
    parmoon_db["problem_type"] = 6;
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);

  //open OUTFILE, this is where all output is written to (additionally to console)
  if(my_rank==0)
  {
    Output::set_outfile(parmoon_db["outfile"]);
  }
  Output::setVerbosity(parmoon_db["verbosity"]);

  if(my_rank==0) //Only one process should do that.
  {
    Database.WriteParamDB(argv[0]);
    Database.WriteTimeDB();
  }
  // Do initial domain refinement
  size_t n_ref = domain.get_n_initial_refinement_steps();
  // set the coordinates and make the boundary parameters
  // conistent on the coarse grid
  std::vector<TCollection *> collArray(n_ref+1);
  TCollection *coll;
  if(parmoon_db["example"].is(7))
  {
    coll = domain.GetCollection(It_Finest,0);
    collArray[0] = coll;
    ChannelTau180::setZCoordinates(coll, 0);
    domain.MakeBdParamsConsistent(coll);
    ChannelTau180::checkZCoordinates(coll, 0);
  }
  for(size_t i = 0; i < n_ref; i++)
  {
    domain.RegRefineAll();
    // set the coordinates and make the boundary parameters
    // conistent on the finer grids
    if(parmoon_db["example"].is(7))
    {
      coll = domain.GetCollection(It_Finest,0);
      collArray[i+1] = coll;
      // set the z Coordinates
      ChannelTau180::setZCoordinates(coll, i+1);
      domain.MakeBdParamsConsistent(coll);
      ChannelTau180::checkZCoordinates(coll, i+1);
    }
  }
  // set the refinement descriptor and the periodic joints
  // for particular examples
  if(parmoon_db["example"].is(7))
  {
    n_ref +=1;
    for(size_t i=0; i<n_ref; i++)
    {
      ChannelTau180::setRefineDesc(collArray[i]);
      ChannelTau180::setPeriodicFaceJoints(collArray[i]);
    }
  }
  
  coll = domain.GetCollection(It_Finest,0);
  TOutput3D output(0,0,0,0,std::addressof(domain),coll);
  
  output.WriteVtk(parmoon_db["output_basename"]);

  // write grid into an Postscript file (before partitioning)
  if(parmoon_db["output_write_ps"] && my_rank==0)
    domain.PS("Domain.ps", It_Finest, 0);

#ifdef _MPI
  // Partition the by now finest grid using Metis and distribute among processes.

  // 1st step: Analyse interfaces and create edge objects,.
  domain.GenerateEdgeInfo();

  // 2nd step: Call the mesh partitioning.

  int maxCellsPerVertex;
  //do the actual partitioning, and examine the return value
  if ( Partition_Mesh3D(comm, &domain, maxCellsPerVertex) == 1)
  {
    /** \todo It is a known issue that Metis does not operate entirely
     * deterministic here. On a coarse grid (as if doing the partitioning on
     * the coarsest grid of a multgrid hierarchy) it can happen, that
     * one process goes entirely without own cells to work on.
     * The workarounds which were used so far (setting another metis type,
     * doing more uniform steps than the user requested ) are unsatisfactoy
     * imo. So this is a FIXME
     *
     * One can reproduce the problem when using the cd3d multigrid test program
     * in MPI and setting LEVELS to 3 and UNIFORM_STEPS to 1.
     *
     * Of course the same issue occurs if one calls this upon too small
     * a grid with too many processes.
     *
     */
    ErrThrow("Partitioning did not succeed.");
  }

  // 3rd step: Generate edge info anew
  //(since distributing changed the domain).
  domain.GenerateEdgeInfo();

  // calculate largest possible number of processes which share one dof
  int maxSubDomainPerDof = MIN(maxCellsPerVertex, size);

#endif

  //print information on the mesh partition on the finest grid
  domain.print_info("TNSE3D domain");
  // set some parameters for time stepping
  SetTimeDiscParameters(0);
  // Construct an object of the Time_NSE3D-problem type.
#ifdef _MPI
  Time_NSE3D tnse3d(domain, parmoon_db, example, maxSubDomainPerDof);
#else
  Time_NSE3D tnse3d(domain, parmoon_db, example);
#endif
  // assemble all matrices and right hand side at start time
  // it assembles A's, B's and M's blocks. Nonlinear blocks are
  // added in the loops thanks to assemble_nonlinear_term()
  if(parmoon_db["example"].is(7))
  {
    ChannelTau180::GetCoordinatesOfDof(tnse3d);
  }
  tnse3d.assemble_initial_time();
  
  double end_time = TDatabase::TimeDB->ENDTIME;
  tnse3d.current_step_ = 0;

  int n_substeps = GetN_SubSteps();

  int image = 0;

  chrono_parts.print_time(std::string("setting up spaces, matrices and initial assembling"));
  chrono_parts.reset();
  TDatabase::TimeDB->CURRENTTIME = 0.0;
  if(parmoon_db["example"].is(7))
  {
    TDatabase::ParamDB->INTERNAL_MEAN_COMPUTATION = 1;
    ChannelTau180::computeMeanVelocity(tnse3d);
  }
  //======================================================================
  // time iteration
  //======================================================================
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    Chrono chrono_timeit;

    tnse3d.current_step_++;

    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    for(int j = 0; j < n_substeps; ++j) // loop over substeps in one time iteration
    {
      // setting the time discretization parameters
      SetTimeDiscParameters(1);
//      if( tnse3d.current_step_ == 1 && my_rank==0) // a few output, not very necessary
//      {
//        Output::print<1>("Theta1: ", TDatabase::TimeDB->THETA1);
//        Output::print<1>("Theta2: ", TDatabase::TimeDB->THETA2);
//        Output::print<1>("Theta3: ", TDatabase::TimeDB->THETA3);
//        Output::print<1>("Theta4: ", TDatabase::TimeDB->THETA4);
//      }
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
          break;

        Chrono chrono_nonlinit;

        tnse3d.solve();

        if(tnse3d.imex_scheme(1))
          continue;

        tnse3d.assemble_nonlinear_term();

        tnse3d.assemble_system();

        chrono_nonlinit.print_time(std::string("solving and reassembling in the nonlinear iteration ") + std::to_string(k));
      }  // end of nonlinear loop

      chrono_timeit.print_time(std::string("solving the time iteration ") + std::to_string(tau));

      tnse3d.output(tnse3d.current_step_,image);
      
      if(parmoon_db["example"].is(7))
      {
        if (TDatabase::TimeDB->CURRENTTIME>=TDatabase::TimeDB->T0)
        {
          if (TDatabase::ParamDB->INTERNAL_MEAN_COMPUTATION == 0)
          {
            TDatabase::ParamDB->INTERNAL_MEAN_COMPUTATION = 1;
            TDatabase::TimeDB->T0 = TDatabase::TimeDB->CURRENTTIME;
          }
        }
        ChannelTau180::computeMeanVelocity(tnse3d);
      }

    } // end of subtime loop
  } // end of time loop

  chrono_parts.print_time(std::string("whole solving procedure "));

  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================

  Output::close_file();
#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;
}
