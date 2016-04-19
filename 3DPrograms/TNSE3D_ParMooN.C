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
  double t_start=GetTime();
  Output::print("<<<<< Running ParMooN: TNSE3D Main Program >>>>>");

  // Construct the ParMooN Databases.
  TDatabase Database;
  TFEDatabase3D feDatabase;

#ifdef _MPI
  //Construct and initialise the default MPI communicator and store it.
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  TDatabase::ParamDB->Comm = comm;

  // Hold mpi rank and size ready, check whether the current processor
  // is responsible for output (usually root, 0).
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  int my_rank = 0;
#endif

  // =====================================================================
  // set the database values and generate mesh
  // =====================================================================
  // Construct domain, thereby read in controls from the input file.
  TDomain domain(argv[1]);

  //open OUTFILE, this is where all output is written to (additionally to console)
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 6;
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);

  if(my_rank==0) //Only one process should do that.
    {
    Database.WriteParamDB(argv[0]);
    Database.WriteTimeDB();
    }

  // Do a makeshift parameter check and the old parameter check of the Database.
  // TODO Adapt the check_parameters() method to the class TNSE3D
  // NSE3D::check_parameters();
  Database.CheckParameterConsistencyNSE();

  // Read in geometry and initialize the mesh.
  domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);

  // Do initial regular grid refinement.
  for(int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
  {
    domain.RegRefineAll();
  }

  // Write grid into a postscript file (before partitioning)
  if(TDatabase::ParamDB->WRITE_PS && my_rank == 0)
  {
    domain.PS("Domain.ps", It_Finest, 0);
  }

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

  //print information on the mesh partitioning
  Output::print("Process ", my_rank, ". N_OwnCells: ",
                domain.GetN_OwnCells(),
                ". N_HaloCells: ",
                domain.GetN_HaloCells());
#endif

  // Create output directory, if not already existing.
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  }

  // set some parameters for time stepping
  SetTimeDiscParameters(0);

  // Choose example according to the value of
  // TDatabase::ParamDB->EXAMPLE, e.g. 101, and construct it.
  Example_NSE3D example;

  // Construct an object of the Time_NSE3D-problem type.
//#ifdef _MPI
//  Time_NSE3D nse3d(domain, example, maxSubDomainPerDof);
//#else
  Time_NSE3D tnse3d(domain, example);
//#endif

  // assemble all matrices and right hand side at start time
  // it assembles A's, B's and M's blocks. Nonlinear blocks are
  // added in the loops thanks to assemble_nonlinear()
  tnse3d.assemble_initial_time();

  double end_time = TDatabase::TimeDB->ENDTIME;
  int step = 0;
  int n_substeps = GetN_SubSteps();

  int image = 0;

//  TODO Check if this is necessary. Done in NSE3D but not in TNSE2D
//  nse3d.stop_it(0);
//  // check initial residuals

  //======================================================================
  // time iteration
  //======================================================================
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    step++;
    // Output::print("memory before ":, GetMemory());
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    for(int j = 0; j < n_substeps; ++j) // loop over substeps in one time iteration
    {
      // setting the time discretization parameters
      SetTimeDiscParameters(1);
      if( step == 1) // a few output, not very necessary
      {
        Output::print<1>("Theta1: ", TDatabase::TimeDB->THETA1);
        Output::print<1>("Theta2: ", TDatabase::TimeDB->THETA2);
        Output::print<1>("Theta3: ", TDatabase::TimeDB->THETA3);
        Output::print<1>("Theta4: ", TDatabase::TimeDB->THETA4);
      }
      // tau may change depending on the time discretization (adaptive time)
      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;

      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      // prepare the right hand side vector - needed only once per time step
      tnse3d.assemble_rhs();

      // assemble the nonlinear matrices
      tnse3d.assemble_nonlinear_term();

      // prepare the matrices for defect computations and solvers
      tnse3d.assemble_system();

      for(unsigned int k=0; ; k++)
      {
        // checking residuals
        if(tnse3d.stop_it(k))
        {
          break;
        }

        Output::print<1>("\nNONLINEAR ITERATION :", setw(3), k);
        Output::print<1>("Residuals :", tnse3d.get_residuals());

        tnse3d.solve();

        tnse3d.assemble_nonlinear_term();

        tnse3d.assemble_system();
      }  // end of nonlinear loop

      tnse3d.output();

    } // end of subtime loop
  } // end of time loop

  Output::close_file();
  /*#ifdef _MPI
  //  MPI_Finalize();
  #endif*/
  return 0;
}
