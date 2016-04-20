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
#include <NSE3D.h>
#include <MeshPartition.h>
#include <Chrono.h>

#include <sys/stat.h>

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
#ifdef _MPI
  //Construct and initialise the default MPI communicator.
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
#endif

  Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
  Chrono chrono_parts;
  chrono_parts.print_time(std::string("program start"));
  // Construct the ParMooN Databases.
  TDatabase Database;

#ifdef _MPI
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

  // Construct domain, thereby read in controls from the input file.
  TDomain domain(argv[1]);

  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);

  if(my_rank==0) //Only one process should do that.
    Database.WriteParamDB(argv[0]);

  // Do a makeshift parameter check and the old parameter check of the Database.
  NSE3D::check_parameters();
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

  // Choose example according to the value of
  // TDatabase::ParamDB->EXAMPLE and construct it.
  Example_NSE3D example;

  // Construct an object of the NSE3D-problem type.
#ifdef _MPI
  NSE3D nse3d(domain, example, maxSubDomainPerDof);
#else
  NSE3D nse3d(domain, example);
#endif

  // assemble all matrices and right hand side
  nse3d.assemble_linear_terms();
  nse3d.stop_it(0);  // check initial residuals

  chrono_parts.print_time(std::string("setting up spaces, matrices, linear assemble"));
  chrono_parts.reset();

  //======================================================================
  for(unsigned int k=1;; k++)
  {
    Chrono chrono_nonlinit;

    if(my_rank==0)
    {
     Output::print("\nNONLINEAR ITERATION ", setw(3), k-1);
     Output::print(" residuals ",nse3d.get_residuals());
    }

    // solve the system
    nse3d.solve();

    nse3d.assemble_non_linear_term();

    chrono_nonlinit.print_time(std::string("nonlinear iteration ") + std::to_string(k-1));

    // checking residuals
    if(nse3d.stop_it(k))
      break;

  }

  chrono_parts.print_time(std::string("solving procedure "));

  nse3d.output();


  Output::close_file();


#ifdef _MPI
  MPI_Finalize();
#endif
  
  return 0;
}
