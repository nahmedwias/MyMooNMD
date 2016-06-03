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

//project specific
#include <CoiledPipe.h>

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

//project specific declaration
struct derived_properties;
int main(int argc, char* argv[])
{
#ifdef _MPI
  //Construct and initialise the default MPI communicator.
  MPI_Init(&argc, &argv);
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
  Chrono chrono_parts;

  // Construct the ParMooN Databases.
  TDatabase Database;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

#ifdef _MPI
  TDatabase::ParamDB->Comm = comm;
#endif

  TFEDatabase3D feDatabase;

  // Construct domain, thereby read in controls from the input file.
  TDomain domain(argv[1], parmoon_db);

  //open OUTFILE, this is where all output is written to (additionally to console)
  if(my_rank==0)
  {
    Output::set_outfile(parmoon_db["outfile"]);
  }
  Output::setVerbosity(parmoon_db["verbosity"]);

  if(my_rank==0) //Only one process should do that.
    Database.WriteParamDB(argv[0]);

  // Do a makeshift parameter check and the old parameter check of the Database.
  NSE3D::check_parameters();
  Database.CheckParameterConsistencyNSE();

  // project specific: prepare the coiled geometry
  size_t n_twists                 = parmoon_db["twisted_pipe_n_twists"];
  size_t n_segments_per_twist     = parmoon_db["twisted_pipe_n_segments_per_twist"];
  double l_inflow                 = parmoon_db["twisted_pipe_l_inflow"];
  size_t n_segments_inflow        = parmoon_db["twisted_pipe_n_segments_inflow"];
  double l_outflow                = parmoon_db["twisted_pipe_l_outflow"];
  size_t n_segments_outflow       = parmoon_db["twisted_pipe_n_segments_outflow"];
  double tube_radius              = parmoon_db["twisted_pipe_tube_radius"];
  double twist_radius             = parmoon_db["twisted_pipe_twist_radius"];
  double space_between_twists     = parmoon_db["twisted_pipe_space_between_twists"];

  CoiledPipe::set_up_geoconsts(
      n_twists,
      n_segments_per_twist,
      l_inflow,
      n_segments_inflow,
      l_outflow,
      n_segments_outflow,
      tube_radius,
      twist_radius,
      space_between_twists
  );
  double drift_x = 0;
  double drift_y = 0;
  double drift_z =CoiledPipe::GeoConsts::l_tube;
  
  // Choose and construct example - in this project, this has to be done before
  // initiing the Domain, as the example itself inluences the geometry by setting
  // a global parameter (which is awful).
  Example_NSE3D example(parmoon_db["example"]);

  // Read in geometry and initialize the mesh.
  domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"],
              drift_x, drift_y, drift_z,
              CoiledPipe::GeoConsts::segment_marks);

  // Initial domain refinement
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
    domain.RegRefineAll();

  // write grid into an Postscript file
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

  //print information on the mesh partitioning
  domain.print_info(std::string("coiled tube"));

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // Construct an object of the NSE3D-problem type.
#ifdef _MPI
  NSE3D nse3d(domain, parmoon_db, example, maxSubDomainPerDof);
#else
  NSE3D nse3d(domain, parmoon_db, example);
#endif

  // assemble all matrices and right hand side
  nse3d.assemble_linear_terms();
  nse3d.stop_it(0);  // check initial residuals

  chrono_parts.print_time(std::string("setting up spaces, matrices, linear assemble"));
  chrono_parts.reset();

  //======================================================================
  for(unsigned int k=1;; k++)
  {
    nse3d.output(k);

    Chrono chrono_nonlinit;

    if(my_rank==0)
    {
     Output::info("NONLINEAR ITERATION ", setw(3), k-1);
     Output::dash(" residuals ",nse3d.get_residuals());
    }

    // solve the system
    nse3d.solve();

    //no nonlinear iteration for Stokes problem
    if(parmoon_db["problem_type"].is(3))
      break;
    
    nse3d.assemble_non_linear_term();

    chrono_nonlinit.print_time(std::string("nonlinear iteration ") + std::to_string(k-1));

    // checking residuals
    if(nse3d.stop_it(k))
      break;

  }

  chrono_parts.print_time(std::string("solving procedure "));

  nse3d.output();

  if(my_rank==0)
    Output::print("<<<<< ParMooN Finished: NSE3D Main Program >>>>>");

  if(my_rank == 0)
    Output::close_file();


#ifdef _MPI
  MPI_Finalize();
#endif
  
  return 0;
}
