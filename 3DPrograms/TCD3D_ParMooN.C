/**
 * @brief Main program for solving a 3D time dependent convection diffusion reaction problem
 * @author Naveed Ahmed
 * @history 14.06.16
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <Time_CD3D.h>
#include <MeshPartition.h>
#include <TimeDiscRout.h>

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
  // Construct the ParMooN Databases.
  TDatabase Database;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  
#ifdef _MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  Output::print("<<<<< Running ParMooN: Time_CD3D Main Program >>>>>");
  Output::info("Time_CD3D", "MPI, using ", size, " processes");
#else
  int my_rank = 0;
  Output::print("<<<<< Running ParMooN: Time_CD3D Main Program >>>>>");
  Output::info("Time_CD3D", "SEQUENTIAL");
#endif
  
  TFEDatabase3D feDatabase;
  TDomain domain(argv[1], parmoon_db);
  
  // open outfile, this is where all output is written (additionally to console)
  if(my_rank==0)
    Output::set_outfile(parmoon_db["outfile"]);
  
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  if(my_rank==0)
    Database.WriteParamDB(argv[1]);
  // reading mesh
  domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"]);
  
  // split the number of refinement steps - some have to be done before,
  // some after the domain partitioning
  int n_ref_total = domain.get_n_initial_refinement_steps();
  size_t n_ref_after =  parmoon_db["multigrid_n_levels"];
  n_ref_after -= 1;
  int n_ref_before =  n_ref_total - n_ref_after;
  if(n_ref_before < 0)
  {
    ErrThrow("Number of multigrid levels is greater than number of refinement "
        "levels. Garbage in, garbage out.")
  }

  for(int i = 0; i < n_ref_before; i++)
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
    /** \todo It is a knwon issue, that Metis does not operate entirely
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
  int maxSubDomainPerDof = MIN(maxCellsPerVertex, mpiSize);
  
#endif
  
  
  // Collect those Collections which will be used in multigrid.
  // ("Collection" in ParMooN means a set of grid cells which form a
  // specific computational domain).
  // This must be done here instead of deep inside CD3D, because the Domain_Crop
  // method disables the use of sensible Collection Iterators. The only possibility
  // to get a certain level of cells is to grab it the moment when it's the finest...
  std::list<TCollection* > gridCollections;
  gridCollections.push_front(domain.GetCollection(It_Finest, 0));

  for(size_t level=0; level <  n_ref_after; level++)
  {
    domain.RegRefineAll();
#ifdef _MPI
    domain.GenerateEdgeInfo();  // has to be called anew after every refinement step
    Domain_Crop(comm, &domain); // remove unwanted cells in the halo after refinement
#endif
    // Grab collection.
    gridCollections.push_front(domain.GetCollection(It_Finest, 0));
  }
  
  // set some parameters for time stepping
  SetTimeDiscParameters(0);
  
  //print information on the mesh partition on the finest grid
  domain.print_info("TCD3D domain");

  // Choose example according to the value of
  Example_TimeCD3D example(parmoon_db["example"],parmoon_db);
  
  // create an object of the class Time_CD3D
#ifdef _MPI
  Time_CD3D tcd3d(gridCollections, parmoon_db, example, maxSubDomainPerDof);
#else
  Time_CD3D tcd3d(gridCollections, parmoon_db, example);
#endif

  // assemble the matrices and right hand side at the start time
  tcd3d.assemble_initial_time();
  int step = 0, image=0;
  
  while(TDatabase::TimeDB->CURRENTTIME < TDatabase::TimeDB->ENDTIME - 1E-10)
  {
    step ++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    SetTimeDiscParameters(1);

    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    TDatabase::TimeDB->CURRENTTIME += tau;
    
    Output::print<1>("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    tcd3d.assemble();
    tcd3d.solve();
    tcd3d.descale_stiffness();
    
    tcd3d.output(step,image);
  }
  
  if(my_rank==0)
    Output::print("<<<<< ParMooN Finished: NSE3D Main Program >>>>>");

  Output::close_file();
  
  #ifdef _MPI
  MPI_Finalize();
#endif
 
  return 0;
}
