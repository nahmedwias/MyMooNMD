/**
 * @brief Main program for solving a 3D stationary scalar equation using ParMooN.
 * @author Sashikumaar Ganesan, Clemens Bartsch
 *
 * Implementation started on 2015/01/23. Rework since 2015/10/19.
 *
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <CD3D.h>
#include <MeshPartition.h>

#include <sys/stat.h>

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

/**
 * @brief Get the current time, when profiling is activated.
 *
 * @return The current value of time as a float, when
 * TDatabase::ParamDB->timeprofiling is true, 0 otherwise.
 */
double getTime()
{
  if(TDatabase::ParamDB->timeprofiling)
  {
  #ifdef _MPI
    return MPI_Wtime();
  #else
    return GetTime();
  #endif
  }
  return 0.;
}

// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // Construct the ParMooN Databases.
  TDatabase Database;

#ifdef _MPI
  //Construct and initialise the default MPI communicator and store it.
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  TDatabase::ParamDB->Comm = comm;

  // Hold mpi rank and size ready, check whether the current processor
  // is responsible for output (usually root, 0).
  int mpiRank, mpiSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  bool iAmOutRank= (mpiRank == TDatabase::ParamDB->Par_P0);
#endif

  TFEDatabase3D feDatabase;

  // Construct domain, thereby read in controls from the input file.
  TDomain domain(argv[1]);

  //do a makeshift parameter check
  CD3D::checkParameters();

  Output::set_outfile(TDatabase::ParamDB->OUTFILE);

  #ifdef _MPI
  if(iAmOutRank) //Only one process should do that.
  #endif
    Database.WriteParamDB(argv[0]);

  // Read in geometry and initialize the mesh. (See code of domain.Init
  // for usage of hard-coded example meshes)
  domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);

  // Do initial regular grid refinement.
  for(int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
  {
	  domain.RegRefineAll();
  }

  // Write grid into a postscript file
  if(TDatabase::ParamDB->WRITE_PS)
  {
#ifdef _MPI
    if(iAmOutRank)
#endif
      domain.PS("Domain.ps", It_Finest, 0);
  }


  #ifdef _MPI
  // Partition the by now finest grid using Metis and distribute among processes.
  //analyse interfaces and create edge objects
  domain.GenerateEdgeInfo();
  int maxCellsPerVertex; //out parameter of mesh partitioning

  Partition_Mesh3D(comm, &domain, maxCellsPerVertex); //do the actual partitioning

  //domain has been reduced, thus generate edge info anew.
  domain.GenerateEdgeInfo();

  //largest possible number of processes which share one dof
  int maxSubDomainPerDof = MIN(maxCellsPerVertex, mpiSize);

  if(iAmOutRank)
  {
	  Output::print("MaxSubDomainPerDof: ", maxSubDomainPerDof);
  }
                domain.GetCollection(It_Finest, 0)->GetN_OwnCells(),
                ". N_HaloCells: ",
                domain.GetCollection(It_Finest, 0)->GetN_HaloCells());
  #endif

  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  }

  // choose example according to the value of TDatabase::ParamDB->EXAMPLE and construct it
  Example_CD3D example;

  // construct the cd3d problem object
#ifdef _MPI
  CD3D cd3d(domain, example, maxSubDomainPerDof);
#else
  CD3D cd3d(domain, example);
#endif

  //=========================================================================
  //Start the actual computations.
  //=========================================================================

  cd3d.assemble(); // assemble matrix and rhs

  cd3d.solve();    // solve the system

  cd3d.output();   // produce nice output

  //=========================================================================

  Output::close_file();

#ifdef _MPI
  MPI_Finalize();
#endif

  return 0;
} // end main


