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
 * @return The current value of time as a float, when TDatabase::ParamDB->timeprofiling is true, 0 otherwise.
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

  //FYI
  OutPut("Hello there. I am process " << mpiRank << " of " << mpiSize << endl);
#endif

  TFEDatabase3D feDatabase;

  // Construct domain, thereby read in controls from the input file.
  TDomain domain(argv[1]);


  //set PROBLEM_TYPE to CD (corresponds to 1)
  TDatabase::ParamDB->PROBLEM_TYPE = 1;

  //TODO CB Write a method checking the parameter consistency
  //(and taking care for setting Problem_type to 1).

  //Write all Parameters to the OUTFILE for later reference
  OpenFiles();
  #ifdef _MPI
	if(iAmOutRank) //Only one process should do that.
  #endif
  Database.WriteParamDB(argv[0]);

  // Read in geometry and initialize the mesh. (See code of domain.Init
  // for usage of hard-coded example meshes)
  domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); // call mesh generator

  // Do initial grid refinement. //CB FIXME Figure out a good usage of control parameters UNIFORM_STEPS and LEVELS
  for(int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
  {
	  domain.RegRefineAll(); //Regular refinement
  }

  // Write grid into a postscript file
  #ifdef _MPI
  if(iAmOutRank && TDatabase::ParamDB->WRITE_PS) //Only one process should do that.
  {
    domain.PS("Domain.ps", It_Finest, 0);
  }
  #else
  if(TDatabase::ParamDB->WRITE_PS)
  {
    domain.PS("Domain.ps", It_Finest, 0);
  }
  #endif

  #ifdef _MPI
  // Partition the by now finest grid using Metis and distribute among processes.
  domain.GenerateEdgeInfo(); //analyse interfaces and create edge objects
  int maxCellsPerVertex;
  int t1 = MPI_Wtime(); //measure time
  Partition_Mesh3D(comm, &domain, maxCellsPerVertex); //do the actual partitioning
  int t2 = MPI_Wtime(); //measure time
  domain.GenerateEdgeInfo(); //domain has been reduced, thus generate edge info anew.

  int maxSubDomainPerDof = MIN(maxCellsPerVertex, mpiSize); //largest possible number of processes which share one dof.

  if(iAmOutRank)
  {
	  int maxSubDomainPerDof = MIN(maxCellsPerVertex, mpiSize);
	  printf("Time taken for Domain Decomposition is %e\n", (t2-t1));
	  OutPut("MaxSubDomainPerDof: " << maxSubDomainPerDof << endl);
  }

  OutPut("Process " << mpiRank << ". N_Cells: " << domain.GetCollection(It_Finest, 0)->GetN_Cells() <<
		". N_OwnCells: " << domain.GetCollection(It_Finest, 0)->GetN_OwnCells() <<
		". N_HaloCells: " << domain.GetCollection(It_Finest, 0)->GetN_HaloCells() << endl);

  // FIXME CB Loop over levels and call to Domain_Crop similiar to below is only necessary when using multigrid.
  //  for(int i=TDatabase::ParamDB->UNIFORM_STEPS;i<TDatabase::ParamDB->LEVELS;i++)
  //  {
  //    printf("************************************LEVEL %d****************************************\n",i);
  //    domain.RegRefineAll();
  //    domain.GenerateEdgeInfo();
  //    Domain_Crop(Comm, &domain);
  //  }
  //    OutPut("Process " << rank << ". N_Cells: " << domain.GetCollection(It_Finest, 0)->GetN_Cells() <<
  //    		". N_OwnCells: " << domain.GetCollection(It_Finest, 0)->GetN_OwnCells() <<
  //			". N_HaloCells: " << domain.GetCollection(It_Finest, 0)->GetN_HaloCells() << endl);
  #endif

  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

  // choose example according to the value of TDatabase::ParamDB->EXAMPLE and construct it
  Example_CD3D example;


#ifdef _MPI
  CD3D cd3d(domain, example, maxSubDomainPerDof);
#else
  CD3D cd3d(domain, example);
#endif
  //Start the actual computations.
  //=========================================================================


  cd3d.assemble();
  cd3d.solve();
  cd3d.output();
  //=========================================================================

  CloseFiles();

#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;
} // end main


