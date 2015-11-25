/**
 * @brief A test program for the solving of CD3D problems.
 *
 * This serves as a test for the solving of CD3D problems. It is intended to
 * perform CD3D calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
 * So far only one such test is implemented.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes wrong)
 * the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then this tests
 * shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 * @date 2015/20/11
 * @author Clemens Bartsch
 *
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <CD3D.h>
#include <MeshPartition.h>

#include <sys/stat.h>

#include <MainUtilities.h> // for measuring of errors

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

int main(int argc, char* argv[])
{
  /** Program 1
   *  This program tests
   */
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
  TDomain domain;

  // Set Database values (this is what is usually done by the input-file)

  TDatabase::ParamDB->PROBLEM_TYPE = 1; // CDR problem type
  TDatabase::ParamDB->EXAMPLE = 0; //simple sine laplace example

  TDatabase::ParamDB->UNIFORM_STEPS = 2; // 2 uniform refinement steps
  TDatabase::ParamDB->DRIFT_Z = 1;
  TDatabase::ParamDB->ANSATZ_ORDER=2; // ANSATZ_ORDER 1 is not working yet
  TDatabase::ParamDB->DISCTYPE = 1; //Galerkin discretization, nothing else implemented

  TDatabase::ParamDB->RE_NR = 0.5; // enters the chosen example as diffusion coefficient: 1/RE_NR

  TDatabase::ParamDB->Par_P0 = 0; // process responsible for the output
  TDatabase::ParamDB->Par_P3 = 1; // use mesh partitioning with halo cells

  TDatabase::ParamDB->SOLVER_TYPE = 1; //use iterative solver

  TDatabase::ParamDB->SC_SOLVER_SCALAR = 11; // Richardson Iteraiton (FixedPointIte)
  TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 1; // Jacobi preconditioner
  TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = 100; // max number of fixed point iterations - reached!
  TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR = 0.0; //relative error breaking criterion (not reached)
  TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = 1e-12; //absolute error breaking criterion (not reached)

  // choose unit cube as test domain and a corresponding initial mesh
  TDatabase::ParamDB->BNDFILE = "Default_UnitCube";
  TDatabase::ParamDB->GEOFILE = "Default_UnitCube_Geo";

  // Initialize geometry and initialize the mesh.
  domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);

  // Do initial regular grid refinement.
  for(int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
  {
    domain.RegRefineAll();
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
  #endif

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

  //=========================================================================

  // instead of calling CD3D::output() - compare errors against reference errors
  // TODO (when there is more than one programming running, this will be encapsulated
  // and out-sourced)

    double errors[4];
    TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
    const TFEFunction3D& function = cd3d.getFunction();
    const TFESpace3D* space = function.GetFESpace3D();

    function.GetErrors(cd3d.getExample().get_exact(0), 4, AllDerivatives,
                   2, L2H1Errors, cd3d.getExample().get_coeffs(),
                   &aux, 1, &space, errors);
    #ifdef _MPI
//    // usual code block to gather information about this
//    // processes role in the mpi communicator
//    MPI_Comm globalComm = TDatabase::ParamDB->Comm;
//    int mpiRank, mpiSize;
//    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
//    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
//    bool iAmOutRank= (mpiRank == TDatabase::ParamDB->Par_P0);
//
//    double errorsReduced[4]; //memory for global (across all processes) error
//
//    MPI_Allreduce(errors, errorsReduced, 2, MPI_DOUBLE, MPI_SUM, globalComm);
//    for(i=0;i<2;i++)
//      errors[i] = sqrt(errorsReduced[i]);
//    if(iAmOutRank){
//      OutPut(endl);
//      OutPut( "L2: " << sqrt(errorsReduced[0]) << endl);
//      OutPut( "H1-semi: " << sqrt(errorsReduced[1]) << endl);
//    }
    #else
      // check L2 error
      if( errors[0] - 0.00166608 > 1e-8 )
      {
        ErrThrow("Program 1: L2 norm not correct.");
      }
      // check H1-semi error
      if( errors[1] - 0.0444529 > 1e-7 )
      {
        ErrThrow("Program 1: H1-semi norm not correct.");
      }
    #endif

#ifdef _MPI
  MPI_Finalize();
#endif
  } //end program 1

  return 0;
} // end main
