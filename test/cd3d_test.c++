/**
 * @brief Test programs for the solving of CD3D problems.
 *
 * This serves as a test for the solving of CD3D problems. It is intended to
 * perform CD3D calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
 * So far only two such tests are implemented.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or some method in
 * between throws an uncaught exception) the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the program's functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then this tests
 * shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 * @note There is a known issue with multiple test programs within a single
 * main. This has to do with the static nature of the Database and FEDatabase
 * classes and is described in the forum.
 *
 *
 * @date 2015/20/11
 * @author Clemens Bartsch
 *
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <CD3D.h>
#include <Example_CD3D.h>
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
   *  This program tests the iterative solving method with the Jacobi preconditioner
   *  both for MPI run (with 4 processes) and sequential run.
   */
  {
  // Construct the ParMooN Databases.
  TDatabase Database;

#ifdef _MPI
  //Construct and initialise the default MPI communicator and store it.
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  TDatabase::ParamDB->Comm = comm;
  // Hold mpi rank and size ready.
  int mpiRank, mpiSize;
  MPI_Comm_rank(comm, &mpiRank);
  MPI_Comm_size(comm, &mpiSize);
#endif

  TFEDatabase3D feDatabase;


  // Set Database values (this is what is usually done by the input-file)
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db["problem_type"] = 1;
  db["example"] = 0; //simple sine laplace example
  db.add("refinement_n_initial_steps", (size_t) 2, "");
  db.add("solver_type", std::string("iterative"), "");
  db.add("preconditioner", std::string("jacobi"), "");
  db["boundary_file"] = "Default_UnitCube";
  db["geo_file"] = "Default_UnitCube_Hexa";

  // Construct domain.
  TDomain domain(db);

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

  // Initialize geometry and initialize the mesh.
  domain.Init(db["boundary_file"], db["geo_file"]);

  // Do initial regular grid refinement.
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
    domain.RegRefineAll();

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

  Example_CD3D example(db["example"],db);


  // extract the finest level collection from the domain
  std::list<TCollection * > cellCollections;
  cellCollections.push_front(domain.GetCollection(It_Finest, 0));

  // construct the cd3d problem object
#ifdef _MPI
  CD3D cd3d(cellCollections, db, example, maxSubDomainPerDof);
#else
  CD3D cd3d(cellCollections, db, example);
#endif

  //=========================================================================
  //Start the actual computations.
  //=========================================================================

  cd3d.assemble(); // assemble matrix and rhs

  cd3d.solve();    // solve the system

  //=========================================================================

  // instead of calling CD3D::output() - compare errors against reference errors

  double errors[5];
  TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
  MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
  const TFEFunction3D& function = cd3d.getFunction();
  const TFESpace3D* space = function.GetFESpace3D();

  function.GetErrors(cd3d.getExample().get_exact(0), 4, AllDerivatives,
                     2, L2H1Errors, cd3d.getExample().get_coeffs(),
                     &aux, 1, &space, errors);
#ifdef _MPI

  // The hard-coded errors here are for an mpirun with 4 processes
  // and are validated by comparison with Sashis program.
  // For some reason (?) iterative and mpi solve do not produce the same
  // output (also there are differences between different numbers of processes)

  double errorsReduced[4]; //memory for global (across all processes) error

  // calculate global error
  MPI_Allreduce(errors, errorsReduced, 2, MPI_DOUBLE, MPI_SUM, comm);
  for(int i=0;i<2;i++)
  {
    errors[i] = sqrt(errorsReduced[i]);
  }

  // check L2 error
  if( fabs(errors[0] -0.00410749) > 1e-8 )
  {
     ErrThrow("Program 1: L2 norm not correct.");
  }
  // check H1-semi error
  if( fabs(errors[1] - 0.0875751) > 1e-7 )
  {
    ErrThrow("Program 1: H1-semi norm not correct.");
  }

  //Do not finalize here, because there are more programs to follow.
  //MPI_Finalize();

#else
  // Hard-coded error values for a sequential run, validated
  // by comparison with Sashis program

  // check L2 error
  if( fabs(errors[0] - 0.00166608) > 1e-8 )
  {
    Output::print(errors[0]);
    ErrThrow("Program 1: L2 norm not correct.");
  }
  // check H1-semi error
  if( fabs(errors[1] - 0.0444529) > 1e-7 )
  {
    ErrThrow("Program 1: H1-semi norm not correct.");
  }

#endif

  Output::print("Program 1 finished. \n");

  } //end program 1


  /** Program 2
   *  This program tests the iterative solving method with
   *  the multgrid preconditioner both for MPI run (with
   *  4 processes) and sequential run.
   */
  {
    // Construct the ParMooN Databases.
    TDatabase Database;
    // Reset the Database (this method is unchecked and unsafe!)
    TDatabase::SetDefaultParameters();

#ifdef _MPI
    // Get the MPI communicator from the Database,
    // it was initialized during run of first program.
    // (...and got reset to MPI_COMM_WORLD by SetDefaultParameters())
    MPI_Comm comm = TDatabase::ParamDB->Comm;
    // Hold mpi rank and size ready.
    int mpiRank, mpiSize;
    MPI_Comm_rank(comm, &mpiRank);
    MPI_Comm_size(comm, &mpiSize);
#endif

    // Set Database values (this is what is usually done by the input-file)
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db["problem_type"] = 1;
    db["example"] = 0;
    db.add("refinement_n_initial_steps", (size_t) 3, " ");
    db.add("solver_type", std::string("iterative"), "");
    db.add("preconditioner", std::string("multigrid"), "");
    db.add("multigrid_n_levels", (size_t)2, "", (size_t)0, (size_t)10);
    db["boundary_file"] = "Default_UnitCube";
    db["geo_file"] = "Default_UnitCube_Hexa";
    
    // Construct domain.
    TDomain domain(db);

    TDatabase::ParamDB->DRIFT_Z = 1;
    TDatabase::ParamDB->ANSATZ_ORDER = 2;
    TDatabase::ParamDB->DISCTYPE = 1; //Galerkin discretization, nothing else implemented

    TDatabase::ParamDB->RE_NR = 0.5; // enters the chosen example as diffusion coefficient: 1/RE_NR

    TDatabase::ParamDB->Par_P0 = 0; // process responsible for the output
    TDatabase::ParamDB->Par_P3 = 1; // use mesh partitioning with halo cells

    TDatabase::ParamDB->SOLVER_TYPE = 1; //use iterative solver

    TDatabase::ParamDB->SC_SOLVER_SCALAR = 11; // Richardson Iteration (FixedPointIte)
    TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 5; // multigrid preconditioner
    TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = 10; // max number of fixed point iterations (not reached)
    TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR = 0.0; //relative error breaking criterion (not reached)
    TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = 1e-12; //absolute error breaking criterion (reached!)

    // specific multigrid parameters
    TDatabase::ParamDB->SC_MG_TYPE_SCALAR = 0; // geometric multigrid TODO so far no MDML implemented (and no distinction)
    TDatabase::ParamDB->SC_MG_CYCLE_SCALAR = 1;
    TDatabase::ParamDB->SC_SMOOTHER_SCALAR = 3; // SSOR smoother on all finer grids
    TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR = 3;
    TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR = 3;
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR = 1.0;
    TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR = 3; // SSOR smoother on coarsest grid, too
    TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR = 10;
    TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR = 0.1;
    TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR = 0; // no slc
    TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR = 0;  // no slc

    
    // Initialize geometry and initialize the mesh.
    domain.Init(db["boundary_file"], db["geo_file"]);
    
    // split the number of refinement steps - some have to be done before,
    // some after the domain partitioning
    int n_ref_total = domain.get_n_initial_refinement_steps();
    size_t mg_levels = db["multigrid_n_levels"];
    size_t n_ref_after = mg_levels > 1 ? mg_levels - 1: 0;
    size_t n_ref_before =  n_ref_total - n_ref_after;
    if(n_ref_before < 0)
    {
      ErrThrow("Number of multigrid levels is greater than number of refinement "
          "levels. Garbage in, garbage out.")
    }
    for(size_t i = 0; i < n_ref_before; i++)
      domain.RegRefineAll();
    
#ifdef _MPI
  // Partition the by now finest grid using Metis and distribute among processes.

  // 1st step: Analyse interfaces and create edge objects,.
  domain.GenerateEdgeInfo();

  // 2nd step: Call the mesh partitioning.
  int maxCellsPerVertex;
  Partition_Mesh3D(comm, &domain, maxCellsPerVertex);

  // 3rd step: Generate edge info anew
  domain.GenerateEdgeInfo();

  // calculate largest possible number of processes which share one dof
  int maxSubDomainPerDof = MIN(maxCellsPerVertex, mpiSize);

#endif

  // Collect those Collections which will be used in multigrid.
  std::list<TCollection* > gridCollections;
  gridCollections.push_front(domain.GetCollection(It_Finest, 0));

  // Further mesh refinement and grabbing of collections.
  for(size_t level=0; level < n_ref_after;level++)
  {
    domain.RegRefineAll();
#ifdef _MPI
    domain.GenerateEdgeInfo();  // has to be called anew after every refinement step
    Domain_Crop(comm, &domain); // remove unwanted cells in the halo after refinement
#endif
    // Grab collection.
    gridCollections.push_front(domain.GetCollection(It_Finest, 0));
  }

  // Choose example.
  Example_CD3D example(db["example"],db);

  // Construct the cd3d problem object.
#ifdef _MPI
  CD3D cd3d(gridCollections, db, example, maxSubDomainPerDof);
#else
  CD3D cd3d(gridCollections, db, example);
#endif

  //=========================================================================
  //Start the actual computations.
  //=========================================================================

  cd3d.assemble(); // assemble matrix and rhs

  cd3d.solve();    // solve the system

    //=========================================================================

    // instead of calling CD3D::output() - compare errors against reference errors

    double errors[5];
    TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
    const TFEFunction3D& function = cd3d.getFunction();
    const TFESpace3D* space = function.GetFESpace3D();

    function.GetErrors(cd3d.getExample().get_exact(0), 4, AllDerivatives,
                       2, L2H1Errors, cd3d.getExample().get_coeffs(),
                       &aux, 1, &space, errors);
#ifdef _MPI

    // The hard-coded errors here are for an mpirun with 4 processes
    // and are validated by comparison with Sashis program.
    // For some reason (?) iterative and mpi solve do not produce the same
    // output (also there are differences between different numbers of processes)

    double errorsReduced[4]; //memory for global (across all processes) error

    // calculate global error
    MPI_Allreduce(errors, errorsReduced, 2, MPI_DOUBLE, MPI_SUM, comm);
    for(int i=0;i<2;i++)
    {
      errors[i] = sqrt(errorsReduced[i]);
    }

    // check L2 error
    if( fabs(errors[0] - 0.000212104) > 1e-8 )
    {
      ErrThrow("Program 1: L2 norm not correct.");
    }
    // check H1-semi error
    if( fabs(errors[1] - 0.0110723) > 1e-7 )
    {
      ErrThrow("Program 1: H1-semi norm not correct.");
    }

    MPI_Finalize();

#else
    // Hard-coded error values for a sequential run, validated
    // by comparison with Sashis program

    // check L2 error
    if( fabs(errors[0] - 0.000212104) > 1e-8 )
    {
      ErrThrow("Program 1: L2 norm not correct.");
    }
    // check H1-semi error
    if( fabs(errors[1] - 0.0110723) > 1e-7 )
    {
      ErrThrow("Program 1: H1-semi norm not correct.");
    }

#endif

  } //end program 2

  return 0;
} // end main
