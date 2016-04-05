/**
 * @brief A test program for the solving of CD3D problems.
 *
 * This test program is intended to check whether different solvers for CD3D
 * are able to solve three artificial CDR problems with analytic solution,
 * where the solution lies in the ansatz space.
 * It can be easily adapted to include more solvers, the control which solver
 * to use should come from the outside.
 *
 * @author Clemens Bartsch (heavily inspired by Najib's NSE2D Test program)
 *
 * @date 2016/03/31
 */
#include <CD3D.h>

#include <Database.h>
#include <FEDatabase3D.h>

#ifdef _MPI
#include <mpi.h>
#include <MeshPartition.h>
double bound = 0;
double timeC = 0;
#endif

void compare(const CD3D& cd3d, std::array<double, int(2)> errors, double tol)
{
  std::array<double, int(2)> computed_errors;
  computed_errors = cd3d.get_errors();

  // check the L2-error
  if( fabs(computed_errors.at(0)-errors.at(0)) > tol )
  {
    ErrThrow("L2 norm: ", computed_errors.at(0), "  ", errors.at(0));
  }
  // check the H1-error
  if( fabs(computed_errors.at(1) - errors.at(1)) > tol )
  {
    ErrThrow("H1-semi norm: ", computed_errors.at(1), "  ", errors.at(1));
  }
}

void check(int example, int ansatz_order, int geo,
           std::array<double, int(2)> errors, double tol)
{
  TDatabase::ParamDB->ANSATZ_ORDER = ansatz_order;
  TDatabase::ParamDB->EXAMPLE = example;
#ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  //int my_rank = 0;
#endif

  // fresh domain object
  TDomain domain;

  if(geo == 6)
    domain.Init(std::string("Default_UnitCube").c_str(),
                std::string("Default_UnitCube_Hexa").c_str());
  else if(geo == 4)
    domain.Init(std::string("Default_UnitCube").c_str(),
                std::string("Default_UnitCube_Tetra").c_str());
  else
    ErrThrow("Chose geo = 6 (hexahedra) or 4 (tets) for this test.");

  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
  {
    domain.RegRefineAll();
  }

#ifdef _MPI
  // Partition the by now finest grid using Metis and distribute among processes.

  // 1st step: Analyse interfaces and create edge objects,.
  domain.GenerateEdgeInfo();

  // 2nd step: Call the mesh partitioning.
  int maxCellsPerVertex;
  Partition_Mesh3D(MPI_COMM_WORLD, &domain, maxCellsPerVertex);

  // 3rd step: Generate edge info anew
  domain.GenerateEdgeInfo();

  // calculate largest possible number of processes which share one dof
  int maxSubDomainPerDof = MIN(maxCellsPerVertex, size);

#endif

  // Collect those Collections which will be used in multigrid.
  std::list<TCollection* > gridCollections;
  gridCollections.push_front(domain.GetCollection(It_Finest, 0));

  // Further mesh refinement and grabbing of collections.
  for(int level=1;level<TDatabase::ParamDB->LEVELS;level++)
  {
    domain.RegRefineAll();
#ifdef _MPI
    domain.GenerateEdgeInfo();  // has to be called anew after every refinement step
    Domain_Crop(MPI_COMM_WORLD, &domain); // remove unwanted cells in the halo after refinement
#endif
    // Grab collection.
    gridCollections.push_front(domain.GetCollection(It_Finest, 0));
  }

  // Choose example according to the value of
  // TDatabase::ParamDB->EXAMPLE and construct it.
  Example_CD3D example_obj;

  // Construct the cd3d problem object.
#ifdef _MPI
  CD3D cd3d(gridCollections, example_obj, maxSubDomainPerDof);
#else
  CD3D cd3d(gridCollections, example_obj);
#endif

  Output::print("Made it here.");

  cd3d.assemble();
  cd3d.solve();
  cd3d.output();

  compare(cd3d, errors, tol);
}

// Choose the solver according to the input string and set global database
// entries accordingly.
void set_solver_globals(std::string solver_name)
{

  if (solver_name.compare("jacobi") == 0)
  {
    TDatabase::ParamDB->SOLVER_TYPE = 1;
    TDatabase::ParamDB->SC_SOLVER_SCALAR = 11;
    TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 1;
    TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = 500;  //Jacobi needs a lot of iterations
    TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR = 0.0;
    TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = 1e-10;
  }
  else if (solver_name.compare("multigrid") == 0)
  {
    TDatabase::ParamDB->SOLVER_TYPE = 1;
    TDatabase::ParamDB->SC_SOLVER_SCALAR = 11;
    TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 5;
    TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = 100; //mg needs fewer iterations
    TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR = 0.0;
    TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = 1e-12; //measures L2only

    // specific multigrid parameters
    TDatabase::ParamDB->SC_MG_TYPE_SCALAR = 0; // geometric multigrid
    TDatabase::ParamDB->SC_MG_CYCLE_SCALAR = 1;
    TDatabase::ParamDB->SC_SMOOTHER_SCALAR = 3; // SSOR smoother
    TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR = 3;
    TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR = 3;
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR = 1.0;
    TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR = 3; // SSOR smoother
    TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR = 10;
    TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR = 0.1;
    TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR = 0; // no slc
    TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR = 0;  // no slc

    Output::setVerbosity(2);

  }
#ifdef _SEQ
  else if(solver_name.compare("umfpack") == 0)
  {
    TDatabase::ParamDB->SOLVER_TYPE = 2;
  }
#endif
#ifdef _MPI
  else if (solver_name.compare("mumps") == 0)
  {
    TDatabase::ParamDB->SOLVER_TYPE = 2;
  }
#endif
  else
  {
    throw std::runtime_error("Unknown solver for CD3D problem!");
  }

}

double get_tolerance(std::string solver_name)
{//solver dependent tolerance?
  if (solver_name.compare("jacobi") == 0)
    return 1e-10;

  else if (solver_name.compare("multigrid") == 0)
    return 1e-10;

#ifdef _SEQ
  else if(solver_name.compare("umfpack") == 0)
    return 1e-13 ;
#endif

#ifdef _MPI
  else if (solver_name.compare("mumps") == 0)
    return 1e-13;
#endif

  else
    throw std::runtime_error("Unknown solver for CD3D problem!");
}


int main(int argc, char* argv[])
{
  //declaration of databases
  TDatabase Database;

#ifdef _MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  TDatabase::ParamDB->Comm = comm;

  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  int my_rank = 0;
#endif

  TFEDatabase3D FEDatabase;

  TDatabase::ParamDB->PROBLEM_TYPE = 1; // CDR problem type

  TDatabase::ParamDB->UNIFORM_STEPS = 2; // 2 uniform refinement steps
  TDatabase::ParamDB->LEVELS = 1;
  TDatabase::ParamDB->DRIFT_Z = 1;
  TDatabase::ParamDB->DISCTYPE = 1; //Galerkin discretization, nothing else implemented

  TDatabase::ParamDB->Par_P0 = 0; // process responsible for the output
  TDatabase::ParamDB->Par_P3 = 1; // use mesh partitioning with halo cells

  TDatabase::ParamDB->MEASURE_ERRORS = 1;



  set_solver_globals(std::string(argv[1]));

  double tol = get_tolerance(std::string(argv[1]));

  //===========================================================
  if(my_rank==0)
    Output::print<1>("Starting computations with solver: "
        , std::string(argv[1]), ".");
  //===========================================================
  std::array<double, int(2)> errors;
  errors = {{0.0, 0.0}};
  size_t geo;
  if(my_rank==0)
    Output::print<1>("Hexahedra grid.");
  geo = 6;
  // Example -1: constant solution, order 1 elements
  check(-1, 1, geo, errors, tol);
  // Example -2: linear solution, order 1 elements
  check(-2, 1, geo, errors, tol);
  // Example -3: quadratic solution, order 2 elements
  check(-3, 2, geo, errors, tol);

  if (std::string(argv[1]).compare("jacobi") != 0)//Jacobi simply fails on this grid (mpi and sequential).
  {
    if(my_rank==0)
      Output::print<1>("Tetrahedra grid.");
    geo = 4;
    // Example -1: constant solution, order 1 elements
    check(-1, 1, geo, errors, tol);
    // Example -2: linear solution, order 1 elements
    check(-2, 1, geo, errors, tol);
    // Example -3: quadratic solution, order 2 elements
    check(-3, 2, geo, errors, tol);
  }
#ifdef _MPI
  MPI_Finalize();
#endif
  exit(0);

}