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
#include <Solver.h>
#include <Multigrid.h>

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

void check(ParameterDatabase& db, int ansatz_order,
           std::array<double, int(2)> errors, double tol)
{
  TDatabase::ParamDB->ANSATZ_ORDER = ansatz_order;
#ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  //int my_rank = 0;
#endif

  // fresh domain object
  TDomain domain(db);

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
  for(size_t level= 0; level < n_ref_after; level++)
  {
    domain.RegRefineAll();
#ifdef _MPI
    domain.GenerateEdgeInfo();  // has to be called anew after every refinement step
    Domain_Crop(MPI_COMM_WORLD, &domain); // remove unwanted cells in the halo after refinement
#endif
    // Grab collection.
    gridCollections.push_front(domain.GetCollection(It_Finest, 0));
  }

  // Choose and construct example.
  Example_CD3D example_obj(db["example"]);

  // Construct the cd3d problem object.
#ifdef _MPI
  CD3D cd3d(gridCollections, db, example_obj, maxSubDomainPerDof);
#else
  CD3D cd3d(gridCollections, db, example_obj);
#endif

  cd3d.assemble();
  cd3d.solve();
  cd3d.output();

  compare(cd3d, errors, tol);
}

// Choose the solver according to the input string and set global database
// entries accordingly.
void set_solver_globals(std::string solver_name, ParameterDatabase& db)
{

  if (solver_name.compare("jacobi") == 0)
  {
    db["solver_type"] = "iterative";
    db["iterative_solver_type"] = "richardson";
    db["preconditioner"] = "jacobi";
    db["residual_tolerance"] = 1.0e-10;
    db["residual_reduction"] =  0.0;
    db["max_n_iterations"] =  1000;
    db["min_n_iterations"] =  5;
  }
  else if (solver_name.compare("multigrid") == 0) //multigrid test will fail in mpi
  {
    db.merge(Multigrid::default_multigrid_database() ,true);

    db["solver_type"] = "iterative";
    db["iterative_solver_type"] = "richardson";
    db["preconditioner"] = "multigrid";
    db["refinement_n_initial_steps"] = 2;
    db["multigrid_n_levels"] = 2;
    db["max_n_iterations"] =  100;
    db["residual_tolerance"] = 1.0e-14;
    db["residual_reduction"] =  0.0;
    // Multigrid parameters
    db["multigrid_n_levels"] = 2;
    db["multigrid_cycle_type"] = "V";
    db["multigrid_smoother"] = "jacobi";
    db["multigrid_smoother_coarse"] = "direct_solve";
    db["multigrid_correction_damp_factor"] = 0.8;
    db["multigrid_n_pre_smooth"] = 2;
    db["multigrid_n_post_smooth"] = 2;

    Output::setVerbosity(2);

  }
#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
  {
    db["solver_type"] = "direct";
    TDatabase::ParamDB->SOLVER_TYPE = 2;
  }
#endif
#ifdef _MPI
  else if (solver_name.compare("mumps") == 0)
  {
    db["solver_type"] = "direct";
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
    return 1e-8;

  else if (solver_name.compare("multigrid") == 0)
    return 1e-12;

#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
    return 1e-13 ;
#else
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
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(Solver<>::default_solver_database() ,true);
  db["problem_type"] = 1;
  db.add("refinement_n_initial_steps",(size_t) 2,"",(size_t) 0, (size_t) 2);
  db.add("multigrid_n_levels", (size_t) 0, "",(size_t) 0, (size_t) 2);
  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Hexa", "",
	 {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
  
  TDatabase::ParamDB->DRIFT_Z = 1;
  TDatabase::ParamDB->DISCTYPE = 1; //Galerkin discretization, nothing else implemented

  TDatabase::ParamDB->Par_P0 = 0; // process responsible for the output
  TDatabase::ParamDB->Par_P3 = 1; // use mesh partitioning with halo cells

  set_solver_globals(std::string(argv[1]), db);

  double tol = get_tolerance(std::string(argv[1]));

  //===========================================================
  if(my_rank==0)
    Output::print<1>("Starting computations with solver: ",
                     std::string(argv[1]), ".");
  //===========================================================
  std::array<double, int(2)> errors;
  errors = {{0.0, 0.0}};
  
  if(my_rank==0)
    Output::print<1>("Hexahedra grid.");
  db["geo_file"] = "Default_UnitCube_Hexa";
  
  db["example"] = -1; // Example -1: constant solution, order 1 elements
  check(db, 1, errors, tol);
  
  db["example"] = -2; // Example -2: linear solution, order 1 elements
  check(db, 1, errors, tol);
  
  db["example"] = -3; // Example -3: quadratic solution, order 2 elements
  check(db, 2, errors, tol);

    if(my_rank==0)
      Output::print<1>("Tetrahedra grid.");
    db["geo_file"] = "Default_UnitCube_Tetra";
    
    db["example"] = -1; // Example -1: constant solution, order 1 elements
    check(db, 1, errors, tol);
    
    db["example"] = -2; // Example -2: linear solution, order 1 elements
    check(db, 1, errors, tol);
    
    db["example"] = -3; // Example -3: quadratic solution, order 2 elements
    check(db, 2, errors, tol);

#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;

}
