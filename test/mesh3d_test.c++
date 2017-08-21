/**
 * @brief A test program for the initialization from .mesh files
 *
 * This test program check whether the conversion from .geo
 * (original mesh format of ParMooN) to .mesh (new mesh format) works correctly.
 * The test is based on the cd3d_solvers_test (by Clemens). It consists in running
 * the cd3d_solvers_test once, writing the domain into a mesh, reading the 
 * written file and running the test again.
 *
 *
 * @author Alfonso Caiazzo
 *
 * @date 2016/08/16
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
  
  // write domain on a mesh file
  TCollection *coll = domain.GetCollection(It_Finest, 0);
  coll->writeMesh("tmp_unicube.mesh");

  // Intial refinement and grabbing of grids for multigrid.
#ifdef _MPI
  int maxSubDomainPerDof = 0;
#endif
  std::list<TCollection* > gridCollections
  = domain.refine_and_get_hierarchy_of_collections( db
  #ifdef _MPI
      , maxSubDomainPerDof
  #endif
      );

  // Choose and construct example.
  Example_CD3D example_obj(db);

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
  else if (solver_name.compare("multigrid") == 0)
  {
    db.merge(Multigrid::default_multigrid_database() ,true);

    db["solver_type"] = "iterative";
    db["iterative_solver_type"] = "richardson";
    db["preconditioner"] = "multigrid";
    db["refinement_n_initial_steps"] = 2;
    db["multigrid_n_levels"] = 2;
    db["max_n_iterations"] =  100;
    db["residual_tolerance"] = 1.0e-15;
    db["residual_reduction"] =  0.0;
    // Multigrid parameters
    db["multigrid_cycle_type"] = "W";
    db["multigrid_smoother"] = "jacobi";
    db["multigrid_smoother_coarse"] = "direct_solve";
    db["multigrid_correction_damp_factor"] = 0.8;
    db["multigrid_n_pre_smooth"] = 3;
    db["multigrid_n_post_smooth"] = 3;

    Output::setVerbosity(2);

  }
  else if(solver_name.compare("petsc") == 0)
  {
    db["solver_type"] = "petsc";
    db["max_n_iterations"] = 1000;
    db["residual_tolerance"] = 1.0e-12;
    db["residual_reduction"] =  0.0;
  }
#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
  {
    db["solver_type"] = "direct";
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
    ErrThrow("Unknown solver for CD3D problem! ", solver_name);
  }

}

double get_tolerance(std::string solver_name)
{//solver dependent tolerance?
  if (solver_name.compare("jacobi") == 0)
    return 1e-8;

  else if (solver_name.compare("multigrid") == 0)
    return 1e-12;
  else if(solver_name.compare("petsc") == 0)
    return 1e-12;
  
#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
    return 1e-13 ;
#else
  else if (solver_name.compare("mumps") == 0)
    return 1e-13;
#endif

  else
    ErrThrow("Unknown solver for CD3D problem!");
  return 0;
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
  db.merge(Example3D::default_example_database());
  db["problem_type"] = 1;
  db.add("refinement_n_initial_steps",(size_t) 2,"",(size_t) 0, (size_t) 3);
  db.add("multigrid_n_levels", (size_t) 0, "",(size_t) 0, (size_t) 3);
  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Hexa", "",
	 {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra",
	     "tmp_unicube.mesh"});
  
  TDatabase::ParamDB->DRIFT_Z = 1;
  db["space_discretization_type"] = "galerkin"; //Galerkin discretization, nothing else implemented

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
    Output::print<1>("Tetrahedra grid.");

  // part 1: solve with the previous format
  db["geo_file"] = "Default_UnitCube_Tetra";
  
  db["example"] = -1; // Example -1: constant solution, order 1 elements
  check(db, 1, errors, tol);

  // part 1: solve with the .mesh input
  // change the database to read the mesh file 
  db["geo_file"] = "tmp_unicube.mesh";
  check(db, 1, errors, tol);
  system("rm tmp_unicube.mesh");
  
#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;

}
