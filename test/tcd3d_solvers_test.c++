/**
 * @brief A test program for the solution of Time dependent convection diffusion reaction equations
 * 
 * This test program is intended to check whether different solvers for CD3D
 * are able to solve three artificial CDR problems with analytic solution,
 * where the solution lies in the ansatz space.
 * It can be easily adapted to include more solvers, the control which solver
 * to use should come from the outside.
 * 
 * @author 
 * @history 14.05.2016
 * 
 */

#include <Time_CD3D.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <Domain.h>
#include <TimeDiscRout.h>
#include <Multigrid.h>

#ifdef _MPI
#include <mpi.h>
#include <MeshPartition.h>
double bound = 0;
double timeC = 0;
#endif

// function to compare the errors
void compare(const Time_CD3D& tcd3d, std::array<double, int(3)>errors, double tol)
{
  std::array<double, int(3)> computed_errors;
  computed_errors = tcd3d.get_errors();
  //check the L2-error
  if(fabs(computed_errors.at(0)-errors.at(0)) > tol )
  {
    ErrThrow("L2(0,t,L2): " , computed_errors.at(0), "  ", errors.at(0));
  }
  if(fabs(computed_errors.at(1) - errors.at(1)) > tol )
  {
    ErrThrow("L2(0,t,H1): " , computed_errors.at(1), "  ", errors.at(1));
  }
}

void check(ParameterDatabase& db, int ansatz_order, int time_disc,
           std::array<double, int(3)> errors, double tol)
{
  TDatabase::ParamDB->ANSATZ_ORDER = ansatz_order;
  TDatabase::ParamDB->EXAMPLE = db["example"];
  TDatabase::TimeDB->TIME_DISC = time_disc;
  TDatabase::TimeDB->STARTTIME=0;
  TDatabase::TimeDB->TIMESTEPLENGTH = 0.1;
  TDatabase::TimeDB->ENDTIME = 1;
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  
#ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
#endif
  TDomain domain(db);
  
  // set the time discretization 
  SetTimeDiscParameters(0);
  
  domain.Init(db["boundary_file"], db["geo_file"]);
  
  int n_ref_total = domain.get_n_initial_refinement_steps();
  size_t mg_levels = db["multigrid_n_levels"];
  size_t n_ref_after = mg_levels > 1 ? mg_levels - 1: 0;
  size_t n_ref_before =  n_ref_total - n_ref_after;
  
  if(n_ref_before < 0)
  {
    ErrThrow("Number of multigrid levels is greater than number of refinement "
        "levels. Garbage in, garbage out.")
  }
  for(size_t i=0; i<n_ref_before; i++)
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
  
  // example object
  Example_CD3D example_obj;
#ifdef _MPI
  Time_CD3D tcd3d(gridCollections, db, example_obj, maxSubDomainPerDof);
#else
  Time_CD3D tcd3d(gridCollections, db, example_obj);
#endif
  
  // assemble the matrices and right hand side at the start time
  tcd3d.assemble_initial_time();
  
  int step = 0;
  int imag=0;

  while(TDatabase::TimeDB->CURRENTTIME < TDatabase::TimeDB->ENDTIME-1e-10)
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
    
    tcd3d.output(step,imag);
    
    compare(tcd3d, errors, tol);
  }
}

// Choose the solver according to the input string and set global database
// entries accordingly.
void set_solver_globals(std::string solver_name, ParameterDatabase& db)
{
  if(solver_name.compare("jacobi")==0)
  {
    db["solver_type"] = "iterative";
    db["preconditioner"] = "jacobi";
    db["residual_tolerance"] = 1e-13;
    // Output::setVerbosity(2);
  }
  else if(solver_name.compare("multigrid") ==0)
  {
    db["solver_type"] = "iterative";
    db["preconditioner"] = "multigrid";
    db["refinement_n_initial_steps"] = 2;
    db["multigrid_n_levels"] = 2;
    db["multigrid_cycle_type"] = "V";
    db["multigrid_smoother"] = "jacobi";
    db["residual_tolerance"] = 1e-13;
    db["multigrid_coarse_max_n_iterations"] = 5;
    db["damping_factor"] = 0.7;
    db["damping_factor_finest_grid"] = 0.7;
    

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
  else if(solver_name.compare("mumps") == 0)
  {
    db["solver_type"] = "direct";
    TDatabase::ParamDB->SOLVER_TYPE = 2;  }
#endif
  else
  {
    throw std::runtime_error("Unknown solver for Time_CD3D problem!");
  }
}

double get_tolerance(std::string solver_name)
{//solver dependent tolerance?
  if (solver_name.compare("jacobi") == 0)
    return 1e-10;

  else if (solver_name.compare("multigrid") == 0)
    return 1e-10;

#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
    return 1e-13 ;
#else
  else if (solver_name.compare("mumps") == 0)
    return 1e-13;
#endif
  else
    throw std::runtime_error("Unknown solver for Time_CD3D problem!");
}

int main(int argc, char* argv[])
{
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
  db.merge(Solver<>::default_solver_database(), true);
  db.merge(Multigrid::default_multigrid_database());
  db["problem_type"] = 1;
  db.add("refinement_n_initial_steps",(size_t) 2,"",(size_t) 0, (size_t) 2);
  db["multigrid_n_levels"] = 1;
  // db.add("n_multigrid_levels", (size_t) 0, "",(size_t) 0, (size_t) 2);
  // db.add("solver_type", std::string("direct"), "", {"direct", "iterative"});
  // db.add("preconditioner", std::string("multigrid"), "",{"jacobi", "multigrid"});
  db["boundary_file"] = "Default_UnitCube";
  // db["output_write_vtk"] = false;
  
  TDatabase::ParamDB->PROBLEM_TYPE = 1; // CDR problem type

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
  std::array<double, int(3)> errors;
  errors = {{0.0, 0.0, 0.0}};
  
  if(my_rank==0)
    Output::print<1>("Hexahedra grid.");
  db["geo_file"] = "Default_UnitCube_Hexa";
    
  db["example"] = -4; // Example -4: linear space and time 
  check(db, 1, 1, errors, tol); // time discretization is also included in the function
  
  db["example"] = -5; // Example -4: quadratic space time
  check(db, 2, 2, errors, tol); // time discretization is also included in the function
  
  if (std::string(argv[1]).compare("jacobi") != 0)//Jacobi simply fails on this grid (mpi and sequential).
  {
    if(my_rank==0)
      Output::print<1>("Tetrahedra grid.");
    db["geo_file"] = "Default_UnitCube_Tetra";
    
    db["example"] = -4; // Example -4: linear space and time 
    check(db, 1, 1, errors, tol); // time discretization is also included in the function
    
    db["example"] = -5; // Example -4: quadratic space time
    check(db, 2, 2, errors, tol); // time discretization is also included in the function
  }
#ifdef _MPI
  MPI_Finalize();
#endif
  
  return 0;
}