/**
 * @brief A test program for the solving of NSE3D problems.
 *
 * This test program is intended to check whether different solvers for NSE3D
 * are able to solve three artificial Navier--Stokes problems with known
 * analytic solution, where the solutions lie in the ansatz spaces.
 * It can be easily adapted to include more solvers, the control which solver
 * to use should come from the outside.
 *
 * The program tests a selection of combinations of
 *  - the polynomial examples -1 to -3,
 *  - the stable finite element pairs P2/P1, P3/P2, Q2/Q1, Q2/P1_disc, Q3/Q2
 *  - the NSTYPES 1,2,3,4
 * on the default unit cube geometry.
 * We're only testing examples whose analytic solution is in the ansatz space,
 * thus we expect very small errors every time.
 * Fixed are DISCTYPE 1, LAPLACETYPE 0, NSE_NONLINEAR_FORM 0
 * and SC_NONLIN_ITE_TYPE_SADDLE = 0.
 * Note that for these discretizations it is futile to choose any other NSTYPE
 * than 1 (see e.g. MooNMD documentation p. 35), but we vary them anyway, just
 * for the sake of testing the different types.
 *
 * So far the test is adapted to:
 *  - testing the umfpack solver when compiled SEQUENTIAL
 *  - testing lsc preconditioned fgmres SEQUENTIAL
 *  - testing multigrid preconditioned fgmres SEQUENTIAL
 *  - testing the mumps solver when compiled MPI
 *
 * The MPI Mumps test contains one example for the three combinations
 * 1 - 2 (hexa), 2 - 12 (hexa), 4 - 2 (hexa) of nstype and velocity space.
 * The tests for 3rd order elements cannot be run, because these elements are
 * not fitted for mpi yet.
 *
 * @author Clemens Bartsch (heavily inspired by Najib's NSE2D Test program)
 *
 * @date 2016/04/04
 */

#include <NSE3D.h>

#include <Database.h>
#include <FEDatabase3D.h>
#include <Multigrid.h>

#ifdef _MPI
#include <mpi.h>
#include <MeshPartition.h>
double bound = 0;
double timeC = 0;
#endif

void compare(const NSE3D& nse3d, std::array<double, int(4)> errors, double tol)
{
  std::array<double, int(4)> computed_errors;
  computed_errors = nse3d.get_errors();

  // check the L2-error of the velcoity
  if( fabs(computed_errors[0]-errors[0]) > tol ||
      computed_errors[0] != computed_errors[0]) //check for nan!
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the H1-error of the velcoity
  if( fabs(computed_errors[1] - errors[1]) > tol )
  {
    ErrThrow("H1 norm of velocity: ", computed_errors[1], "  ", errors[1]);
  }
  // check the L2-error of the pressure
  if( fabs(computed_errors[2] - errors[2]) > tol)
  {
    ErrThrow("L2 norm of pressure: ", computed_errors[2], "  ", errors[2]);
  }
  // check the H1-error of the pressure
  if(fabs(computed_errors[3] - errors[3]) > tol )
  {
    ErrThrow("H1 norm of pressure: ", computed_errors[3], "  ", errors[3]);
  }
}
#ifndef _MPI
void check(ParameterDatabase& db, int example, const TDomain& domain,
           int velocity_order, int pressure_order,
           int nstype, std::array<double, int(4)> errors, double tol)
#else
void check(ParameterDatabase& db, int example, const TDomain& domain,
           int maxSubDomainPerDof, int velocity_order, int pressure_order,
           int nstype, std::array<double, int(4)> errors, double tol)
#endif
{
#ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  int my_rank = 0;
#endif

  if (my_rank ==0)
  {
    Output::print("******* Check ",example, " ",velocity_order, " ",
                  pressure_order, " ",nstype," *******");
  }

  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  TDatabase::ParamDB->NSTYPE = nstype;

  Example_NSE3D example_obj(example);

  //Perform usual checks on the parameter consistency
  NSE3D::check_parameters(); //makeshift check
  TDatabase::CheckParameterConsistencyNSE(); //old check

  // Construct the nse3d problem object.
#ifndef _MPI
  NSE3D nse3d(domain, db, example_obj);
#else
  NSE3D nse3d(domain, db, example_obj, maxSubDomainPerDof);
#endif

  nse3d.assemble_linear_terms();

  // check stopping criterion
  nse3d.stop_it(0);
  for(unsigned int k=1;; k++)
  {
    Output::print<1>("nonlinear step " , setw(3), k-1, "\t",
                     nse3d.get_residuals());
    nse3d.solve();

    // checking the first nonlinear iteration
    nse3d.assemble_non_linear_term();;
    if(nse3d.stop_it(k))
      break;
  }
  nse3d.output();
  // now compare the errors
  compare(nse3d, errors, tol);
}

// Choose the solver according to the input string and set global database
// entries accordingly.
void set_solver_globals(std::string solver_name, ParameterDatabase& db)
{
  db["solver_type"] = std::string("iterative");
  db["direct_solver_type"] = std::string("umfpack");
  db["iterative_solver_type"] = std::string("fgmres");
  db["preconditioner"] = std::string("no_preconditioner");
  db["residual_tolerance"] = 1.0e-13;
  
  if (solver_name.compare("lsc") == 0)
  {
    db["preconditioner"] = "least_squares_commutator";
    db["nonlinloop_epsilon"] = 1e-12;
    // just to not distract 'NSE3D::check_parameters'
    TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE = 20;
  }
  else if (solver_name.compare("multigrid") == 0)
  {
    db.merge(Multigrid::default_multigrid_database());
    db["preconditioner"] = "multigrid";
    db["refinement_n_initial_steps"] = 1;
    //control nonlinear loop
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_maxit"] = 5;
    // New multigrid parameters
    db["multigrid_n_levels"] = 2;
    db["multigrid_cycle_type"] = "V";
    db["multigrid_smoother"] = "batch_vanka";
    db["multigrid_smoother_coarse"] = "nodal_vanka";
    db["multigrid_correction_damp_factor"] = 1.0;
    db["multigrid_n_pre_smooth"] = 1;
    db["multigrid_n_post_smooth"] = 1;
    db["multigrid_coarse_residual"] = 1.0e-1;
    db["multigrid_coarse_max_n_iterations"] = 5;
    db["multigrid_vanka_damp_factor"]=1.0;

  }
#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
  {
    db["solver_type"] = "direct";
    db["direct_solver_type"] = "umfpack";
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_maxit"] = 5;
  }
  else if(solver_name.compare("pardiso") == 0)
  {
    db["solver_type"] = "direct";
    db["direct_solver_type"] = "pardiso";
    ErrThrow("pardiso not yet set!");
  }
#else
  else if (solver_name.compare("mumps") == 0)
  {
    db["solver_type"] = "direct";
    db["direct_solver_type"] = "mumps";
    db["nonlinloop_epsilon"] = 1e-15;
    db["nonlinloop_maxit"] = 5;
  }
#endif
  else
  {
    throw std::runtime_error("Unknown solver for NSE3D problem!");
  }

}

double get_tolerance(std::string solver_name)
{//solver dependent tolerance?

#ifndef _MPI
  if(solver_name.compare("umfpack") == 0)
    return 1e-9;
  if(solver_name.compare("lsc") == 0)
    return 1e-9;
  if(solver_name.compare("multigrid") == 0)
    return 1e-8;
#else
  if(solver_name.compare("mumps") == 0)
    return 1e-9 ;
#endif
    throw std::runtime_error("Unknown solver for NSE3D problem!");

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
  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_nonlinit_database());

  db["problem_type"].set<size_t>(5);
  
  db.add("refinement_n_initial_steps", (size_t) 1,"", (size_t) 0, (size_t) 2);

  TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 5; // flow problem type
  TDatabase::ParamDB->PROBLEM_TYPE = 5; // to be on the safe side...

  TDatabase::ParamDB->DRIFT_Z = 1;

  TDatabase::ParamDB->DISCTYPE = 1; //Galerkin discretization, nothing else implemented
  TDatabase::ParamDB->LAPLACETYPE = 0;
  TDatabase::ParamDB->NSE_NONLINEAR_FORM = 0;
  TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE = 0;

  TDatabase::ParamDB->Par_P0 = 0; // process responsible for the output
  TDatabase::ParamDB->Par_P3 = 1; // use mesh partitioning with halo cells

  set_solver_globals(std::string(argv[1]), db);

  double tol = get_tolerance(std::string(argv[1]));

  //===========================================================
  if(my_rank==0)
    Output::print<1>(">>>>> Starting computations with solver: <<<<<"
        , std::string(argv[1]), ".");
  //===========================================================
  std::array<double, int(4)> errors;
  errors = {{0.0, 0.0, 0.0, 0.0}};

  //============= Tests on hexa grid ==========================
  if(my_rank==0)
    Output::print<1>(">>>>> Hexahedra grid. <<<<<");
  //===========================================================
  {
    //do the domain thingy
    TDomain domain_hex(db);
    domain_hex.Init("Default_UnitCube", "Default_UnitCube_Hexa");

    size_t n_ref = domain_hex.get_n_initial_refinement_steps();
    for(size_t i=0; i< n_ref ; i++)
    {
      domain_hex.RegRefineAll();
    }

#ifdef _MPI
    // Partition the by now finest grid using Metis and distribute among processes.

    // 1st step: Analyse interfaces and create edge objects,.
    domain_hex.GenerateEdgeInfo();

    // 2nd step: Call the mesh partitioning.
    int maxCellsPerVertex;
    Partition_Mesh3D(MPI_COMM_WORLD, &domain_hex, maxCellsPerVertex);

    // 3rd step: Generate edge info anew
    domain_hex.GenerateEdgeInfo();

    // calculate largest possible number of processes which share one dof
    int maxSubDomainPerDof = MIN(maxCellsPerVertex, size);

#endif

    {
      if(my_rank==0)
        Output::print<1>("\n>>>>> Q2/Q1 element on hexahedral grid. <<<<<");
      size_t exmpl = -3;
      size_t nstype = 1;
#ifndef _MPI
      check(db, exmpl, domain_hex, 2, -4711, nstype, errors, tol);
#else
      check(db, exmpl, domain_hex, maxSubDomainPerDof, 2, -4711, nstype,
            errors, tol);
#endif
    }
    {
      if(my_rank==0)
        Output::print<1>("\n>>>>> Q2/P1^disc element on hexahedral grid. <<<<<");
      size_t exmpl = -3;
      size_t nstype = 2;
#ifndef _MPI
      check(db, exmpl, domain_hex, 12, -4711, nstype, errors, tol);
#else
      check(db, exmpl, domain_hex, maxSubDomainPerDof, 12, -4711, nstype,
            errors, tol);
#endif
    }
#ifndef _MPI//only for seq, 3rd order elements are not yet adapted for parallel
    {
      if(my_rank==0)
        Output::print<1>("\n>>>>> Q3/Q2 element on hexahedral grid. <<<<<");
      size_t exmpl = -4;
      size_t nstype = 3;
      check(db, exmpl, domain_hex, 3, -4711, nstype, errors, tol);
    }
#endif
  }

  //============= Tests on tetra grid =========================
  if(my_rank==0)
    Output::print<1>(">>>>> Tetrahedral grid. <<<<<");
  //===========================================================
  {
    //do the domain thingy
    TDomain domain_tet(db);
    domain_tet.Init("Default_UnitCube", "Default_UnitCube_Tetra");
    for(size_t i=0; i< domain_tet.get_n_initial_refinement_steps(); i++)
    {
      domain_tet.RegRefineAll();
    }

#ifdef _MPI
    // Partition the by now finest grid using Metis and distribute among processes.

    // 1st step: Analyse interfaces and create edge objects,.
    domain_tet.GenerateEdgeInfo();

    // 2nd step: Call the mesh partitioning.
    int maxCellsPerVertex;
    Partition_Mesh3D(MPI_COMM_WORLD, &domain_tet, maxCellsPerVertex);

    // 3rd step: Generate edge info anew
    domain_tet.GenerateEdgeInfo();

    // calculate largest possible number of processes which share one dof
    int maxSubDomainPerDof = MIN(maxCellsPerVertex, size);

#endif

    {
      size_t exmpl = -3;
      size_t nstype = 4;
      if(my_rank==0)
        Output::print<1>("\n>>>>> P2/P1 element on tetrahedral grid. <<<<<");
#ifndef _MPI
      check(db, exmpl, domain_tet, 2,-4711, nstype, errors, tol);
#else
      check(db, exmpl, domain_tet, maxSubDomainPerDof, 2,-4711, nstype, errors, 
            tol);
#endif
    }
#ifndef _MPI
    {
      //FIXME This test does currently not converge for multigrid! Investigate!
      // Update: It seems in the old version, it did not run either - the old
      // multigrid test just used an fgmres preconditioned with a direct solver
      // - switching the smoothers to actual nodal vanka in the old implementation
      // leads to no convergence!
      if(std::string(argv[1]) == std::string("multigrid"))
        return 0;
      if(my_rank==0)
        Output::print<1>("\n>>>>> P3/P2 element on tetrahedral grid. <<<<<");
      size_t exmpl = -4;
      size_t nstype = 4; //TODO 14
      check(db, exmpl, domain_tet, 3,-4711, nstype, errors, tol);
    }
#endif
  }

#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;

}
