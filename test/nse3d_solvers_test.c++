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
 * So far the test is only adapted for testing the umfpack solver in sequential.
 *
 * @todo TODO Fix the example -4 (something wrong with error computation) and let the
 * order 3 elements deal with that one instead.
 * @todo TODO Enable this test for: mumps (mpi only), pardiso,
 * fgmres with multigrid, fgmres with lsc preconditioner
 *
 * @author Clemens Bartsch (heavily inspired by Najib's NSE2D Test program)
 *
 * @date 2016/04/04
 */

#include <NSE3D.h>

#include <Database.h>
#include <FEDatabase3D.h>

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
void check(int example, const TDomain& domain,
           int velocity_order, int pressure_order,
           int nstype,
           std::array<double, int(4)> errors, double tol)
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



  TDatabase::ParamDB->EXAMPLE = example;
  Example_NSE3D example_obj;

  //Perform usual checks on the parameter consistency
  NSE3D::check_parameters(); //makeshift check
  TDatabase::CheckParameterConsistencyNSE(); //old check

  // Construct the nse3d problem object.
#ifdef _MPI
  NSE3D nse3d(domain, example_obj, maxSubDomainPerDof);
#else
  NSE3D nse3d(domain, example_obj);
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
void set_solver_globals(std::string solver_name)
{

  if (solver_name.compare("lsc") == 0)
  {
    ErrThrow("LSC not yet set!");
  }
  else if (solver_name.compare("multigrid") == 0)
  {
    ErrThrow("Multigrid not yet set!");
  }
#ifdef _SEQ
  else if(solver_name.compare("umfpack") == 0)
  {
    TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE = 1e-10;
    TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 5;
    TDatabase::ParamDB->SOLVER_TYPE = 2;
  }
  else if(solver_name.compare("pardiso") == 0)
  {
    ErrThrow("pardiso not yet set!");
  }
#endif
#ifdef _MPI
  else if (solver_name.compare("mumps") == 0)
  {
    ErrThrow("mumps not yet set!");
  }
#endif
  else
  {
    throw std::runtime_error("Unknown solver for NSE3D problem!");
  }

}

double get_tolerance(std::string solver_name)
{//solver dependent tolerance?

  #ifdef _SEQ
  if(solver_name.compare("umfpack") == 0)
    return 1e-9 ;
#endif

  else
    throw std::runtime_error("Unknown solver for CD3D problem!");

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

  TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 5; // flow problem type
  TDatabase::ParamDB->PROBLEM_TYPE = 5; // to be on the safe side...

  TDatabase::ParamDB->UNIFORM_STEPS = 1; // 1 uniform refinement step
  TDatabase::ParamDB->LEVELS = 1;
  TDatabase::ParamDB->DRIFT_Z = 1;

  TDatabase::ParamDB->DISCTYPE = 1; //Galerkin discretization, nothing else implemented
  TDatabase::ParamDB->LAPLACETYPE = 0;
  TDatabase::ParamDB->NSE_NONLINEAR_FORM = 0;
  TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE = 0;

  TDatabase::ParamDB->Par_P0 = 0; // process responsible for the output
  TDatabase::ParamDB->Par_P3 = 1; // use mesh partitioning with halo cells

  TDatabase::ParamDB->MEASURE_ERRORS = 1;

  set_solver_globals(std::string(argv[1]));

  double tol = get_tolerance(std::string(argv[1]));

  TDatabase::ParamDB->BNDFILE = "Default_UnitCube";
  TDatabase::ParamDB->GEOFILE = "not_specified_globally";


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
    TDomain domain_hex;
    domain_hex.Init(TDatabase::ParamDB->BNDFILE,
                    "Default_UnitCube_Hexa");
    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    {
      domain_hex.RegRefineAll();
    }

//#ifdef _MPI
//  // Partition the by now finest grid using Metis and distribute among processes.
//
//  // 1st step: Analyse interfaces and create edge objects,.
//  domain_hex.GenerateEdgeInfo();
//
//  // 2nd step: Call the mesh partitioning.
//  int maxCellsPerVertex;
//  Partition_Mesh3D(MPI_COMM_WORLD, &domain_hex, maxCellsPerVertex);
//
//  // 3rd step: Generate edge info anew
//  domain_hex.GenerateEdgeInfo();
//
//  // calculate largest possible number of processes which share one dof
//  int maxSubDomainPerDof = MIN(maxCellsPerVertex, size);
//
//#endif

  {
    if(my_rank==0)
      Output::print<1>("\n>>>>> Q2/Q1 element on hexahedral grid. <<<<<");
    size_t exmpl = -3;
    size_t nstype = 1;
    check(exmpl, domain_hex, 2, -4711, nstype, errors, tol);
  }
  {
    if(my_rank==0)
      Output::print<1>("\n>>>>> Q2/P1^disc element on hexahedral grid. <<<<<");
    size_t exmpl = -3;
    size_t nstype = 2;
    check(exmpl, domain_hex, 12, -4711, nstype, errors, tol);
  }
  {
    if(my_rank==0)
      Output::print<1>("\n>>>>> Q3/Q2 element on hexahedral grid. <<<<<");
    size_t exmpl = -3; //TODO -4
    size_t nstype = 3;
    check(exmpl, domain_hex, 3, -4711, nstype, errors, tol);
  }
  }

  //============= Tests on tetra grid =========================
  if(my_rank==0)
    Output::print<1>(">>>>> Tetrahedral grid. <<<<<");
  //===========================================================
  {
    //do the domain thingy
    TDomain domain_tet;
    domain_tet.Init(TDatabase::ParamDB->BNDFILE,
                    "Default_UnitCube_Tetra");
    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    {
      domain_tet.RegRefineAll();
    }

    //#ifdef _MPI
    //  // Partition the by now finest grid using Metis and distribute among processes.
    //
    //  // 1st step: Analyse interfaces and create edge objects,.
    //  domain_tet.GenerateEdgeInfo();
    //
    //  // 2nd step: Call the mesh partitioning.
    //  int maxCellsPerVertex;
    //  Partition_Mesh3D(MPI_COMM_WORLD, &domain_tet, maxCellsPerVertex);
    //
    //  // 3rd step: Generate edge info anew
    //  domain_tet.GenerateEdgeInfo();
    //
    //  // calculate largest possible number of processes which share one dof
    //  int maxSubDomainPerDof = MIN(maxCellsPerVertex, size);
    //
    //#endif

  {
    size_t exmpl = -3;
    size_t nstype = 4;
    if(my_rank==0)
      Output::print<1>("\n>>>>> P2/P1 element on tetrahedral grid. <<<<<");

    check(exmpl, domain_tet, 2,-4711, nstype, errors, tol);
  }
  {
    if(my_rank==0)
      Output::print<1>("\n>>>>> P3/P2 element on tetrahedral grid. <<<<<");
    size_t exmpl = -3; //TODO -4
    size_t geo = 4;
    size_t nstype = 4; //TODO 14
    check(exmpl, domain_tet, 3,-4711, nstype, errors, tol);
  }
  }

#ifdef _MPI
  MPI_Finalize();
#endif
  exit(0);

}
