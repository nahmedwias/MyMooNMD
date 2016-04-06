/**
 * @brief A test program for the solving of NSE3D problems.
 *
 * This test program is intended to check whether different solvers for NSE3D
 * are able to solve three artificial Navier--Stokes problems with known
 * analytic solution, where the solutions lie in the ansatz spaces.
 * It can be easily adapted to include more solvers, the control which solver
 * to use should come from the outside.
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
  if( fabs(computed_errors[0]-errors[0]) > tol )
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the H1-error of the velcoity
  if( fabs(computed_errors[1] - errors[1]) > tol )
  {
    ErrThrow("H1 norm of velocity: ", computed_errors[1], "  ", errors[1]);
  }
//  // check the L2-error of the pressure
//  if( fabs(computed_errors[2] - errors[2]) > tol)
//  {
//    ErrThrow("L2 norm of pressure: ", computed_errors[2], "  ", errors[2]);
//  }
  // check the H1-error of the pressure
  if(fabs(computed_errors[3] - errors[3]) > tol )
  {
    ErrThrow("H1 norm of pressure: ", computed_errors[3], "  ", errors[3]);
  }
}
void check(int example, int geo,
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
    Output::print();
    Output::print("******* Check ",example, " ",geo, " ",velocity_order, " ",
                  pressure_order, " ",nstype," *******");
  }

  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  TDatabase::ParamDB->NSTYPE = nstype;
  // form of the Laplacian and the non-linear term set fixed
  TDatabase::ParamDB->LAPLACETYPE = 0;
  TDatabase::ParamDB->NSE_NONLINEAR_FORM = 0;

  TDatabase::ParamDB->EXAMPLE = example;
  Example_NSE3D example_obj;

  /* *****************Start domain creation. ************************ */
  // I'd love to put this in some method, but that's a bit out of reach for
  // the domain is not safely movable so far.
  // FIMXE Maybe it would be better if NSE3D was equipped with a constructor
  // taking a Domain - it does not have an mpi multigrid...

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
  /* *****************End domain creation. ************************ */


  // Construct the nse3d problem object.
#ifdef _MPI
  NSE3D nse3d(gridCollections, example_obj, maxSubDomainPerDof);
#else
  NSE3D nse3d(gridCollections, example_obj);
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
    TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 50;
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
    return 1e-8 ;
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
    Output::print<1>(">>>>> Starting computations with solver: <<<<<"
        , std::string(argv[1]), ".");
  //===========================================================
  std::array<double, int(4)> errors;
  errors = {{0.0, 0.0, 0.0, 0.0}};

  //============= Tests on hexa grid ==========================
  if(my_rank==0)
    Output::print<1>(">>>>> Hexahedra grid. <<<<<");
  //===========================================================
  if(my_rank==0)
    Output::print<1>(">>>>> Q2/Q1 element on hexahedral grid. <<<<<");
  for(int exmpl = -1; exmpl > -4 ; --exmpl)
  {//four test examples
    size_t geo = 6;
    for(int nstype = 1; nstype < 5; ++nstype)
    {// all nstypes but 14
      check(exmpl, geo, 2,-4711, nstype, errors, tol);
    }
  }
  if(my_rank==0)
    Output::print<1>(">>>>> Q2/P1^disc element on hexahedral grid. <<<<<");
  for(int exmpl = -1; exmpl > -4 ; --exmpl)
  {//four test examples
    size_t geo = 6;
    for(int nstype = 1; nstype < 5; ++nstype)
    {// all nstypes but 14
      check(exmpl, geo, 12,-4711, nstype, errors, tol);
    }
  }

  //============= Tests on tetra grid =========================
  if(my_rank==0)
    Output::print<1>(">>>>> Tetrahedral grid. <<<<<");
  //===========================================================
  if(my_rank==0)
    Output::print<1>(">>>>> P2/P1 element on tetrahedral grid. <<<<<");
  for(int exmpl = -1; exmpl > -4 ; --exmpl)
  {//four test examples
    size_t geo = 4;
    for(int nstype = 1; nstype < 5; ++nstype)
    {// all nstypes but 14
      check(exmpl, geo, 2,-4711, nstype, errors, tol);
    }
  }
  if(my_rank==0)
    Output::print<1>(">>>>> P3/P2 element on tetrahedral grid. <<<<<");
  for(int exmpl = -1; exmpl > - 4 ; --exmpl)
  {//four test examples
    size_t geo = 4;
    for(int nstype = 1; nstype < 5; ++nstype)
    {// all nstypes but 14
      check(exmpl, geo, 3,-4711, nstype, errors, tol);
    }
  }

#ifdef _MPI
  MPI_Finalize();
#endif
  exit(0);

}
