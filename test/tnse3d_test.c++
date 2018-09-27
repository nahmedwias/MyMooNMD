/**
 * @brief A test program for the solving of TNSE3D problems.
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
 * Fixed are discretization_type 1 (galerkin), LAPLACETYPE 0,
 * NSE_NONLINEAR_FORM 0 and SC_NONLIN_ITE_TYPE_SADDLE = 0.
 * Note that for these discretizations it is futile to choose any other NSTYPE
 * than 1 (see e.g. MooNMD documentation p. 35), but we vary them anyway, just
 * for the sake of testing the different types.
 *
 * So far the test is adapted to:
 *  - testing the umfpack solver when compiled SEQUENTIAL
 *  - testing lsc preconditioned fgmres SEQUENTIAL/MPI
 *  - testing multigrid preconditioned fgmres SEQUENTIAL
 *  - testing the mumps solver when compiled MPI
 *
 * The MPI Mumps test contains one example for the three combinations
 * 1 - 2 (hexa), 2 - 12 (hexa), 4 - 2 (hexa) of nstype and velocity space.
 * The tests for 3rd order elements cannot be run, because these elements are
 * not fitted for mpi yet.
 *
 * @author Najib Alia
 *
 * @date 2016/05/18
 */

#include <Time_NSE3D.h>

#include <Database.h>
#include <FEDatabase3D.h>
#include <Multigrid.h>
#include <TimeDiscRout.h>

#ifdef _MPI
#include <mpi.h>
#include <MeshPartition.h>
double bound = 0;
double timeC = 0;
#endif

void compare(const Time_NSE3D& tnse3d, std::array<double, int(4)> errors, double tol)
{
  std::array<double, int(6)> computed_errors;
  computed_errors = tnse3d.get_errors();

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

void compute(std::list<TCollection* > grid_collections, ParameterDatabase& db,
             std::array<std::array<double, int(4)>,3> errors, double tol
#ifdef _MPI
             , int maxSubDomainPerDof
#endif
             )
{
  #ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  int my_rank = 0;
#endif

  // set some parameters for time stepping
  TDatabase::TimeDB->STARTTIME=0;
  TDatabase::TimeDB->TIMESTEPLENGTH=0.05;
  TDatabase::TimeDB->ENDTIME=1;
  TDatabase::TimeDB->CURRENTTIME=  TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters(0);

  // Construct example object
  Example_TimeNSE3D example(db);
  // Construct Time_NSE3D object
#ifdef _MPI
  Time_NSE3D tnse3d(grid_collections, db, example, maxSubDomainPerDof);
#else
  Time_NSE3D tnse3d(grid_collections, db, example);
#endif
  
  tnse3d.get_time_stepping_scheme().current_step_ =0;
  tnse3d.get_time_stepping_scheme().set_time_disc_parameters();

  int image = 0;

  tnse3d.assemble_initial_time();
  //======================================================================
  // time iteration
  //======================================================================
  while(TDatabase::TimeDB->CURRENTTIME < TDatabase::TimeDB->ENDTIME-1e-10)
  {
    tnse3d.get_time_stepping_scheme().current_step_++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    // SetTimeDiscParameters(1);
    tnse3d.get_time_stepping_scheme().set_time_disc_parameters();
    double tau = tnse3d.get_time_stepping_scheme().get_step_length();
    //TODO: set also the current_time in the main class
    TDatabase::TimeDB->CURRENTTIME += tau;
    if (my_rank==0)
      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    tnse3d.assemble_rhs();
    tnse3d.assemble_nonlinear_term();
    tnse3d.assemble_system();
    for(unsigned int k=0; ; k++)
    {
      tnse3d.compute_residuals();
      // checking residuals
      if(tnse3d.stop_it(k))
        break;
      tnse3d.solve();
      tnse3d.assemble_nonlinear_term();
      tnse3d.assemble_system();
    }  // end of nonlinear loop
    cout<<" current step : " << tnse3d.get_time_stepping_scheme().current_step_<<endl;
    
    tnse3d.output(tnse3d.get_time_stepping_scheme().current_step_,image);
    cout<<" current step : " << tnse3d.get_time_stepping_scheme().current_step_<<endl;
    // check the errors
    if(tnse3d.get_time_stepping_scheme().current_step_==1)
      compare(tnse3d, errors[0],tol);
    else if(tnse3d.get_time_stepping_scheme().current_step_ ==2)
      compare(tnse3d, errors[1],tol);
    else if(tnse3d.get_time_stepping_scheme().current_step_ ==20)
      compare(tnse3d, errors[2],tol);
  } // end of time loop
}

void check(ParameterDatabase& db, std::list<TCollection* > grid_collections,
           int velocity_order, int pressure_order,
           int nstype, int laplacetype, int nonlineartype, int time_discretizationtype,
           std::array<std::array<double, int(4)>,3> errors, double tol
#ifdef _MPI
           , int maxSubDomainPerDof
#endif
          )
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  TDatabase::ParamDB->NSTYPE = nstype;
  TDatabase::ParamDB->NSE_NONLINEAR_FORM = nonlineartype;
  TDatabase::ParamDB->LAPLACETYPE = laplacetype;
  TDatabase::TimeDB->TIME_DISC = time_discretizationtype;
  
  // rough check
  if(time_discretizationtype == 1)
    db["time_discretization"] = "backward_euler";
  else if(time_discretizationtype == 2)
    db["time_discretization"] = "crank_nicolson";
  else if(time_discretizationtype == 3)
    db["time_discretization"]= "bdf_two";    
  
  db["time_start"]= 0.;
  db["time_end"]= 1.;
  db["time_step_length"]= 0.05;
  
  db["imex_scheme_"] = false;
  
  Output::print("============Time Discretization===============", db["time_discretization"]);

#ifdef _MPI
  compute(grid_collections,db,errors,tol,maxSubDomainPerDof);
#else
  compute(grid_collections,db,errors,tol);
#endif
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
    db["multigrid_smoother"] = "patch_vanka_store";
    db["multigrid_smoother_coarse"] = "nodal_vanka_store";
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

  if(solver_name.compare("lsc") == 0)
    return 1e-7;

#ifndef _MPI
  if(solver_name.compare("umfpack") == 0)
    return 1e-7;
  if(solver_name.compare("multigrid") == 0)
    return 1e-7;
#else
  if(solver_name.compare("mumps") == 0)
    return 1e-7;
#endif
    throw std::runtime_error("Unknown solver for NSE3D problem!");

  return 0;
}

void set_errors(int example, int velocity_order, int nstype,
                int timediscretizationtype, std::string solver_name,
                bool istetra,
                std::array<std::array<double, int(4)>,3>& errors)
{
  // Note that these errors remain the same between NSTypes and between SEQ AND MPI
  // If it is not the case => THERE IS A PROBLEM!
  // Errors[0] are in the first time step (t=0.05)
  // Errors[1] are in the second time step (t=0.1)
  // Errors[2] are in the last time step (t=1)

  if (example == 0) // Errors for the example Linear_space_time.h
  {
    errors[0] = {{0.0, 0.0, 0, 0}};
    errors[1] = {{0.0, 0.0, 0, 0}};
    errors[2] = {{0.0, 0.0, 0, 0}};
  }
  else if (example == 1) // Example AnsatzLinConst
  {
    errors[0] = {{0.0, 0.0, 0.0, 0.0}};
    errors[1] = {{0.0, 0.0, 0.0, 0.0}};
    errors[2] = {{0.0, 0.0, 0.0, 0.0}};
  }
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
  db.merge(Example3D::default_example_database());
  db.merge(TimeDiscretization::default_TimeDiscretization_database());
  db.merge(LocalAssembling3D::default_local_assembling_database());
  db["problem_type"].set<size_t>(6);
  db["nonlinloop_slowfactor"]=1.;
  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Hexa", "", 
         {"Default_UnitCube_Hexa","Default_UnitCube_Tetra"});
  db.add("refinement_n_initial_steps", (size_t) 1,"", (size_t) 0, (size_t) 2);
  TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 6; // flow problem type
  db["space_discretization_type"] = "galerkin"; //Galerkin discretization, nothing else implemented
  TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE = 0;
  TDatabase::ParamDB->Par_P0 = 0; // process responsible for the output
  TDatabase::ParamDB->Par_P3 = 1; // use mesh partitioning with halo cells

  set_solver_globals(std::string(argv[1]), db);
  double tol = get_tolerance(std::string(argv[1]));
  std::array<std::array<double, int(4)>,3> errors;

  //=======================================================================
  //============= PROGRAM 1 : HEXAHEDRA GRID ==============================
  //=======================================================================
  {
    // Construct domain and refine (default = once)
    TDomain domain_hex(db);
//     size_t n_ref = domain_hex.get_n_initial_refinement_steps();
//     for(size_t i=0; i< n_ref ; i++)
//     {
//       domain_hex.RegRefineAll();
//     }

#ifdef _MPI
    int maxSubDomainPerDof = 0;
#endif
    std::list<TCollection* > grid_collections
    = domain_hex.refine_and_get_hierarchy_of_collections(
        db
#ifdef _MPI
        , maxSubDomainPerDof
#endif
       );

    //=============================================================================
    // examples ... (0 to 5)
    db["example"] = 0;
    int laplacetype = 0; int nonlineartype = 0;
    //=============================================================================
    // CRANK-NICHOLSON TIME STEPPING SCHEME========================================
    int timediscretizationtype = 2;
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q2/P1-disc elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
      set_errors(db["example"], 12, 1, timediscretizationtype, 
                 std::string(argv[1]),0, errors);
      check(db, grid_collections, 12, -4711, 1, laplacetype, nonlineartype,
            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 12, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 12, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 12, -4711, 4, laplacetype, nonlineartype,
//          timediscretizationtype, errors, tol);
//
#else
      set_errors(db["example"], 12, 1, timediscretizationtype, 
                 std::string(argv[1]),0, errors);
//      check(db, domain_hex, 12, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      set_errors(db["example"], 12, 2, timediscretizationtype, 
//                 std::string(argv[1]),0, errors);
//      check(db, domain_hex, 12, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_hex, 12, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
      check(db, grid_collections, 12, -4711, 4, laplacetype, nonlineartype,
            timediscretizationtype, errors, tol,maxSubDomainPerDof);
#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q3/P2-disc elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
//      set_errors(db["example"], 12, 1, timediscretizationtype, 
//                 std::string(argv[1]),0, errors);
//      check(db, domain_hex, 13, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 13, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 13, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 13, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
#else
      // Q3/P2-disc elements are not implemented yet in MPI
#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q2/Q1 elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
//      set_errors(db["example"], 12, 1, timediscretizationtype, 
//                 std::string(argv[1]),0, errors);
//      check(db, domain_hex, 2, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 2, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 2, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 2, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
#else
      check(db, grid_collections, 2, -4711, 1, laplacetype, nonlineartype,
            timediscretizationtype, errors, tol,maxSubDomainPerDof);

      check(db, grid_collections, 2, -4711, 2, laplacetype, nonlineartype,
            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_hex, 2, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_hex, 2, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q3/Q2 elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
//      set_errors(db["example"], 12, 1, timediscretizationtype, 
//                 std::string(argv[1]),0, errors);
//      check(db, domain_hex, 3, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 3, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 3, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 3, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
#else
      // Q3/Q2 elements are not implemented yet in MPI
#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q4/Q3 elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
//      set_errors(db["example"], 12, 1, timediscretizationtype, 
//                 std::string(argv[1]),0, errors);
//      check(db, domain_hex, 4, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 4, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 4, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 4, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
#else
      // Q4/Q3 elements are not implemented yet in MPI
#endif
    //=============================================================================
    // BACKWARD-EULER TIME STEPPING SCHEME=========================================
    timediscretizationtype = 1;
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q2/P1-disc elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
      set_errors(db["example"], 12, 1, timediscretizationtype,
                 std::string(argv[1]),0, errors);
      check(db, grid_collections, 12, -4711, 1, laplacetype, nonlineartype,
            timediscretizationtype, errors, tol);

     check(db, grid_collections, 12, -4711, 2, laplacetype, nonlineartype,
           timediscretizationtype, errors, tol);

     check(db, grid_collections, 12, -4711, 3, laplacetype, nonlineartype,
           timediscretizationtype, errors, tol);

     check(db, grid_collections, 12, -4711, 4, laplacetype, nonlineartype,
           timediscretizationtype, errors, tol);
#else
      set_errors(db["example"], 12, 1, timediscretizationtype,
                 std::string(argv[1]),0, errors);
      check(db, grid_collections, 12, -4711, 1, laplacetype, nonlineartype,
            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_hex, 12, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_hex, 12, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_hex, 12, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q3/P2-disc elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
//      set_errors(db["example"], 12, 1, timediscretizationtype, 
//                 std::string(argv[1]),0, errors);
//      check(db, domain_hex, 13, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 13, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 13, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 13, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
#else
      // Q3/P2-disc not implemented yet in MPI
#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q2/Q1 elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
      set_errors(db["example"], 12, 1, timediscretizationtype,
                 std::string(argv[1]),0, errors);
     check(db, grid_collections, 2, -4711, 1, laplacetype, nonlineartype,
           timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 2, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);

     check(db, grid_collections, 2, -4711, 3, laplacetype, nonlineartype,
           timediscretizationtype, errors, tol);
//
//      check(db, domain_hex, 2, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
#else
//      check(db, domain_hex, 2, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
      check(db, grid_collections, 2, -4711, 2, laplacetype, nonlineartype,
            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_hex, 2, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_hex, 2, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q3/Q2 elements for several NSTypes");
    //=============================================================================
//#ifndef _MPI // solve with umfpack in SEQ case
//      check(db, domain_hex, 3, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//      check(db, domain_hex, 3, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//      check(db, domain_hex, 3, -4711, 3, laplacetype, nonlineartype,
//           timediscretizationtype, errors, tol);
//      check(db, domain_hex, 3, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//#else
//      // Q3/Q2 not implemented yet in MPI
//#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing Q4/Q3 elements for several NSTypes");
    //=============================================================================
//#ifndef _MPI // solve with umfpack in SEQ case
//      check(db, domain_hex, 4, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//      check(db, domain_hex, 4, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//      check(db, domain_hex, 4, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//      check(db, domain_hex, 4, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//#else
//      // Q4/Q3 not implemented yet in MPI
//#endif
  }

  //=======================================================================
  //============= PROGRAM 2 : TETRAHEDRA GRID ==============================
  //=======================================================================
  {

    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Solver<>::default_solver_database());
    db.merge(ParameterDatabase::default_nonlinit_database());
    db.merge(Example3D::default_example_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database());
    db["problem_type"].set<size_t>(6);
    db["nonlinloop_slowfactor"]=1.;
    db.add("boundary_file", "Default_UnitCube", "");
    db.add("geo_file", "Default_UnitCube_Tetra", "",
           {"Default_UnitCube_Hexa","Default_UnitCube_Tetra"});
    db.add("refinement_n_initial_steps", (size_t) 1,"", (size_t) 0, (size_t) 2);

    TDomain domain_tet(db);
//    size_t n_ref = domain_tet.get_n_initial_refinement_steps();
//     for(size_t i=0; i< n_ref ; i++)
//     {
//       domain_tet.RegRefineAll();
//     }

#ifdef _MPI
    int maxSubDomainPerDof = 0;
#endif
    std::list<TCollection* > grid_collections
    = domain_tet.refine_and_get_hierarchy_of_collections(
        db
#ifdef _MPI
        , maxSubDomainPerDof
#endif
       );
    //=============================================================================
    // examples ... (0 to 5)
    db["example"] = 0;
    int laplacetype = 0; int nonlineartype = 0;
    //=============================================================================
    // CRANK-NICHOLSON TIME STEPPING SCHEME========================================
    int timediscretizationtype = 2;
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing P2/P1 elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
      set_errors(db["example"], 2, 1, timediscretizationtype,
                 std::string(argv[1]),1,errors);
//      check(db, domain_tet, 2, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
      check(db, grid_collections, 2, -4711, 2, laplacetype, nonlineartype,
            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 2, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 2, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
#else
      set_errors(db["example"], 2, 1, timediscretizationtype,
                 std::string(argv[1]),1,errors);
      check(db, grid_collections, 2, -4711, 1, laplacetype, nonlineartype,
            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_tet, 2, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_tet, 2, -4711, 3, laplacetype, nonlineartype,
//                 timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_tet, 2, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing P3/P2 elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
//      set_errors(db["example"], 2, 1, timediscretizationtype, 
//                 std::string(argv[1]),1,errors);
//      check(db, domain_tet, 3, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 3, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 3, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 3, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
#else
      // P3/P2 not implemented yet in MPI
#endif
    //=============================================================================
    // BACKWARD EULER TIME STEPPING SCHEME=========================================
    timediscretizationtype = 1;
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing P2/P1 elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
//      set_errors(db["example"], 2, 1, timediscretizationtype, 
//                 std::string(argv[1]),1,errors);
//      check(db, domain_tet, 2, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 2, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 2, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 2, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
#else
//      set_errors(db["example"], 2, 1, timediscretizationtype, 
//                 std::string(argv[1]),1,errors);
//      check(db, domain_tet, 2, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_tet, 2, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_tet, 2, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
//
//      check(db, domain_tet, 2, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol,maxSubDomainPerDof);
#endif
    //=============================================================================
    if (my_rank == 0)
      Output::print<1>("Testing P3/P2 elements for several NSTypes");
    //=============================================================================
#ifndef _MPI // solve with umfpack in SEQ case
//      set_errors(db["example"], 2, 1, timediscretizationtype, 
//                 std::string(argv[1]),1,errors);
//      check(db, domain_tet, 3, -4711, 1, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 3, -4711, 2, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 3, -4711, 3, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
//
//      check(db, domain_tet, 3, -4711, 4, laplacetype, nonlineartype,
//            timediscretizationtype, errors, tol);
#else
      // P3/P2 not implemented yet in MPI
#endif
      }

#ifdef _MPI
  MPI_Finalize();
#endif

  return 0;
}
