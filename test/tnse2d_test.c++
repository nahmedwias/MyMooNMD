/**
 * @brief A test program to test Time-dependent Navier--Stokes
 * 
 * Residuals are checked at different time steps
 *
 * @history Reworked 22.10.2018 by Najib -
 * added solvers as arguments + activated lsc test
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Time_NSE2D.h>
#include <TimeDiscRout.h>
#include <MainUtilities.h>
#include <FEFunction2D.h>
#include <Multigrid.h> // newly implemented by clemens

void compare(Time_NSE2D& tnse2d, std::array<double, int(4)> errors, double tol)
{
  std::array<double, int(6)> computed_errors;
  computed_errors = tnse2d.get_errors();

  // check the L2-error of the velocity
  if( fabs(computed_errors[0]-errors[0]) > tol )
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the H1-error of the velocity
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

void compute(TDomain& domain, ParameterDatabase& db,
             std::array<std::array<double, int(4)>,4> errors, double tol,
             int time_disc)
{
  TDatabase::TimeDB->TIME_DISC=time_disc;
  TDatabase::TimeDB->CURRENTTIME = db["time_start"];

  Time_NSE2D tnse2d(domain, db);

  tnse2d.get_time_stepping_scheme().current_step_ = 0;
  tnse2d.get_time_stepping_scheme().set_time_disc_parameters();

  tnse2d.assemble_initial_time();
  double end_time = db["time_end"];
  while(TDatabase::TimeDB->CURRENTTIME < end_time-1e-10)
  {
    tnse2d.get_time_stepping_scheme().current_step_++;
    TDatabase::TimeDB->INTERNAL_STARTTIME
       = TDatabase::TimeDB->CURRENTTIME;

    tnse2d.get_time_stepping_scheme().set_time_disc_parameters();
    double tau = db["time_step_length"];
    TDatabase::TimeDB->CURRENTTIME += tau;

    Output::print<1>("\nCURRENT TIME: ",
                     TDatabase::TimeDB->CURRENTTIME);

    // prepare the right hand side vector
    tnse2d.assemble_matrices_rhs(0);
    // nonlinear iteration
    for(unsigned int i=0;; i++)
    {
      if(tnse2d.stopIte(i))
        break;
      tnse2d.solve();

      tnse2d.assemble_matrices_rhs(i+1);
    }
    // post processing: error computations
    Output::print("TIME DISC FOR ME: ", time_disc);
    tnse2d.output(tnse2d.get_time_stepping_scheme().current_step_);
    // check the errors
    if(tnse2d.get_time_stepping_scheme().current_step_==1)
      compare(tnse2d, errors[0], tol);
    else if(tnse2d.get_time_stepping_scheme().current_step_ ==2)
      compare(tnse2d, errors[1], tol);
    else if(tnse2d.get_time_stepping_scheme().current_step_ ==20)
      compare(tnse2d, errors[2], tol);
  }
}

void check(TDomain& domain, ParameterDatabase& db, int velocity_order, int pressure_order,
           int nstype, int laplace_type, int nonlinear_form, int time_disc,
           std::array<std::array<double, int(4)>,4> errors, double tol)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->NSTYPE = nstype;
  TDatabase::ParamDB->LAPLACETYPE = laplace_type;
  db["nse_nonlinear_form"] = "convective";
  if(nonlinear_form != 0)
    ErrThrow("other nonlinear forms are not yet tested");

  // rough check
  if(time_disc == 1)
    db["time_discretization"] = "backward_euler";
  else if(time_disc == 2)
    db["time_discretization"] = "crank_nicolson";
  else if(time_disc == 3){
    db["time_discretization"]= "bdf_two";
  }

  db["time_start"]= 0.;
  db["time_end"]= 1.;
  db["time_step_length"]= 0.05;
//  db["current_time_step_length"]= 0.05;  //not defined yet in the default TimeDiscdb

  db["imex_scheme_"] = false;
  // not defined yet in the default TimeDiscdb
//  db["extrapolate_velocity"] = false;
//  db.add("extrapolation_type", "no_extrapolation" , " ",
//         {"no_extrapolation", "constant_extrapolate", "linear_extrapolate", "quadratic_extrapolate"});

  Output::print("============Time Discretization===============", db["time_discretization"]);
  compute(domain,db,errors,tol,time_disc);
}

void set_solver_globals(std::string solver_name, ParameterDatabase& db)
{
  db["solver_type"] = std::string("iterative");
  db["direct_solver_type"] = std::string("umfpack");
  db["iterative_solver_type"] = std::string("fgmres");
  db["preconditioner"] = std::string("no_preconditioner");
  db["residual_tolerance"] = 1.0e-12;

  if (solver_name.compare("lsc") == 0)
  {
    db["preconditioner"] = "least_squares_commutator";
    db["nonlinloop_epsilon"] = 1e-10;
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
    db["multigrid_smoother"] = "nodal_vanka";
    db["multigrid_smoother_coarse"] = "nodal_vanka";
    db["multigrid_correction_damp_factor"] = 0.8;
    db["multigrid_n_pre_smooth"] = 2;
    db["multigrid_n_post_smooth"] = 2;
    db["multigrid_coarse_residual"] = 1.0e-1;
    db["multigrid_coarse_max_n_iterations"] = 5;
    db["multigrid_vanka_damp_factor"]=0.7;

  }
  else if(solver_name.compare("umfpack") == 0)
  {
    db["solver_type"] = "direct";
    db["direct_solver_type"] = "umfpack";
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_maxit"] = 2;
  }
  else
  {
    throw std::runtime_error("Unknown solver for Time_NSE2D problem!");
  }
}

double get_tolerance(std::string solver_name)
{
  double val;
  if(solver_name.compare("umfpack") == 0)
    val = 1e-9;
  if(solver_name.compare("lsc") == 0)
    val = 1e-5;
  if(solver_name.compare("multigrid") == 0)
    val = 1e-6;
  return val;
}

int main(int argc, char* argv[])
{
  
  // test for quads
  {
    TDatabase Database;
    TFEDatabase2D FEDatabase;

    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Solver<>::default_solver_database());
    db.merge(ParameterDatabase::default_nonlinit_database());
    db.merge(Example2D::default_example_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database());
    db.merge(LocalAssembling<2>::default_local_assembling_database());
    db.merge(ParameterDatabase::default_output_database());
    db["problem_type"] = 5;
    db["example"] = 0;
    db["reynolds_number"] = 1;

    db["output_compute_errors"] = true;

    db.add("refinement_n_initial_steps", (size_t) 1,"");
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});

    db["space_discretization_type"] = "galerkin";
    set_solver_globals(std::string(argv[1]), db);

    TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 5;
    db["time_end"]=1;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.05;

    //  declaration of databases
    TDomain domain(db);
    SetTimeDiscParameters(0);
    
    size_t n_ref = domain.get_n_initial_refinement_steps();
    for(unsigned int i=0; i< n_ref; ++i)
      domain.RegRefineAll();


    // test here
    std::array<std::array<double, int(4)>, 4> errors;
    double tol = get_tolerance(std::string(argv[1]));
    // errors[0], errors[1] are at first two time steps
    // and errors[2] at end time for comparison
    //domain, velocity_order, pressure_order,  nstype=1,2,3,4, 
    // laplace_type=0, nonlinear_form, time_disc, 
    //errors)
    //=============================================================================
    // CRANK-NICOLSON TIME STEPPING SCHEME
    int time_disc=2; int laplace_type; int nl_form=0;
    //=============================================================================
    Output::print<1>("Testing the Q2/P1-disc elements");
    //=============================================================================
    Output::print<1>("LAPLACETYPE: ", 0, " NSTYPE's: ",1,", ",2,", ",3,", ",4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.001112092818, 0.01455279725, 0.02027240218, 0.09770369527}};
    errors[1] = {{0.002261228092, 0.02900169754, 0.02290308963, 0.1548609258}};
    errors[2] = {{0.01917368068, 0.2444306449, 0.1019443925, 1.213868823}};
    laplace_type = 0;
    check(domain, db, 12, -4711, 1, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 2, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    
    Output::print<1>("LAPLACETYPE: ", 1, " NSTYPE's: ", 3, ", " , 4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.001100290828, 0.01461603719, 0.02009574353, 0.09725227236}};
    errors[1] = {{0.002224014719, 0.02909117943, 0.0220984703, 0.1527635601}};
    errors[2] = {{0.01881760343, 0.245035241, 0.08781391363, 1.185368926}};
    laplace_type = 1;
    check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    //=============================================================================
    Output::print<1>("Testing the Q3/P2-disc elements");
    //=============================================================================
    Output::print<1>("LAPLACETYPE: ", 0, " NSTYPE's: ",1,", ",2,", ",3,", ",4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{9.921682115e-05, 0.001899852996, 0.01927905366, 0.07644264853}};
    errors[1] = {{0.0001999247874, 0.003792810533, 0.01926073977, 0.07440290132}};
    errors[2] = {{0.001688306366, 0.03196972361, 0.01508981038, 0.2468378611}};
    laplace_type=0;
    check(domain, db, 13, -4711, 1, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 13, -4711, 2, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 13, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 13, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    
    Output::print<1>("LAPLACETYPE: ", 1, " NSTYPE's: ", 3, ", " , 4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.0001001227494, 0.001903725617, 0.01927900725, 0.07644570017}};
    errors[1] = {{0.0002016432211, 0.003804972955, 0.01926093964, 0.07442742209}};
    errors[2] = {{0.001701216631, 0.03207217596, 0.01515796884, 0.2477427189}};
    laplace_type=1; 
    check(domain, db, 13, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 13, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    //=============================================================================
    // more elements 
    //=============================================================================
    //=============================================================================
    // BACKWARD-EULER TIME STEPPING SCHEME
    time_disc = 1;
    //=============================================================================
    //=============================================================================
    Output::print<1>("Testing the Q2/P1-disc elements");
    //=============================================================================
    Output::print<1>("LAPLACETYPE: ", 0, " NSTYPE's: ",1,", ",2,", ",3,", ",4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.001123315065, 0.01452798821, 0.008151783113, 0.07430027476}};
    errors[1] = {{0.002257854029, 0.02900402091, 0.01381823763,0.1459569354}};
    errors[2] = {{0.0191888938, 0.2444304958, 0.1000448675, 1.21617188}};
    laplace_type = 0;
    check(domain, db, 12, -4711, 1, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 2, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    
    Output::print<1>("LAPLACETYPE: ", 1, " NSTYPE's: ", 3, ", " , 4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.001107694207, 0.01458218406, 0.006922680518, 0.07227851039}};
    errors[1] = {{0.002221829805, 0.02909724362, 0.01167668865, 0.1422431078}};
    errors[2] = {{0.01883002947, 0.2450131733, 0.0861689347, 1.187591355}};
    laplace_type = 1;
    check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    
    //=============================================================================
    // BACKWARD-EULER TIME STEPPING SCHEME
    time_disc = 3;
    //=============================================================================
    Output::print<1>("LAPLACETYPE: ", 0, " NSTYPE's: ",1,", ",2,", ",3,", ",4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.001123315065, 0.01452798821, 0.008151783113, 0.07430027476}};
    errors[1] = {{0.00225621577, 0.02900533136, 0.01401455865,0.1462434268}};
    errors[2] = {{0.01917397891, 0.2444306295, 0.1023554738, 1.219801748}};
    laplace_type = 0;
    check(domain, db, 12, -4711, 1, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 2, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
  }
  
  Output::print<1>("TEST SUCCESFULL: ");
  return 0;
}
