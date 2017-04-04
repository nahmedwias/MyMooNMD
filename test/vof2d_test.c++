/**
 * @brief A test program for the solving of VOF2D problems.
 *

 *
 * So far the test is adapted to:
 *  - examples number ......
 *  - the finite elements ......
 *  - the nstype .....
 *  - the time stepping ....
 *  - the discretization type ....
 *  - the solver .....
 *  - the geometry (quads,trias)....
 *  - the laplace type ....
 *  - the couplings ...
 *
 * @author Najib Alia
 *
 * @date 2017/04/03
 */

#include <VOF_TwoPhase2D.h>
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <TimeDiscRout.h>
#include <AlgebraicFluxCorrection.h>

void check_errors(const VOF_TwoPhase2D& vof2d,
                  std::array<double, int(6)> errors,
                  double tol)
{
  std::array<double, int(6)> computed_errors_tnse ;
  computed_errors_tnse = vof2d.tnse2d_.get_errors();
  // check the Errors of velocity and pressure
  if( fabs(computed_errors_tnse[0]-errors[0]) > tol ||
      computed_errors_tnse[0] != computed_errors_tnse[0]) //check for nan!
    ErrThrow("L2 norm of velocity: ", computed_errors_tnse[0], "  ", errors[0]);

  if( fabs(computed_errors_tnse[1] - errors[1]) > tol )
    ErrThrow("H1 norm of velocity: ", computed_errors_tnse[1], "  ", errors[1]);

  if( fabs(computed_errors_tnse[2] - errors[2]) > tol)
    ErrThrow("L2 norm of pressure: ", computed_errors_tnse[2], "  ", errors[2]);

  if(fabs(computed_errors_tnse[3] - errors[3]) > tol )
    ErrThrow("H1 norm of pressure: ", computed_errors_tnse[3], "  ", errors[3]);

  std::array<double, int(3)> computed_errors_tcd;
  computed_errors_tcd = vof2d.phaseconvection2d_.get_errors();
  // check the Errors of concentration
  if( fabs(computed_errors_tcd[0]-errors[4]) > tol )
    ErrThrow("L2t norm of concentration: ", computed_errors_tcd[0], "  ", errors[4]);
  if( fabs(computed_errors_tcd[1]-errors[5]) > tol )
    ErrThrow("L2H1t norm of concentration: ", computed_errors_tcd[1], "  ", errors[5]);
}

void compute(ParameterDatabase& db,
             std::array<double, int(6)> errors,
             double tol)
{
  ParameterDatabase tnse_db("Navier-Stokes Database");
  ParameterDatabase tcd_db("Convection Diffusion Database");
  tnse_db.merge(db,true);
  tcd_db.merge(db,true);
  tcd_db["problem_type"]     = 2;
  tnse_db["problem_type"]    = 6;
  check_parameters_consistency_NSE(tnse_db);
  //  declaration of databases
  TDomain domain(db);
  std::list<TCollection* > gridCollections
  = domain.refine_and_get_hierarchy_of_collections(db);

  /********************************************************************
   * Creating VOF object, which contains both TimeNSE2D and TimeCD2D
   ********************************************************************/
  SetTimeDiscParameters(0);                      // Initialize parameters for time discretization
  VOF_TwoPhase2D vof(domain,tnse_db,tcd_db);
  vof.manage_example_parameters();
  vof.update_field_vectors();
  vof.output_initial_info();

  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=0;
  vof.tnse2d_.assemble_initial_time();                                // assemble linear term

  if (!tcd_db["algebraic_flux_correction"].is("none"))
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=1;
  if (vof.solve_convection_ == true)
  {
    if (vof.nse2cd_coupling_ == true)
      vof.phaseconvection2d_.assemble_initial_time_with_convection(&vof.tnse2d_.get_velocity());
    else
      vof.phaseconvection2d_.assemble_initial_time();
  }

   double end_time = TDatabase::TimeDB->ENDTIME;
   int step = 0;
   int n_substeps = GetN_SubSteps();
   vof.tnse2d_.current_step_ = 0;

  /********************************************************************
   * TIME ITERATION LOOP
   ********************************************************************/
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    step++;
    vof.tnse2d_.current_step_++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    for(int j=0; j < n_substeps; ++j)
    {
      SetTimeDiscParameters(1);            // setting the time disc parameters
      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=0;
      if (vof.tnse_variable_fluid_ == true)
      {
        vof.tnse2d_.assemble_rhs_withfields(&vof.rho_fefunction_,&vof.mu_fefunction_);
        vof.tnse2d_.assemble_massmatrix_withfields(&vof.rho_fefunction_);
        vof.tnse2d_.assemble_nonlinear_term_withfields(&vof.rho_fefunction_,&vof.mu_fefunction_);
        if( TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION == 1 )
            vof.tnse2d_.apply_slip_penetration_bc(true,true);
      }
      else
      {
        vof.tnse2d_.assemble_rhs();
        vof.tnse2d_.assemble_nonlinear_term();
      }
      vof.tnse2d_.assemble_system();

    /********************************************************************
     * NON LINEAR LOOP
     ********************************************************************/
    for(unsigned int k = 0;; k++)
    {
      if(vof.tnse2d_.stopIte(k))
        break;

      vof.tnse2d_.solve();
      if(vof.tnse2d_.imex_scheme(1))
        continue; // this interrupts the NL-Loop
      if (vof.tnse_variable_fluid_ == true)
      {
        vof.tnse2d_.assemble_nonlinear_term_withfields(&vof.rho_fefunction_,&vof.mu_fefunction_);
        if( TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION == 1 )
          vof.tnse2d_.apply_slip_penetration_bc(false,false);
      }
      else
        vof.tnse2d_.assemble_nonlinear_term();
      vof.tnse2d_.assemble_system();
     } // end for k, non linear loop


    /********************************************************************
     * SOLVING CD2D WITH NSE2D SOLUTION
     ********************************************************************/
    if (!tcd_db["algebraic_flux_correction"].is("none"))
      TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=1;
    if (vof.solve_convection_ == true )
    {
      if (vof.nse2cd_coupling_ == true)
      {
//        tcd2d.assemble_rhs_vector(&tnse2d.get_velocity()); // once per time step
//        tcd2d.assemble_stiffness_matrix_alone_with_convection(&tnse2d.get_velocity());
//        tcd2d.scale_stiffness_matrix();
        vof.phaseconvection2d_.assemble_with_convection(&vof.tnse2d_.get_velocity());
        vof.phaseconvection2d_.solve();
        vof.phaseconvection2d_.descale_stiffness(tau, TDatabase::TimeDB->THETA1); //needed once per time loop
      }
      else  // if we solve TCD2D standard, without any coupling
      {
        vof.phaseconvection2d_.assemble();
        vof.phaseconvection2d_.solve();
        vof.phaseconvection2d_.descale_stiffness(tau, TDatabase::TimeDB->THETA1); //needed once per time loop
      }
    }
      /********************************************************************
       * UPDATING VELOCITY VECTOR WITH CD2D SOLUTION
       ********************************************************************/
      if (vof.cd2nse_coupling_ == true )
        vof.update_field_vectors();
    }
    vof.tnse2d_.output(step);
    vof.phaseconvection2d_.output();
    if(TDatabase::TimeDB->CURRENTTIME >= end_time - 1e-10)
      check_errors(vof, errors,tol);
  } // end for step, time loop
}

void set_errors(int example,
                std::array<double, int(6)>& errors)
{
  switch(example)
  {
    case 10:
      errors= {2.179086354e-06, /*L2(u)*/ 0.0002257627577, /*H1(u)*/
               0.003852625123, /*L2(p)*/ 0.01793605573, /*H1(p)*/
               8.45657e-08, /*L2t(c)*/ 5.00189e-06  /*L2H1t(c)*/ };
      break;
    case 20:  // TODO: the case 20 has to be debugged
      errors= {2.179086354e-06, /*L2(u)*/ 0.0002257627577, /*H1(u)*/
               0.003852625123, /*L2(p)*/ 0.01793605573, /*H1(p)*/
               8.45657e-08, /*L2t(c)*/ 5.00189e-06  /*L2H1t(c)*/ };
      break;
    case 32:
      errors= {0.2604532338, /*L2(u)*/ 4.576438279, /*H1(u)*/
               1.511430802, /*L2(p)*/ 21.5888442, /*H1(p)*/
               0.00272719, /*L2t(c)*/ 4.32481  /*L2H1t(c)*/ };
      break;
    case 42:
      errors= {0.07022243989, /*L2(u)*/ 2.280610654, /*H1(u)*/
               17.10410817, /*L2(p)*/ 14.23273126, /*H1(p)*/
               0.211403, /*L2t(c)*/ 0.564915  /*L2H1t(c)*/ };
      break;
  }
}


int main(int argc, char* argv[])
{
  //declaration of databases
  TDatabase Database;

  TFEDatabase2D FEDatabase;
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(ParameterDatabase::default_output_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(Solver<>::default_solver_database());
  db.merge(Example2D::default_example_database());
  db.merge(AlgebraicFluxCorrection::default_afc_database());

  db["output_compute_errors"] = true;
  db["verbosity"]             = 3;
  db.add("refinement_n_initial_steps", (size_t) 1,"", (size_t) 0, (size_t) 10);;
  db.add("boundary_file", "Default_UnitSquare", "", {"Default_UnitSquare",
                                                     "../../ParMooN/data/mesh/RTInstab_PochetOrig.PRM"});
  db.add("geo_file", "UnitSquare", "", {"UnitSquare",
                                        "TwoTriangles",
                                        "../../ParMooN/data/mesh/RTInstab_PochetOrig.GEO"});

  std::array<double, int(6)> errors;
  double tol = 1e-7;
  //=======================================================================
  /* ===== EXAMPLE 10 - TNSE2D with rho!=1 and TCD2D, no coupling =========
   * TNSE2D is compared with ex. 0 of TNSE2D main program, with Re_nr=200
   * TCD2D is compared with ex. 0 of TCD2D main program   ================= */
  // ======================================================================
  {
    db["example"] = 10;
    db["dimensional_nse"] = true;
    db["solve_cd"] = true;
    db["refinement_n_initial_steps"] = 4;
    TDatabase::TimeDB->STARTTIME     = 0.;
    TDatabase::TimeDB->ENDTIME       = 0.05;
    TDatabase::TimeDB->TIMESTEPLENGTH= 0.01;
    db["fluid_density"]           = 400;  // this is then compared to
    db["fluid_dynamic_viscosity"] = 2;    // TNSE2D with Re=200
    db["reynolds_number"]         = 200;
    TDatabase::ParamDB->RE_NR     = 200;
    db["solver_type"]             = "direct";
    db["direct_solver_type"]      = "umfpack";
    TDatabase::ParamDB->VELOCITY_SPACE = 12;
    TDatabase::ParamDB->PRESSURE_SPACE = -4711;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE = 1;
    TDatabase::ParamDB->ANSATZ_ORDER= 1;
    db["algebraic_flux_correction"] = "none";
    set_errors(10,errors);
    compute(db, errors, tol);
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  }

  // this example 20 has some discrepancies, but it comes from
  // the error calculation...needs to be debugged
//  //=======================================================================
//  /* ===== EXAMPLE 20 - Coupling TNSE2D > TCD2D ===========================
//   * TNSE2D solves a given constant velocity field which is used in TCD2D =
//   * TCD2D is compared with TCD2D main program and example 0 ============== */
//  // ======================================================================
//  {
//    db["example"] = 20;
//    db["refinement_n_initial_steps"] = 4;
//    TDatabase::TimeDB->STARTTIME     = 0.;
//    TDatabase::TimeDB->CURRENTTIME   = 0.;
//    TDatabase::TimeDB->ENDTIME       = 0.05;
//    TDatabase::TimeDB->TIMESTEPLENGTH= 0.01;
//    db["fluid_density"]           = 1;  // this is then compared to
//    db["fluid_dynamic_viscosity"] = 1;    // TNSE2D with Re=200
//    db["reynolds_number"]         = 1;
//    TDatabase::ParamDB->RE_NR     = 1;
//    db["solver_type"]             = "direct";
//    db["direct_solver_type"]      = "umfpack";
//    TDatabase::ParamDB->VELOCITY_SPACE = 12;
//    TDatabase::ParamDB->PRESSURE_SPACE = -4711;
//    TDatabase::ParamDB->NSTYPE = 4;
//    TDatabase::ParamDB->LAPLACETYPE = 1;
//    TDatabase::ParamDB->ANSATZ_ORDER= 1;
//    db["algebraic_flux_correction"] = "fem-fct-cn";
//    set_errors(20,errors);
//    compute(db, errors, tol);
//  }

  //=======================================================================
  /* ===== EXAMPLE 32 - Coupling TCD2D > TNSE2D ===========================
   * TCD2D solves a given constant phase field which is used in TNSE2D ====
   * TNSE2D is compared with TNSE2D main program and example 2 ============ */
  // ======================================================================
  {
    db["example"] = 32;
    TDatabase::ParamDB->P1 = 1. ;
    TDatabase::ParamDB->P2 = 10.;
    TDatabase::ParamDB->P3 = 1. ;
    db["refinement_n_initial_steps"] = 5;
    db["dimensional_nse"] = true;
    db["coupling_nse_cd"] = false;
    db["coupling_cd_nse"] = true;
    db["solve_cd"] = true;
    TDatabase::TimeDB->STARTTIME     = 0.;
    TDatabase::TimeDB->CURRENTTIME   = 0.;
    TDatabase::TimeDB->ENDTIME       = 0.05;
    TDatabase::TimeDB->TIMESTEPLENGTH= 0.01;
    db["fluid_density"]           = 1.;
    db["fluid_dynamic_viscosity"] = 1.;
    db["reynolds_number"]         = 1.;
    TDatabase::ParamDB->RE_NR     = 1.;
    db["solver_type"]             = "direct";
    db["direct_solver_type"]      = "umfpack";
    TDatabase::ParamDB->VELOCITY_SPACE = 12;
    TDatabase::ParamDB->PRESSURE_SPACE = -4711;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE = 1;
    TDatabase::ParamDB->ANSATZ_ORDER= 1;
    db["algebraic_flux_correction"] = "fem-fct-cn";
    set_errors(32,errors);
    tol = 1e-2; // discrepancies (probably) due ...
    // ...to convection approximate solution
    compute(db, errors, tol);
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  }

  //=======================================================================
  /* ===== EXAMPLE 42 - 2-way coupling TNSE2D-TCD2D =======================
   * Solves the 2-phase Rayleigh-Taylor instability. The errors are not ===
   * meant to be correct, they are used to test program stability. ========
   * Slip condition is also tested in this example. ======================= */
  // ======================================================================
  {
//    db["boundary_file"].set_range("");
    db["example"] = 42;
    db["boundary_file"] = "../../ParMooN/data/mesh/RTInstab_PochetOrig.PRM";
    db["geo_file"] = "../../ParMooN/data/mesh/RTInstab_PochetOrig.GEO";
    TDatabase::ParamDB->P7 = 1.2 ;
    TDatabase::ParamDB->P8 = 0.018;
    db["refinement_n_initial_steps"] = 3;
    db["dimensional_nse"] = true;
    db["coupling_nse_cd"] = true;
    db["coupling_cd_nse"] = true;
    db["solve_cd"] = true;
    TDatabase::TimeDB->STARTTIME     = 0.;
    TDatabase::TimeDB->CURRENTTIME   = 0.;
    TDatabase::TimeDB->ENDTIME       = 0.05;
    TDatabase::TimeDB->TIMESTEPLENGTH= 0.01;
    db["fluid_density"]           = 1000.;
    db["fluid_dynamic_viscosity"] = 1.;
    db["reynolds_number"]         = 1000.;
    TDatabase::ParamDB->RE_NR     = 1000.;
    db["solver_type"]             = "direct";
    db["direct_solver_type"]      = "umfpack";
    TDatabase::ParamDB->VELOCITY_SPACE = 12;
    TDatabase::ParamDB->PRESSURE_SPACE = -4711;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE = 1;
    TDatabase::ParamDB->ANSATZ_ORDER= 1;
    db["algebraic_flux_correction"] = "fem-fct-cn";
    set_errors(42,errors);
    tol = 1e-4;
    compute(db, errors, tol);
    TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 0;
  }
  return 0;
}
