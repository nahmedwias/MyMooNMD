// =======================================================================
// Purpose:
//
// Author:      NA
//
// History:     Start 17.10.2016
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Time_NSE2D.h>
#include <Time_CD2D.h>
#include <Example_TimeNSE2D.h>
#include <Chrono.h>
#include <LoopInfo.h>
#include <ParameterDatabase.h>
#include <TimeDiscRout.h>
#include <FEFunctionInterpolator.h>
#include <VOF_TwoPhase2D.h>


// ***** MAIN PROGRAM ***** //
int main(int argc, char* argv[])
{
  Chrono        stopwatch;        // Start a stopwatch for time measurement during execution

  TDatabase     Database;         // Initialize User Input Databases.
  TFEDatabase2D FEDatabase;       // Initialize FE Database.

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  ParameterDatabase tnse_db("Navier-Stokes Database");
  ParameterDatabase tcd_db("Convection Diffusion Database");

  std::ifstream fs(argv[1]); parmoon_db.read(fs);  fs.close();
  tnse_db.merge(parmoon_db,true);
  tcd_db.merge(parmoon_db,true);

//  tcd_db["example"]          = 10;
  tcd_db["problem_type"]     = 2;
  tcd_db["output_basename"]  = "multiphase_tconvection_output";

//  tnse_db["example"]         = 18;
  tnse_db["problem_type"]    = 6;
  tnse_db["output_basename"] = "multiphase_tnse_output";


  TDomain domain(argv[1], parmoon_db);            // Initialize geometry

  /********************************************************************
   * WRITE PARAMETERS TO OUTFILE
   ********************************************************************/
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  Database.CheckParameterConsistencyNSE();
  Database.WriteParamDB(argv[0]);

  std::list<TCollection* > gridCollections
  = domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  if(parmoon_db["output_write_ps"]) domain.PS("Domain.ps", It_Finest, 0);

  domain.print_info("Multiphase2D domain");      // Output domain info


  /********************************************************************
   * Creating VOF object, which contains both TimeNSE2D and TimeCD2D
   ********************************************************************/
  SetTimeDiscParameters(0);                      // Initialize parameters for time discretization
  VOF_TwoPhase2D vof(domain,tnse_db,tcd_db);
  vof.manage_example_parameters();

   /* SOME OUTPUT AND INFORMATION SET */
  vof.output_initial_info();

  /* This calculates rho and mu vectors depending on example number
   * Check that the vectors are as expected using "output_vectors(..)" */
  vof.update_field_vectors();
  vof.output_vectors("vector_phi_init","vector_rho_init","vector_mu_init");



  /********************************************************************
   * START ASSEMBLING TimeNSE2D WITH GIVEN FIELDS
   ********************************************************************/
  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=0;
  if (vof.tnse_variable_fluid_ == true)
    vof.tnse2d_.assemble_initial_time_withfields(&vof.rho_fefunction_,&vof.mu_fefunction_); // assemble linear term
  else
    vof.tnse2d_.assemble_initial_time();                                // assemble linear term
  vof.tnse2d_.apply_slip_penetration_bc(true,true);
//  cout << "SUCCESS" << endl;
//  exit(0);
  if (!tcd_db["algebraic_flux_correction"].is("none"))
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=1;
  if (vof.solve_convection_ == true)
  {
    if (vof.nse2cd_coupling_ == true)
    {
      vof.phaseconvection2d_.assemble_initial_time_with_convection(&vof.tnse2d_.get_velocity());
    }
    else
      vof.phaseconvection2d_.assemble_initial_time();
  }

   double end_time = TDatabase::TimeDB->ENDTIME;
   int step = 0;
   int n_substeps = GetN_SubSteps();
   vof.tnse2d_.current_step_ = 0;

  /********************************************************************
   * SOME OUTPUT AND INFORMATION SET FOR THE LOOP
   ********************************************************************/
  LoopInfo  loop_info("nonlinear");
  loop_info.print_time_every_step = true;
  loop_info.verbosity_threshold   = 1;            // full verbosity
//  loop_info.print(0, tnse2d.getFullResidual());

  stopwatch.print_total_time("setting up spaces, matrices, linear assemble");
  stopwatch.reset();
  stopwatch.start();

  Chrono nse_nl_stopwatch;
  Chrono nse_timeit_stopwatch;






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
      nse_timeit_stopwatch.restart_and_print("preparation of NSE iterations");
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
      nse_nl_stopwatch.restart_and_print("solving and reassembling NL iter. "
                                          + std::to_string(k));
    } // end for k, non linear loop

    nse_timeit_stopwatch.restart_and_print("total NSE time iter. "
                                  +std::to_string(TDatabase::TimeDB->CURRENTTIME));


    /********************************************************************
     * SOLVING CD2D WITH NSE2D SOLUTION
     ********************************************************************/
    if (!tcd_db["algebraic_flux_correction"].is("none"))
      TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=1;
    if (vof.solve_convection_ == true )
    {
      Output::print<1>("<<<<<<<<<<<<<<<<<< NOW SOLVING CONVECTION  >>>>>>>>>>>>>");
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
      Output::print<1>("<<<<<<<<<<<<<<<<<< END SOLVING CONVECTION >>>>>>>>>>>>>>");


      /********************************************************************
       * UPDATING VELOCITY VECTOR WITH CD2D SOLUTION
       ********************************************************************/
      if (vof.cd2nse_coupling_ == true )
      {
        vof.update_field_vectors();
        vof.output_vectors("vector_phi_updated","vector_rho_updated","vector_mu_updated");
      }
    }

    stopwatch.restart_and_print("total whole iter. " +
                                std::to_string(TDatabase::TimeDB->CURRENTTIME));

    vof.tnse2d_.output(step);
    if(vof.solve_convection_ == true)
    {
      if((step-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
        vof.phaseconvection2d_.output();
    }
    //    vof.tnse2d_.get_solution().write("solution_velocity");
    }
  } // end for step, time loop

  stopwatch.print_total_time("total solving duration: ");
  Output::close_file();

  return 0;
}
// end main
