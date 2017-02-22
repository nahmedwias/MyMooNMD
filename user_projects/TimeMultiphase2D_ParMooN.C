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


// ***** LIST OF FUNCTIONS USED IN MAIN PROGRAMM ***** //
TFEFunction2D update_fieldfunction(const TFESpace2D* feSpace_, BlockVector& vector_, char* name_)
{
  TFEFunction2D result_fieldfunction(feSpace_, name_, name_,
                                     vector_.block(0),
                                     vector_.length(0));

  return result_fieldfunction;
}

BlockVector   update_fieldvector(double property_field1, double property_field2,
                               BlockVector phasefraction_vector_, std::string text_)
{
  BlockVector result_fieldvector = phasefraction_vector_;
  BlockVector unity              = phasefraction_vector_; unity = 1;
  result_fieldvector.scale(property_field1 - property_field2);
  result_fieldvector.add_scaled(unity, property_field2);
//  result_fieldvector.write(text_);
  return result_fieldvector;
}






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
  int tcd_example_number     = tcd_db["example"]; // just for convenience

//  tnse_db["example"]         = 18;
  tnse_db["problem_type"]    = 6;
  tnse_db["output_basename"] = "multiphase_tnse_output";
//  int tnse_example_number     = tnse_db["example"]; // just for convenience

  /********************************************************************
   * MANAGE PARAMETERS FOR BENCHMARK PROBLEMS!
   ********************************************************************/
  switch(tcd_example_number)
  {
    case 10:
      tcd_db["coupling_nse_cd"]  = false;
      tnse_db["coupling_cd_nse"] = false;
      break;
    case 20: case 21: case 22:
      tnse_db["dimensional_nse"] = true;
      tcd_db["coupling_nse_cd"]  = true;
      tnse_db["coupling_cd_nse"] = false;
      break;
    case 30: case 31: case 32:
//      tnse_db["dimensional_nse"] = true;
      tcd_db["coupling_nse_cd"]  = false;
//      tnse_db["coupling_cd_nse"] = true;
//      tcd_db["solve_cd"] = true;
      break;
    case 40: case 42:
//      tnse_db["dimensional_nse"] = true;
//      tcd_db["coupling_nse_cd"]  = true;
//      tnse_db["coupling_cd_nse"] = true;
      break;
  }


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
   * DECLARING OBJECTS FOR TimeNSE2D AND TimeCD2D
   ********************************************************************/
  SetTimeDiscParameters(0);                                   // Initialize parameters for time discretization

  Example_TimeNSE2D example_tnse2d(tnse_db);                  // Construct Example for NSE
  Time_NSE2D        tnse2d(domain, tnse_db, example_tnse2d);  // Construct NSE system
  Example_TimeCD2D  example_tcd2d(tcd_db);                    // Construct Example for CD
  Time_CD2D         tcd2d(domain, tcd_db, example_tcd2d);     // Construct CD system

  double rho1 = tnse_db["fluid_density"];   // density constant of fluid1, eg 1000
  double rho2 = 0;                          // density constant of fluid2, eg 0
  double mu1  = tnse_db["fluid_dynamic_viscosity"];   // mu constant of fluid1, eg 1e-3
  double mu2  = 0;                                    // mu constant of fluid2, eg 0


  /********************************************************************
   * SOME OUTPUT AND INFORMATION SET
   ********************************************************************/
  Output::print<1>("The velocity space is ", TDatabase::ParamDB->VELOCITY_SPACE);
  Output::print<1>("The pressure space is ", TDatabase::ParamDB->PRESSURE_SPACE);
  Output::print<1>("The ansatz   space is ", TDatabase::ParamDB->ANSATZ_ORDER);
  Output::print<1>("Convection_example number     ", tcd_db["example"]);
  Output::print<1>("TimeNSE_example number        ", tnse_db["example"]);



  /********************************************************************
   * INITIALIZING OBJECTS FOR MULTIPHASE
   ********************************************************************/
  BlockVector phase_field = tcd2d.get_solution(); // copy vector structure
  BlockVector rho_vector  = update_fieldvector(rho1,rho2,phase_field,"rho_vector");
  BlockVector mu_vector   = update_fieldvector(mu1, mu2, phase_field,"mu_vector" );

  // Set phase field = 1 everywhere when we dont use cd>nse coupling.
  // If we use it, then do nothing, this will take the initial solution of tcd2d.
  if (tnse_db["coupling_cd_nse"].is(false))
  {
    phase_field = 1;
    rho_vector = rho1;
    mu_vector  = mu1;
  }
  else  // for the case rho = constant, and only viscosity depends on TCD2D
  {
    mu_vector = mu1; // uncomment in case mu must stay constant
//    rho_vector = rho1;  // uncomment in case rho must stay constant
  }

  /** @brief Finite Element function for density and viscosity field */
  TFEFunction2D rho_field = update_fieldfunction(&tcd2d.get_space(),rho_vector,(char*)"r");
  TFEFunction2D mu_field  = update_fieldfunction(&tcd2d.get_space(),mu_vector,(char*)"m");


  /********************************************************************
   * START ASSEMBLING TimeNSE2D WITH GIVEN FIELDS
   ********************************************************************/
  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=0;
  if (tnse_db["dimensional_nse"].is(true))
    tnse2d.assemble_initial_time_withfields(&rho_field,&mu_field); // assemble linear term
  else
    tnse2d.assemble_initial_time();                                // assemble linear term

  if (!tcd_db["algebraic_flux_correction"].is("none"))
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=1;
  if (tcd_db["solve_cd"].is(true))
  {
    if (tcd_db["coupling_nse_cd"].is(true))
    {
      tcd2d.assemble_initial_time_with_convection(&tnse2d.get_velocity());
    }
    else
      tcd2d.assemble_initial_time();
  }

   double end_time = TDatabase::TimeDB->ENDTIME;
   int step = 0;
   int n_substeps = GetN_SubSteps();
   tnse2d.current_step_ = 0;

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
    tnse2d.current_step_++;

    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    for(int j=0; j < n_substeps; ++j)
    {
      SetTimeDiscParameters(1);            // setting the time disc parameters
      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=0;
      if (tnse_db["dimensional_nse"].is(true))
      {
        tnse2d.assemble_rhs_withfields(&rho_field,&mu_field);
        tnse2d.assemble_nonlinear_term_withfields(&rho_field,&mu_field);
        tnse2d.assemble_massmatrix_withfields(&rho_field);
      }
      else
      {
        tnse2d.assemble_rhs();
        tnse2d.assemble_nonlinear_term();
      }
      tnse2d.assemble_system();

    /********************************************************************
     * NON LINEAR LOOP
     ********************************************************************/
      nse_timeit_stopwatch.restart_and_print("preparation of NSE iterations");
    for(unsigned int k = 0;; k++)
    {
      if(tnse2d.stopIte(k))
        break;

      tnse2d.solve();

      if(tnse2d.imex_scheme(1))
        continue; // this interrupts the NL-Loop

      if (tnse_db["dimensional_nse"].is(true))
        tnse2d.assemble_nonlinear_term_withfields(&rho_field,&mu_field);
      else
        tnse2d.assemble_nonlinear_term();

      tnse2d.assemble_system();
      nse_nl_stopwatch.restart_and_print("solving+reassembling NL iter. "
                                          + std::to_string(k));
    } // end for k, non linear loop

    nse_timeit_stopwatch.restart_and_print("total NSE time iter. "
                                  +std::to_string(TDatabase::TimeDB->CURRENTTIME));

    /********************************************************************
     * SOLVING CD2D WITH NSE2D SOLUTION
     ********************************************************************/
    if (!tcd_db["algebraic_flux_correction"].is("none"))
      TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=1;
    if (tcd_db["solve_cd"].is(true))
    {
      Output::print<1>("<<<<<<<<<<<<<<<<<< NOW SOLVING CONVECTION  >>>>>>>>>>>>>");
      if (tcd_db["coupling_nse_cd"].is(true))
      {
//        tcd2d.assemble_rhs_vector(&tnse2d.get_velocity()); // once per time step
//        tcd2d.assemble_stiffness_matrix_alone_with_convection(&tnse2d.get_velocity());
//        tcd2d.scale_stiffness_matrix();
        tcd2d.assemble_with_convection(&tnse2d.get_velocity());
        tcd2d.solve();
        tcd2d.descale_stiffness(tau, TDatabase::TimeDB->THETA1); //needed once per time loop
      }
      else  // if we solve TCD2D standard, without any coupling
      {
        tcd2d.assemble();
        tcd2d.solve();
        tcd2d.descale_stiffness(tau, TDatabase::TimeDB->THETA1); //needed once per time loop
      }
      Output::print<1>("<<<<<<<<<<<<<<<<<< END SOLVING CONVECTION >>>>>>>>>>>>>>");

      /********************************************************************
       * UPDATING VELOCITY VECTOR WITH CD2D SOLUTION
       ********************************************************************/
      if (tnse_db["coupling_cd_nse"].is(true))
      {
        BlockVector new_phase_field = tcd2d.get_solution();

        /** @brief Finite Element function for density field */
        BlockVector   new_rho_vector = update_fieldvector(rho1,rho2,new_phase_field,"rho_vector");
//        new_rho_vector=1; // for the case where rho = constant and only viscosity depends on TCD2D
        TFEFunction2D new_rho_field  = update_fieldfunction(&tcd2d.get_space(),new_rho_vector,(char*) "q");
        rho_field = new_rho_field;

        /** @brief Finite Element function for dynamic viscosity field */
        BlockVector   new_mu_vector = update_fieldvector(mu1, mu2, new_phase_field,"mu_vector" );
        new_mu_vector = 1; // for the case mu=constant and only density depends on TCD2D
        TFEFunction2D new_mu_field  = update_fieldfunction(&tcd2d.get_space(),new_mu_vector,(char*) "s");
        mu_field = new_mu_field;
      }
    }

    stopwatch.restart_and_print("total whole iter. " +
                                std::to_string(TDatabase::TimeDB->CURRENTTIME));

    tnse2d.output(step);
    if(tcd_db["solve_cd"].is(true))
    {
      if((step-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
        tcd2d.output();
    }
    //    tnse2d.get_solution().write("solution_velocity");
    }
  } // end for step, time loop

  stopwatch.print_total_time("total solving duration: ");
  Output::close_file();
  return 0;
}
// end main
