// =======================================================================
// Purpose:
//
// Author:      NA
//
// History:     Start 13.09.2016
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <NSE2D.h>
#include <CD2D.h>
#include <Example_NSE2D.h>
#include <Chrono.h>
#include <LoopInfo.h>
#include <ParameterDatabase.h>



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
  result_fieldvector.scale(property_field1 - property_field2);
  result_fieldvector.add_scaled(phasefraction_vector_, property_field2);
  result_fieldvector.write(text_);
  return result_fieldvector;
}






// ***** MAIN PROGRAM ***** //
int main(int argc, char* argv[])
{
  Chrono        stopwatch;        // Start a stopwatch for time measurement during execution

  TDatabase     Database;         // Initialize User Input Databases.
  TFEDatabase2D FEDatabase;       // Initialize FE Database.

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  ParameterDatabase nse_db("Navier-Stokes Database");
  ParameterDatabase cd_db("Convection Diffusion Database");

  std::ifstream fs(argv[1]); parmoon_db.read(fs);  fs.close();
  nse_db.merge(parmoon_db,true);
  cd_db.merge(parmoon_db,true);

  cd_db["example"]          = 1;    // 1 = TwoInteriorLayer, 4 = Multiphase2D
  cd_db["problem_type"]     = 1;
  cd_db["output_basename"]  = "multiphase_convection_output";

//  nse_db["example"]         = 0;
  nse_db["problem_type"]    = 5;
  nse_db["output_basename"] = "multiphase_nse_output";

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
   * DECLARING OBJECTS FOR NSE2D AND CD2D
   ********************************************************************/
  Example_NSE2D example_nse2d(nse_db);                  // Construct Example for NSE
  NSE2D         nse2d(domain, nse_db, example_nse2d);   // Construct NSE system
  CD2D          cd2d(domain, cd_db);                    // Construct CD system

  double rho1 = nse_db["fluid_density"];  // density constant of fluid1, eg 1000
  double rho2 = 0;                        // density constant of fluid2, eg 0
  double mu1  = nse_db["fluid_dynamic_viscosity"];    // mu constant of fluid1, eg 1e-3
  double mu2  = 0;                                    // mu constant of fluid2, eg 0




  /********************************************************************
   * INITIALIZING OBJECTS FOR MULTIPHASE
   ********************************************************************/
  BlockVector phase_field = cd2d.get_solution(); // copy vector structure
  phase_field = 1;     // set phase field = 1 everywhere

  /** @brief Finite Element function for density field */
  BlockVector   rho_vector = update_fieldvector(rho1,rho2,phase_field,"rho_vector");
  TFEFunction2D rho_field  = update_fieldfunction(&cd2d.get_space(),rho_vector,(char*)"r");

  /** @brief Finite Element function for dynamic viscosity field */
  BlockVector   mu_vector = update_fieldvector(mu1, mu2, phase_field,"mu_vector" );
  TFEFunction2D mu_field  = update_fieldfunction(&cd2d.get_space(),mu_vector,(char*)"m");



  /********************************************************************
   * START ASSEMBLING NSE2D WITH GIVEN FIELDS
   ********************************************************************/
  if (nse_db["dimensional_nse"].is(true))
    nse2d.assemble_withfields(&rho_field,&mu_field); // assemble linear term
  else
    nse2d.assemble();                                // assemble linear term
  nse2d.stopIt(0);                                   // check initial residual



  /********************************************************************
   * SOME OUTPUT AND INFORMATION SET FOR THE LOOP
   ********************************************************************/
  Output::print<1>("The velocity space is ", TDatabase::ParamDB->VELOCITY_SPACE);
  Output::print<1>("The pressure space is ", TDatabase::ParamDB->PRESSURE_SPACE);
  Output::print<1>("The ansatz space is   ", TDatabase::ParamDB->ANSATZ_ORDER);
  Output::print<1>("Convection_example number ", cd_db["example"]);
  Output::print<1>("NSE_example number ", nse_db["example"]);

  LoopInfo  loop_info("nonlinear");
  loop_info.print_time_every_step = true;
  loop_info.verbosity_threshold   = 1;            // full verbosity
  //loop_info.print(0, nse2d.getFullResidual());

  stopwatch.print_total_time("setting up spaces, matrices, linear assemble");
  stopwatch.reset();



  /********************************************************************
   * NON LINEAR LOOP
   ********************************************************************/
  for(unsigned int k = 1;; k++)
  {
    nse2d.solve();


    /********************************************************************
     * SOLVING CD2D WITH NSE2D SOLUTION
     ********************************************************************/
    if (cd_db["coupling_nse_cd"].is(true))
    {
      Output::print<1>("<<<<<<<<<<<<<<<<<< NOW SOLVING CONVECTION  >>>>>>>>>>>>>");
      Output::print<1>("================== JE COMMENCE A ASSEMBLER =============");
      cd2d.assemble(); // this line is outcommented when you want to make hand tests
//      cd2d.assemble_with_convection(&nse2d.get_velocity());
      Output::print<1>("================== JE COMMENCE A RESOUDRE  =============");
      cd2d.solve();
      Output::print<1>("<<<<<<<<<<<<<<<<<< END SOLVING CONVECTION >>>>>>>>>>>>>>");



      /********************************************************************
       * UPDATING VELOCITY VECTOR WITH CD2D SOLUTION
       ********************************************************************/
      if (nse_db["coupling_cd_nse"].is(true))
      {
        BlockVector new_phase_field = cd2d.get_solution();
        // THIS LOOP HAS TO BE RECONSIDERED
//        for (int indice=0; indice < phase_fraction.length(); indice++)
//        {
//          if (phase_fraction.at(indice) >= 1) phase_fraction.at(indice)=1;
//          if (phase_fraction.at(indice) <= 0) phase_fraction.at(indice)=0;
//        }

        /** @brief Finite Element function for density field */
        BlockVector   new_rho_vector = update_fieldvector(rho1,rho2,new_phase_field,"rho_vector");
        TFEFunction2D new_rho_field  = update_fieldfunction(&cd2d.get_space(),rho_vector,(char*) "q");

        /** @brief Finite Element function for dynamic viscosity field */
        BlockVector   new_mu_vector = update_fieldvector(mu1, mu2, new_phase_field,"mu_vector" );
        TFEFunction2D new_mu_field  = update_fieldfunction(&cd2d.get_space(),mu_vector,(char*) "s");



        /********************************************************************
         * REASSEMBLE AND CALCULATE RESIDUALS FOR NSE
         ********************************************************************/
        nse2d.assemble_nonlinear_term_withfields(&new_rho_field,&new_mu_field);
      }
      else if (nse_db["dimensional_nse"].is(true))  // if 2way coupling is deactivated but 1way is active
      { nse2d.assemble_nonlinear_term_withfields(&rho_field,&mu_field); }
      else
      { nse2d.assemble_nonlinear_term(); }
    }
    else if (nse_db["dimensional_nse"].is(true))  // if 1way coupling is deactivated but dimensional is active
    { nse2d.assemble_nonlinear_term_withfields(&rho_field,&mu_field); cout << "I AM HERE, IN MAIN PROGRAM, IF CONDITION " << endl; }
    else
    { nse2d.assemble_nonlinear_term(); }


    if(nse2d.stopIt(k)) // Check residuals
    {
      loop_info.finish(k, nse2d.getFullResidual());
      break;
    }
    else loop_info.print(k, nse2d.getFullResidual());

  } // end for k, non linear loop

  stopwatch.print_total_time("total solving duration: ");

  nse2d.output();
  cd2d.output();
//  nse2d.get_solution().write("solution_velocity");

  Output::close_file();
  return 0;
}
// end main
