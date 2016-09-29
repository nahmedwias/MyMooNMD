// =======================================================================
//
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
TFEFunction2D update_fieldfunction(const TFESpace2D* feSpace_, BlockVector vector_)
{
  TFEFunction2D result_fieldfunction(feSpace_, (char*) "c", (char*)"c",
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

  cd_db["example"]          = 10;    // 1 = TwoInteriorLayer, 4 = Multiphase2D
  cd_db["problem_type"]     = 1;
  cd_db["output_basename"]  = "multiphase_convection_output";

  nse_db["example"]         = 0;
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

  double rho1 = 1;      // density constant of fluid1
  double rho2 = 1000;   // density constant of fluid2
  double mu1  = 1e-7;   // mu      constant of fluid1
  double mu2  = 1e-3;   // mu      constant of fluid2




  /********************************************************************
   * INITIALIZING OBJECTS FOR MULTIPHASE
   ********************************************************************/
  BlockVector phase_field = cd2d.get_solution(); // copy vector structure
  phase_field = 1;     // set phase field = 1 everywhere

  /** @brief Finite Element function for density field */
  BlockVector   rho_vector = update_fieldvector(rho1,rho2,phase_field,"rho_vector");
  TFEFunction2D rho_field  = update_fieldfunction(&cd2d.get_space(),rho_vector);

  /** @brief Finite Element function for dynamic viscosity field */
  BlockVector   mu_vector = update_fieldvector(mu1, mu2, phase_field,"mu_vector" );
  TFEFunction2D mu_field  = update_fieldfunction(&cd2d.get_space(),mu_vector);




  /********************************************************************
   * START ASSEMBLING NSE2D WITH GIVEN FIELDS
   ********************************************************************/
  nse2d.assemble_withfields(&rho_field,&mu_field); // assemble linear term
  nse2d.stopIt(0);                                 // check initial residual



  /********************************************************************
   * SOME OUTPUT AND INFORMATION SET FOR THE LOOP
   ********************************************************************/
  Output::print<1>("The velocity space is ", TDatabase::ParamDB->VELOCITY_SPACE);
  Output::print<1>("The pressure space is ", TDatabase::ParamDB->PRESSURE_SPACE);
  Output::print<1>("The ansatz space is ", TDatabase::ParamDB->ANSATZ_ORDER);
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
    Output::print<1>("<<<<<<<<<<<<<<<<<< NOW SOLVING CONVECTION  >>>>>>>>>>>>>");
    Output::print<1>("================== JE COMMENCE A ASSEMBLER =============");
    //cd2d.assemble();
    cd2d.assemble_with_convection(&nse2d.get_velocity());
    Output::print<1>("================== JE COMMENCE A RESOUDRE  =============");
    cd2d.solve();
    Output::print<1>("<<<<<<<<<<<<<<<<<< END SOLVING CONVECTION >>>>>>>>>>>>>>");


    /********************************************************************
     * UPDATING VELOCITY VECTOR WITH CD2D SOLUTION
     ********************************************************************/
    //nse2d.update_solution(cd2d.get_solution());
    BlockVector phase_fraction = cd2d.get_solution();
    int longueur               = phase_fraction.length();
    for (int indice=0; indice < longueur; indice++)
    {
      if (phase_fraction.at(indice) >= 1) phase_fraction.at(indice)=1;
      if (phase_fraction.at(indice) <= 0) phase_fraction.at(indice)=0;
    }
    BlockVector rho_vect = phase_fraction;
    BlockVector mu_vect  = phase_fraction;
    rho_vect.scale(rho1-rho2);
    mu_vect .scale(mu1-mu2);
    rho_vect.add_scaled(phase_field,rho2);
    mu_vect .add_scaled(phase_field,mu2);
    rho_vect.write("rho_vector");
    mu_vect .write("mu_vector");

    /** @brief Finite Element function for density field */
    TFEFunction2D rho_fiel(&cd2d.get_space(), (char*) "c", (char*)"c",
                            rho_vect.block(0),
                            rho_vect.length(0));

    /** @brief Finite Element function for dynamic viscosity field */
    TFEFunction2D mu_fiel(&cd2d.get_space(), (char*) "c", (char*)"c",
                           mu_vect.block(0),
                           mu_vect.length(0));

    /********************************************************************
     * REASSEMBLE AND CALCULATE RESIDUALS FOR NSE
     ********************************************************************/
    nse2d.assemble_nonlinear_term_withfields(&rho_fiel,&mu_fiel);

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
  nse2d.get_solution().write("solution_velocity");

  Output::close_file();
  return 0;
}
// end main
