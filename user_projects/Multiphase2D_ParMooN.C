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

  cd_db["example"]          = 4;    // 1 = TwoInteriorLayer, 4 = Multiphase2D
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



  /********************************************************************
   * INITIALIZING OBJECTS FOR MULTIPHASE
   ********************************************************************/
  TCollection   *collection = domain.GetCollection(It_Finest, 0, -4711);
  Example_CD2D  ghost_example(parmoon_db);
  TFESpace2D    ghost_space(collection,(char*)"space", (char*)"cd2d fe_space",
                         ghost_example.get_bc(0),
                         2, nullptr);
  BlockFEMatrix ghost_matrix({&ghost_space});
  BlockVector   convection_vector(ghost_matrix,false);

  convection_vector.write("convection_vector0");

  int longueur = convection_vector.length();

  for (int indice=0; indice < longueur; indice++)
  {    convection_vector.at(indice) = 1 ;  }

  convection_vector.write("convection_vector");

  TFEFunction2D phase_fraction(&ghost_space, (char*) "c", (char*)"c",
                               convection_vector.block(0),
                               convection_vector.length(0));

  /* density constant of the fluids */
  double rho1 = 1;
  double rho2 = 1000;
  double mu1 = 1e-7;
  double mu2 = 1e-3 ;

  BlockVector rho_vector = convection_vector;
  BlockVector mu_vector  = convection_vector;
  rho_vector.scale(rho1-rho2);
  mu_vector .scale(mu1-mu2);
  rho_vector.add_scaled(convection_vector,rho2);
  mu_vector .add_scaled(convection_vector,mu2);
  rho_vector.write("rho_vector");
  mu_vector .write("mu_vector");

  /** @brief Finite Element function for density field */
  TFEFunction2D rho_field(&ghost_space, (char*) "c", (char*)"c",
                          rho_vector.block(0),
                          rho_vector.length(0));
  /** @brief Finite Element function for dynamic viscosity field */
  TFEFunction2D mu_field(&ghost_space, (char*) "c", (char*)"c",
                         mu_vector.block(0),
                         mu_vector.length(0));



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
    rho_vect.add_scaled(convection_vector,rho2);
    mu_vect .add_scaled(convection_vector,mu2);
    rho_vect.write("rho_vector");
    mu_vect .write("mu_vector");

    /** @brief Finite Element function for density field */
    TFEFunction2D rho_fiel(&ghost_space, (char*) "c", (char*)"c",
                            rho_vect.block(0),
                            rho_vect.length(0));

    /** @brief Finite Element function for dynamic viscosity field */
    TFEFunction2D mu_fiel(&ghost_space, (char*) "c", (char*)"c",
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
