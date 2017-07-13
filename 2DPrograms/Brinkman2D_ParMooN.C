// =======================================================================
//
// Purpose:     main program for solving a stationary Brinkman equation in ParMooN
//
// Author:      Laura Blank
//
// History:     Implementation started on  14.03.2016

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <Brinkman2D.h>
#include <MainUtilities.h>
#include <LocalAssembling2D.h>
#include <Example_Brinkman2D.h>
#include <Chrono.h>

#include <ParameterDatabase.h>

#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
    Output::print(" ");
    Output::print("################################################################################################################");
    Output::print("################################################################################################################");
    //for(refinement_n_initial_steps=1; refinement_n_initial_steps <= 6;++refinement_n_initial_steps)
    //{
    
    Chrono timer;
    
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain Domain(argv[1], parmoon_db);
    
  //open OUTFILE, this is where all output is written to (addionally to console)

  Output::set_outfile(parmoon_db["outfile"]);

  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);

  // refine grid
  size_t n_ref =  Domain.get_n_initial_refinement_steps();
  //-----------WITHOUT REFINEMENT (for n_initial_refinement_steps: 0)
    // write grid into a Postscript file
    if(parmoon_db["output_write_ps"])
        Domain.PS("Domain.ps", It_Finest, 0);
    Example_Brinkman2D example(parmoon_db);
    timer.restart_and_print("the setup of the domain, database, and example: )");
    
    //=========================================================================
    // create an object of the Brinkman class
    Brinkman2D brinkman2d(Domain, parmoon_db, example);
    timer.restart_and_print("constructing the Brinkman2D object: ");
    Output::print<>("Database info: ");
    parmoon_db.info();
    brinkman2d.assemble();
    timer.restart_and_print("assembling: ");
    brinkman2d.solve();
    timer.restart_and_print("solving: ");
    brinkman2d.output();
    timer.restart_and_print("creating the output: ");

    
    //=========================================================================
    Output::print("<<<<< ParMooN Finished: Brinkman3D Main Program >>>>>");
    
    timer.print_total_time("this Level (total): ");
    timer.reset();
    
    //      Output::print("geo_file: ",TDatabase::ParamDB->geo_file);
    //      Output::print("example: ",TDatabase::ParamDB->example);
    //      Output::print(": ",TDatabase::ParamDB->VELOCITY_SPACE);
    //      Output::print(": ",TDatabase::ParamDB->PRESSURE_SPACE);
    //      Output::print("PERMEABILITY: ",TDatabase::ParamDB->PERMEABILITY);
    //      Output::print("n_neumann_boundary: ",TDatabase::ParamDB->n_neumann_boundary);
    //      Output::print("n_g_v_boundary: ",TDatabase::ParamDB->n_g_v_boundary);
    //      Output::print("n_u_v_boundary: ",TDatabase::ParamDB->n_u_v_boundary);
    //      Output::print("n_nitsche_boundary: ",TDatabase::ParamDB->n_nitsche_boundary);
    //      Output::print("PkPk_stab: ",TDatabase::ParamDB->PkPk_stab);
    //      Output::print("equal_order_stab_weight_PkPk: ",TDatabase::ParamDB->equal_order_stab_weight_PkPk);
    //------------------------------------
//-----------WITH REFINEMENT
  for(int i=0; i< n_ref; i++)
  {
      Domain.RegRefineAll();
  } //end for(int i=0; i< n_ref; i++), Compute on highest level only
//      Output::print("========================================================================================================");
//      Output::print("Level: ", i+1);
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);
  Example_Brinkman2D example(parmoon_db);
        timer.restart_and_print("the setup of the domain, database, and example: )");

  //=========================================================================
  // create an object of the Brinkman class
  Brinkman2D brinkman2d(Domain, parmoon_db, example);
       timer.restart_and_print("constructing the Brinkman2D object: ");
    Output::print<>("Database info: ");
    parmoon_db.info();
  brinkman2d.assemble();
      timer.restart_and_print("assembling: ");
  brinkman2d.solve();
      timer.restart_and_print("solving: ");
  brinkman2d.output();
      timer.restart_and_print("creating the output: ");

  //=========================================================================
      Output::print("<<<<< ParMooN Finished: Brinkman3D Main Program >>>>>");
      
      timer.print_total_time("this Level (total): ");
      timer.reset();
      
//      Output::print("geo_file: ",TDatabase::ParamDB->geo_file);
//      Output::print("example: ",TDatabase::ParamDB->example);
//      Output::print(": ",TDatabase::ParamDB->VELOCITY_SPACE);
//      Output::print(": ",TDatabase::ParamDB->PRESSURE_SPACE);
//      Output::print("PERMEABILITY: ",TDatabase::ParamDB->PERMEABILITY);
//      Output::print("n_neumann_boundary: ",TDatabase::ParamDB->n_neumann_boundary);
//      Output::print("n_g_v_boundary: ",TDatabase::ParamDB->n_g_v_boundary);
//      Output::print("n_u_v_boundary: ",TDatabase::ParamDB->n_u_v_boundary);
//      Output::print("n_nitsche_boundary: ",TDatabase::ParamDB->n_nitsche_boundary);
//      Output::print("PkPk_stab: ",TDatabase::ParamDB->PkPk_stab);
//      Output::print("equal_order_stab_weight_PkPk: ",TDatabase::ParamDB->equal_order_stab_weight_PkPk);

      
      
  //} //end for(int i=0; i< n_ref; i++), Compute on all levels
        
    
  Output::close_file();
    return 0;
         // }
} // end main
