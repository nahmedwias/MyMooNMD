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
#include <ParameterDatabase.h>
#include <LinAlg.h>
#include <Brinkman2D.h>
#include <MainUtilities.h>
#include <LocalAssembling2D.h>
#include <Example_Brinkman2D.h>
#include <Chrono.h>
#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>

 void analytic_coefficient_function(double x, double y, double * values)
    {
      if ( (x < 0.5 && y<0.5) || (x > 0.5 && y>0.5) )
        values[0] = 2.0;
      else 
        values[0] = 200.0;
    }


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

  // Start a stopwatch which measures the time spent in different parts of the program
  Chrono timer;

  // Declaration of the ParMooN Database (ParamDB) and FE2D Database (basis functions etc.), this is obligatory in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]); // the .dat file is transformed into a stream
  parmoon_db.read(fs); // all parameters identified (according to read()) in the stream(.dat-file) are saved in parmoon_db
  fs.close();

  // Set each variables' value in TDatabase using argv[1] (*.dat file)
  TDomain Domain(parmoon_db, argv[1]);

  // Open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  // Write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);

  // Refine the grid
  size_t n_ref =  Domain.get_n_initial_refinement_steps();

  Output::print("========================================================================================================");
  Output::print("Level: ", 0);



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
  // LB NEW 17.04.18 start
  // TO DO define Brinkman or Example parameter in database spacially_varying_coefficient_function default false
  //bool use_spatially_varying_coefficient_function = true;
  if (parmoon_db["coefficient_function_type"] == 1)
  {
    // create a FEFunction which will serve as the coefficient_function:
    auto coll = brinkman2d.get_pressure_space().GetCollection();
    // fe space of piecewise constant functions
    TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
        BoundConditionNoBoundCondition, 0, nullptr);
    BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());
    TFEFunction2D coefficient_function(
        &coefficient_function_FEspace, "coefficient_function", "coefficient_function",
        &coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());

   
    coefficient_function.Interpolate(analytic_coefficient_function);
    //coefficient_function_vector.print("coefficient_function");
TFEFunction2D *coefficient_function_ptr;
coefficient_function_ptr = &coefficient_function;
    brinkman2d.assemble(coefficient_function_ptr);
  }
  else 
    brinkman2d.assemble( );
  // LB NEW 17.04.18 end
/*
  // LB OLD 17.04.18 start
  brinkman2d.assemble();
  // LB OLD 17.04.18 end
  */
  timer.restart_and_print("assembling: ");
  brinkman2d.solve();
  timer.restart_and_print("solving: ");
  brinkman2d.output(0);
  timer.restart_and_print("creating the output: ");


  //=========================================================================
  Output::print("<<<<< ParMooN Finished: Brinkman2D Main Program >>>>>");

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
  for(size_t i = 0; i < n_ref; i++)
  {
    Domain.RegRefineAll();
    //  } //end for(int i=0; i< n_ref; i++), Compute on highest level only
    Output::print("========================================================================================================");
    Output::print("Level: ", i+1);

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
    // LB NEW 16.04.18 start
    // TO DO define Brinkman or Example parameter in database spacially_varying_coefficient_function default false
    if (parmoon_db["coefficient_function_type"] == 1)
    {
      // create a FEFunction which will serve as the coefficient_function:
      auto coll = brinkman2d.get_pressure_space().GetCollection();
      // fe space of piecewise constant functions
      TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
          BoundConditionNoBoundCondition, 0, nullptr);
      BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());
      TFEFunction2D coefficient_function(
          &coefficient_function_FEspace, "coefficient_function", "coefficient_function",
          &coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());

      coefficient_function.Interpolate(analytic_coefficient_function);
      //coefficient_function_vector.print("coefficient_function");

TFEFunction2D *coefficient_function_ptr;
coefficient_function_ptr = &coefficient_function;
    brinkman2d.assemble(coefficient_function_ptr);
     }
    else 
      brinkman2d.assemble( );
    // LB NEW 16.04.18 end

    /*
    // LB OLD 16.04.18 start
    brinkman2d.assemble();
    // LB OLD 16.04.18 start
     */
    timer.restart_and_print("assembling: ");
    brinkman2d.solve();
    timer.restart_and_print("solving: ");
    brinkman2d.output(i+1);
    timer.restart_and_print("creating the output: ");
    //=========================================================================
    Output::print("<<<<< ParMooN Finished: Brinkman2D Main Program >>>>>");

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

} //end for(int i = 0; i< n_ref; i++), Compute on all levels


Output::close_file();
return 0;
// }
} // end main
