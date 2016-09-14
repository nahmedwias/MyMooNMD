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
  Chrono stopwatch; // Start a stopwatch for time measurement during execution

  TDatabase Database;  // Initialize User Input Databases.
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  TFEDatabase2D FEDatabase; // Initialize FE Database.

  TDomain domain(argv[1], parmoon_db); // Initialize geometry

  // ***** Write parameters to output file ***** //
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  Database.CheckParameterConsistencyNSE();
  Database.WriteParamDB(argv[0]);

  std::list<TCollection* > gridCollections
  = domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  if(parmoon_db["output_write_ps"]) domain.PS("Domain.ps", It_Finest, 0);

  domain.print_info("Multiphase2D domain");  // Output domain info

  // ***** Declaring objects for NSE2D ***** //
//  Example_NSE2D example_nse2d(parmoon_db); // Construct Example for NSE
//  NSE2D nse2d(domain, parmoon_db, example_nse2d); // Construct NSE system
//  nse2d.assemble(); // assemble linear term
//  nse2d.stopIt(0); // check initial residual

  // ***** Initialize object for iterations residual output ***** //
  LoopInfo loop_info("nonlinear");
  loop_info.print_time_every_step = true;
  loop_info.verbosity_threshold = 1; // full verbosity
//  loop_info.print(0, nse2d.getFullResidual());

  stopwatch.print_time("setting up spaces, matrices, linear assemble");
  stopwatch.reset();
//
//  //================================================================
//  // ****** NON-LINEAR LOOP ****** //
//  for(unsigned int k = 1;; k++)
//  {
//    nse2d.solve();
//    if(parmoon_db["problem_type"].is(3)) break; // Stokes
//    nse2d.assemble_nonlinear_term();
//
//    if(nse2d.stopIt(k)) // Check residuals
//    {
//      loop_info.finish(k, nse2d.getFullResidual());
//      break;
//    }
//    else loop_info.print(k, nse2d.getFullResidual());
//  } // end for k, non linear loop
//  stopwatch.print_time("total solving duration: ");
//  //================================================================
//  nse2d.output();




  Output::print<1>("<<<<<<<<<<<<<<< NOW SOLVING CONVECTION >>>>>>>>>>>>");
  //================================================================
  CD2D cd2d(domain, parmoon_db);

  TCollection *collection = domain.GetCollection(It_Finest, 0, -4711);
  Example_CD2D ghost_example(parmoon_db);
  TFESpace2D ghost_space(collection,(char*)"space", (char*)"cd2d fe_space",
                         ghost_example.get_bc(0),
                         2, nullptr);
  BlockFEMatrix ghost_matrix({&ghost_space});
  BlockVector convection_vector(ghost_matrix,false);

  convection_vector.write("convection_vector0");

  int longueur = convection_vector.length();
  for (int indice=0; indice < longueur; indice++)
  {
    convection_vector.at(indice) = 0.5 ;
  }
  convection_vector.write("convection_vector");

  TFEFunction2D convecting_function(&ghost_space, (char*) "c", (char*)"c",
                                    convection_vector.get_entries(),
                                    convection_vector.length());

  double resultats[3];

  convecting_function.FindGradient(0.2, 0.7, resultats);
  cout << resultats[0] << endl;
  cout << resultats[1] << endl;
  cout << resultats[2] << endl;



  Output::print<1>("================== JE COMMENCE A ASSEMBLER =============");
  cd2d.assemble();
  Output::print<1>("================== JE COMMENCE A RESOUDRE =============");
  cd2d.solve();
  cd2d.output();
  //================================================================




  Output::close_file();
  return 0;
}
// end main
