// =======================================================================
//
// Purpose:     main program for solving a stationary scalar equation using ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Darcy.h>
#include <Example_Darcy2D.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  double t_start = GetTime();
  //  declaration of database, you need this in every program
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(std::string(argv[1]));
  
  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  TDomain domain(parmoon_db);
 
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  
  // refine grid
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
    domain.RegRefineAll();
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    domain.PS("Domain.ps", It_Finest, 0);
  
  // choose example according to the value of db["example"]
  Example_Darcy2D example(parmoon_db);
  
  //=========================================================================
  Darcy<2> darcy2d(domain, parmoon_db, example);
  darcy2d.assemble();
  darcy2d.solve();
  darcy2d.output();
  //=========================================================================
  
  Output::print<1>("used time: ", GetTime() - t_start, " seconds");
  Output::close_file();
  return 0;
} // end main
