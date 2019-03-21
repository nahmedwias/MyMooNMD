// =======================================================================
//
// Purpose:  main program for solving a stationary scalar equation using ParMooN
//
// Author:   Sashikumaar Ganesan
//
// History:  Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <CD2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_CD2D.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();  
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(parmoon_db, argv[1]);

  //open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
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
   
  //=========================================================================
  CD2D cd2d(domain, parmoon_db);
  cd2d.assemble();
  cd2d.solve();
  
//   if( cd2d.get_db()["algebraic_flux_correction"].is("afc") )
//   {//nonlinear loop necessary
//     size_t Max_It = cd2d.get_db()["afc_nonlinloop_maxit"];
//     for(unsigned int k = 1;; k++)
//     {
//       bool converged;
//       
//       converged = cd2d.solve(k);
// 
//       if ((converged)||(k>= Max_It))
// 	break;
//     }
//   }
  
  cd2d.output();  
  
  //=========================================================================
  //db.info(false);
  Output::close_file();
  return 0;
} // end main
