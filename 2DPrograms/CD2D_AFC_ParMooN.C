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
#include "ConvectionDiffusion_AFC.h"

#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// main program
// =======================================================================
int main(int, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);
  
  //open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  TDomain domain(parmoon_db);
  
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    domain.PS("Domain.ps", It_Finest, 0);
   
  //=========================================================================
  ConvectionDiffusion_AFC<2> cd2d(domain, parmoon_db);
  cd2d.assemble(0);
  cd2d.solve(0);
  
  if( cd2d.get_db()["algebraic_flux_correction"].is("afc") )
  {//nonlinear loop necessary
    size_t Max_It = cd2d.get_db()["afc_nonlinloop_maxit"];
    for(unsigned int k = 1;; k++)
    {
      bool converged;
      converged = cd2d.solve(k);

      if ((converged)||(k>= Max_It))
        break;
    }
  }
  
  cd2d.output();
  //=========================================================================
  //db.info(false);
  Output::close_file();
  return 0;
} // end main
