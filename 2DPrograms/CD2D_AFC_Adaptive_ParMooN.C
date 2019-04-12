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
#include "CDErrorEstimator_AFC.h"
#include "RefinementStrategy.h"
#include "LoopInfo.h"
#include "Chrono.h"
#include "LocalAssembling.h"
#include "Multigrid.h"

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
   
  std::string base_name = parmoon_db["output_basename"];
  
  //=========================================================================
  int estimator_type = 1;
  bool conform_grid = true;
  CDErrorEstimator_AFC<2> estimator(estimator_type, conform_grid,
                                Parameter(parmoon_db["space_discretization_type"]));
  RefinementStrategy<2> refinementStrategy(parmoon_db);
  LoopInfo loop_info("adaptive", true, true, 1);
  BlockVector values;
  size_t n_adaptive_steps = domain.get_database()["refinement_max_n_adaptive_steps"];
  bool adaptive_converged;
  for(size_t curr_level = 0; ; ++curr_level)
  {
    Output::print("\nadaptive loop ", curr_level);
    std::ostringstream ostr;
    ostr << base_name << "_" << std::setw(2) << std::setfill('0') << curr_level;
    parmoon_db["output_basename"].set(ostr.str(), false);
    
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
    estimator.estimate(cd2d.get_example(), cd2d.get_function(), 
                       cd2d.get_afc_D_entries(), cd2d.get_afc_alphas());
    estimator.info();
    
    
    ConvectionDiffusion<2>::FEFunction estimated_error(estimator, values);
    cd2d.add_to_output(&estimated_error);
    cd2d.output();
    
    if(estimator.get_estimated_global_error()[estimator_type]<1.0e-3 ||
      curr_level == n_adaptive_steps)
      adaptive_converged = true;
    if(adaptive_converged)
    {
      loop_info.finish(curr_level,
                       estimator.get_estimated_global_error()[estimator_type]);
      break;
    }
    else
    {
      loop_info.print(curr_level,
                      estimator.get_estimated_global_error()[estimator_type]);
      refinementStrategy.apply_estimator(estimator);
      domain.RefineByRefinementStrategy(refinementStrategy, true);
    }
  }
  //=========================================================================
  //db.info(false);
  Output::close_file();
  return 0;
} // end main
