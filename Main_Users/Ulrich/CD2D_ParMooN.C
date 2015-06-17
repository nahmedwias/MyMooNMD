// =============================================================================
//
// Purpose:  main program for solving a stationary scalar equation using ParMooN
//
// Author:   Sashikumaar Ganesan
//
// History:  Implementation started on 22.08.2014
//           rewrite using CD2D class, Ulrich, March, 2015
// =============================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <CD2D.h>
#include <CD2DErrorEstimator.h>

// =============================================================================
// main program
// =============================================================================
int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain Domain(argv[1]);  

  //set PROBLEM_TYPE to CD (convection diffusion) if not yet set
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 1;
  //open OUTFILE, where all output is written to (in addition to console)
  OpenFiles(); // call CloseFiles() at the end of the program
 
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
  
  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */
  Domain.Init(NULL, TDatabase::ParamDB->GEOFILE); // call mesh generator
   
  // refine grid up to the coarsest level
  for(int i=0; i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR; i++)
    Domain.RegRefineAll();  
  
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain.PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  // for adaptive grid refinement:
  double maximal_local_error, estimated_global_error[7], *eta_K = NULL;
  int current_estimator = TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION;
  // flag for either conforming closures or hanging nodes
  bool conforming = TDatabase::ParamDB->GRID_TYPE != 0;
  
  // choose example according to the value of TDatabase::ParamDB->EXAMPLE
  Example_CD2D example;
  
  TCollection* coll = Domain.GetCollection(It_Finest, 0);
  int LEVELS = TDatabase::ParamDB->LEVELS;
  for(int curr_level = 0; curr_level < LEVELS; curr_level++)
  {
    double t_start_level = GetTime(); // time for level
    OutPut("\n\n***********************************************************\n");
    OutPut("GEOMETRY  LEVEL " << curr_level << endl);
    
   
    // refine grid if level is greater than 0
    if((curr_level > 0))
    {
      // regular refinement if
      // adaptive procedure should not yet start
      // or no error estimation
      if((curr_level <= TDatabase::ParamDB->UNIFORM_STEPS)
         ||(! TDatabase::ParamDB->ESTIMATE_ERRORS))
      {
        Domain.RegRefineAll();
      }
      else
        Domain.RefineByErrorEstimator(coll, eta_K, maximal_local_error,
                                      estimated_global_error[current_estimator],
                                      conforming);
      delete [] eta_K; // eta on this level no longer needed
      eta_K = NULL;
      delete coll; // old collection no longer needed
      coll = Domain.GetCollection(It_Finest, 0);
    }
    
    CD2D cd2d(&Domain, &example);
    cd2d.assemble();
    cd2d.solve();
    cd2d.output(curr_level);
    
    if(TDatabase::ParamDB->ESTIMATE_ERRORS)
    {
      eta_K = new double[coll->GetN_Cells()];
      MultiIndex2D Derivatives_All[5] = { D10, D01, D00, D20, D02 };
      TCD2DErrorEstimator CDErrorEstimator(
        TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION,
        cd2d.get_function(), TDatabase::ParamDB->ERROR_CONTROL);
      TFESpace2D* space = cd2d.getSpace();
      TAuxParam2D aux;
      CDErrorEstimator.GetErrorEstimate(5, Derivatives_All, 
                                        example.get_coeffs(), example.get_bc(),
                                        example.get_bd(), &aux,
                                        1, &space, eta_K, &maximal_local_error,
                                        &estimated_global_error[0]);
      OutPut("estimated global error " << setw(10)
             << estimated_global_error[current_estimator] << endl);
    }
    OutPut("Time for level " << curr_level << ": " << GetTime()-t_start_level
           << endl);
  }
  
  CloseFiles();
  return 0;
} // end main
