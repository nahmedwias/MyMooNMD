// =======================================================================
//
// Purpose:     main program for scalar equations with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 08.08.2014

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>

#include <Time_CD2D.h>


#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_CD2D.h>
#include <TimeDiscRout.h>


using namespace std;

int main(int argc, char* argv[])
{
  double t_start = GetTime();
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  TDomain Domain(argv[1]);  
  
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 2;
  OpenFiles();

  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();
  
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh  
  Domain.ReadGeo(TDatabase::ParamDB->GEOFILE);

  // refine grid up to the coarsest level
  for(int i=0; i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR
         + TDatabase::ParamDB->LEVELS; i++)
    Domain.RegRefineAll();  
  
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain.PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  
  Example_CD2D example;
  Time_CD2D tcd(Domain, example);
  // ======================================================================
  // assemble matrices and right hand side at start time  
  tcd.assemble_initial_time();
  // ======================================================================
  
  double end_time = TDatabase::TimeDB->ENDTIME; 
  int step = 0;
  int n_substeps = GetN_SubSteps();
    
  int image=0;
  double errors[5];
  for(int i=0; i<5; i++)
    errors[i] = 0.;
  tcd.output(0,image, errors);
  // ======================================================================
  // time iteration
  // ======================================================================
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    step++;
    // OutPut("mem before: " << GetMemory()<<endl);
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    for(int j=0; j < n_substeps; ++j)
    {
      SetTimeDiscParameters(1);
      if(step==1)
      {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
      }
      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
      
      OutPut("\nCURRENT TIME: " << TDatabase::TimeDB->CURRENTTIME << endl);
      
      tcd.assemble();
      
      tcd.solve();
      
      tcd.output(step, image, errors);
    }
    // OutPut("mem after: " << GetMemory()<<endl);
  }
  // ======================================================================
  OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0) << " MB" << endl);
  OutPut("used time: " << GetTime() - t_start << "s" << endl);
  // ======================================================================
  CloseFiles();
  return 0;
} // end main









