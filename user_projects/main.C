#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <PrRobustNSE2D.h>
#include <Example_NSE2D.h>
#include <PrRobustTime_NSE2D.h>
#include <TimeDiscRout.h>

#include <MooNMD_Io.h>
#include <sys/stat.h>


int main(int argc, char *argv[])
{
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  TDomain Domain(argv[1]);
  // set PROBLEM_TYPE to NSE or STOKES
  // PROBLEM_TYPE = 3, STOKES
  // PROBLEM_TYPE = 5, NSE
  if(TDatabase::ParamDB->PROBLEM_TYPE != 3 &&
     TDatabase::ParamDB->PROBLEM_TYPE !=5)
    TDatabase::ParamDB->PROBLEM_TYPE = 5;
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);
  
  Database.CheckParameterConsistencyNSE();
  Database.WriteParamDB(argv[0]);
  
  Domain.Init(TDatabase::ParamDB->BNDFILE,
              TDatabase::ParamDB->GEOFILE);
  
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain.RegRefineAll();
  
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR,0777);
  
  Example_NSE2D example;
  
  
  for(int i=0; i<TDatabase::ParamDB->LEVELS;i++)
  {
    Output::print("*******************************************************" );
    Output::print("******           GEOMETRY  LEVEL ", i                    );
    Output::print("*******************************************************" );    
    Domain.RegRefineAll();
    PrRobustNSE2D stokes_prbst(Domain, example);
    
    if(TDatabase::ParamDB->PROBLEM_TYPE == 3)
    {
      stokes_prbst.assembleMatrixRhs();
      stokes_prbst.solve();
      stokes_prbst.output(i);
    }    
  }
  
  Output::print<1>('\n');
  Output::print('\n','\n',"*******************************************************" );  
  Output::print<1>("STARTED TIME DEPENDENT STOKES PROBLEM ");
  Output::print("*******************************************************" );
  
  if(TDatabase::ParamDB->PROBLEM_TYPE != 3 &&
     TDatabase::ParamDB->PROBLEM_TYPE !=5)
    TDatabase::ParamDB->PROBLEM_TYPE = 5;
  Output::set_outfile("outtimense.out");
  
  Database.CheckParameterConsistencyNSE();
  Database.WriteParamDB(argv[0]);
  
  TDatabase::ParamDB->EXAMPLE = 108;
  
  // TDatabase::ParamDB->DISCTYPE = 300;
  // TDatabase::ParamDB->RE_NR = 1.0;
  // TDatabase::ParamDB->LAPLACETYPE = 1;
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain.PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  SetTimeDiscParameters(0);
  Domain.Init(TDatabase::ParamDB->BNDFILE,
              TDatabase::ParamDB->GEOFILE);
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain.RegRefineAll();
  
  Example_NSE2D example_time;
  PrRobustTime_NSE2D stokes_prbst_time(Domain, example_time);
  stokes_prbst_time.assemble_initial_timeNS();
  
  
  int step = 0;
  int image=0;
  int n_substeps = GetN_SubSteps();
  double end_time = TDatabase::TimeDB->ENDTIME;
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    step++;
     // Output::print("mem before: ", GetMemory());
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
     for(int j=0; j < n_substeps; ++j)
     {
       SetTimeDiscParameters(1);
       if(step==1)
       {
         Output::print<1>("Theta1: ", TDatabase::TimeDB->THETA1);
         Output::print<1>("Theta2: ", TDatabase::TimeDB->THETA2);
         Output::print<1>("Theta3: ", TDatabase::TimeDB->THETA3);
         Output::print<1>("Theta4: ", TDatabase::TimeDB->THETA4);
       }
       double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
       TDatabase::TimeDB->CURRENTTIME += tau;
       
       Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

       stokes_prbst_time.assemble_rhsNS();
       
       stokes_prbst_time.assemble_system_matrix();
       
       stokes_prbst_time.solve();
       
       stokes_prbst_time.output(step, image);
       
     }    
   // exit(0);
  }
  
  
}