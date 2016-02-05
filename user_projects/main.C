#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <PrRobustNSE2D.h>
#include<Example_NSE2D.h>

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
  return 0;
}