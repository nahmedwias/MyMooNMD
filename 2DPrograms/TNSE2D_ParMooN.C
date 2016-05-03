#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>


#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>

#include <Example_NSE2D.h>
#include <Time_NSE2D.h>
#include <TimeDiscRout.h>

using namespace std;

int main(int argc, char* argv[])
{
  double t_start=GetTime();
  TDatabase Database;
  TFEDatabase2D FEDatabase2D;
  
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  TDomain Domain(argv[1]);  
  
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 6;

  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh  
  Domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"]);

  // refine grid up to the coarsest level
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    Domain.RegRefineAll();  

  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain.PS("Domain.ps", It_Finest, 0);

  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  // set some parameters for time stepping 
  //TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters(0);

  Example_NSE2D example;
  // create an object of Time_NSE2D class
  Time_NSE2D tnse2d(Domain, parmoon_db, example);
  // assemble everything at the start time
  // this includes assembling of all A's, B's
  // and M's blocks that are necessary 
  tnse2d.assemble_initial_time();
  
  double end_time = TDatabase::TimeDB->ENDTIME; 
  int step = 0;
  int n_substeps = GetN_SubSteps();
    
  int image=0;

  // ======================================================================
  // time iteration
  // ======================================================================
   while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
   {
     step++;
     // Output::print("mem before: ", GetMemory());
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
     for(int j=0; j < n_substeps; ++j)
     {
       // setting the time disc parameters
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
       // prepare the right hand side vector
       // only needed once per time step
       tnse2d.assemble_rhs();
       // assemble the nonlinear matrices
       tnse2d.assemble_nonlinear_term();
       // prepare the matrices for defect computations
       // and solvers
       tnse2d.assemble_system();
       // nonlinear iteration
       for(unsigned int k=0;; k++)
       {
         if(tnse2d.stopIte(k))
           break;
         tnse2d.solve();
         // assemble the nonlinear matrices 
         tnse2d.assemble_nonlinear_term();
         // prepare the matrices for next nonlinear iteration
         tnse2d.assemble_system();
       }
       // post processing: error computations
       // and solutions for visualization
       tnse2d.output(step,image);
     }
   }
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();
  return 0;
}
