#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Example_TimeNSE2D.h>
#include <Time_NSE2D.h>
#include <TimeDiscRout.h>
#include <LoopInfo.h>

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
  TDomain Domain(argv[1], parmoon_db);

  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  parmoon_db.write(Output::get_outfile());
  check_parameters_consistency_NSE(parmoon_db);
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();
  
  // refine grid up to the coarsest level
  size_t n_ref = Domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
    Domain.RegRefineAll();

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);
  
  // set some parameters for time stepping 
  //TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters(0);

  Example_TimeNSE2D example( parmoon_db );
  // create an object of Time_NSE2D class
  Time_NSE2D tnse2d(Domain, parmoon_db, example);
  // assemble everything at the start time
  // this includes assembling of all A's, B's
  // and M's blocks that are necessary 
  tnse2d.assemble_initial_time();
  
  LoopInfo loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1; // full verbosity
 
  double end_time = TDatabase::TimeDB->ENDTIME; 
  int step = 0;
  int n_substeps = GetN_SubSteps();
  int linear_iterations = 0;
    
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
       LoopInfo loop_info("nonlinear");
       loop_info.print_time_every_step = true;
       loop_info.verbosity_threshold = 1; // full verbosity
       for(unsigned int k=0;; k++)
       {
         if(tnse2d.stopIte(k))
         {
           loop_info.finish(k, tnse2d.getFullResidual());
	   linear_iterations+=k;
	   /// @todo provide all parts of the residual 
	   /// @todo loop_info restricted to the solver only
           loop_info_time.print(linear_iterations, tnse2d.getFullResidual());
	   break;
         }
         else
           loop_info.print(k, tnse2d.getFullResidual());
         tnse2d.solve();
         // assemble the nonlinear matrices 
         tnse2d.assemble_nonlinear_term();
         // prepare the matrices for next nonlinear iteration
         tnse2d.assemble_system();
       }
       // post processing: error computations
       // and solutions for visualization
       // FIXME CB: The call to output depends on 'step' but is in the loop
       // over 'substeps' - shouldn't it be outside of the loop?
       tnse2d.output(step);
     }
   }
   
  loop_info_time.finish(linear_iterations, tnse2d.getFullResidual());
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();
  return 0;
}
