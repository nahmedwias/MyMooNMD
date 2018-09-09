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

#include <TimeDiscRout.h>
#include <SNAPS.h>



using namespace std;

int main(int argc, char* argv[])
{
  double t_start = GetTime();
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  TDomain Domain(parmoon_db, argv[1]);
  
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();
  
  
  // Initialize class for storing snapshots
  SNAPS snaps( parmoon_db );
  // refine grid up to the coarsest level
  size_t n_ref = Domain.get_n_initial_refinement_steps();
  for(unsigned int i=0; i<n_ref; i++){
    Domain.RegRefineAll();  
  }
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);
  
  Example_TimeCD2D example( parmoon_db );
  Time_CD2D tcd(Domain, parmoon_db, example);
  // ======================================================================
  // assemble matrices and right hand side at start time  
  tcd.assemble_initial_time();
  // ======================================================================
  
  double end_time = TDatabase::TimeDB->ENDTIME; 
  int step = 0;
  int n_substeps = GetN_SubSteps();
    
  tcd.output();
  // store initial condition as snapshot
  snaps.write_data(tcd.get_solution());
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
      
      Output::print<1>("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
      
      tcd.assemble();
      
      tcd.solve();
      // write the snap shots
      snaps.write_data(tcd.get_solution(), step);
      tcd.descale_stiffness(tau, TDatabase::TimeDB->THETA1);

      if((step-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
        tcd.output();
    }
    // OutPut("mem after: " << GetMemory()<<endl);
  }
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();
  return 0;
} // end main
