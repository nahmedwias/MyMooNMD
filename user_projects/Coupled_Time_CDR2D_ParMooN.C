/*!
 * @brief Main program to run an example of a system of coupled
 * time dependent CDR equations.
 *
 * @author Clemens Bartsch
 *
 * @date Feb 10, 2016
 */


#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Output2D.h>
#include <Coupled_Time_CDR_2D.h>
#include <Example_TimeCoupledCDR2D.h>

#include <TimeDiscRout.h>

int main(int argc, char* argv[])
{
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
  TDomain Domain(argv[1], parmoon_db);

  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

  //Domain creation
  Domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"]);

  // refine grid up to the coarsest level
  size_t n_ref = Domain.get_n_initial_refinement_steps();
  for(size_t i=0; i<n_ref; i++){
    Domain.RegRefineAll();
  }
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);

  // ======================================================================
  // Here the calls to Coupled_Time_CDR_2D start.
  // ======================================================================

  // Construct the CoupledCDR_2D object
  Example_TimeCoupledCDR2D example(parmoon_db["example"]);
  Coupled_Time_CDR_2D tcdr_system(Domain, parmoon_db, example);

  // ======================================================================
  // assemble matrices and right hand side at start time
  tcdr_system.assemble_initial_time();
  // ======================================================================

  double end_time = TDatabase::TimeDB->ENDTIME;
  int step = 0;
  int n_substeps = GetN_SubSteps();

  tcdr_system.output();
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

      tcdr_system.assemble_uncoupled_part();

      tcdr_system.couple_and_solve();

      tcdr_system.output();
    }

  }

  Output::close_file();

  return 0;
}
