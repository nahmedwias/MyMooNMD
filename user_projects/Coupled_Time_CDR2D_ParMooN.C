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
#include <Example_CoupledCDR2D.h>

#include <TimeDiscRout.h>

int main(int argc, char* argv[])
{
  //take the starting time.
  double t_start = GetTime();

  // Put up domain and databases

  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain Domain(argv[1]);

  //===========================================================================
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);
  OutFile.setf(std::ios::scientific);
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

  //===========================================================================
  /* include the mesh from a mesh generator, for a standard mesh use the
   * build-in function. The GEOFILE describes the boundary of the domain. */
   Domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);

  //===========================================================================
  // do initial refinements of the domain
    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
      Domain.RegRefineAll();

    // ======================================================================
    // Here the calls to Coupled_Time_CDR_2D start.
    // ======================================================================

    // Construct the CoupledCDR_2D object
    // Which example gets constructed is determined by the input file.
    // 100 - test_time dependent example
    // else - Error
    Example_CoupledCDR2D example;
    Coupled_Time_CDR_2D tcdr_system(Domain, example);

    // ======================================================================
    // assemble matrices and right hand side at start time
    tcdr_system.assemble_initial_time();
    // ======================================================================

    double end_time = TDatabase::TimeDB->ENDTIME;
    int step = 0;
    int n_substeps = GetN_SubSteps();

    int image=0;

    tcdr_system.output(image);
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

      tcdr_system.output(step);
    }
    // OutPut("mem after: " << GetMemory()<<endl);
  }
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();

  return 0;
}
