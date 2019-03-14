// =======================================================================
//
// Purpose:     main program for ROM of scalar 2D time-dep. conv.-diffusion
//				problems (at the moment for problems with stationary Dirichlet
//				boundary conditions and problem-specific coefficients)
//
// Author:      Swetlana Giere &  
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Example_TimeCD2D.h>
#include <TimeDiscRout.h>
#include <TimeConvectionDiffusionROM.h>

using namespace std;

int main(int argc, char* argv[])
{
  double t_start = GetTime();
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;
  
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);
  
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  TDomain Domain(parmoon_db);
  
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();
  
  // refine grid
  Domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
   Domain.PS("Domain.ps", It_Finest, 0);
   
  Example_TimeCD2D example( parmoon_db );

  TimeConvectionDiffusionROM<2> rom_2d_rep(parmoon_db, example);
  
  TimeDiscretization& tss = rom_2d_rep.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  
  double end_time = tss.get_end_time();
  double tau = tss.get_step_length();

  rom_2d_rep.assemble_matrices_rhs(true);
  
  rom_2d_rep.compute_initial_solution();
  rom_2d_rep.output(0);

  
  rom_2d_rep.set_system_matrix();

  // ======================================================================
  // time iteration
  // ======================================================================
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    tss.current_step_++;
    TDatabase::TimeDB->CURRENTTIME += tau;
    Output::print<1>("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

    rom_2d_rep.assemble_matrices_rhs(false);
    /* set rom_2d.set_system_rhs(false) if the source term should not be reassembled */
    rom_2d_rep.set_system_rhs();
    rom_2d_rep.solve();
    rom_2d_rep.output(tss.current_step_);
  }
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();
  return 0;
} // end main
