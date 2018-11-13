// =======================================================================
//
// Purpose:     main program for ROM of scalar 2D time-dep. conv.-diffusion
//				problems (at the moment for problems with stationary Dirichlet
//				boundary conditions and problem-specific coefficients)
//
// Author:      Swetlana Giere
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>

#include <ROM_TCDR2D.h>
#include <ROM_TCDR2D_Rep.h>

#include <Example_TimeCD2D.h>
#include <TimeDiscRout.h>


using namespace std;

int main(int argc, char* argv[])
{
  double t_start = GetTime();
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  for (int i=1; i<argc; ++i)
  {
    std::ifstream fs(argv[i]);
    parmoon_db.read(fs);
    fs.close();
  }

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

  // refine grid up to the coarsest level
  size_t n_ref = Domain.get_n_initial_refinement_steps();
  for(unsigned int i=0; i<n_ref; i++){
    Domain.RegRefineAll();  
  }
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);
  
  Example_TimeCD2D example( parmoon_db );

  ROM_TCDR2D_Rep rom_2d_rep(Domain, parmoon_db, example);
  
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
