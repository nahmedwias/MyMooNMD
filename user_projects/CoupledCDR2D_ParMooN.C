/*!
 * @brief Main program to run an example of a system of coupled CDR equations.
 *
 * @author Clemens Bartsch
 *
 * @date May 8, 2015
 *
 * @note This was first implemented in MooNMD and then imported to ParMooN. Start of "import":
 * May 26, 2015
 */


#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Output2D.h>
#include <CoupledCDR_2D.h>
#include <Example_CoupledCDR2D.h>

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

  //Domain creation
  Domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"]);

  // refine grid up to the coarsest level
  size_t n_ref = Domain.get_n_initial_refinement_steps();
  for(int i=0; i<n_ref; i++){
    Domain.RegRefineAll();
  }
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);

	/**********************************
	 * Here the calls to CoupledCDR_2D start.
	 **********************************/

	// Construct the CoupledCDR_2D object
	// Which example gets constructed is determined by the input file.
	// 0 - constant_function example
	// else - Error
	Example_CoupledCDR2D example(parmoon_db["example"]);
	CoupledCDR_2D cdrsysObject(Domain, parmoon_db,
	                           example, CoupledCDR_2D::SolvingStrategy::linear_decoupled);
	// Assemble whatever matrices and vectors are necessary.
	cdrsysObject.assembleCDPart();
	// Solve the system.
	cdrsysObject.solve();
	// Call the output method.
	cdrsysObject.output();

  Output::close_file();

	return 0;
}
