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

#include <FEFunction2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <Domain.h>
#include <Database.h>
#include <Output2D.h>
#include <LinAlg.h> //IntoL20FEFunction

#include <DirectSolver.h>
#include <PardisoSolver.h>
#include <ConvDiff2D.h>

#include <CoupledCDR_2D.h>
#include <Example_CoupledCDR2D.h>

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
	OpenFiles();
	OutFile.setf(std::ios::scientific);
	Database.WriteParamDB(argv[0]);

	//===========================================================================
	/* include the mesh from a mesh generator, for a standard mesh use the
	 * build-in function. The GEOFILE describes the boundary of the domain. */
	 Domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); // call mesh generator

	//===========================================================================
	// do initial refinements of the domain
	  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
	    Domain.RegRefineAll();

	/**********************************
	 * Here the calls to CoupledCDR_2D start.
	 **********************************/

	// Construct the CoupledCDR_2D object
	// Which example gets constructed is determined by the input file.
	// 0 - constant_function example
	// else - Error
	Example_CoupledCDR2D example;
	CoupledCDR_2D cdrsysObject(Domain, example, CoupledCDR_2D::SolvingStrategy::linear_decoupled);
	// Assemble whatever matrices and vectors are necessary.
	cdrsysObject.assembleCDPart();
	// Solve the system.
	cdrsysObject.solve();
	// Call the output method.
	cdrsysObject.output();

	/**********************************
	 * Print information on used time and memory (Unchanged code from Stokes.C)
	 **********************************/
	OutPut("MEMORY   : " << setw(8) << GetMemory()/(1048576.0) << " MB" << endl);
	OutPut("used time: " << GetTime() - t_start << "s" << endl);
	CloseFiles();

	return 0;
}
