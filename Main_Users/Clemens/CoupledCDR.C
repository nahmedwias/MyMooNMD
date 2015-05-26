/*!
 * @brief Main program to run an example of a system of coupled CDR equations.
 *
 * @author Clemens Bartsch
 *
 * @date May 8, 2015
 */



#include <FEFunction2D.h>
#include <NSE2DMainRoutines.h>
//#include <LocalIntegrals2D.h>
#include <LocalForms2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <Domain.h>
#include <Database.h>
#include <Output2D.h>
#include <LinAlg.h> //IntoL20FEFunction

#include <DirectSolver.h>
#include <PardisoSolver.h>
#include <ConvDiff2D.h>
#include <BlockVector2D.h>

#include <CDR_2D_System.h>
#include <Example_CoupledCDR2D.h>
#include <CDRDatabase.h>

int main(int argc, char* argv[])
{
	// Put up domain and databases
	double t_start = GetTime();
	TDomain *Domain = new TDomain();
	TDatabase *Database = new TDatabase();
	TFEDatabase2D *FEDatabase = new TFEDatabase2D();

	//Create a CDRDatabase, which is a hack needed for the assembling (see CDRDatabase.h).
	CDRDatabase* auxiliaryDataBaseForAssembling = new CDRDatabase();

	//===========================================================================
	// read input
	Domain->ReadParam(argv[1]);
	OpenFiles();
	OutFile.setf(std::ios::scientific);
	Database->WriteParamDB(argv[0]);

	//===========================================================================
	// create domain, do initial refinements
	Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);
	for(int i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE
	+ TDatabase::ParamDB->LEVELS; i++)
		Domain->RegRefineAll();
	//===========================================================================
	TCollection *coll = Domain->GetCollection(It_Finest,0);
	if(TDatabase::ParamDB->WRITE_PS) Domain->PS("grid.ps", coll);
	//===========================================================================

	/**********************************
	 * Here the calls to CDR_2D_Systems start.
	 **********************************/

	// Construct the object CDR _2D_System object
	// Which example gets constructed is determined by the input file.
	// 0 - constant_function example
	// else - Error
	Example_CoupledCDR2D* example = new Example_CoupledCDR2D();
	CDR_2D_System cdrsysObject(coll, example, CDR_2D_System::SolvingStrategy::linear_decoupled);
	// Assemble whatever matrices and vectors are necessary.
	cdrsysObject.assemble();
	// Solve the system.
	cdrsysObject.solve();
	// Call the output method.
	cdrsysObject.output();
	/**********************************
	 * Unchanged code from Stokes.C
	 **********************************/
	OutPut("MEMORY   : " << setw(8) << GetMemory()/(1048576.0) << " MB" << endl);
	OutPut("used time: " << GetTime() - t_start << "s" << endl);
	CloseFiles();

	return 0;
}
