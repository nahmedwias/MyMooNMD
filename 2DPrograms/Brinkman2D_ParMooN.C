// =======================================================================
//
// Purpose:     main program for solving a stationary Brinkman equation in ParMooN
// Author:      Laura Blank
// History:     Implementation started on  14.03.2016

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <ParameterDatabase.h>
#include <LinAlg.h>
#include <Brinkman2D.h>
#include <MainUtilities.h>
#include <LocalAssembling.h>
#include <Example_Brinkman2D.h>
#include <Chrono.h>
#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>


/* ************************************************************************** */
// This function refers to the case parmoon_db["coefficient_function_type"].is(1)
// It describes, e.g., a permeability field in an analytic way
// A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function) for each refinement step and then it is used in the assemble routine

// for Discacciati_Flow.h : permeability
/*
void analytic_coefficient_function(double x, double y, double * values)
{
	if ( ( ( (x >=7  ) && (x <= 8  ) )  && ( (y >= 4) && (y <= 6) ) ) )
	  {
	  	 values[0] =  0.000001;
	  }
	else if ( ( ( (x >= 9) && (x <= 10) )  && ( (y>=2) && (y <= 4.2) ) ) )
	  {
	  	 values[0] =  0.000001;
	  }
  else if ( ( (x >=5  ) && ( (y>=2) && (y <= 6) ) ) )
  {
    values[0] =   0.01;
  }
  else
  {
    values[0] =  10000000000;
  }
}*/

/* ************************************************************************** */
// This function refers to the case parmoon_db["coefficient_function_type"].is(1)
// It describes, e.g., a permeability field in an analytic way
// A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function) for each refinement step and then it is used in the assemble routine
/*
void analytic_coefficient_function(double x, double y, double * values)
{
  if ( ( (x < 0.5) && (y<0.5) ) || ( (x > 0.5) && (y>0.5) ) )
  {
    values[0] = 0.01;
  }
  //   else if ((x >0.79) && (x<0.8) && (y < 0.3) && (y > 0.29) )
  //     values[0]=0.00001;
  else
  {
    values[0] = 1000.0;
  }
}
*/

/* ************************************************************************** */
// This function refers to the case parmoon_db["coefficient_function_type"].is(1)
// It describes, e.g., a certain permeability field in an analytic way
// A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function) for each refinement step and then it is used in the assemble routine
/*
void analytic_coefficient_function(double x, double y, double * values)
{
  if (  (y-x >= 0.2)  || (y-0.8/0.7 * x <= - 0.8/0.7 * 0.3)  )
  {
    values[0] = 0.01;
  }
  else
  {
    values[0] = 1000.0;
  }
}
*/

/* ************************************************************************** */
/* This function refers to the case parmoon_db["coefficient_function_type"].is(1)
* and the riverbed example
* It describes, e.g., a certain permeability field in an analytic way
* A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function)
*  for each refinement step and then it is used in the assemble routine
*  */
void analytic_coefficient_function(double x, double y, double * values)
{
  if (y <= 1.5)
  {
    values[0] = 0.01;
  }
  else
  {
    values[0] = 1000.0;
  }
}

/* ************************************************************************** */
// This function refers to the case parmoon_db["coefficient_function_type"].is(1)
// It describes sources and sinks via an approximate delta distribution.
// A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function) for each refinement step and then it is used in the assemble routine
/* void analytic_coefficient_function(double x, double y, double * values)
  {
  //values[0] = ( 1/sqrt((Pi*0.01*0.01) )* exp(-(((x-0.4)*(x-0.4))+((y-0.4)*(y-0.4)))/(0.01*0.01)) - 1/sqrt((Pi*0.01*0.01) )* exp(-(((x-0.6)*(x-0.6))+((y-0.6)*(y-0.6)))/(0.01*0.01)) )*10000;

  //another approx.
  	double a = 0.1;
  	std::vector<double> center_source = {0.4,0.4};
  	std::vector<double > center_sink = {0.6,0.4};
  	if ( ( ((x-center_source[0])/a) < 1 && ((x-center_source[0])/a) > -1) && (((y-center_source[1])/a) < 1 && ((y-center_source[1])/a) > -1))
  	{
  		 values[0] = 1/(a*a)* ( (cos(Pi*(x-center_source[0])/a)+1)/2 )* ( (cos(Pi*(y-center_source[1])/a)+1)/2 ); // source
  	}
  	else
  	{
  		values[0] = 0;
  	}

  	if ( ( ((x-center_sink[0])/a) < 1 && ((x-center_sink[0])/a) > -1) && (((y-center_sink[1])/a) < 1 && ((y-center_sink[1])/a) > -1))
  	{
  		 values[0] +=  - 1/(a*a)* ( (cos(Pi*(x-center_sink[0])/a)+1)/2 )* ( (cos(Pi*(y-center_sink[1])/a)+1)/2 ); // sink
  	}
  	else
  	{
  		values[0] = 0;
  	}
  }
*/

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
	Output::print(" ");
	Output::print("###################################################################");
	Output::print("###################################################################");

	// Start a stopwatch which measures the time spent in different parts of the program
	Chrono timer;

	// Declaration of the ParMooN Database (ParamDB) and FE2D Database
	// (basis functions etc.), this is obligatory in every program
	TDatabase Database;
	TFEDatabase2D FEDatabase;

	ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
	std::ifstream fs(argv[1]); // the .dat file is transformed into a stream
	parmoon_db.read(fs); // all parameters identified (according to read()) in the stream(.dat-file) are saved in parmoon_db
	fs.close();

	// Set each variables' value in TDatabase using argv[1] (*.dat file)
	TDomain Domain(parmoon_db, argv[1]);

	// Test the mesh:
	//TCollection *coll = Domain.GetCollection(It_Finest, 0);
	//coll->writeMesh("test_mesh_file.mesh");

	// Open OUTFILE, this is where all output is written to (addionally to console)
	Output::set_outfile(parmoon_db["outfile"]);
	Output::setVerbosity(parmoon_db["verbosity"]);

	// Write all Parameters to the OUTFILE (not to console) for later reference
	parmoon_db.write(Output::get_outfile());
	Database.WriteParamDB(argv[0]);

	// Refine the grid
	size_t n_ref =  Domain.get_n_initial_refinement_steps();

	//******************   Refinement steps   ********************//
	switch((int) parmoon_db["coefficient_function_type"])
	{
	case 3:
	{
		for (size_t k = 0; k < n_ref; k++)
		{
			Domain.RegRefineAll();
		}

		Output::print("===============================================================");
		Output::print("Level: ", n_ref);

		// write grid into an Postscript file
		if(parmoon_db["output_write_ps"])
			Domain.PS("Domain.ps", It_Finest, 0);

		Example_Brinkman2D example(parmoon_db);

		timer.restart_and_print("the setup of the domain, database, and example: )");

		// create an object of the Brinkman class
		Brinkman2D brinkman2d(Domain, parmoon_db, example);
		timer.restart_and_print("constructing the Brinkman2D object: ");
		Output::print<>("Database info: ");
		parmoon_db.info();

		/////////////// Routines for periodic boundary conditions /////////////////
		if(parmoon_db["example"].is(11))
		{
			brinkman2d.findPeriodicDOFs();
			brinkman2d.checkPeriodicDOFs();
			brinkman2d.makePeriodicBoundary();
		}

		timer.restart_and_print("assembling: ");
		brinkman2d.solve();

		timer.restart_and_print("solving: ");
		brinkman2d.output(n_ref);
		timer.restart_and_print("creating the output: ");

		//=========================================================================
		Output::print("<<<<< ParMooN Finished: Brinkman2D Main Program >>>>>");

		timer.print_total_time("this Level (total): ");
		timer.reset();
	}
	break;
	case 0: case 1: case 2:
	{
		for(size_t i = 0; i < n_ref+1; i++)
		{
			if (i != 0)
				Domain.RegRefineAll();

			Output::print("==============================================================");
			Output::print("Level: ", n_ref);

			// write grid into an Postscript file
			if(parmoon_db["output_write_ps"])
				Domain.PS("Domain.ps", It_Finest, 0);

			Example_Brinkman2D example(parmoon_db);

			timer.restart_and_print("the setup of the domain, database, and example: )");

			// create an object of the Brinkman class
			Brinkman2D brinkman2d(Domain, parmoon_db, example);
			timer.restart_and_print("constructing the Brinkman2D object: ");
			Output::print<>("Database info: ");
			parmoon_db.info();

			/* ********************************************************************* */
			// use an analytic coefficient function (defined on top of this file)
			// In the .dat-file, permeability has to be set equal to -1
			switch((int) parmoon_db["coefficient_function_type"])
			{
			case 1:
			{
				cout <<" *****coefficient_function_type 1 detected **** "<<endl;
				Output::print("A spatially varying coefficient function is detected and used in the assemble routine.");
				// create a FEFunction which will serve as the coefficient_function:
				auto coll = brinkman2d.get_pressure_space().GetCollection();
				// fe space of piecewise constant functions
				TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
						BoundConditionNoBoundCondition, 0);
				BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());
				TFEFunction2D coefficient_function(
						&coefficient_function_FEspace, "coefficient_function", "coefficient_function",
						&coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());
				coefficient_function.Interpolate(analytic_coefficient_function);
				//coefficient_function_vector.print("coefficient_function");

				////coefficient_function.WriteSol("/Home/flow/blank/PARMOON/Tests/Brinkman2D/Poiseuille_Hannukainen/Tria_mesh/Varying_Coefficient_Function", "written_coefficient_function_new");

				//NEW LB Debug test:
				// double tmp_values[3];
				// coefficient_function.FindGradient(0.25,0.25,tmp_values);
				// Output::print(" value = ",tmp_values[0]);

				TFEFunction2D *coefficient_function_ptr;
				coefficient_function_ptr = &coefficient_function;
				brinkman2d.assemble(i, coefficient_function_ptr);
			}
			break;
			/* ********************************************************************* */
			// use an external file containing information on, e.g., the permeability field in the format of an FEFunction2D (see ReadSol() and WriteSol())
			// in combination with an appropriate mesh-file as input.
			// Note that the input is then geo_file: ....mesh and in 2D a boundary_file: ...PRM
			// In the .dat-file, permeability has to be set equal to -1
			case 2:
			{
				cout <<" ***** coefficient_function_type 2 detected **** "<< endl;
				Output::print("It is assumed that a mesh and a fitting TFEFunction2D ReadSol() file are provided.");

				// fe space of piecewise constant functions
				TCollection *coll = brinkman2d.get_pressure_space().GetCollection();

				//Domain.GetCollection();
				//   cout << "*********** number of cells:" << coll->GetN_Cells()<<endl;
				/// TCollection tmp_collection = TCollection(brinkman2d.get_pressure_space().GetCollection()->GetN_Cells(),
				/// brinkman2d.get_pressure_space().GetCollection()->GetCells());
				/// read_coll = *tmp_collection;

				TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
						BoundConditionNoBoundCondition, 0);
				BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());

				TFEFunction2D coefficient_function(&coefficient_function_FEspace, "coefficient_function", "coefficient_function",
						&coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());

				coefficient_function.ReadSol(parmoon_db["read_coefficient_function_directory"]);

				////coefficient_function.Interpolate(&read_coefficient_function);

				TFEFunction2D *coefficient_function_ptr;
				coefficient_function_ptr = &coefficient_function;

				brinkman2d.assemble(i, coefficient_function_ptr);
			}
			break;
			/* ********************************************************************* */
			// just hand in a constant permeability
			case 0:
			{
				brinkman2d.assemble( i );
			}
			break;
			default:
				ErrMsg("unknown coefficicent_function_type for Brinkman2D");
				throw(std::runtime_error("unknown coefficicent_function_type for Brinkman2D"));
				break;
			}

			/////////////// Routines for periodic boundary conditions /////////////////
			if(parmoon_db["example"].is(11))
			{
				brinkman2d.findPeriodicDOFs();
				brinkman2d.checkPeriodicDOFs();
				brinkman2d.makePeriodicBoundary();
			}

			///////////////// ///////////////// ///////////////// ////////////////
			/*
   {
    	     auto blocks = brinkman2d.get_matrix().get_blocks_uniquely();
    			 int const * const rowPtr = blocks[8]->GetRowPtr();
						int const * const colPtr = blocks[8]->GetKCol();
						double const * const entries = blocks[8]->GetEntries();
						const int n_rows = blocks[8]->GetN_Rows();
    }
			 */

			//New LB 28.08.18
			// brinkman2d.get_rhs().write("righthandside_Periodic_P2P1_Observations_Poiseuille");
			// brinkman2d.get_matrix().get_combined_matrix()->write("matrix_Periodic_P2P1_Observations_Poiseuille");

			// brinkman2d.get_rhs().write("righthandside_periodic_P2P1");
			// brinkman2d.get_matrix().get_combined_matrix()->write("matrix_periodic_P2P1");
			// brinkman2d.get_rhs().write("righthandside_Neumann_P2P1");
			// brinkman2d.get_matrix().get_combined_matrix()->write("matrix_Neumann_P2P1");

			timer.restart_and_print("assembling: ");
			brinkman2d.solve();

			// brinkman2d.get_solution().write("solution_Neumann_test");

			timer.restart_and_print("solving: ");
			brinkman2d.output(i);
			timer.restart_and_print("creating the output: ");

			//=========================================================================
			Output::print("<<<<< ParMooN Finished: Brinkman2D Main Program >>>>>");

			timer.print_total_time("this Level (total): ");
			timer.reset();
		}
		break;
	default:
		ErrMsg("unknown coefficicent_function_type for Brinkman2D");
		throw(std::runtime_error("unknown coefficicent_function_type for Brinkman2D"));
		break;

	}
	}
	Output::close_file();
	return 0;
} // end main

