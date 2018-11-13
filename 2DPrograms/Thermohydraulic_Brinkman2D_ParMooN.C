// =======================================================================
//
// Purpose:     main program for solving a 2D stationary thermohydraulic Brinkman equation in ParMooN
// Author:      Laura Blank
// History:     Implementation started on  11.05.2018

/** @brief This is a 2D program, which uses existing problem classes (CD2D and Brinkman2D).
 * The basic functionality is to solve a coupled flow-porous problem by firstly solving a
 * Brinkman2D problem (with possibly specific permeability/viscosity/effective_viscosity functions
 * and sources/sinks of mass which are specified by the parameters 'coefficient_function_type' and
 * 'use_source_sink_penalty'), whose solution might be used as convection
 * (as a coefficient function via the example and param) in a conv.-diff.-react.
 * problem equipped with a ssource/ink of species/temperature in the domain.
 */

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Brinkman2D.h>
#include <Example_Brinkman2D.h>
#include <Chrono.h>
#include <sys/stat.h>
#include <sys/types.h>


#include <CD2D.h>
#include <Example_CD2D.h>


// ================================================================================
// Functions given as, e.g., permeability, sources/sinks of mass,... for Brinkman2d
// ================================================================================

/* ************************************************************************** */
// This function refers to the case parmoon_db["coefficient_function_type"].is(1)
// It describes, e.g., a permeability field in an analytic way
// A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function)
// for each refinement step and then it is used in the assemble routine

/*void analytic_coefficient_function(double x, double y, double * values)
{
  if ((y < 0.6) || (y > 0.8))
  {
    values[0] = 0.00001;
  }
  else
  {
    values[0] = 1000.0;
  }
}
 */

/*void analytic_coefficient_function(double x, double y, double * values)
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

/*void analytic_coefficient_function(double x, double y, double * values)
{
	//if ( y> 0.8)
	//	values[0] = 0.001;
	//else
	//	values[0] = 0.001;//0.001;
	//
	if (  (y-x >= 0.2)  || (y-0.8/0.7 * x <= - 0.8/0.7 * 0.3)  )
  {
    values[0] = 0.0001;
  }
  else
  {
	  if (    ( (( x >= 0.3) && (x <= 0.4) ) && ( (y >= 0.3) && ( y <= 0.35) ))
		  ||  ( (( x >= 0.5) && (x <= 0.6) ) && ( (y >= 0.5) && ( y <= 0.55) ))
		  ||  ( (( x >= 0.7) && (x <= 0.8) ) && ( (y >= 0.7) && ( y <= 0.75) ))    )
	  {
		  values[0] = 0.0001;
	  }
	  else
	  {
    values[0] = 100000.0;
	  }
	  }
}
 */

/* ************************************************************************** */
// This function refers to the case parmoon_db["coefficient_function_type"].is(1)
// It describes sources and sinks via an approximate delta distribution.
// A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function) for each refinement step and then it is used in the assemble routine

void analytic_coefficient_function(double x, double y, double * values)
{
	//diagonal
	//values[0] = (  1/sqrt((Pi*0.01*0.01) )* exp(-(((x-0.4)*(x-0.4))+((y-0.4)*(y-0.4)))/(0.01*0.01))
	//		         - 1/sqrt((Pi*0.01*0.01) )* exp(-(((x-0.6)*(x-0.6))+((y-0.6)*(y-0.6)))/(0.01*0.01)) )*10000;

	//horizontal
	//values[0] = (  1/sqrt((Pi*0.01*0.01) )* exp(-(((x-0.4)*(x-0.4))+((y-0.4)*(y-0.4)))/(0.01*0.01))
	//		         - 1/sqrt((Pi*0.01*0.01) )* exp(-(((x-0.6)*(x-0.6))+((y-0.4)*(y-0.4)))/(0.01*0.01)) )*10000;

	/*// an approximation to the delta distribution with domain [0,1]^2
  	double a = 0.05;
  	std::vector<double> center_source = {0.4, 15./40.};
  	std::vector<double > center_sink = {0.6, 15./40.};
  	values[0] = 0.;
  	if ( (((x-center_source[0])/a) < 1 && ((x-center_source[0])/a) > -1) &&
  	     (((y-center_source[1])/a) < 1 && ((y-center_source[1])/a) > -1)   )
  	{
  		 values[0] = 1/(a*a)* ( (cos(Pi*(x-center_source[0])/a)+1)/2 )* ( (cos(Pi*(y-center_source[1])/a)+1)/2 )*1; // source
  	}

	if ( (((x-center_sink[0])/a) < 1 && ((x-center_sink[0])/a) > -1) &&
	     (((y-center_sink[1])/a) < 1 && ((y-center_sink[1])/a) > -1)   )
  	{
  		 values[0] -=  1/(a*a)* ( (cos(Pi*(x-center_sink[0])/a)+1)/2 )* ( (cos(Pi*(y-center_sink[1])/a)+1)/2 )*1; // sink
  	}
	 */

	// an approximation to the delta distribution with domain [0, 10] x [0, 6]
	double a = 0.05;
	double r_Wells = 5*0.001;
	double Q_in = 0.03*Pi*r_Wells*r_Wells * 1000;
	std::vector<double> center_source = {4.5, 3};
	std::vector<double > center_sink = {5.5, 3};
	values[0] = 0.;
	if (  (((x-center_source[0])/a) < 1 && ((x-center_source[0])/a) > -1) &&
			(((y-center_source[1])/a) < 1 && ((y-center_source[1])/a) > -1)    )
	{
		values[0] = 1/(a*a)* ( (cos(Pi*(x-center_source[0])/a)+1)/2 )* ( (cos(Pi*(y-center_source[1])/a)+1)/2 ) * Q_in; // source
	}

	if (  (((x-center_sink[0])/a) < 1 && ((x-center_sink[0])/a) > -1) &&
			(((y-center_sink[1])/a) < 1 && ((y-center_sink[1])/a) > -1)   )
	{
		values[0] -=  1/(a*a)* ( (cos(Pi*(x-center_sink[0])/a)+1)/2 )* ( (cos(Pi*(y-center_sink[1])/a)+1)/2 ) * Q_in; // sink
	}
}

// ================================================================================
// Functions given as, e.g., convection,... for CD2D
// ================================================================================

void analytic_convection_function1(double x, double y, double * values)
{
	values[0] = cos(Pi/18);
}


void analytic_convection_function2(double x, double y, double * values)
{
	values[0] = sin(Pi/18);
}

// ================================================================================
// main program
// ================================================================================
int main(int argc, char* argv[])
{

	Output::print(" ");
	Output::print("################################################################");
	Output::print("<<< ParMooN Started: Thermohydraulic Brinkman2D Main Program >>>");
	Output::print("################################################################");
	//for(refinement_n_initial_steps=1; refinement_n_initial_steps <= 6;++refinement_n_initial_steps)
	//{

	// Start a stopwatch which measures the time spent in different parts of the program
	Chrono timer;

	// Declaration of the ParMooN Database (ParamDB) and FE2D Database (basis functions etc.),
	// this is obligatory in every program
	TDatabase Database(argv[1]);
	TFEDatabase2D FEDatabase;

	ParameterDatabase parmoon_db_brinkman = ParameterDatabase::parmoon_default_database();
	std::ifstream fs(argv[1]); // the .dat file is transformed into a stream
	parmoon_db_brinkman.read(fs); // all parameters identified (according to read()) in the
	// stream(.dat-file) are saved in parmoon_db
	fs.close();

	TDomain Domain(parmoon_db_brinkman);

	// Test the mesh:
	//TCollection *coll = Domain.GetCollection(It_Finest, 0);
	//coll->writeMesh("test_mesh_file.mesh");

	// Open OUTFILE, this is where all output is written to (additionally to console)
	Output::set_outfile(parmoon_db_brinkman["outfile"]);
	Output::setVerbosity(parmoon_db_brinkman["verbosity"]);

	// Write all Parameters to the OUTFILE (not to console) for later reference
	parmoon_db_brinkman.write(Output::get_outfile());
	Database.WriteParamDB(argv[0]);

	TFEFunction2D *u1_ptr;
	TFEFunction2D *u2_ptr;

	size_t n_ref =  Domain.get_n_initial_refinement_steps();

	// It is important tu have one single collection for the brinkman2d and the cd2d problems below,
	// such that the brinkman velocity can be used as convection without 'problems'.
	// Therefore, it is declared here in the main scope.
	TCollection *coll;

	// make sure that the brinkman2d object and all information it contains is not deleted
	// before the whole main() has finished
	std::shared_ptr<Brinkman2D> brinkman2d(nullptr);


	// ################################################################################
	// ############################   BRINKMAN2D    ###################################
	// ################################################################################

	switch((int) parmoon_db_brinkman["coefficient_function_type"])
	{
	case 3:
	{
		cout <<" ***** coefficient_function_type 3 detected **** "<< endl;
		for(size_t i = 0; i < n_ref; i++)
		{
			Domain.RegRefineAll();
		}

		Output::print("===============================================================");
		Output::print("Level: ", n_ref);

		// write grid into an Postscript file
		if(parmoon_db_brinkman["output_write_ps"])
			Domain.PS("Domain.ps", It_Finest, 0);

		Example_Brinkman2D example_brinkman(parmoon_db_brinkman);

		timer.restart_and_print("the setup of the domain, database, and example: )");

		// create an object of the Brinkman class
		brinkman2d = std::make_shared<Brinkman2D>(Domain, parmoon_db_brinkman, example_brinkman);

		timer.restart_and_print("constructing the Brinkman2D object: ");
		Output::print<>("Database info: ");
		parmoon_db_brinkman.info();

		Output::print("It is assumed that a mesh and a fitting TFEFunction2D ReadSol() file "
				"are provided and only the finest grid is of interest.");

		// get collection w.r.t. pressure space
		// (if velocity and pressure grid do not coincide, something else has to be done here)
		coll = brinkman2d->get_pressure_space().GetCollection();

		TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
				BoundConditionNoBoundCondition, 0);
		BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());

		TFEFunction2D coefficient_function(&coefficient_function_FEspace, "coefficient_function", "coefficient_function",
				&coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());

		coefficient_function.ReadSol(parmoon_db_brinkman["read_coefficient_function_directory"]);

		TFEFunction2D *coefficient_function_ptr;
		////coefficient_function.Interpolate(&read_coefficient_function);
		coefficient_function_ptr = &coefficient_function;
		brinkman2d->assemble(n_ref, coefficient_function_ptr);
		////brinkman2d->assemble(n_ref);

		timer.restart_and_print("assembling: ");
		brinkman2d->solve();
		timer.restart_and_print("solving: ");
		brinkman2d->output(n_ref);

		// save velocity solution for convection input in CD2D
		u1_ptr = brinkman2d->get_velocity_component(0);
		u2_ptr = brinkman2d->get_velocity_component(1);

		bool write_streamed_velocity = false;
		if (write_streamed_velocity)
		{
			u1_ptr->WriteSol(parmoon_db_brinkman["write_velocity1_directory"], "written_brinkman2d_u1_");
			u2_ptr->WriteSol(parmoon_db_brinkman["write_velocity2_directory"], "written_brinkman2d_u2_");
		}

		timer.restart_and_print("creating the output: ");

		//=========================================================================
		Output::print("<<<<< ParMooN Finished: Thermohydraulic Brinkman2D Main Program >>>>>");

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

			Output::print("======================================================================");
			Output::print("Level: ", i);

			// write grid into an Postscript file
			if(parmoon_db_brinkman["output_write_ps"])
				Domain.PS("Domain.ps", It_Finest, 0);

			Example_Brinkman2D example_brinkman(parmoon_db_brinkman);

			timer.restart_and_print("the setup of the domain, database, and example: )");

			// create an object of the Brinkman2D class
			brinkman2d = std::make_shared<Brinkman2D>(Domain, parmoon_db_brinkman, example_brinkman);

			timer.restart_and_print("constructing the Brinkman2D object: ");
			Output::print<>("Database info: ");
			parmoon_db_brinkman.info();
			coll = brinkman2d->get_pressure_space().GetCollection();

			/* ********************************************************************* */
			// use an analytic coefficient function (defined on top of this file)
			// In the .dat-file, depending on the example, it might be the case that the
			// permeability has to be set equal to -1

			switch((int) parmoon_db_brinkman["coefficient_function_type"])
			{
			case 1:
			{
				cout <<" *****coefficient_function_type 1 detected **** "<<endl;

				Output::print("A spatially varying coefficient function is detected and used in the assemble routine.");
				// create a FEFunction which will serve as the coefficient_function:
				// fe space of piecewise constant functions
				TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
						BoundConditionNoBoundCondition, 0);
				BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());
				TFEFunction2D coefficient_function(
						&coefficient_function_FEspace, "coefficient_function", "coefficient_function",
						&coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());
				coefficient_function.Interpolate(analytic_coefficient_function);
				//coefficient_function_vector.print("coefficient_function");

				bool write_streamed_velocity = false;
				if (write_streamed_velocity)
				{
					coefficient_function.WriteSol(parmoon_db_brinkman["write_coefficient_function_directory"], "written_coeff_fct");
				}

				////NEW LB Debug test:
				// double tmp_values[3];
				// coefficient_function.FindGradient(0.25,0.25,tmp_values);
				// Output::print(" value = ",tmp_values[0]);

				TFEFunction2D *coefficient_function_ptr;
				coefficient_function_ptr = &coefficient_function;
				brinkman2d->assemble(i, coefficient_function_ptr);
			}
			break;
			/* ********************************************************************* */
			// use an external file containing information on, e.g., the permeability field in the format of an FEFunction2D (see ReadSol() and WriteSol())
			// in combination with an appropriate mesh-file as input.
			// Note that the input is then geo_file: '...'.mesh and in 2D a boundary_file: '...'.PRM
			// In the .dat-file, permeability maybe has to be set equal to -1 (see example)

			case 2:
			{
				cout <<" *****coefficient_function_type 2 detected **** "<<endl;
				Output::print("It is assumed that a mesh and a fitting TFEFunction2D ReadSol() file are provided.");

				TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
						BoundConditionNoBoundCondition, 0);
				BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());

				TFEFunction2D coefficient_function(&coefficient_function_FEspace, "coefficient_function", "coefficient_function",
						&coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());

				coefficient_function.ReadSol(parmoon_db_brinkman["read_coefficient_function_directory"]);

				////coefficient_function.Interpolate(&read_coefficient_function);

				TFEFunction2D *coefficient_function_ptr;
				coefficient_function_ptr = &coefficient_function;

				brinkman2d->assemble(i, coefficient_function_ptr);

				break;
			}
			/* ********************************************************************* */
			// default: no coefficient_function used as param
			case 0:
			{
				cout <<" *****coefficient_function_type 0 detected **** "<<endl;

				brinkman2d->assemble( i );

			}
			break;
			default:
				ErrMsg("unknown coefficicent_function_type for Brinkman2D");
				throw(std::runtime_error("unknown coefficicent_function_type for Brinkman2D"));
				break;
			}

			timer.restart_and_print("assembling: ");
			brinkman2d->solve();
			timer.restart_and_print("solving: ");
			brinkman2d->output(i);
			timer.restart_and_print("creating the output: ");


			// save velocity solution for convection input of cd2d
			u1_ptr = brinkman2d->get_velocity_component(0);
			u2_ptr = brinkman2d->get_velocity_component(1);

			if (i == n_ref)
			{
				u1_ptr->WriteSol(parmoon_db_brinkman["write_velocity1_directory"], "written_brinkman2d_u1_");
				u2_ptr->WriteSol(parmoon_db_brinkman["write_velocity2_directory"], "written_brinkman2d_u2_");
			}

			//=========================================================================
			Output::print("<<<<< ParMooN Finished: Thermohydraulic Brinkman2D Main Program >>>>>");

			timer.print_total_time("this Level (total): ");
			timer.reset();

		}

		Output::close_file();
		cout <<" Saved the velocity solution obtained by the Brinkman routine --> into /Home/flow/blank/PARMOON/Tests/Thermohydraulic_Brinkman2D/Brinkman_u1 and .../Brinkman_u2 "<<endl;
	}
	break;
	default:
		ErrMsg("unknown coefficicent_function_type for Brinkman2D");
		throw(std::runtime_error("unknown coefficicent_function_type for Brinkman2D"));
		break;
	}


	// ################################################################################
	// ###############################   CD2D    ######################################
	// ################################################################################

	//  declaration of database, you need this in every program
	/*TDatabase Database;
  TFEDatabase2D FEDatabase;*/

	ParameterDatabase parmoon_db_cd2d = ParameterDatabase::parmoon_default_database();

	std::ifstream fs2(argv[2]);

	parmoon_db_cd2d.read(fs2);
	fs2.close();

	/** set variables' value in TDatabase using argv[1] (*.dat file) */
	TDomain domain(parmoon_db_cd2d, argv[2]);

	//open OUTFILE, this is where all output is written to (additionally to console)
	Output::set_outfile(parmoon_db_cd2d["outfile"]);
	Output::setVerbosity(parmoon_db_cd2d["verbosity"]);

	// write all Parameters to the OUTFILE (not to console) for later reference
	parmoon_db_cd2d.write(Output::get_outfile());
	Database.WriteParamDB(argv[0]);

	// refine grid. Only one input, namely that in Brinkman2d_db is used here.
	for(size_t i = 0; i < n_ref; i++)
	{
		if( i != 0)
			domain.RegRefineAll();
	}
	// write grid into an Postscript file
	if(parmoon_db_cd2d["output_write_ps"])
		domain.PS("Domain.ps", It_Finest, 0);


	//=========================================================================
	CD2D cd2d(domain, parmoon_db_cd2d);

	/* ********************************************************************* */
	// use an analytic coefficient function (defined on top of this file) for, e.g., convection
	// In the .dat-file, permeability has to be set equal to -1

	switch((int) parmoon_db_cd2d["coefficient_function_type"])
	{
	case 3:
	{
		Output::print("A spatially varying coefficient function is detected and used in the assemble routine.");
		// create a FEFunction which will serve as the coefficient_function:
		// OLD 21.06.18: auto coll = cd2d.get_space().GetCollection();
		// fe space of piecewise constant functions
		TFESpace2D convection_function1_FEspace(coll, "convection_function1_FEspace", "s",
				BoundConditionNoBoundCondition, 0);
		TFESpace2D convection_function2_FEspace(coll, "convection_function2_FEspace", "s",
				BoundConditionNoBoundCondition, 0);
		BlockVector convection_function1_vector(convection_function1_FEspace.GetN_DegreesOfFreedom());
		BlockVector convection_function2_vector(convection_function2_FEspace.GetN_DegreesOfFreedom());
		TFEFunction2D convection_function1(
				&convection_function1_FEspace, "convection_function1", "convection_function1",
				&convection_function1_vector.at(0), convection_function1_FEspace.GetN_DegreesOfFreedom());
		TFEFunction2D convection_function2(
				&convection_function2_FEspace, "convection_function2", "convection_function2",
				&convection_function2_vector.at(0), convection_function2_FEspace.GetN_DegreesOfFreedom());
		convection_function1.Interpolate(analytic_convection_function1);
		convection_function2.Interpolate(analytic_convection_function2);

		bool write_streamed_velocity = false;
		if (write_streamed_velocity)
		{
			convection_function1.WriteSol(parmoon_db_cd2d["write_velocity1_directory"], "written_convection1_");
			convection_function2.WriteSol(parmoon_db_cd2d["write_velocity2_directory"], "written_convection2_");
		}

		//NEW LB Debug test:
		// double tmp_values[3];
		// coefficient_function1.FindGradient(0.25,0.25,tmp_values);
		// Output::print(" value = ",tmp_values[0]);

		TFEFunction2D *convection_function1_ptr, *convection_function2_ptr;
		convection_function1_ptr = &convection_function1;
		convection_function2_ptr = &convection_function2;
		cd2d.assemble(convection_function1_ptr, convection_function2_ptr);
		cout <<" *****coefficient_function_type 3 detected **** "<<endl;
	}
	break;
	case 4:
	{
		cout <<" *****coefficient_function_type 4 detected **** "<<endl;
		Output::print("Two spatially varying convection functions from prior Brinkman problem will be used in the assemble routine.");
		// create a FEFunction which will serve as the coefficient_function:
		// OLD 21.06.18:  auto coll = cd2d.get_space().GetCollection();
		// fe space of piecewise constant functions

		//LB Debug 11.09.18	TFESpace2D convection_function1_FEspace(coll, "convection_function1_FEspace", "s",
		//LB Debug 11.09.18			BoundConditionNoBoundCondition, 0, nullptr);
		//LB Debug 11.09.18	TFESpace2D convection_function2_FEspace(coll, "convection_function2_FEspace", "s",
		//LB Debug 11.09.18			BoundConditionNoBoundCondition, 0, nullptr);
		//LB Debug 11.09.18	BlockVector convection_function1_vector(convection_function1_FEspace.GetN_DegreesOfFreedom());
		//LB Debug 11.09.18	BlockVector convection_function2_vector(convection_function2_FEspace.GetN_DegreesOfFreedom());
		//LB Debug 11.09.18	TFEFunction2D convection_function1(
		//LB Debug 11.09.18			&convection_function1_FEspace, "convection_function1", "convection_function1",
		//LB Debug 11.09.18			&convection_function1_vector.at(0), convection_function1_FEspace.GetN_DegreesOfFreedom());
		//LB Debug 11.09.18	TFEFunction2D convection_function2(
		//LB Debug 11.09.18			&convection_function2_FEspace, "convection_function2", "convection_function2",
		//LB Debug 11.09.18			&convection_function2_vector.at(0), convection_function2_FEspace.GetN_DegreesOfFreedom());

		bool write_streamed_velocity = false;
		if (write_streamed_velocity)
		{
			u1_ptr->WriteSol(parmoon_db_cd2d["write_velocity1_directory"], "written_convection1_");
			u2_ptr->WriteSol(parmoon_db_cd2d["write_velocity2_directory"], "written_convection2_");
		}


		//LB Debug 11.09.18	  convection_function1.Interpolate(u1_ptr);
		//LB Debug 11.09.18		convection_function2.Interpolate(u2_ptr);

		////////// convection_function1.WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "written_convection_function1");
		////////// convection_function2.WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "written_convection_function2");

		//NEW LB Debug test:
		// double tmp_values[3];
		// coefficient_function1.FindGradient(0.25,0.25,tmp_values);
		// Output::print(" value = ",tmp_values[0]);

		/*	u1_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u1_incd2d");
		u2_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u2_incd2d");
		 */
		//LB Debug 11.09.18	TFEFunction2D *convection_function1_ptr, *convection_function2_ptr;
		//LB Debug 11.09.18		convection_function1_ptr = &convection_function1;
		//LB Debug 11.09.18	convection_function2_ptr = &convection_function2;
		cd2d.assemble(u1_ptr, u2_ptr);
		//LB Debug 11.09.18	cd2d.assemble(convection_function1_ptr, convection_function2_ptr);

	}
	break;
	/* ********************************************************************* */
	// use an analytic coefficient function (defined on top of this file)
	case 1:
	{
		Output::print("A spatially varying coefficient function is detected and used in the assemble routine.");
		// create a FEFunction which will serve as the coefficient_function:
		// OLD 21.06.18: auto coll = cd2d.get_space().GetCollection();
		// fe space of piecewise constant functions
		TFESpace2D convection_function1_FEspace(coll, "convection_function1_FEspace", "s",
				BoundConditionNoBoundCondition, 0);
		BlockVector convection_function1_vector(convection_function1_FEspace.GetN_DegreesOfFreedom());
		TFEFunction2D convection_function1(
				&convection_function1_FEspace, "convection_function1", "convection_function1",
				&convection_function1_vector.at(0), convection_function1_FEspace.GetN_DegreesOfFreedom());
		convection_function1.Interpolate(analytic_convection_function1);
		//coefficient_function_vector.print("coefficient_function");

		bool write_streamed_velocity = false;
		if (write_streamed_velocity)
		{
			convection_function1.WriteSol(parmoon_db_cd2d["write_coefficient_function_directory"], "written_conv_fct");
		}

		TFEFunction2D *convection_function1_ptr;
		convection_function1_ptr = &convection_function1;
		cd2d.assemble(convection_function1_ptr);
	}
	break;
	/* ********************************************************************* */
	// use an external file containing information on, e.g., the convection field in the format of an FEFunction2D (see ReadSol() and WriteSol())
	// in combination with an appropriate mesh-file as input.
	// Note that the input is then geo_file: ....mesh and in 2D a boundary_file: ...PRM
	case 2:
	{
		Output::print("It is assumed that a mesh and a fitting TFEFunction2D ReadSol() file are provided.");

		// fe space of piecewise constant functions
		// OLD 21.06.18: TCollection *coll = cd2d.get_space().GetCollection();

		//Domain.GetCollection();
		//   cout << "*********** number of cells:" << coll->GetN_Cells()<<endl;
		/// TCollection tmp_collection = TCollection(brinkman2d.get_pressure_space().GetCollection()->GetN_Cells(),
		/// brinkman2d.get_pressure_space().GetCollection()->GetCells());
		/// read_coll = *tmp_collection;

		TFESpace2D coefficient_function_u1_FEspace(coll, "coefficient_function_u1_FEspace", "s",
				BoundConditionNoBoundCondition, 0);
		TFESpace2D coefficient_function_u2_FEspace(coll, "coefficient_function_u2_FEspace", "s",
				BoundConditionNoBoundCondition, 0);
		BlockVector coefficient_function_u1_vector(coefficient_function_u1_FEspace.GetN_DegreesOfFreedom());
		BlockVector coefficient_function_u2_vector(coefficient_function_u2_FEspace.GetN_DegreesOfFreedom());

		TFEFunction2D coefficient_function_u1(&coefficient_function_u1_FEspace, "coefficient_function_u1", "coefficient_function_u1",
				&coefficient_function_u1_vector.at(0), coefficient_function_u1_FEspace.GetN_DegreesOfFreedom());
		TFEFunction2D coefficient_function_u2(&coefficient_function_u2_FEspace, "coefficient_function_u2", "coefficient_function_u2",
				&coefficient_function_u2_vector.at(0), coefficient_function_u2_FEspace.GetN_DegreesOfFreedom());

		//coefficient_function.ReadSol("/Home/flow/blank/PARMOON/Tests/Brinkman2D/Poiseuille_Hannukainen/Tria_mesh/Varying_Coefficient_Function/written_coefficient_function0.Sol");
		//// coefficient_function.ReadSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D/Brinkman_ux0.Sol");

		coefficient_function_u1.ReadSol(parmoon_db_cd2d["read_velocity1_directory"]);
		TFEFunction2D *coefficient_function_u1_ptr = &coefficient_function_u1;

		coefficient_function_u2.ReadSol(parmoon_db_cd2d["read_velocity2_directory"]);
		TFEFunction2D *coefficient_function_u2_ptr = &coefficient_function_u2;

		cd2d.assemble(coefficient_function_u1_ptr, coefficient_function_u2_ptr);
	}
	break;
	case 0:
	{
		cd2d.assemble( );
	}
	break;
	default:
		ErrMsg("unknown coefficicent_function_type for CD2D");
		throw(std::runtime_error("unknown coefficicent_function_type for CD2D"));
		break;
	}

	cd2d.solve();
	cd2d.output();
	//=========================================================================
	//db.info(false);
	Output::close_file();
	return 0;

}

