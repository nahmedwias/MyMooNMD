// =======================================================================
//
// Purpose:     main program for solving a 2D stationary thermohydraulic Brinkman equation in ParMooN
// Author:      Laura Blank
// History:     Implementation started on  11.05.2018

// A Thermohydraulic Brinkman problem is solved, whose solution is used as convection (as a coefficient function via the example) in a conv-diff problem for the temperature

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


/* ************************************************************************** */
// This function refers to the case parmoon_db["coefficient_function_type"].is(1)
// It describes, e.g., a permeability field in an analytic way
// A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function) for each refinement step and then it is used in the assemble routine

void analytic_coefficient_function(double x, double y, double * values)
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


/* ************************************************************************** */
// This function refers to the case parmoon_db["coefficient_function_type"].is(1)
// It describes, e.g., a permeability field in an analytic way
// A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function) for each refinement step and then it is used in the assemble routine

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



/* ************************************************************************** */
// This function refers to the case parmoon_db["coefficient_function_type"].is(1)
// It describes, e.g., a certain permeability field in an analytic way
// A corresponding FEFunction2D is constructed via interpolate(analytic_coefficient_function) for each refinement step and then it is used in the assemble routine
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

/*
  void analytic_coefficient_function(double x, double y, double * values)
  {
  values[0] = ( 1/sqrt((Pi*0.01*0.01) )* exp(-(((x-0.4)*(x-0.4))+((y-0.4)*(y-0.4)))/(0.01*0.01)) - 1/sqrt((Pi*0.01*0.01) )* exp(-(((x-0.6)*(x-0.6))+((y-0.6)*(y-0.6)))/(0.01*0.01)) )*100;
  }
 */

void analytic_convection_function1(double x, double y, double * values)
{
	//	if ( y> 0.8)
	values[0] = cos(Pi/18);
	//	else
	//	values[0] = 0.001;//0.001;
}

void analytic_convection_function2(double x, double y, double * values)
{
	//	if ( y> 0.8)
	values[0] = sin(Pi/18);
	//	else
	//	values[0] = 0.001;//0.001;
}

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{

	Output::print(" ");
	Output::print("################################################################################################################");
	Output::print("################################################################################################################");
	//for(refinement_n_initial_steps=1; refinement_n_initial_steps <= 6;++refinement_n_initial_steps)
	//{

	// Start a stopwatch which measures the time spent in different parts of the program
	Chrono timer;

	// Declaration of the ParMooN Database (ParamDB) and FE2D Database (basis functions etc.), this is obligatory in every program
	TDatabase Database;
	TFEDatabase2D FEDatabase;

	ParameterDatabase parmoon_db_brinkman = ParameterDatabase::parmoon_default_database();
	std::ifstream fs(argv[1]); // the .dat file is transformed into a stream
	parmoon_db_brinkman.read(fs); // all parameters identified (according to read()) in the stream(.dat-file) are saved in parmoon_db
	fs.close();

	// Set each variables' value in TDatabase using argv[1] (*.dat file)
	TDomain Domain(parmoon_db_brinkman, argv[1]);

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
	// Refine the grid
	size_t n_ref =  Domain.get_n_initial_refinement_steps();
	TCollection *coll; // New 21.06.18

	// todo shared pointer brinkman objekt NEW 21.06.18
	////std::shared_ptr<Brinkman2D> brinkman2d_global(nullptr);
	std::shared_ptr<Brinkman2D> brinkman2d(nullptr);

	// ############################################################################################################### //
	// ###########################################    BRINKMAN2D    ################################################## //

	if  (n_ref-1 > 0 && parmoon_db_brinkman["coefficient_function_type"].is(3))
	{
		for(size_t i = 0; i < n_ref; i++)
		{
			Domain.RegRefineAll();
		}
		Output::print("========================================================================================================");
		Output::print("Level: ", n_ref-1);

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

		Output::print("It is assumed that a mesh and a fitting TFEFunction2D ReadSol() file are provided.");

		// fe space of piecewise constant functions
		//OLD 21.06.18: TCollection *coll = brinkman2d.get_pressure_space().GetCollection();
		coll = brinkman2d->get_pressure_space().GetCollection(); // New 21.06.18

		//Domain.GetCollection();
		//   cout << "*********** number of cells:" << coll->GetN_Cells()<<endl;
		/// TCollection tmp_collection = TCollection(brinkman2d.get_pressure_space().GetCollection()->GetN_Cells(),
		/// brinkman2d.get_pressure_space().GetCollection()->GetCells());
		/// read_coll = *tmp_collection;

		TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
				BoundConditionNoBoundCondition, 0, nullptr);
		BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());

		TFEFunction2D coefficient_function(&coefficient_function_FEspace, "coefficient_function", "coefficient_function",
				&coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());

		////////// coefficient_function.ReadSol("/Home/flow/blank/PARMOON/Tests/Thermohydraulic_Brinkman2D/written_coefficient_function0.Sol");

		////coefficient_function.Interpolate(&read_coefficient_function);

		TFEFunction2D *coefficient_function_ptr;
		coefficient_function_ptr = &coefficient_function;

		////////// brinkman2d.assemble(coefficient_function_ptr);
		brinkman2d->assemble();

		timer.restart_and_print("assembling: ");
		brinkman2d->solve();
		timer.restart_and_print("solving: ");
		brinkman2d->output(n_ref);

		// save velocity solution for convection input of cd2d
		u1_ptr = brinkman2d->get_velocity_component(0);
		u2_ptr = brinkman2d->get_velocity_component(1);

		u1_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u1");
		u2_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u2");


		timer.restart_and_print("creating the output: ");

		//=========================================================================
		Output::print("<<<<< ParMooN Finished: Thermohydraulic Brinkman2D Main Program >>>>>");

		timer.print_total_time("this Level (total): ");
		timer.reset();

	}
	else
	{
		for(size_t i = 0; i < n_ref; i++)
		{
			Domain.RegRefineAll();

			Output::print("========================================================================================================");
			Output::print("Level: ", i);

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
			coll = brinkman2d->get_pressure_space().GetCollection();

			/* ********************************************************************* */
			// use an analytic coefficient function (defined on top of this file)
			// In the .dat-file, permeability has to be set equal to -1
			if (parmoon_db_brinkman["coefficient_function_type"].is(1))
			{
				Output::print("A spatially varying coefficient function is detected and used in the assemble routine.");
				// create a FEFunction which will serve as the coefficient_function:
				// Old 21.06.18: auto coll = brinkman2d.get_pressure_space().GetCollection();
				// fe space of piecewise constant functions
				TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
						BoundConditionNoBoundCondition, 0, nullptr);
				BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());
				TFEFunction2D coefficient_function(
						&coefficient_function_FEspace, "coefficient_function", "coefficient_function",
						&coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());
				coefficient_function.Interpolate(analytic_coefficient_function);
				//coefficient_function_vector.print("coefficient_function");

				coefficient_function.WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "written_coefficient_function_new");

				//NEW LB Debug test:
				// double tmp_values[3];
				// coefficient_function.FindGradient(0.25,0.25,tmp_values);
				// Output::print(" value = ",tmp_values[0]);

				TFEFunction2D *coefficient_function_ptr;
				coefficient_function_ptr = &coefficient_function;
				brinkman2d->assemble(coefficient_function_ptr);


			}

			/* ********************************************************************* */
			// use an external file containing information on, e.g., the permeability field in the format of an FEFunction2D (see ReadSol() and WriteSol())
			// in combination with an appropriate mesh-file as input.
			// Note that the input is then geo_file: ....mesh and in 2D a boundary_file: ...PRM
			// In the .dat-file, permeability has to be set equal to -1
			else if (parmoon_db_brinkman["coefficient_function_type"].is(2))
			{
				Output::print("It is assumed that a mesh and a fitting TFEFunction2D ReadSol() file are provided.");

				// fe space of piecewise constant functions
				//Old 21.06.18: TCollection *coll = brinkman2d.get_pressure_space().GetCollection();


				//Domain.GetCollection();
				//   cout << "*********** number of cells:" << coll->GetN_Cells()<<endl;
				/// TCollection tmp_collection = TCollection(brinkman2d.get_pressure_space().GetCollection()->GetN_Cells(),
				/// brinkman2d.get_pressure_space().GetCollection()->GetCells());
				/// read_coll = *tmp_collection;

				TFESpace2D coefficient_function_FEspace(coll, "coefficient_function_FEspace", "s",
						BoundConditionNoBoundCondition, 0, nullptr);
				BlockVector coefficient_function_vector(coefficient_function_FEspace.GetN_DegreesOfFreedom());

				TFEFunction2D coefficient_function(&coefficient_function_FEspace, "coefficient_function", "coefficient_function",
						&coefficient_function_vector.at(0), coefficient_function_FEspace.GetN_DegreesOfFreedom());

				coefficient_function.ReadSol("/Home/flow/blank/PARMOON/Tests/Thermohydraulic_Brinkman2D/written_coefficient_function0.Sol");

				////coefficient_function.Interpolate(&read_coefficient_function);

				TFEFunction2D *coefficient_function_ptr;
				coefficient_function_ptr = &coefficient_function;

				brinkman2d->assemble(coefficient_function_ptr);
			}

			/* ********************************************************************* */
			// just hand in a constant permeability
			else
			{
				brinkman2d->assemble( );
			}

			timer.restart_and_print("assembling: ");
			brinkman2d->solve();
			timer.restart_and_print("solving: ");
			brinkman2d->output(i+1);
			timer.restart_and_print("creating the output: ");


			// save velocity solution for convection input of cd2d
			u1_ptr = brinkman2d->get_velocity_component(0);
			u2_ptr = brinkman2d->get_velocity_component(1);

			if (i == n_ref-1)
			{
			u1_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u1");
			u2_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u2");
			}

			//=========================================================================
			Output::print("<<<<< ParMooN Finished: Thermohydraulic Brinkman2D Main Program >>>>>");

			timer.print_total_time("this Level (total): ");
			timer.reset();

		}

		Output::close_file();



		cout <<" Wrote out the solution obtained by the Brinkman routine --> into /Home/flow/blank/PARMOON/Tests/Thermohydraulic_Brinkman2D/Brinkman_ui "<<endl;

	}


	// ############################################################################################################### //
	// ##############################################    CD2D    ##################################################### //

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

	// refine grid
	//size_t n_ref_cd2d = domain.get_n_initial_refinement_steps();
	for(size_t i = 0; i < n_ref; i++)
	{
		domain.RegRefineAll();
	}
	// write grid into an Postscript file
	if(parmoon_db_cd2d["output_write_ps"])
		domain.PS("Domain.ps", It_Finest, 0);


	//=========================================================================
	CD2D cd2d(domain, parmoon_db_cd2d);

	/* ********************************************************************* */
	// use an analytic coefficient function (defined on top of this file)
	// In the .dat-file, permeability has to be set equal to -1
	if (parmoon_db_cd2d["coefficient_function_type"].is(3))
	{
		Output::print("A spatially varying coefficient function is detected and used in the assemble routine.");
		// create a FEFunction which will serve as the coefficient_function:
		// OLD 21.06.18: auto coll = cd2d.get_space().GetCollection();
		// fe space of piecewise constant functions
		TFESpace2D convection_function1_FEspace(coll, "convection_function1_FEspace", "s",
				BoundConditionNoBoundCondition, 0, nullptr);
		TFESpace2D convection_function2_FEspace(coll, "convection_function2_FEspace", "s",
				BoundConditionNoBoundCondition, 0, nullptr);
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


		convection_function1.WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "written_convection_function1");
		convection_function2.WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "written_convection_function2");

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
	else if (parmoon_db_cd2d["coefficient_function_type"].is(4))
	{
		Output::print("Two spatially varying convection functions from prior Brinkman problem will be used in the assemble routine.");
		// create a FEFunction which will serve as the coefficient_function:
		// OLD 21.06.18:  auto coll = cd2d.get_space().GetCollection();
		// fe space of piecewise constant functions
		TFESpace2D convection_function1_FEspace(coll, "convection_function1_FEspace", "s",
				BoundConditionNoBoundCondition, 0, nullptr);
		TFESpace2D convection_function2_FEspace(coll, "convection_function2_FEspace", "s",
				BoundConditionNoBoundCondition, 0, nullptr);
		BlockVector convection_function1_vector(convection_function1_FEspace.GetN_DegreesOfFreedom());
		BlockVector convection_function2_vector(convection_function2_FEspace.GetN_DegreesOfFreedom());
		TFEFunction2D convection_function1(
				&convection_function1_FEspace, "convection_function1", "convection_function1",
				&convection_function1_vector.at(0), convection_function1_FEspace.GetN_DegreesOfFreedom());
		TFEFunction2D convection_function2(
				&convection_function2_FEspace, "convection_function2", "convection_function2",
				&convection_function2_vector.at(0), convection_function2_FEspace.GetN_DegreesOfFreedom());

		u1_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u1_incd2d");
		u2_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u2_incd2d");


		convection_function1.Interpolate(u1_ptr);
		convection_function2.Interpolate(u2_ptr);

		////////// convection_function1.WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "written_convection_function1");
		////////// convection_function2.WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "written_convection_function2");

		//NEW LB Debug test:
		// double tmp_values[3];
		// coefficient_function1.FindGradient(0.25,0.25,tmp_values);
		// Output::print(" value = ",tmp_values[0]);

	/*	u1_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u1_incd2d");
		u2_ptr->WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "u2_incd2d");
*/
		TFEFunction2D *convection_function1_ptr, *convection_function2_ptr;
		convection_function1_ptr = &convection_function1;
		convection_function2_ptr = &convection_function2;

		cd2d.assemble(convection_function1_ptr, convection_function2_ptr);
		cout <<" *****coefficient_function_type 4 detected **** "<<endl;
	}
	/* ********************************************************************* */
	// use an analytic coefficient function (defined on top of this file)
	// In the .dat-file, permeability has to be set equal to -1
	else if (parmoon_db_cd2d["coefficient_function_type"].is(1))
	{
		Output::print("A spatially varying coefficient function is detected and used in the assemble routine.");
		// create a FEFunction which will serve as the coefficient_function:
		// OLD 21.06.18: auto coll = cd2d.get_space().GetCollection();
		// fe space of piecewise constant functions
		TFESpace2D convection_function1_FEspace(coll, "convection_function1_FEspace", "s",
				BoundConditionNoBoundCondition, 0, nullptr);
		BlockVector convection_function1_vector(convection_function1_FEspace.GetN_DegreesOfFreedom());
		TFEFunction2D convection_function1(
				&convection_function1_FEspace, "convection_function1", "convection_function1",
				&convection_function1_vector.at(0), convection_function1_FEspace.GetN_DegreesOfFreedom());
		convection_function1.Interpolate(analytic_convection_function1);
		//coefficient_function_vector.print("coefficient_function");


		convection_function1.WriteSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D", "written_convection_function");

		//NEW LB Debug test:
		// double tmp_values[3];
		// coefficient_function.FindGradient(0.25,0.25,tmp_values);
		// Output::print(" value = ",tmp_values[0]);


		TFEFunction2D *convection_function1_ptr;
		convection_function1_ptr = &convection_function1;
		cd2d.assemble(convection_function1_ptr);

	}
	/* ********************************************************************* */
	// use an external file containing information on, e.g., the permeability field in the format of an FEFunction2D (see ReadSol() and WriteSol())
	// in combination with an appropriate mesh-file as input.
	// Note that the input is then geo_file: ....mesh and in 2D a boundary_file: ...PRM
	// In the .dat-file, permeability has to be set equal to -1
	else if (parmoon_db_cd2d["coefficient_function_type"].is(2))
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
				BoundConditionNoBoundCondition, 0, nullptr);
		TFESpace2D coefficient_function_u2_FEspace(coll, "coefficient_function_u2_FEspace", "s",
				BoundConditionNoBoundCondition, 0, nullptr);
		BlockVector coefficient_function_u1_vector(coefficient_function_u1_FEspace.GetN_DegreesOfFreedom());
		BlockVector coefficient_function_u2_vector(coefficient_function_u2_FEspace.GetN_DegreesOfFreedom());

		TFEFunction2D coefficient_function_u1(&coefficient_function_u1_FEspace, "coefficient_function_u1", "coefficient_function_u1",
				&coefficient_function_u1_vector.at(0), coefficient_function_u1_FEspace.GetN_DegreesOfFreedom());
		TFEFunction2D coefficient_function_u2(&coefficient_function_u2_FEspace, "coefficient_function_u2", "coefficient_function_u2",
				&coefficient_function_u2_vector.at(0), coefficient_function_u2_FEspace.GetN_DegreesOfFreedom());

		//coefficient_function.ReadSol("/Home/flow/blank/PARMOON/Tests/Brinkman2D/Poiseuille_Hannukainen/Tria_mesh/Varying_Coefficient_Function/written_coefficient_function0.Sol");
		//// coefficient_function.ReadSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D/Brinkman_ux0.Sol");

		coefficient_function_u1.ReadSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D/Brinkman_ux0.Sol");
		TFEFunction2D *coefficient_function_u1_ptr = &coefficient_function_u1;


		coefficient_function_u2.ReadSol("/Users/blank/ParMooN/Tests/Thermohydraulic_Brinkman2D/Brinkman_uy1.Sol");
		TFEFunction2D *coefficient_function_u2_ptr = &coefficient_function_u2;

		cd2d.assemble(coefficient_function_u1_ptr, coefficient_function_u2_ptr);
	}

	// just hand in a constant convection
	else
	{
		cd2d.assemble( );
	}

	cd2d.solve();
	cd2d.output();
	//=========================================================================
	//db.info(false);
	Output::close_file();
	return 0;



}

