/*!
 * @brief Main program for a system of coupled time dependent convection-
 * diffusion reaction equations, which are influenced by a time-dependent flow
 * and coupled to a particle distribution.
 *
 * This will be the main program for the first simulation for my PhD Thesis
 * and shall be used to simulate a flow crystallizer as reported by Eder et al.
 * 2010.
 *
 * We aim to use
 *    - an analytical or precomputed stationary flow field (e.g. Hagen--Poiseuille)
 * and within a time loop
 *    - the module for coupled Time_CD2D equations (couple temperature and one species)
 *    - the Brush interface module for stochastic calculation of the particles
 *
 * @author Clemens Bartsch
 * @date April 22, 2016
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <TimeDiscRout.h>

//for stationary flow
#include <NSE2D_axisymmetric.h>
#include <Example_NSE2D.h>

//for the time cd2d object
#include <Coupled_Time_CDR_2D.h>
#include <Example_TimeCD2D.h>

//for the particles object
#include <BrushWrapper.h>

//Timer fun
#include <Chrono.h>

//CB DEBUG
#include <GridTransferTool.h>
void DirichletBoundaryConditions_TEST(int BdComp, double t, BoundCond &cond)
{
      cond = DIRICHLET;
}
//END DEBUG

int main(int argc, char* argv[])
{

  bool axisymmetric = true;

  TDatabase Database;
  TFEDatabase2D FEDatabase;


  // read in the database and distribute them
  std::ifstream fs(argv[1]);
  ParameterDatabase general_database("General Database");
  general_database.read(fs);
  fs.close();
  ParameterDatabase flow_database = general_database.get_nested_database("Flow Database");
  ParameterDatabase conc_database = general_database.get_nested_database("Conc Database");
  ParameterDatabase particle_database = general_database.get_nested_database("Particle Database");

  // put some needed parameters into the respective DBs
  flow_database.merge(general_database, true);
  conc_database.merge(general_database, true);
  particle_database.merge(general_database, true);

  Output::print(">>>>>>", "Starting ParMooN Program: Flow, Concentration, Particles (2D).");
  Chrono total_chrono;

  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(general_database, argv[1]);

  //open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(general_database["outfile"]);
  Output::setVerbosity(general_database["verbosity"]);


  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

  std::string brush_grid_control = particle_database["brush_grid_control"];

  //===========================================================================
  // do initial refinements of the domain
  size_t n_ref = domain.get_n_initial_refinement_steps();
  TCollection* brush_grid;
  for(size_t i=0; i<n_ref; i++)
  {
    if(brush_grid_control == "coarser" && i == n_ref - 1)
    {
      Output::info("GRID", "Picking one level coarser grid for Brush.");
      brush_grid = domain.GetCollection(It_Finest, 0); //brush collection: one level less
    }
    domain.RegRefineAll();
  }
  if(brush_grid_control == "same")
  {
   Output::info("GRID", "Picking same grid for Brush.");
   brush_grid = domain.GetCollection(It_Finest, 0);
  }

  // PART: PRECOMPUTE STATIONARY FLOW FIELD //////////////////////////////////
  Output::info("PROGRAM PART", "Precomputing velocity.");
  Example_NSE2D example_flow(flow_database);
  NSE2D_axisymmetric flow_object({domain.GetCollection(It_Finest, 0)}, flow_database, example_flow);

  flow_object.assemble_linear_terms();
  flow_object.stop_iteration(0);
  //STOKES CASE
  flow_object.solve();
  flow_object.output();

  TFESpace2D new_space(brush_grid,(char*)"new_space", (char*)"Space to "
                 "project the velocity on.",
                 DirichletBoundaryConditions_TEST, 0,
                 nullptr);
  GridTransferTool to_new_space_tool(&new_space, GridTransferType::MultiGrid);
  std::vector<double> new_function_vals(new_space.GetN_DegreesOfFreedom(), 0.0);
  TFEFunction2D new_function(&new_space,(char*)"new_function",
                    (char*) "Velocity projected to Q0.",
                    &new_function_vals.at(0), new_function_vals.size());

  to_new_space_tool.transfer(flow_object.get_axial_velocity(), new_function, new_function_vals);


//  //NAVIER--STOKES CASE
//  // nonlinear loop
//  // in function 'stopIt' termination condition is checked
//  for(unsigned int k = 1;; k++)
//  {
//    flow_object.solve();
//
//    //flow_object.assemble_nonlinear_term();
//
//    if(flow_object.stop_iteration(k))
//      break;
//  } // end for k
//
//  flow_object.output();



  // PART: SET UP PARTICLE- AND CDRE OBJECT //////////////////////////////////
  Output::info("PROGRAM PART", "Constructing particles- and cdre object.");

  Example_TimeCoupledCDR2D example_conc(conc_database);
  Coupled_Time_CDR_2D conc_object(domain, conc_database, example_conc);

  if(brush_grid_control == "finer")
  {//one additional refinement for Brush
	  domain.RegRefineAll();
	  Output::info("GRID", "Picking one level finer grid for Brush.");
	  brush_grid = domain.GetCollection(It_Finest, 0);
  }
  // the particles object which wraps up Brush
  BrushWrapper part_object(brush_grid, domain.GetCollection(It_Finest, 0), particle_database, axisymmetric);

  // PARTS: SET UP INITIAL STATES /////////////////////////////////////////////
  Output::info("PROGRAM PART", "Setting up initial states.");

  //get the velocity field of the precomputed velo object
  //TFEFunction2D& uz = flow_object.get_axial_velocity();
  TFEFunction2D& uz = new_function;
  TFEFunction2D& ur = flow_object.get_radial_velocity();

  const TFEFunction2D& pressure = flow_object.get_pressure();

  // concentrations: assemble matrices and right hand side at start time
  conc_object.assemble_initial_time(&uz, &ur);
  conc_object.output();

  //set fluid phase in the particle object
  std::vector<TFEFunction2D*> fcts = conc_object.get_fe_functions();
  part_object.reset_fluid_phase(uz,ur, pressure, fcts);
  part_object.output(TDatabase::TimeDB->CURRENTTIME);

  // PART: SOLVE THE SYSTEM IN A TIME LOOP ///////////////////////////////////
  Output::info("PROGRAM PART", "Solving the coupled system.");

  double end_time = TDatabase::TimeDB->ENDTIME;
  int step = 0;

  int output_steps_parts = particle_database["output_all_k_steps"];
  int output_steps_concs = conc_database["output_all_k_steps"];

  // time iteration
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    step++;

    SetTimeDiscParameters(1);

    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    TDatabase::TimeDB->CURRENTTIME += tau;

    Output::print<1>("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

    //update and solve particles
    std::vector<TFEFunction2D*> fcts = conc_object.get_fe_functions();
    part_object.reset_fluid_phase(uz,ur, pressure, fcts);
    part_object.solve(TDatabase::TimeDB->CURRENTTIME, TDatabase::TimeDB->CURRENTTIME + tau);

    //solve cdr system
    conc_object.assemble_uncoupled_part(&uz, &ur, part_object.sources_and_sinks());
    conc_object.couple_and_solve();

    if(step %  output_steps_parts == 0)
      part_object.output(TDatabase::TimeDB->CURRENTTIME);
    if(step % output_steps_concs == 0)
      conc_object.output();
  }

  total_chrono.print_total_time("The whole shebang took: ");

  Output::print(">>>>>> Program completed <<<<<<");

  Output::close_file();

  //clean up the collections required for the particle object
  delete brush_grid;

  return 0;
}