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
 * Tasks:
 *  - interface Brush
 *
 * @author Clemens Bartsch
 * @date April 22, 2016
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <TimeDiscRout.h>

//for the time nse 2d flow object
#include <Time_NSE2D.h>
#include <NSE2D.h> //for stationary flow
#include <Example_NSE2D.h>

//for the time cd2d object
#include <Coupled_Time_CDR_2D.h>
#include <Example_TimeCD2D.h>

//for the particles object
#include <BrushWrapper.h>


int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  // read all input parameters into four different databases
  ParameterDatabase general_database("General Database");
  ParameterDatabase flow_database("Flow Database");
  ParameterDatabase conc_database("Conc Database");
  ParameterDatabase particle_database("Particle Database");

  std::ifstream fs(argv[1]);
  general_database.read(fs);
  flow_database.read(fs);
  conc_database.read(fs);
  particle_database.read(fs);
  fs.close();

  flow_database.merge(general_database, true);
  conc_database.merge(general_database, true);

  Output::print(">>>>>>", "Starting ParMooN Program: Flow, Concentration, Particles (2D).");

  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(argv[1], general_database);

  //open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(general_database["outfile"]);
  Output::setVerbosity(general_database["verbosity"]);


  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

  /* include the mesh from a mesh generator, for a standard mesh use the
   * build-in function. The GEOFILE describes the boundary of the domain. */
  domain.InitFromMesh(general_database["boundary_file"], general_database["mesh_file"]);

  //===========================================================================
  // do initial refinements of the domain
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(size_t i=0; i<n_ref; i++)
    domain.RegRefineAll();

  // PART 2: SET UP SPECIFIC PROBLEM OBJECTS ///////////////////////////////////

  //set global parameters which are to be used for NSE construction here...
  TDatabase::ParamDB->VELOCITY_SPACE = 2;
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;

  Example_NSE2D example_flow(flow_database);
  NSE2D flow_object(domain, flow_database, example_flow);

  // PART 3: PRECOMPUTE STATIONARY FLOW FIELD //////////////////////////////////
  Output::info("PROGRAM PART", "Precomputing velocity.");
  flow_object.assemble();
  flow_object.stopIt(0);

  // nonlinear loop
  // in function 'stopIt' termination condition is checked
  for(unsigned int k = 1;; k++)
  {
    flow_object.solve();

    flow_object.assemble_nonlinear_term();

    if(flow_object.stopIt(k))
      break;
  } // end for k

  flow_object.output();


  // PART 4: SET UP PARTICLE- AND CDRE OBJECT //////////////////////////////////
  Output::info("PROGRAM PART", "Constructing particles- and cdre object.");

  //set global parameters which are to be used for Time_CD2D construction here...
  TDatabase::ParamDB->ANSATZ_ORDER = 1;
  Example_TimeCoupledCDR2D example_conc(conc_database);
  Coupled_Time_CDR_2D conc_object(domain, conc_database, example_conc);

  // the particles object which wraps up Brush
  TCollection* coll = domain.GetCollection(It_Finest, 0); //take finest grid collection for now
  BrushWrapper part_object(coll, particle_database);

  // PART 5: SET UP INITIAL STATES /////////////////////////////////////////////
  Output::info("PROGRAM PART", "Setting up initial states.");

  //get the velocity field of the precomputed velo object
  const TFEVectFunct2D& velo_field = flow_object.get_velocity();
  const TFEFunction2D& pressure = flow_object.get_pressure();

  // concentrations: assemble matrices and right hand side at start time
  conc_object.assemble_initial_time(&velo_field);
  conc_object.output();

  //set fluid phase in the particle object
  std::vector<const TFEFunction2D*> fcts = conc_object.get_fe_functions();
  part_object.reset_fluid_phase(velo_field, pressure, fcts);
  part_object.output(TDatabase::TimeDB->CURRENTTIME);

  // PART 6: SOLVE THE SYSTEM IN A TIME LOOP ///////////////////////////////////
  Output::info("PROGRAM PART", "Solving the coupled system.");

  double end_time = TDatabase::TimeDB->ENDTIME;
  int step = 0;
  int n_substeps = GetN_SubSteps();
  // time iteration
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    step++;

    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    for(int j=0; j < n_substeps; ++j)
    {
      SetTimeDiscParameters(1);

      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;

      Output::print<1>("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      //solve cdr system
      conc_object.assemble_uncoupled_part(&velo_field);
      conc_object.couple_and_solve(&velo_field); //TODO sources/sinks due to particles must enter here!
      conc_object.output();

      //update and solve particles
      std::vector<const TFEFunction2D*> fcts = conc_object.get_fe_functions();
      part_object.reset_fluid_phase(velo_field, pressure, fcts);
      part_object.solve(TDatabase::TimeDB->CURRENTTIME, TDatabase::TimeDB->CURRENTTIME + tau);
      part_object.output(TDatabase::TimeDB->CURRENTTIME);
    }
  }


  Output::close_file();

  delete coll; //the collection which was required for the particle object

  return 0;
}
