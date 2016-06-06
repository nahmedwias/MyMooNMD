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
 *  - write example file for Eder example
 *  - enable flow and psd to get into the coupling module
 *  - interface Brush
 *
 * Questions:
 *  - which part of the program will take care of particle growth??
 *
 * @author Clemens Bartsch
 * @date April 22, 2016
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Output2D.h>
#include <TimeDiscRout.h>

//for the time nse 2d flow object
#include <Time_NSE2D.h>
#include <NSE2D.h> //for stationary flow
#include <Example_NSE2D.h>

//for the time cd2d object
#include <Time_CD2D.h>
#include <Example_CD2D.h>


//the interpolator
#include <FEFunctionInterpolator.h>


int main(int argc, char* argv[])
{
  // TODO it might be a problem that we have no global "problem type" defined
  // - on the other hand this is good, because that urges us to find out
  // whether we can remove that global parameter entirely

  // PART 1: USUAL PARMOON SETUP //////////////////////////////////////////////

  // Put up domain and databases.
  // These are the usual ParMooN-initialisation steps.
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  // set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain Domain(argv[1]);

  //===========================================================================
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);
  OutFile.setf(std::ios::scientific);
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

  //===========================================================================
  /* include the mesh from a mesh generator, for a standard mesh use the
   * build-in function. The GEOFILE describes the boundary of the domain. */
  Domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);

  //===========================================================================
  // do initial refinements of the domain
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    Domain.RegRefineAll();

  // PART 2: SET UP SPECIFIC PROBLEM OBJECTS ///////////////////////////////////

  // TODO Here is where the work begins.

  //set global parameters which are to be used for NSE construction here...
  TDatabase::ParamDB->VELOCITY_SPACE = 2;
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;

  TDatabase::ParamDB->EXAMPLE = 0; //global example variable for nse
  NSE2D flow_object(Domain, Example_NSE2D());

  //set global parameters which are to be used for Time_CD2D construction here...
  TDatabase::ParamDB->ANSATZ_ORDER = 1;

  TDatabase::ParamDB->EXAMPLE = 104; //global example variable for tcd2d
  Time_CD2D conc_object(Domain, Example_CD2D());

  // PART 3: PRECOMPUTE STATIONARY FLOW FIELD //////////////////////////////////

  flow_object.assemble();
  flow_object.stopIt(0);

  // nonlinear loop
  // in function 'stopIt' termination condition is checked
  for(unsigned int k = 1;; k++)
  {
    Output::print<1>("nonlinear iteration step ", setw(3), k-1, "\t",
                     flow_object.getResiduals());
    flow_object.solve();

    //no nonlinear iteration for Stokes problem
    if(TDatabase::ParamDB->PROBLEM_TYPE == 3)
      break;

    flow_object.assemble_nonlinear_term();

    if(flow_object.stopIt(k))
      break;
  } // end for k

    // Have a look at the precomputed velocity field - do we like it?
    TDatabase::ParamDB->WRITE_VTK = 1; //write a VTK
    TDatabase::ParamDB->MEASURE_ERRORS = 0; //error measuring is evil (i.e. unsafe methods, incorrect aux object...)
    flow_object.output();


  // PART 5: SOLVING THE CDRE IN A TIME LOOP //////////////////////////////////
  //change problem type...globally...
  TDatabase::ParamDB->PROBLEM_TYPE = 2;


  //change the outfile...globally...
  TDatabase::ParamDB->BASENAME = (char*)"conc";
  TDatabase::ParamDB->VTKBASENAME = (char*)"conc";

  //get the velocity field of the precomputed velo object
  const TFEVectFunct2D& velo_field = flow_object.get_velocity();

  // assemble matrices and right hand side at start time
  //conc_object.assemble_initial_time(&velo_field);
  conc_object.assemble_initial_time(&velo_field);

  double end_time = TDatabase::TimeDB->ENDTIME;
  int step = 0;
  int n_substeps = GetN_SubSteps();

  int image=0;

  conc_object.output(0,image);

  // time iteration
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    step++;
    // Output::print("mem before: ", GetMemory());
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    for(int j=0; j < n_substeps; ++j)
    {
      SetTimeDiscParameters(1);
      if(step==1)
      {
        Output::print<1>("Theta1: ", TDatabase::TimeDB->THETA1);
        Output::print<1>("Theta2: ", TDatabase::TimeDB->THETA2);
        Output::print<1>("Theta3: ", TDatabase::TimeDB->THETA3);
        Output::print<1>("Theta4: ", TDatabase::TimeDB->THETA4);
      }
      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;

      Output::print<1>("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      conc_object.assemble(&velo_field);

      conc_object.solve();

      conc_object.output(step, image);
    }
  }



  Output::close_file();
  return 0;
}
