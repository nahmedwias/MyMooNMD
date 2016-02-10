/*!
 * @brief Main program for a system of coupled time dependent convection-
 * diffusion reaction equations, which are influenced by a time-dependent flow.
 *
 * This main program starts is life as one which adminsters a one-way coupling
 * between a (simple) Time_NSE2D obejct (flow) and a Time_CD2D object
 * (concentration). Later the Time_CD2D object will be replaced by an
 * instance of the coupled module.
 *
 * @author Clemens Bartsch
 * @date February 3, 2016
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Output2D.h>

//for the time nse 2d flow object
#include <Time_NSE2D.h>
#include <Example_NSE2D.h>

//for the time cd2d object
#include <Time_CD2D.h>
#include <Example_CD2D.h>


int main(int argc, char* argv[])
{
  // TODO it might be a problem that we have no global "problem type" defined
  // - on the other hand this is good, because that urges us to find out
  // whether we can remove that global parameter entirely

  // //////////////////////////////////////////////////////////////////////////
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

  //===========================================================================
  /* include the mesh from a mesh generator, for a standard mesh use the
   * build-in function. The GEOFILE describes the boundary of the domain. */
  Domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);

  //===========================================================================
  // do initial refinements of the domain
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    Domain.RegRefineAll();
  // //////////////////////////////////////////////////////////////////////////

  // TODO Here is where the work begins.

  Output::close_file();

  return 0;
}
