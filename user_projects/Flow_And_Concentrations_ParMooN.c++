/*!
 * @brief Main program for a system of coupled time dependent convection-
 * diffusion reaction equations, which are influenced by a time-dependent flow.
 *
 * This main program starts is life as one which adminsters a one-way coupling
 * between a (simple) Time_NSE2D obejct (flow) and a Time_CD2D object
 * (concentration). Later the Time_CD2D object will be replaced by an
 * instance of the coupled module.
 *
 * TODO Might be this will need some big changes as soon as the default branch
 * is merged here.
 *
 * @author Clemens Bartsch
 * @date February 3, 2016
 */

int main(int argc, char* argv[])
{
  // Put up domain and databases.
  // These are the usual ParMooN-initialisation steps.
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  // set variables' value in TDatabase using argv[1] (*.dat file) */
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

  return 0;
}
