/*************************************************************************
 *
 * A ParMooN main program. Solves a CD2D problem on meshes succesively
 * gained by applying Franco Dassi's mesh adaptation algorithm.
 *
 *************************************************************************/
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <CD2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <Example_CD2D.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(argv[1]);

  //set PROBLEM_TYPE to CD if not yet set
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 1;
  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);

  // Creation of an initial grid - ideally this is done by SurgGen:
  // then ParMooN does not ever have to communicate a grid to SurfGen
  // (which it is utterly unfit for)

  /* include the mesh from a mesh generator, for a standard mesh use the
   * build-in function. The GEOFILE describes the boundary of the domain. */
  domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); // call mesh generator

  // refine grid up to the coarsest level
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS+
    TDatabase::ParamDB->LEVELS; i++)
    domain.RegRefineAll();

  // write grid into a postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    domain.PS("Domain.ps", It_Finest, 0);

  // end creation of initial mesh


 // The solving procedure
  //=========================================================================
  for (size_t step = 0; step < 3 ; ++step) // replace by useful stopping criterion
  {
    //build and solve CD2D problem
    CD2D cd2d(domain);
    cd2d.assemble();
    cd2d.solve();

    // // temporary output of the grid after each solution step
    // Output::print("Output after step ", step);
    // cd2d.output();


    // check some stopping criterion

    // 1) vectors with function values and gradient recoveries at the mesh vertices
    // use:
    //    const TFEFunction& solution_function = cd2d.get_function();
    //    for index i ueber alle gitterpunkte (aus SurfGen)
    //      double x = x-Wert des Gitterpunkts
    //      double y = y-Wert des Giterpunkts
    //      double[3] values;
    //      solution_function.FindGradient(x,y,values);
    //      u[i] = values[0];
    //      ux[i] = values[1];
    //      uy[i] = values[2];
    //    endfor

    // 2) feed surfGen with u, ux, uy and call it

    // 3) re-initialize the domain by calling MODIFIED (TODO,
    //    depends on exact output of Francos code) domain.MakeGrid()


  }
  //=========================================================================
  Output::close_file();
  return 0;
} // end main
