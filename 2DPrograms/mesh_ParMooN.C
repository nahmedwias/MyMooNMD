// =======================================================================
//
// Purpose:  main program for test the different functionalities
//           of mesh-conversion routines
//
// Author:   Alfonso Caiazzo
//
// History:  Implementation started on 22.03.2016
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Mesh.h>

#include <sys/stat.h>
#include <sys/types.h>

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

  //open OUTFILE, this is where all output is written to (addionally to console)
  Output::set_outfile(TDatabase::ParamDB->OUTFILE);
 
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
  
  /* initialization using GEO file (grid) and PRM file (boundary description) 
   */
  domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); // call mesh generator
   
  // refine grid up to the coarsest level
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS+
    TDatabase::ParamDB->LEVELS; i++)
    domain.RegRefineAll();
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  //=========================================================================
  // write the current mesh (using the finest collection)
  TCollection *coll = domain.GetCollection(It_Finest, 0);
  coll->writeMesh("test.mesh");
  //=========================================================================
  
  //=========================================================================
  // read the .mesh from file and initialize a Mesh object
  Mesh m("test.mesh");
  // display mesh properties
  m.info();
  m.writeToMesh("a.mesh");
  //=========================================================================

  
  //=========================================================================
  // read PRM file into a Boundary class
  m.setBoundary(TDatabase::ParamDB->BNDFILE);
  m.boundary.info();
  //=========================================================================

  //=========================================================================
  //write to (x)GEO
  m.writeToGEO("a.xGEO");
  //=========================================================================
  
  //=========================================================================
  Output::close_file();
  return 0;
} // end main
