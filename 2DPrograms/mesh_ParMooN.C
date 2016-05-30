// =======================================================================
//
// Purpose:  main program for testing the Mesh class and the 
//           mesh-conversion routines
//
// Author:   Alfonso Caiazzo
//
// History:  Implementation started on 22.03.2016
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Mesh.h>
#include <ParameterDatabase.h>

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

  // read input file (new version)
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(argv[1], parmoon_db);

  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // write all Parameters to the OUTFILE (not to console) for later reference
  //Database.WriteParamDB(argv[0]);  
  // create output directory, if not already existing
  //if(TDatabase::ParamDB->WRITE_VTK)
  //  mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

  unsigned int testType = 3;
  if (testType==0) {
    Output::print(" Test: ");
    Output::print("   * read .PRM and .GEO file and write the corresponding .mesh");

    // initialize domain
    domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"]);

    // write .mesh file
    TCollection *coll = domain.GetCollection(It_Finest, 0);
    coll->writeMesh(parmoon_db["mesh_file"]);

    // create mesh object and display info
    Mesh m(parmoon_db["mesh_file"]);
    m.info();
    
  } else if (testType==1) {
    //=========================================================================
    Output::print(" Test: ");
    Output::print("   * read .PRM and .GEO file and write the corresponding .mesh");
    Output::print("   * create a Mesh and display info");
    Output::print("   * rewirte a new .GEO file");

    // initialize domain
    domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"]);

    // write .mesh file
    TCollection *coll = domain.GetCollection(It_Finest, 0);
    coll->writeMesh("output/meshFromPRMandGEO.mesh");

    // create mesh object and display info
    Mesh m("output/meshFromPRMandGEO.mesh");
    m.info();
    
    // set boundary description
    m.setBoundary(parmoon_db["boundary_file"]);
    // display boundary info
    m.boundary.info();
    // write new .GEO
    m.writeToGEO("output/GEOFromMesh.xGEO");

    Output::print(" ... done. Compare the new and the old .GEO files. ");
    
    } else if (testType==2) {

    Output::print(" Test: ");
    Output::print("   * read .PRM and the new .GEO file and write the corresponding .mesh");
    domain.Init(parmoon_db["boundary_file"], "output/GEOFromMesh.xGEO");
    TCollection *coll = domain.GetCollection(It_Finest, 0);
    coll->writeMesh("output/meshFromMesh.mesh");

    Output::print(" ... done. Compare the mesh file with the one generated with testType = 1. ");

     } else if (testType==3) {

    Output::print(" Test: ");
    Output::print("   * read .PRM and the .mesh to initialize the domain");
    Output::print("   * write out the mesh file corresponding to the grid");
    domain.InitFromMesh(parmoon_db["boundary_file"],
		"output/meshFromMesh.mesh"); 
    TCollection *coll = domain.GetCollection(It_Finest, 0);
    coll->writeMesh("output/meshFromPRMandMesh.mesh");

    Output::print(" ... done. Compare now meshFromPRMandGEO.mesh and meshFromPRMandMesh.mesh");

  }
  
  
  //=========================================================================
  Output::close_file();
  return 0;
} // end main