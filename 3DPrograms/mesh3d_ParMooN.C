// =======================================================================
//
// Purpose:  main program for testing the Mesh class (3D) and the 
//           initialization of Domain class using different format
//           (mainly: .mesh, .smesh through TetGen)
//
// Author:   Alfonso Caiazzo
//
// History:  Implemetation started on 12.08.2016
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
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
  TFEDatabase3D FEDatabase;

  // read input file 
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  // merge with default tetgen parameters
  parmoon_db.merge(ParameterDatabase::default_tetgen_database(), true);
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  TDomain domain(parmoon_db, argv[1]);

  //Output::print(" ... writing 3d mesh ... ");
  //TCollection *coll = domain.GetCollection(It_Finest, 0);
  //coll->writeMesh("test.mesh");
  

  
  //=========================================================================
  Output::close_file();
  return 0;
} // end main
