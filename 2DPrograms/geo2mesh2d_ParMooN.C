// =======================================================================
//
// Purpose:  main program for converting .(x)GEO files into .mesh
//
// Author:   Alfonso Caiazzo
//
// History:  Implementation started on 10.05.2016
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
int main(int, char* argv[])
{
    
  //  declaration of database, you need this in every program
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;

  // read input file (new version)
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

  TDomain domain(parmoon_db);

  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // write .mesh file
  TCollection *coll = domain.GetCollection(It_Finest, 0);
  coll->writeMesh(parmoon_db["mesh_file"]);
  
  // create mesh object and display info for checking
  Mesh m(parmoon_db["mesh_file"]);
  m.info();
  // set boundary description
  m.setBoundary(parmoon_db["boundary_file"]);
  // display boundary info
  m.boundary.info();
  
  
  //=========================================================================
  Output::close_file();
  return 0;
} // end main
