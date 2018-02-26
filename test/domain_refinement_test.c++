#include "Domain.h"
#include "MooNMD_Io.h"
#ifdef __2D__
  #include "FEDatabase2D.h"
  #include "NSE2D.h"
#else
  #include "FEDatabase3D.h"
  #include "NSE3D.h"
#endif
#include "Database.h"


int main(int , char** )
{
  //  declaration of databases
  TDatabase Database;
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.add("refinement_n_initial_steps", (size_t) 1, "");
  Output::setVerbosity(3);
  
#ifdef __2D__
  TFEDatabase2D FEDatabase;
  db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "TwoTriangles", "", {"UnitSquare", "TwoTriangles"});
#else // 3D
  TFEDatabase3D FEDatabase;
  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Tetra", "", 
         {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
#endif // 2D
  db.add("output_write_vtk", true, "");
  
  TDomain domain(db);
  domain.barycentric_refinement();
  
#ifdef __2D__
  NSE2D nse(domain, db);
  Output::print("created nse2d object");
  nse.assemble();
#else // 3D
  auto gridCollections = domain.refine_and_get_hierarchy_of_collections(db
  #ifdef _MPI
      , maxSubDomainPerDof
  #endif
      );

  // Choose and construct example.
  Example_NSE3D example_obj(db);

  // Construct the nse3d problem object.
#ifdef _MPI
  NSE3D nse(gridCollections, db, example_obj, maxSubDomainPerDof);
#else // no mpi
  NSE3D nse(gridCollections, db, example_obj);
#endif
  nse.assemble_linear_terms();
#endif
  
  nse.solve();
  nse.output();
}

