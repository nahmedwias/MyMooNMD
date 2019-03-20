#include "Domain.h"
#include "Database.h"
#include <ParMooN_repository_info.h>
#ifdef __3D__
#include "FEDatabase3D.h"
using OldDatabase = TFEDatabase3D;
#else
#include "FEDatabase2D.h"
using OldDatabase = TFEDatabase2D;
#endif

const std::string path_to_repo = parmoon::source_directory;
const std::string path_to_meshes = path_to_repo + "/data/mesh/";


struct TestObject
{
    std::string Geo_file;
    std::string Boundary_file;
};

TestObject prepare_test0()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "unit_square/unit_square_quad6.mesh";
  std::string bd_file = path_to_meshes + "unit_square/unit_square.PRM";
#else
  std::string geo_file = "Default_UnitCube_Tetra";
  std::string bd_file = "Default_UnitCube";
#endif
  return TestObject{geo_file, bd_file};
}
TestObject prepare_test1()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "unit_square/unit_square_tria6.mesh";
  std::string bd_file = path_to_meshes + "unit_square/unit_square.PRM";
#else
  std::string geo_file = "Default_UnitCube_Hexa";
  std::string bd_file = "Default_UnitCube";
#endif
  return TestObject{geo_file, bd_file};
}
TestObject prepare_test2()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "Hemker_tria.mesh";
  std::string bd_file = path_to_meshes + "Hemker.PRM";
#else
  std::string geo_file = path_to_meshes + "cylinder.3d.3K.mesh";
  std::string bd_file = "";
#endif
  return TestObject{geo_file, bd_file};
}
TestObject prepare_test3()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "Hemker_quad.mesh";
  std::string bd_file = path_to_meshes + "Hemker.PRM";
#else
  std::string geo_file = path_to_meshes + "channel.3d.mesh";
  std::string bd_file = "";
#endif
  return TestObject{geo_file, bd_file};
}


void check_collection(TCollection& coll)
{
  auto n_cells = coll.GetN_Cells();
  Output::print("check collection with ", n_cells, " cells");
  for(int i = 0; i < n_cells; ++i)
  {
    auto cell = coll.GetCell(i);
    if(i != coll.get_cell_index(cell))
    {
      ErrThrow("wrong cell index, ", i, " ", coll.get_cell_index(cell));
    }
  }
}

int main()
{
  // Construct the ParMooN Databases.
  TDatabase Database;
  OldDatabase FEDatabase;
  std::vector<TestObject> TestObjects = {prepare_test0(), prepare_test1(),
                                         prepare_test2(), prepare_test3()};
  for(auto &Test : TestObjects)
  {
    Output::print("\ntesting geo_file, ", Test.Geo_file);

    // Construct the ParMooN Databases.
    ParameterDatabase parmoon_db =
        ParameterDatabase::parmoon_default_database();
    parmoon_db.merge(TDomain::default_domain_parameters());
    parmoon_db["boundary_file"].set(Test.Boundary_file, false);
    parmoon_db["geo_file"].set(Test.Geo_file, false);
    
    // Construct domain, thereby read in controls from the input file.
    TDomain domain(parmoon_db);
    
    domain.RegRefineAll();

    auto coll = domain.GetCollection(It_Finest, 0);
    check_collection(*coll);
    delete coll;
  }
}
