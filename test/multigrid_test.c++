/**
 * @file Unit testing of methods used in the multigrid framework.
 *
 * Up to now the test is a joke - it only executes the prolongation/restriction
 * of constant functions and allows you to look at the results.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#include <GridTransfer.h>
#include <Domain.h>
#include <FEDatabase2D.h>
#include <FEDatabase3D.h>
#include <FESpace2D.h>
#include <FESpace3D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <MooNMD_Io.h>
#include <ParameterDatabase.h>

#include <vector>
#include <tuple>
#include <MooNMD_Io.h>

int main(int argc, char* argv[])
{

  //  declaration of databases
  TDatabase Database;
  TFEDatabase2D FEDatabase_2d;
#ifdef __3D__
  TFEDatabase3D FEDatabase_3d;
#endif

  ParameterDatabase db = ParameterDatabase::parmoon_default_database();

#ifdef __3D__
  db["boundary_file"] =std::string("Default_UnitCube");
  db["geo_file"] = std::string("Default_UnitCube_Tetra");
#endif

  db.add("refinement_n_initial_steps", (size_t) 1, "");
  db.add("n_multigrid_levels", (size_t) 2, "");
  Output::setVerbosity(db["verbosity"]);

  // default construct a domain object
  TDomain domain(db);
  // the domain is initialised with default description and default
  // initial mesh

  domain.Init(db["boundary_file"], db["geo_file"]);

  // refine grid up to the coarsest level
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(size_t i=0; i < n_ref; i++)
  {
    domain.RegRefineAll();
  }

  //Get two collections.
  TCollection* coll_coarse = domain.GetCollection(It_EQ, 0);
  TCollection* coll_fine = domain.GetCollection(It_EQ, 1);

  Output::print("Coarse cells: ", coll_coarse->GetN_Cells());
  Output::print("Fine cells: ", coll_fine->GetN_Cells());

  // Create FeSpaces to fiddle around with

  size_t order = 2;

#ifdef __2D__
  TFESpace2D space_coarse(coll_coarse, (char*)"space_coarse",
               (char*)"space_coarse",
               BoundConditionNSE, order, nullptr);

  TFESpace2D space_fine(coll_fine, (char*)"space_fine",
                          (char*)"space_fine",
                          BoundConditionNSE, order, nullptr);
#endif
#ifdef __3D__
  TFESpace3D space_coarse(coll_coarse, (char*)"space_coarse",
               (char*)"space_coarse",
               BoundCondition_FEM_FCT, order);

  TFESpace3D space_fine(coll_fine, (char*)"space_fine",
                          (char*)"space_fine",
                          BoundCondition_FEM_FCT, order);
#endif


  std::vector<double> function_coarse(space_coarse.GetN_DegreesOfFreedom(), 1.0);
  std::vector<double> function_fine(space_fine.GetN_DegreesOfFreedom(), 0.0);

  // Start the checks

  //Prolongation
  GridTransfer::Prolongate(
      space_coarse, space_fine,
      &function_coarse[0],  function_coarse.size(),
      &function_fine[0], function_fine.size());

  Output::print("Coarse function");
  for(auto c : function_coarse)
    Output::print(c);

  Output::print("Fine function");
  for(auto f : function_fine)
    Output::print(f);

  //Function restriction
  std::fill(function_coarse.begin(), function_coarse.end(), 0.0);
  std::fill(function_fine.begin(), function_fine.end(), 2.0);

  GridTransfer::RestrictFunction(
      space_coarse, space_fine,
      &function_coarse[0],  function_coarse.size(),
      &function_fine[0], function_fine.size());

  Output::print("Coarse function");
  for(auto c : function_coarse)
    Output::print(c);

  Output::print("Fine function");
  for(auto f : function_fine)
    Output::print(f);

  //Defect Restriction
  std::fill(function_coarse.begin(), function_coarse.end(), 0.0);
  std::fill(function_fine.begin(), function_fine.end(), 2.0);

  GridTransfer::DefectRestriction(
      space_coarse, space_fine,
      &function_coarse[0],  function_coarse.size(),
      &function_fine[0], function_fine.size());

  Output::print("Coarse function");
  for(auto c : function_coarse)
    Output::print(c);

  Output::print("Fine function");
  for(auto f : function_fine)
    Output::print(f);

  //CB DEBUG - Old defect restriction
  double aux[1000];

  std::fill(function_coarse.begin(), function_coarse.end(), 0.0);
  std::fill(function_fine.begin(), function_fine.end(), 2.0);

  DefectRestriction(
      &space_coarse, &space_fine,
      &function_coarse.at(0), &function_fine.at(0),
      aux);

  Output::print("Coarse function");
  for(auto c : function_coarse)
    Output::print(c);

  Output::print("Fine function");
  for(auto f : function_fine)
    Output::print(f);

  //END DEBUG

  Output::print("Test program finished.");
}


