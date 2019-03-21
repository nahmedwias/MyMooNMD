/**
 * @brief A test program for the solving of Algebraic Flux Correction (AFC) in CD3D problems.
 *
 * This serves as a test for the solving of AFC schemes for CD3D problems. It is intended to
 * perform CD3D calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes wrong)
 * the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then this tests
 * shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 *
 * @date 2019/03/13
 * @author Abhinav Jha
 *
 */
#include <AlgebraicFluxCorrection.h>
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <CD3D_AFC.h>
#include <Multigrid.h>

#include <LocalAssembling3D.h>
#include <MainUtilities.h> //for error measuring

// compare the computed errors in the CD3D object with the given ones in 
// the array
void compareErrors(const CD3D_AFC& cd3d, std::array<double, 2> errors)
{
  const double eps = 1e-11;
  std::array<double, int(2)> computed_errors;
  computed_errors = cd3d.get_errors();
  
  // check the errors
  if( fabs(computed_errors.at(0) - errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 error not correct. ",
             computed_errors.at(0) - errors[0]);
  }
  if( fabs(computed_errors.at(1) - errors[1]) > eps )
  {
    ErrThrow("Program 1: H1-semi error not correct. ",
             computed_errors.at(1) - errors[1]);
  }
}

// Here the actual computations take place
void check_cd3d(std::list<TCollection* > gridCollections, ParameterDatabase& db, std::array<double,2> errors)
{  
 
  // Choose and construct example.
  Example_CD3D example_obj(db);
  CD3D_AFC cd3d(gridCollections, db, example_obj);
  cd3d.assemble(0);
  cd3d.solve(0);
  //maximum number of iterations for non lienar loop
  unsigned int max_it=1000;
  for(unsigned int k=1;;k++)
  {
    bool converged;
    converged = cd3d.solve(k);
    if ((converged)||(k>= max_it))
      break;
  }
  cd3d.output();
  //comparison of errors
  compareErrors(cd3d, errors); 
}

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  Output::setVerbosity(2);
  //  declaration of databases
  TDatabase Database;
  TFEDatabase3D FEDatabase;
  
  Output::print("\ntesting with algebraic flux correction");
  
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  
  db.merge(Example3D::default_example_database());
  db["example"] = 3; //Sharp Boundary Layer Example part of Example Database
  
  //only direct solver for CD3D class
  db.add("solver_type", std::string("direct"), "");
  db.add("refinement_n_initial_steps", (size_t) 3,"");
  db.add("multigrid_n_levels", (size_t) 0, "");
  
  //addition for AFC database
  db.merge(AlgebraicFluxCorrection::default_afc_database(), true);
  db["algebraic_flux_correction"].set("afc");
  db["diffusion_coefficient"]=1.0e-1;
  
  // default construct a domain object
  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Tetra", "", {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
  TDomain domain(db);
  
  std::list<TCollection* > gridCollections = domain.refine_and_get_hierarchy_of_collections(db);
  
  //AFC only applicable to P1 elements
  TDatabase::ParamDB->ANSATZ_ORDER = 1; //P1 elements
  db["space_discretization_type"] = "galerkin"; //Galerkin Desicreitzation
  
  //maximum number of iterations for non linear loop
  db["afc_nonlinloop_maxit"]=1000;
 
  
  //Computation for P1 elements
  //P1 elements and Dynamic Damping WITHOUT Anderson Acceleration
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + fixed_point_rhs + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  std::array<double,2> errors = {{0.00149380803286, 0.013230165352384}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + newton + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
 
  errors = {{0.0014938084093367 , 0.013230166506521}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + fixed_point_matrix + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.0014938084639286, 0.013230166672609}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_rhs + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  errors = {{0.000416420573027, 0.010287280962814}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + newton + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  errors = {{0.00041642068192927, 0.010287281981352}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_matrix + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  errors = {{0.00041642061583269, 0.010287285359938}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  
  //P1 elements and Dynamic Damping WITH Anderson Acceleration
  //=========================================================================
  Output::print("\n\n ---------  P1 + kuzmin + fixed_point_rhs + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  db["afc_anderson_damping"].set<>("no");
  errors = {{0.00149380803286, 0.013230165352384}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + newton + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  db["afc_anderson_damping"].set<>("no");
  errors = {{0.0014938084093367, 0.013230166506521}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + fixed_point_matrix + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  db["afc_anderson_damping"].set<>("no");
  errors = {{0.0014938084639286, 0.013230166672609}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_rhs + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  db["afc_anderson_damping"].set<>("no");
  errors = {{0.00041642065861425, 0.010287280925369}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + newton + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  db["afc_anderson_damping"].set<>("no");
  errors = {{0.00041642072608981, 0.010287281031888}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_matrix + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  db["afc_anderson_damping"].set<>("no");
  errors = {{0.0004164208055072, 0.010287281557243}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================

  //Computations for Q1 elements
  //Note: BJK17 limiter not applicable for Q1 elements.
  
  db["geo_file"].set<>("Default_UnitCube_Hexa");
  TDomain domain_quad(db);
  gridCollections = domain_quad.refine_and_get_hierarchy_of_collections(db);
  
  //Q1 elements and Dynamic Damping WITHOUT Anderson Acceleration
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_rhs + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.00039916511745917, 0.0053003916489301}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + newton + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.000399165125613, 0.0053003919795139}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_matrix + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.0003991651006792, 0.0053003919707961}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  
  //Q1 elements and Dynamic Damping WITH Anderson Acceleration
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_rhs + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  db["afc_anderson_damping"].set<>("no");
  errors = {{0.0003991651071879, 0.005300391980387}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + newton + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  db["afc_anderson_damping"].set<>("no");
  errors = {{0.00039916502859844, 0.0053003921620624}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_matrix + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  db["afc_anderson_damping"].set<>("no");
  errors = {{0.00039916491935108, 0.0053003924157625}};
  check_cd3d(gridCollections, db, errors);
  //=========================================================================
  
  return 0;
}
