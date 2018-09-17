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

#ifdef __2D__
void check_2d(ParameterDatabase db, TDomain &domain, int velocity_order)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  // automatically choose inf-sup stable pressure space
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  
  NSE2D nse(domain, db);
  Output::print("created nse2d object");
  nse.assemble();
  nse.solve();
  nse.output();
  
  auto div_error = nse.get_errors()[1];
  if(div_error > 1.e-12)
  {
    ErrThrow("Scott-Vogelius elements should lead to zero divergence, but ",
             div_error, " has been computed");
  }
}
#else // 2D -> 3D
void check_3d(ParameterDatabase db, const std::list<TCollection*>& colls, 
              int velocity_order)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  // automatically choose inf-sup stable pressure space
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  
  Example_NSE3D example_obj(db);
  // Construct the nse3d problem object.
#ifdef _MPI
  NSE3D nse(colls, db, example_obj, maxSubDomainPerDof);
#else // no mpi
  NSE3D nse(colls, db, example_obj);
#endif
  nse.assemble_linear_terms();
  nse.solve();
  nse.output();
  
  auto div_error = nse.get_errors()[1];
  if(div_error > 1.e-12)
  {
    ErrThrow("Scott-Vogelius elements should lead to zero divergence, but ",
             div_error, " has been computed");
  }
}
#endif // 2D

int main(int , char** )
{
  //  declaration of databases
  TDatabase Database;
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(TDomain::default_domain_parameters(), true);
#ifdef __2D__
  db.merge(Example2D::default_example_database(), true);
#else // 3D
  db.merge(Example3D::default_example_database(), true);
#endif // 2D
  db["refinement_n_initial_steps"] = 2;
  db["refinement_final_step_barycentric"] = true;
  db["example"] = 1;
  Output::setVerbosity(3);
  
#ifdef __2D__
  TFEDatabase2D FEDatabase;
  db["boundary_file"] = "Default_UnitSquare";
  db["geo_file"] = "TwoTriangles";
#else // 3D
  TFEDatabase3D FEDatabase;
  db["boundary_file"] = "Default_UnitCube";
  db["geo_file"] = "Default_UnitCube_Tetra";
  db["refinement_n_initial_steps"] = 1;
#endif // 2D
  //db.add("output_write_vtk", true, "");
  
  TDomain domain(db);
  // refinement: in 2D this has to be done here, in 3D this is done in a domain
  // method
#ifdef __2D__
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(unsigned int i=0; i < n_ref; i++)
  {
    domain.RegRefineAll();
  }
  if(domain.get_database()["refinement_final_step_barycentric"])
    domain.barycentric_refinement();
#else // 2D -> 3D
  auto gridCollections = domain.refine_and_get_hierarchy_of_collections(db
  #ifdef _MPI
      , maxSubDomainPerDof
  #endif
      );
#endif
  
#ifdef __2D__
  check_2d(db, domain, 12);
  check_2d(db, domain, 13);
  check_2d(db, domain, 14);
#else // 2D -> 3D
  check_3d(db, gridCollections, 13);
  // check_3d(db, gridCollections, 14); // P4 on tetrahedra not implemented!
#endif // 2D
}

