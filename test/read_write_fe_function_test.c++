#include <MooNMD_Io.h>
#include <FEFunction2D.h>
#include <BlockVector.h>
#include <FEDatabase2D.h>
#include <Database.h>
#include <MainUtilities.h>

// some functions dummy given as an analytical representation
void analytic_function(double x, double y, double * values)
{
  values[0] = x*y;
}


int main(int argc, char **argv)
{
    //  declaration of old databases (unfortunately still needed in ParMooN)
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  Output::setVerbosity(2);

  // new database
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();

  // set a few parameters in the database to construct a unit square with one
  // cell. then refine uniformly twice
  db.add("refinement_n_initial_steps", (size_t) 2, "");
  db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "UnitSquare", "");
  // construct a domain object
  TDomain domain(db);

  // refine grid
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(unsigned int i = 0; i < n_ref; i++)
  {
    domain.RegRefineAll();
  }

  // create finite element space, for that we need a
  // collection of grid cells ( = grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0);
  // type of finite element on each cell (here: Q2)
  size_t ansatz_order = 2;
  TFESpace2D fe_space(coll, (char*) "fe_space", (char*) "fe_space",
                      BoundConditionNoBoundCondition, ansatz_order, nullptr);
  fe_space.info();

  // create some vector which could be seen as a finite element function
  // representation
  BlockVector v1(fe_space.GetN_DegreesOfFreedom());
  {
    // fill the vector with some values, according to some given function
    TFEFunction2D f(&fe_space, (char*)"dummy", (char*) "dummy",
                    v1.get_entries(), v1.length());
    f.Interpolate(analytic_function);
    //Output::print("function norm ", v1.norm());
  }

  // write the vector
  std::stringstream out_in_stream;
  v1.write_to_stream(out_in_stream);

  // create second vector, read into it
  BlockVector v2(fe_space.GetN_DegreesOfFreedom());
  v2.read_from_stream(out_in_stream);

  // substract first vector (difference should be zero now)
  v2.add_scaled(v1, -1.);

  double err = v2.norm();
  if(err != 0.0)
  {
    ErrThrow("writing and then reading a BlockVector failed with error ", err);
  }
  Output::print("test successful");
}
