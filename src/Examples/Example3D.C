#include <Example3D.h>


ParameterDatabase get_default_Example3D_parameters()
{
  Output::print<3>("creating a default Example3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Example3D parameter database");

  //NSE3D requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default example database
  ParameterDatabase example_db = ParameterDatabase::default_example_database();
  db.merge(example_db, true);

  return db;
}

Example3D::Example3D()
 : example_database(get_default_Example3D_parameters()), exact_solution(),
   boundary_conditions(), boundary_data(), problem_coefficients(NULL)
{
}

Example3D::Example3D(std::vector <DoubleFunct3D*> exact,
                     std::vector <BoundCondFunct3D*> bc,
                     std::vector <BoundValueFunct3D*> bd, CoeffFct3D *coeffs)
 : example_database(get_default_Example3D_parameters()), exact_solution(exact),
   boundary_conditions(bc), boundary_data(bd), problem_coefficients(coeffs)
{ 
}
