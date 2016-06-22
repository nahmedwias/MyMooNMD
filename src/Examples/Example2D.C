#include <Example2D.h>


ParameterDatabase get_default_Example2D_parameters()
{
  Output::print<3>("creating a default Example2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Example2D parameter database");

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

Example2D::Example2D() 
 : example_database(get_default_Example2D_parameters()), exact_solution(),
   boundary_conditions(), boundary_data(), problem_coefficients(NULL)
{
}

Example2D::Example2D(std::vector <DoubleFunct2D*> exact,
                     std::vector <BoundCondFunct2D*> bc,
                     std::vector <BoundValueFunct2D*> bd, CoeffFct2D *coeffs)
 : example_database(get_default_Example2D_parameters()), exact_solution(exact),
   boundary_conditions(bc), boundary_data(bd), problem_coefficients(coeffs)
{ 
}
