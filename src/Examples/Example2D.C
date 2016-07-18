#include <Example2D.h>


ParameterDatabase Example2D::default_example_database()
{
  Output::print<3>("creating a default Example2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db("Example2D parameter database");
  
  db.add("example", 0,
         "Choose which example to run. \nNote that depending on the type of "
         "problem you want to solve, different values are meaningful here. See "
         "the derived classes of 'Example2D'.", -5, 200);

  /** TDatabase::ParamDB->RE_NR */
     db.add("reynolds_number", 1.,
            "Reynolds number: dimensionless number which describes how viscous the  "
            "flow is (laminar, turbulent). Re = U.L/nu, where nu is the kinematic "
            "viscosity (=mu/rho)."
            "The higher it is, the more turbulent "
            "the flow is. Reynolds number can often be in the order of "
            "magnitude of millions. Then, the classical NSE is not relevant anymore"
            "to describe the flow. One should in these cases use turbulence Models,"
            "like Smagorinsky or k-eps. Maximum value is 1000,"
            "which already corresponds to a turbulent flow. Default value is 1."
            "Note that this is also equal to 1/eps, where eps is the DIMENSIONLESS "
            "viscosity (sometimes misleadingly named as nu) "
            "which sometimes appears in the dimensionless Navier-Stokes equations"
            "instead of 1/Re.",
            0., 1000.);

     /** TDatabase::ParamDB->PE_NR */
     db.add("peclet_number", 1.,
            "Peclet number: dimensionless number which is used to compare the "
            "convective and diffusive terms in transport processes, also called "
            "diffusion coefficient.",
            0., 1000.);

     /** TDatabase::ParamDB->P0 */
     db.add("variable_parameter0", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P1 */
     db.add("variable_parameter1", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P2 */
     db.add("variable_parameter2", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P3 */
     db.add("variable_parameter3", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P4 */
     db.add("variable_parameter4", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P5 */
     db.add("variable_parameter5", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P6 */
     db.add("variable_parameter6", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P7 */
     db.add("variable_parameter7", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example class, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P8 */
     db.add("variable_parameter8", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P9 */
     db.add("variable_parameter9", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P10 */
     db.add("variable_parameter10", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P11 */
     db.add("variable_parameter11", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P12 */
     db.add("variable_parameter12", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P13 */
     db.add("variable_parameter13", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P14 */
     db.add("variable_parameter14", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

     /** TDatabase::ParamDB->P15 */
     db.add("variable_parameter15", (size_t) 0,
            "A variable parameter which can be used to change some value of the"
            "Example, directly from the .dat file",
            (size_t) -1000., (size_t) 1000.);

  return db;
}

Example2D::Example2D(const ParameterDatabase & db) 
 : example_database(Example2D::default_example_database()), exact_solution(),
   boundary_conditions(), boundary_data(), problem_coefficients(NULL)
{
  this->example_database.merge(db, false);
}

Example2D::Example2D(std::vector <DoubleFunct2D*> exact,
                     std::vector <BoundCondFunct2D*> bc,
                     std::vector <BoundValueFunct2D*> bd, CoeffFct2D *coeffs)
 : example_database(Example2D::default_example_database()), exact_solution(exact),
   boundary_conditions(bc), boundary_data(bd), problem_coefficients(coeffs)
{ 
}
