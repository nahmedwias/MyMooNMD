#include <Example_TimeLinElastic2D.h>
#include <Database.h>
#include <MainUtilities.h>

#include <string>


namespace test_tlinelastic // example number 0
{
 #include <0_Test_TLinElastic2D.h>
}

//=========================================
Example_TimeLinElastic2D::Example_TimeLinElastic2D(
    const ParameterDatabase& user_input_parameter_db)
: Example_NonStationary2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];

  switch(example_code)
    {
      case 0:       // 0_Test_TLinElastic2D.h
        /** output example name */
        test_tlinelastic::ExampleFile();

        /** exact_solution */
        exact_solution.push_back( test_tlinelastic::ExactU1 );
        exact_solution.push_back( test_tlinelastic::ExactU2 );

        /** initial condition */
        initialCOndtion.push_back(test_tlinelastic::InitialU1);
        initialCOndtion.push_back(test_tlinelastic::InitialU2);
        initialCOndtion.push_back(test_tlinelastic::InitialV1);
        initialCOndtion.push_back(test_tlinelastic::InitialV2);

        /** boundary condition */
        boundary_conditions.push_back( test_tlinelastic::BoundCondition );

        /** boundary values */
        boundary_data.push_back( test_tlinelastic::BoundValueU1 );
        boundary_data.push_back( test_tlinelastic::BoundValueU2 );

        /** coefficients */
        problem_coefficients = test_tlinelastic::Coefficients;

        this->timeDependentRhs   =test_tlinelastic::rhs_depends_on_time;
        this->timeDependentCoeffs=test_tlinelastic::coefficients_depend_on_time;

        break;
      default:
        ErrThrow("Unknown example");
    }

  cout << "CONSTRUCTOR OF EXAMPLE_TimeLinElastic2D OK!" << endl;
}
