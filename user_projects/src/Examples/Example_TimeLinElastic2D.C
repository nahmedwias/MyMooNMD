#include <Example_TimeLinElastic2D.h>
#include <Database.h>
#include <MainUtilities.h>

#include <string>


//=========================================
Example_TimeLinElastic2D::Example_TimeLinElastic2D(
    const ParameterDatabase& user_input_parameter_db)
: Example_NonStationary2D(user_input_parameter_db)
{
  cout << "CONSTRUCTOR OF EXAMPLE_TimeLinElastic2D OK!" << endl;
}
