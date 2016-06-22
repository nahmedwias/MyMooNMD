#include <Example_TimeCD3D.h>
#include <Database.h>
#include <string.h>

namespace linear_space_time
{
  #include "TCD_3D/linear_space_time.h"
}

namespace quad_space_time
{
#include "TCD_3D/quadratic_space_time.h"
}

namespace concentration
{
#include "TCD_3D/concentrationOfSpecies_3d.h"
}

Example_TimeCD3D::Example_TimeCD3D(int example_code)
{
  switch(example_code)
  {
    case -2:
    {
      // linear space and time solution example 
      using namespace linear_space_time;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initialCondition.push_back(InitialCondition);
      ExampleFile();
    }
      break;
    case -1:
    {
      using namespace quad_space_time;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initialCondition.push_back(InitialCondition);
      ExampleFile();
    }
      break;
    case 0:
      using namespace concentration;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initialCondition.push_back(InitialCondition);
      ExampleFile();
      break;
    default:
      ErrThrow("Unknown Example_TimeCD3D example!");
  }
}
