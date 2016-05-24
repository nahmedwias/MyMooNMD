#include <Example_CD3D.h>

#include <Database.h>
#include <FEFunction3D.h>
#include <FEDatabase3D.h>
#include <SquareMatrix3D.h>
#include <string.h>

/* examples */

namespace sine_laplace_3D
{
  #include "CD_3D/Laplace.h"
}

namespace test_p0_zero
{
  #include <test_p0.h>
}

namespace test_p1
{
  #include "CD_3D/test_p1.h"
}

namespace test_p2
{
  #include "CD_3D/test_p2.h"
}

//=========================================================================
// time dependent problems
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
//=========================================================================

Example_CD3D::Example_CD3D(int example_code) : Example3D()
{
  switch( example_code )
  {
    //steady-state problems
    case 0:
      /** exact_solution */
      exact_solution.push_back( sine_laplace_3D::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( sine_laplace_3D::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sine_laplace_3D::BoundValue );
      
      /** coefficients */
      problem_coefficients = sine_laplace_3D::BilinearCoeffs;
      
      sine_laplace_3D::ExampleFile();
      break;

    //negative integers are reserved for pure test examples
    case -1:
    {//constant zero solution example
      using namespace test_p0_zero;
      exact_solution.push_back( Exact );
      boundary_conditions.push_back( BoundCondition );
      boundary_data.push_back( BoundValue );
      problem_coefficients = BilinearCoeffs;
      ExampleFile();
      break;
    }
    case -2:
    {//linear solution example
      using namespace test_p1;
      exact_solution.push_back( Exact );
      boundary_conditions.push_back( BoundCondition );
      boundary_data.push_back( BoundValue );
      problem_coefficients = BilinearCoeffs;
      ExampleFile();
      break;
    }
    case -3:
    {//quadratic solution example
      using namespace test_p2;
      exact_solution.push_back( Exact );
      boundary_conditions.push_back( BoundCondition );
      boundary_data.push_back( BoundValue );
      problem_coefficients = BilinearCoeffs;
      ExampleFile();
      break;
    }
    //=====================================================================
    // time dependent problems
    case -4:
    {
      // linear space and time solution example 
      using namespace linear_space_time;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initial_conditions.push_back(InitialCondition);
      ExampleFile();
    }
    break;
    case -5:
    {
      using namespace quad_space_time;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initial_conditions.push_back(InitialCondition);
      ExampleFile();
    }
    break;
    case 101:
      using namespace concentration;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initial_conditions.push_back(InitialCondition);
      ExampleFile();
      break;
    default:
      ErrThrow("Unknown name of the convection-diffusion-reaction example CD3D or Time_CD3D!");
  }
}

      
