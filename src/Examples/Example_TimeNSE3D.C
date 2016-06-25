#include <Example_TimeNSE3D.h>
#include <FEDatabase3D.h>
#include <Database.h>
#include <MainUtilities.h>

namespace lin_space_time
{
  #include "TNSE_3D/linear_space_time.h"  // 0
}
namespace AnsatzLinConst
{
  #include "TNSE_3D/AnsatzLinConst.h"     // 1
}
namespace Bsp0
{
  #include "TNSE_3D/Bsp0.h"   // 2
}
namespace Bsp1
{
  #include "TNSE_3D/Bsp1.h"   // 3
}
namespace Bsp2
{
 #include "TNSE_3D/Bsp2.h"    // 4
}
namespace Bsp3
{
  #include "TNSE_3D/Bsp3.h"   // 5
}

namespace tunnel_tdp1
{
#include "TNSE_3D/Tunnel1.h"
}
//=========================================================
Example_TimeNSE3D::Example_TimeNSE3D(int example_code)
{
  switch(example_code)
  {
    case 0:
    {
      using namespace lin_space_time;
      /** exact_solution */
      exact_solution.push_back( lin_space_time::ExactU1 );
      exact_solution.push_back( lin_space_time::ExactU2 );
      exact_solution.push_back( lin_space_time::ExactU3 );
      exact_solution.push_back( lin_space_time::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( lin_space_time::U1BoundValue );
      boundary_data.push_back( lin_space_time::U2BoundValue );
      boundary_data.push_back( lin_space_time::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = lin_space_time::LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( lin_space_time::InitialU1 );
      initialCondtion.push_back( lin_space_time::InitialU2 );
      initialCondtion.push_back( lin_space_time::InitialU3 );

      ExampleFile();
      break;
    }
    case 1:
    {
      using namespace AnsatzLinConst;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    case 2:
    {
      using namespace Bsp0;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    case 3:
    {
      using namespace Bsp1;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }    
    case 4:
    {
      using namespace Bsp2;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    case 5:
    {
      using namespace Bsp3;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    case 6:
    {
      using namespace tunnel_tdp1;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );

      ExampleFile();
      break;
    }
    default:
      ErrThrow("Unknown Example_TimeNSE3D example!");  
  }
}

void Example_TimeNSE3D::do_post_processing(Time_NSE3D& tnse3d) const
{
  //TODO
}

double Example_TimeNSE3D::get_nu() const
{
  if(nu==-1)
    ErrThrow("Kinematic viscosity is not set in this example!");
  
  return nu;
}

