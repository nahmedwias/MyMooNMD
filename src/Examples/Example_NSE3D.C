#include <Example_NSE3D.h>

#include <FEDatabase3D.h>
#include <Database.h>
#include <MainUtilities.h>


/* examples */


namespace driven_cavity3d
{
#include "NSE_3D/DrivenCavity3D.h"
}
namespace cos_sin_simple
{
#include "NSE_3D/CosSin_simple.h"
}
namespace ansatz_quad_lin
{
#include "NSE_3D/AnsatzQuadLin.h"
}
namespace ansatz_lin_const
{
#include "NSE_3D/AnsatzLinConst.h"
}

//test examples
namespace test_u_0_p_0
{
#include "NSE_3D/test_u_0_p_0.h"
}
namespace test_u_1_p_0_zero
{
#include "NSE_3D/test_u_1_p_0_zero.h"
}
namespace test_u_1_p_0
{
#include "NSE_3D/test_u_1_p_0.h"
}
namespace test_u_1_p_1_x
{
#include "NSE_3D/test_u_1_p_1_x.h"
}
namespace test_u_1_p_1_z
{
#include "NSE_3D/test_u_1_p_1_z.h"
}

Example_NSE3D::Example_NSE3D() : Example3D()
{
  switch( TDatabase::ParamDB->EXAMPLE ) 
  {
    case 0:
      /** exact_solution */
      exact_solution.push_back( ansatz_lin_const::ExactU1 );
      exact_solution.push_back( ansatz_lin_const::ExactU2 );
      exact_solution.push_back( ansatz_lin_const::ExactU3 );
      exact_solution.push_back( ansatz_lin_const::ExactP );

      /* boundary condition */
      boundary_conditions.push_back( ansatz_lin_const::BoundCondition );
      boundary_conditions.push_back( ansatz_lin_const::BoundCondition );
      boundary_conditions.push_back( ansatz_lin_const::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /* boundary values */
      boundary_data.push_back( ansatz_lin_const::U1BoundValue );
      boundary_data.push_back( ansatz_lin_const::U2BoundValue );
      boundary_data.push_back( ansatz_lin_const::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /* coefficients */
      problem_coefficients = ansatz_lin_const::LinCoeffs;

      ansatz_lin_const::ExampleFile();
      break;
    case 1:
      /** exact_solution */
      exact_solution.push_back( ansatz_quad_lin::ExactU1 );
      exact_solution.push_back( ansatz_quad_lin::ExactU2 );
      exact_solution.push_back( ansatz_quad_lin::ExactU3 );
      exact_solution.push_back( ansatz_quad_lin::ExactP );

      /* boundary condition */
      boundary_conditions.push_back( ansatz_quad_lin::BoundCondition );
      boundary_conditions.push_back( ansatz_quad_lin::BoundCondition );
      boundary_conditions.push_back( ansatz_quad_lin::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /* boundary values */
      boundary_data.push_back( ansatz_quad_lin::U1BoundValue );
      boundary_data.push_back( ansatz_quad_lin::U2BoundValue );
      boundary_data.push_back( ansatz_quad_lin::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /* coefficients */
      problem_coefficients = ansatz_quad_lin::LinCoeffs;

      ansatz_quad_lin::ExampleFile();
      break;
    case 2:
      /** exact_solution */
      exact_solution.push_back( cos_sin_simple::ExactU1 );
      exact_solution.push_back( cos_sin_simple::ExactU2 );
      exact_solution.push_back( cos_sin_simple::ExactU3 );
      exact_solution.push_back( cos_sin_simple::ExactP );

      /* boundary condition */
      boundary_conditions.push_back( cos_sin_simple::BoundCondition );
      boundary_conditions.push_back( cos_sin_simple::BoundCondition );
      boundary_conditions.push_back( cos_sin_simple::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /* boundary values */
      boundary_data.push_back( cos_sin_simple::U1BoundValue );
      boundary_data.push_back( cos_sin_simple::U2BoundValue );
      boundary_data.push_back( cos_sin_simple::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /* coefficients */
      problem_coefficients = cos_sin_simple::LinCoeffs;

      cos_sin_simple::ExampleFile();
      break;
    
      
      case 3:
      /** exact_solution */
      exact_solution.push_back( driven_cavity3d::ExactU1 );
      exact_solution.push_back( driven_cavity3d::ExactU2 );
      exact_solution.push_back( driven_cavity3d::ExactU3 );
      exact_solution.push_back( driven_cavity3d::ExactP );
  
      /** boundary condition */
      boundary_conditions.push_back( driven_cavity3d::BoundCondition );
      boundary_conditions.push_back( driven_cavity3d::BoundCondition );
      boundary_conditions.push_back( driven_cavity3d::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
  
      /** boundary values */
      boundary_data.push_back( driven_cavity3d::U1BoundValue );
      boundary_data.push_back( driven_cavity3d::U2BoundValue );
      boundary_data.push_back( driven_cavity3d::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
  
      /** coefficients */
      problem_coefficients = driven_cavity3d::LinCoeffs;
      
      driven_cavity3d::ExampleFile();
      break;

      case -1:
      {
        using namespace test_u_0_p_0;
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

        ExampleFile();
        break;
      }
      case -2:
      {
        using namespace test_u_1_p_0_zero;
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

        ExampleFile();
        break;
      }
      case -3:
      {
        using namespace test_u_1_p_0;
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

        ExampleFile();
        break;
      }
      case -4:
      {
        using namespace test_u_1_p_1_x;
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

        ExampleFile();
        break;
      }
      case -5:
      {
        using namespace test_u_1_p_1_z;
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

        ExampleFile();
        break;
      }
    default:
      ErrThrow("Unknown Navier-Stokes example!");
  }
}

      
