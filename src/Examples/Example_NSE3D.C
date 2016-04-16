#include <Example_NSE3D.h>

#include <FEDatabase3D.h>
#include <Database.h>
#include <MainUtilities.h>


/* examples */


namespace ansatz_lin_const //0
{
#include "NSE_3D/AnsatzLinConst.h"
}
namespace ansatz_quad_lin //1
{
#include "NSE_3D/AnsatzQuadLin.h"
}
namespace cos_sin_simple //2
{
#include "NSE_3D/CosSin_simple.h"
}
namespace driven_cavity3d //3
{
#include "NSE_3D/DrivenCavity3D.h"
}
namespace flow_around_cylinder_stat
{
#include "NSE_3D/FlowAroundCylinder_stat.h"
}

//test examples
namespace test_u_0_p_0 //-1
{
#include "NSE_3D/test_u_0_p_0.h"
}
namespace test_u_1_p_0 //-2
{
#include "NSE_3D/test_u_1_p_0.h"
}
namespace test_u_2_p_1 //-3
{
#include "NSE_3D/test_u_2_p_1.h"
}
namespace test_u_3_p_2 //-4
{
#include "NSE_3D/test_u_3_p_2.h"
}

//========================================
// time dependent cases
namespace lin_space_time
{
  #include "TNSE_3D/linear_space_time.h"  // 101
}
//========================================

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

      case 4:
      {
        using namespace flow_around_cylinder_stat;
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
      case -3:
      {
        using namespace test_u_2_p_1;
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
        using namespace test_u_3_p_2;
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
      case 101:
      {
        using namespace lin_space_time;
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
        initial_conditions.push_back( InitialU1 );
        initial_conditions.push_back( InitialU2 );
        initial_conditions.push_back( InitialU3 );

        ExampleFile();
        break;
      }

    default:
      ErrThrow("Unknown Navier-Stokes example!");
  }
}

      
