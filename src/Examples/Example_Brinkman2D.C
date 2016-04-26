#include <Example_Brinkman2D.h>

#include <Database.h>
#include <FEDatabase2D.h>
#include <SquareMatrix2D.h>
#include <string.h>
#include <MainUtilities.h>


/* examples */

namespace poiseuille
{
  #include "Brinkman_2D/Poiseuille.h"
}
namespace driven_cavity
{
  #include "Brinkman_2D/DrivenCavity.h"
}
namespace sine_cosine
{
  #include "Brinkman_2D/SinCos.h"
}
namespace flow_around_cylinder
{
  #include "Brinkman_2D/flow_around_cylinder.h"
}
namespace poiseuille_Hannukainen
{
#include "Brinkman_2D/Poiseuille_Hannukainen.h"
}


Example_Brinkman2D::Example_Brinkman2D() : Example2D()
{
  switch( TDatabase::ParamDB->EXAMPLE ) 
  {
   case 0:
     /** exact_solution */
     exact_solution.push_back( poiseuille::ExactU1 );
     exact_solution.push_back( poiseuille::ExactU2 );
     exact_solution.push_back( poiseuille::ExactP );
 
     /** boundary condition */
     boundary_conditions.push_back( poiseuille::BoundCondition );
     boundary_conditions.push_back( poiseuille::BoundCondition );
     boundary_conditions.push_back( BoundConditionNoBoundCondition );
 
     /** boundary values */
     boundary_data.push_back( poiseuille::U1BoundValue );
     boundary_data.push_back( poiseuille::U2BoundValue );
     boundary_data.push_back( BoundaryValueHomogenous );
  
      /** coefficients */
      problem_coefficients = poiseuille::LinCoeffs;
  
      poiseuille::ExampleFile();
      break;
    case 1:
      /** exact_solution */
      exact_solution.push_back( driven_cavity::ExactU1 );
      exact_solution.push_back( driven_cavity::ExactU2 );
      exact_solution.push_back( driven_cavity::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( driven_cavity::BoundCondition );
      boundary_conditions.push_back( driven_cavity::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( driven_cavity::U1BoundValue );
      boundary_data.push_back( driven_cavity::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = driven_cavity::LinCoeffs;
      
      driven_cavity::ExampleFile();
      break;
    case 2:
      /** exact_solution */
      exact_solution.push_back( sine_cosine::ExactU1 );
      exact_solution.push_back( sine_cosine::ExactU2 );
      exact_solution.push_back( sine_cosine::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( sine_cosine::BoundCondition );
      boundary_conditions.push_back( sine_cosine::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sine_cosine::U1BoundValue );
      boundary_data.push_back( sine_cosine::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = sine_cosine::LinCoeffs;
      
      sine_cosine::ExampleFile();
      break;
    case 3:
      /** exact_solution */
      exact_solution.push_back( flow_around_cylinder::ExactU1 );
      exact_solution.push_back( flow_around_cylinder::ExactU2 );
      exact_solution.push_back( flow_around_cylinder::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( flow_around_cylinder::BoundCondition );
      boundary_conditions.push_back( flow_around_cylinder::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( flow_around_cylinder::U1BoundValue );
      boundary_data.push_back( flow_around_cylinder::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = flow_around_cylinder::LinCoeffs;
      
      flow_around_cylinder::ExampleFile();
      break;
      case 4:
          /** exact_solution */
          exact_solution.push_back( poiseuille_Hannukainen::ExactU1 );
          exact_solution.push_back( poiseuille_Hannukainen::ExactU2 );
          exact_solution.push_back( poiseuille_Hannukainen::ExactP );
          
          /** boundary condition */
          boundary_conditions.push_back( poiseuille_Hannukainen::BoundCondition );
          boundary_conditions.push_back( poiseuille_Hannukainen::BoundCondition );
          boundary_conditions.push_back( BoundConditionNoBoundCondition );
          
          /** boundary values */
          boundary_data.push_back( poiseuille_Hannukainen::U1BoundValue );
          boundary_data.push_back( poiseuille_Hannukainen::U2BoundValue );
          boundary_data.push_back( BoundaryValueHomogenous );
          
          /** coefficients */
          problem_coefficients = poiseuille_Hannukainen::LinCoeffs;
          
          poiseuille_Hannukainen::ExampleFile();
          break;
    default:
      ErrThrow("Unknown Brinkman example!");
  }
}

      