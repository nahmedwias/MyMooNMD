#include <Example_CD2D.h>

#include <Database.h>
#include <FEFunction2D.h>
#include <FEDatabase2D.h>
#include <SquareMatrix2D.h>
#include <string.h>

//===========================================================
// examples for stationary convection-diffusion-reaction
// problems
//===========================================================

namespace sine_laplace
{
  #include "CD_2D/SineLaplace.h"
}

namespace two_interior_layers
{
  #include "CD_2D/TwoInteriorLayers.h"
}

namespace hemker_1996
{
  #include "CD_2D/Hemker1996.h"
}

namespace sharp_boundary_layer
{
  #include "CD_2D/SharpBoundaryLayer.h"
}

Example_CD2D::Example_CD2D(int example_code) : Example2D()
{
  switch( example_code)
  {
    case 0:
      /** exact_solution */
      exact_solution.push_back( sine_laplace::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( sine_laplace::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sine_laplace::BoundValue );
      
      /** coefficients */
      problem_coefficients = sine_laplace::BilinearCoeffs;
      
      sine_laplace::ExampleFile();
      break;
    case 1:
      /** exact_solution */
      exact_solution.push_back( two_interior_layers::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( two_interior_layers::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( two_interior_layers::BoundValue );
      
      /** coefficients */
      problem_coefficients = two_interior_layers::BilinearCoeffs;
      
      two_interior_layers::ExampleFile();
      break;
    case 2:
      /** exact_solution */
      exact_solution.push_back( hemker_1996::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( hemker_1996::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( hemker_1996::BoundValue );
      
      /** coefficients */
      problem_coefficients = hemker_1996::BilinearCoeffs;
      
      hemker_1996::ExampleFile();
      break;

    case 3:
      /** exact_solution */
      exact_solution.push_back( sharp_boundary_layer::Exact );

      /** boundary condition */
      boundary_conditions.push_back( sharp_boundary_layer::BoundCondition );

      /** boundary values */
      boundary_data.push_back( sharp_boundary_layer::BoundValue );

      /** coefficients */
      problem_coefficients = sharp_boundary_layer::BilinearCoeffs;

      sharp_boundary_layer::ExampleFile();
      break;

    default:
      ErrThrow("Unknown name of the convection-diffusion (CD2D) example!", 
               example_code);
  }
}

      
