#include <Example_CD2D.h>
#include <ConvectionDiffusion.h>

#include <Database.h>
#include <FEFunction2D.h>
#include <FEDatabase2D.h>
#include <SquareMatrix2D.h>
#include "BaseCell.h"
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

namespace smooth_solution
{
  #include "CD_2D/Smooth.h"
}

namespace HMM1986
{
  #include "CD_2D/HMM1986.h"
}

namespace hump
{
  #include "CD_2D/hump.h"
}

namespace boundary_layer_known
{
  #include "CD_2D/boundary_layer_known.h"
}

namespace example_1_ABR17
{
  #include "CD_2D/example_1_ABR17.h"
}

Example_CD2D::Example_CD2D(const ParameterDatabase& user_input_parameter_db) 
 : Example2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
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
      
    case 4:
      /** exact_solution */
      exact_solution.push_back( smooth_solution::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( smooth_solution::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( smooth_solution::BoundValue );
      
      /** coefficients */
      problem_coefficients = smooth_solution::BilinearCoeffs;
      
      smooth_solution::ExampleFile();
      break;
      
    case 5:
      /** exact_solution */
      exact_solution.push_back( HMM1986::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( HMM1986::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( HMM1986::BoundValue );
      
      /** coefficients */
      problem_coefficients = HMM1986::BilinearCoeffs;
      
      HMM1986::ExampleFile();
      break;
    
    case 6:
      /** exact_solution */
      exact_solution.push_back( hump::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( hump::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( hump::BoundValue );
      
      /** coefficients */
      problem_coefficients = hump::BilinearCoeffs;
      
      hump::ExampleFile();
      break;
    
    case 7:
      /** exact_solution */
      exact_solution.push_back( boundary_layer_known::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( boundary_layer_known::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( boundary_layer_known::BoundValue );
      
      /** coefficients */
      problem_coefficients = boundary_layer_known::BilinearCoeffs;
      
      boundary_layer_known::ExampleFile();
      break;
    case 8:
      /** exact_solution */
      exact_solution.push_back( example_1_ABR17::Exact );
      
      /** boundary condition */
      boundary_conditions.push_back( example_1_ABR17::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( example_1_ABR17::BoundValue );
      
      /** coefficients */
      problem_coefficients = example_1_ABR17::BilinearCoeffs;
      
      example_1_ABR17::ExampleFile();
      break;

    default:
      ErrThrow("Unknown name of the convection-diffusion (CD2D) example!", 
               example_code);
  }
}

Example_CD2D::Example_CD2D(const std::vector<DoubleFunct2D*>& exact,
                           const std::vector<BoundCondFunct2D*>& bc,
                           const std::vector<BoundValueFunct2D*>& bd,
                           const CoeffFct2D& coeffs, double nu)
: Example2D(exact, bc, bd, coeffs)
{
  this->example_database["diffusion_coefficient"] = nu;
}


void Example_CD2D::do_post_processing(ConvectionDiffusion<2>& cd2d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(cd2d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_CD2D","No post processing done for the current example.");
  }
}

double Example_CD2D::get_nu() const
{
  double diffusion_coefficient = this->example_database["diffusion_coefficient"];
  return diffusion_coefficient;
}

