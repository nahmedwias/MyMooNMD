#include <Example_Brinkman2D.h>

#include <Database.h>
#include <FEDatabase2D.h>
#include <SquareMatrix2D.h>
#include <string.h>
#include <MainUtilities.h>


/* examples */

namespace poiseuille_brinkman
{
#include "Brinkman_2D/Poiseuille.h"
}

namespace Poiseuille_Hannukainen
{
#include "Brinkman_2D/Poiseuille_Hannukainen.h"
}

namespace Poiseuille_Hannukainen_with_inscribed_physical_sphere
{
#include "Brinkman_2D/Poiseuille_Hannukainen_with_inscribed_physical_sphere.h"
}

namespace sine_cosine_brinkman
{
#include "Brinkman_2D/SinCos.h"
}

namespace sine2_sine2
{
#include "Brinkman_2D/Sin2Sin2.h"
}
//namespace driven_cavity
//{
//#include "Brinkman_2D/DrivenCavity.h"
//}
//
//namespace flow_around_cylinder
//{
//#include "Brinkman_2D/flow_around_cylinder.h"
//}

namespace poiseuille_channel
{
#include "Brinkman_2D/Poiseuille_Channel_1x2.h"
}




Example_Brinkman2D::Example_Brinkman2D(
  const ParameterDatabase& user_input_parameter_db) 
 : Example2D(user_input_parameter_db)
{
    int example_code = this->example_database["example"];
    switch( example_code )
    {
        case 0:
            /** exact_solution */
            exact_solution.push_back( poiseuille_brinkman::ExactU1 );
            exact_solution.push_back( poiseuille_brinkman::ExactU2 );
            exact_solution.push_back( poiseuille_brinkman::ExactP );
            
            /** boundary condition */
            boundary_conditions.push_back( poiseuille_brinkman::BoundCondition );
            boundary_conditions.push_back( poiseuille_brinkman::BoundCondition );
            boundary_conditions.push_back( BoundConditionNoBoundCondition );
            
            /** boundary values */
            boundary_data.push_back( poiseuille_brinkman::U1BoundValue );
            boundary_data.push_back( poiseuille_brinkman::U2BoundValue );
            boundary_data.push_back( BoundaryValueHomogenous );
            
            /** coefficients */
            problem_coefficients = poiseuille_brinkman::LinCoeffs;
            
            //this->example_database["equal_order_stab_weight"] = 500;
            ////poiseuille_brinkman::stab_weight = this->example_database["equal_order_stab_weight"];
            //double equal_order_stab_weight = this->example_database["equal_order_stab_weight"];
            
            poiseuille_brinkman::ExampleFile();

            break;
        case 1:
            /** exact_solution */
            exact_solution.push_back( Poiseuille_Hannukainen::ExactU1 );
            exact_solution.push_back( Poiseuille_Hannukainen::ExactU2 );
            exact_solution.push_back( Poiseuille_Hannukainen::ExactP );
            
            /** boundary condition */
            boundary_conditions.push_back( Poiseuille_Hannukainen::BoundCondition );
            boundary_conditions.push_back( Poiseuille_Hannukainen::BoundCondition );
            boundary_conditions.push_back( BoundConditionNoBoundCondition );
            
            /** boundary values */
            boundary_data.push_back( Poiseuille_Hannukainen::U1BoundValue );
            boundary_data.push_back( Poiseuille_Hannukainen::U2BoundValue );
            boundary_data.push_back( BoundaryValueHomogenous );
            
            /** coefficients */
            problem_coefficients = Poiseuille_Hannukainen::LinCoeffs;
            
            Poiseuille_Hannukainen::ExampleFile();
            break;
        case 2:
            /** exact_solution */
            exact_solution.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::ExactU1 );
            exact_solution.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::ExactU2 );
            exact_solution.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::ExactP );
            
            /** boundary condition */
            boundary_conditions.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::BoundCondition );
            boundary_conditions.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::BoundCondition );
            boundary_conditions.push_back( BoundConditionNoBoundCondition );
            
            /** boundary values */
            boundary_data.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::U1BoundValue );
            boundary_data.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::U2BoundValue );
            boundary_data.push_back( BoundaryValueHomogenous );
            
            /** coefficients */
            problem_coefficients = Poiseuille_Hannukainen_with_inscribed_physical_sphere::LinCoeffs;
            
            Poiseuille_Hannukainen_with_inscribed_physical_sphere::ExampleFile();
            break;
        case 3:
            /** exact_solution */
            exact_solution.push_back( sine_cosine_brinkman::ExactU1 );
            exact_solution.push_back( sine_cosine_brinkman::ExactU2 );
            exact_solution.push_back( sine_cosine_brinkman::ExactP );
            
            /** boundary condition */
            boundary_conditions.push_back(sine_cosine_brinkman::BoundCondition);
            boundary_conditions.push_back(sine_cosine_brinkman::BoundCondition);
            boundary_conditions.push_back(BoundConditionNoBoundCondition);
            
            /** boundary values */
            boundary_data.push_back( sine_cosine_brinkman::U1BoundValue );
            boundary_data.push_back( sine_cosine_brinkman::U2BoundValue );
            boundary_data.push_back( BoundaryValueHomogenous );
            
            /** coefficients */
            problem_coefficients = sine_cosine_brinkman::LinCoeffs;
            
            sine_cosine_brinkman::ExampleFile();
            break;
        case 4:
            /** exact_solution */
            exact_solution.push_back( sine2_sine2::ExactU1 );
            exact_solution.push_back( sine2_sine2::ExactU2 );
            exact_solution.push_back( sine2_sine2::ExactP );
            
            /** boundary condition */
            boundary_conditions.push_back( sine2_sine2::BoundCondition );
            boundary_conditions.push_back( sine2_sine2::BoundCondition );
            boundary_conditions.push_back( BoundConditionNoBoundCondition );
            
            /** boundary values */
            boundary_data.push_back( sine2_sine2::U1BoundValue );
            boundary_data.push_back( sine2_sine2::U2BoundValue );
            boundary_data.push_back( BoundaryValueHomogenous );
            
            /** coefficients */
            problem_coefficients = sine2_sine2::LinCoeffs;
            
            sine2_sine2::ExampleFile();
            break;
            
            //        case 3:
            //            /** exact_solution */
            //            exact_solution.push_back( driven_cavity::ExactU1 );
            //            exact_solution.push_back( driven_cavity::ExactU2 );
            //            exact_solution.push_back( driven_cavity::ExactP );
            //
            //            /** boundary condition */
            //            boundary_conditions.push_back( driven_cavity::BoundCondition );
            //            boundary_conditions.push_back( driven_cavity::BoundCondition );
            //            boundary_conditions.push_back( BoundConditionNoBoundCondition );
            //
            //            /** boundary values */
            //            boundary_data.push_back( driven_cavity::U1BoundValue );
            //            boundary_data.push_back( driven_cavity::U2BoundValue );
            //            boundary_data.push_back( BoundaryValueHomogenous );
            //
            //            /** coefficients */
            //            problem_coefficients = driven_cavity::LinCoeffs;
            //
            //            driven_cavity::ExampleFile();
            //            break;
            //        case 4:
            //            /** exact_solution */
            //            exact_solution.push_back( flow_around_cylinder::ExactU1 );
            //            exact_solution.push_back( flow_around_cylinder::ExactU2 );
            //            exact_solution.push_back( flow_around_cylinder::ExactP );
            //
            //            /** boundary condition */
            //            boundary_conditions.push_back( flow_around_cylinder::BoundCondition );
            //            boundary_conditions.push_back( flow_around_cylinder::BoundCondition );
            //            boundary_conditions.push_back( BoundConditionNoBoundCondition );
            //
            //            /** boundary values */
            //            boundary_data.push_back( flow_around_cylinder::U1BoundValue );
            //            boundary_data.push_back( flow_around_cylinder::U2BoundValue );
            //            boundary_data.push_back( BoundaryValueHomogenous );
            //            
            //            /** coefficients */
            //            problem_coefficients = flow_around_cylinder::LinCoeffs;
            //            
            //            flow_around_cylinder::ExampleFile();
            //            break;
        case 7:
            /** exact_solution */
            exact_solution.push_back( poiseuille_channel::ExactU1 );
            exact_solution.push_back( poiseuille_channel::ExactU2 );
            exact_solution.push_back( poiseuille_channel::ExactP );
            
            /** boundary condition */
            boundary_conditions.push_back( poiseuille_channel::BoundCondition );
            boundary_conditions.push_back( poiseuille_channel::BoundCondition );
            boundary_conditions.push_back( BoundConditionNoBoundCondition );
            
            /** boundary values */
            boundary_data.push_back( poiseuille_channel::U1BoundValue );
            boundary_data.push_back( poiseuille_channel::U2BoundValue );
            boundary_data.push_back( BoundaryValueHomogenous );
            
            /** coefficients */
            problem_coefficients = poiseuille_channel::LinCoeffs;
            
            poiseuille_channel::ExampleFile();
            break;

            
        default:
            ErrThrow("Unknown Brinkman example!");
    }
}

void Example_Brinkman2D::do_post_processing(Brinkman2D& brinkman2d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(brinkman2d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_Brinkman2D","No post processing done for the current example.");
  }
}


//double Example_Brinkman2D::get_stab() const
//{
//  double equal_order_stab_weight = this->example_database["equal_order_stab_weight"];
//  return equal_order_stab_weight;
//}
