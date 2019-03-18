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

namespace poiseuille_channel
{
#include "Brinkman_2D/Poiseuille_Channel_1x2.h"
}

namespace SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D
{
#include "Brinkman_2D/SinCos_BadiaCodina_DarcyFlow_Test.h"
}


namespace Discacciati_Flow
{
#include "Brinkman_2D/Discacciati_Flow.h"
}


namespace Geothermal_Energy_Brinkman2D
{
#include "Brinkman_2D/Geothermal_Energy_Brinkman2D.h"
}

namespace Riverbed_Brinkman2D
{
#include "Brinkman_2D/Riverbed_Brinkman2D.h"
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

      // read parameters from local database
      poiseuille_brinkman::viscosity = get_viscosity();
      poiseuille_brinkman::effective_viscosity = get_effective_viscosity();
      poiseuille_brinkman::permeability = get_permeablity();

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

      // read parameters from local database
      Poiseuille_Hannukainen::viscosity = get_viscosity();
      Poiseuille_Hannukainen::effective_viscosity = get_effective_viscosity();
      Poiseuille_Hannukainen::permeability = get_permeablity();

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

      // read parameters from local database
      Poiseuille_Hannukainen_with_inscribed_physical_sphere::viscosity = get_viscosity();
      Poiseuille_Hannukainen_with_inscribed_physical_sphere::effective_viscosity = get_effective_viscosity();
      Poiseuille_Hannukainen_with_inscribed_physical_sphere::permeability = get_permeablity();

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

      // read parameters from local database
      sine_cosine_brinkman::viscosity = get_viscosity();
      sine_cosine_brinkman::effective_viscosity = get_effective_viscosity();
      sine_cosine_brinkman::permeability = get_permeablity();

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

      // read parameters from local database
      sine2_sine2::viscosity = get_viscosity();
      sine2_sine2::effective_viscosity = get_effective_viscosity();
      sine2_sine2::permeability = get_permeablity();

      sine2_sine2::ExampleFile();
      break;

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

      // read parameters from local database
      poiseuille_channel::viscosity = get_viscosity();
      poiseuille_channel::effective_viscosity = get_effective_viscosity();
      poiseuille_channel::permeability = get_permeablity();

      poiseuille_channel::ExampleFile();
      break;

    case 8:
      /** exact_solution */
      exact_solution.push_back( SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::ExactU1 );
      exact_solution.push_back( SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::ExactU2 );
      exact_solution.push_back( SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::BoundCondition );
      boundary_conditions.push_back( SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::U1BoundValue );
      boundary_data.push_back( SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::LinCoeffs;

      // read parameters from local database
      SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::viscosity = get_viscosity();
      SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::effective_viscosity = get_effective_viscosity();
      SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::permeability = get_permeablity();

      SinCos_BadiaCodina_ForDarcyLimitOfBrinkman2D::ExampleFile();
      break;

    case 9:
      /** exact_solution */
      exact_solution.push_back( Discacciati_Flow::ExactU1 );
      exact_solution.push_back( Discacciati_Flow::ExactU2 );
      exact_solution.push_back( Discacciati_Flow::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( Discacciati_Flow::BoundCondition );
      boundary_conditions.push_back( Discacciati_Flow::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( Discacciati_Flow::U1BoundValue );
      boundary_data.push_back( Discacciati_Flow::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = Discacciati_Flow::LinCoeffs;

      // read parameters from local database
      Discacciati_Flow::viscosity = get_viscosity();
      Discacciati_Flow::effective_viscosity = get_effective_viscosity();
      Discacciati_Flow::permeability = get_permeablity();

      Discacciati_Flow::ExampleFile();
      break;

case 10: // Geothermal_Energy_Brinkman2D
      /** exact_solution */
      exact_solution.push_back( Geothermal_Energy_Brinkman2D::ExactU1 );
      exact_solution.push_back( Geothermal_Energy_Brinkman2D::ExactU2 );
      exact_solution.push_back( Geothermal_Energy_Brinkman2D::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( Geothermal_Energy_Brinkman2D::BoundCondition );
      boundary_conditions.push_back( Geothermal_Energy_Brinkman2D::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( Geothermal_Energy_Brinkman2D::U1BoundValue );
      boundary_data.push_back( Geothermal_Energy_Brinkman2D::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = Geothermal_Energy_Brinkman2D::LinCoeffs;

      // read parameters from local database
      Geothermal_Energy_Brinkman2D::viscosity = get_viscosity();
      Geothermal_Energy_Brinkman2D::effective_viscosity = get_effective_viscosity();
      Geothermal_Energy_Brinkman2D::permeability = get_permeablity();

      Geothermal_Energy_Brinkman2D::ExampleFile();
      break;


case 11: // Riverbed
      /** exact_solution */
      exact_solution.push_back( Riverbed_Brinkman2D::ExactU1 );
      exact_solution.push_back( Riverbed_Brinkman2D::ExactU2 );
      exact_solution.push_back( Riverbed_Brinkman2D::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( Riverbed_Brinkman2D::BoundCondition );
      boundary_conditions.push_back( Riverbed_Brinkman2D::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( Riverbed_Brinkman2D::U1BoundValue );
      boundary_data.push_back( Riverbed_Brinkman2D::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = Riverbed_Brinkman2D::LinCoeffs;

      // read parameters from local database
      Riverbed_Brinkman2D::viscosity = get_viscosity();
      Riverbed_Brinkman2D::effective_viscosity = get_effective_viscosity();
      Riverbed_Brinkman2D::permeability = get_permeablity();

      Riverbed_Brinkman2D::ExampleFile();
      break;

    default:
      ErrThrow("Unknown Brinkman example!");
  }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// functions to avoid global parameters (TDatabase) in the example file
double Example_Brinkman2D::get_viscosity() const
{ 
  double viscosity = this->example_database["viscosity"];
  return viscosity;
}

double Example_Brinkman2D::get_effective_viscosity() const
{ 
  double effective_viscosity = this->example_database["effective_viscosity"];
  return effective_viscosity; 
}

double Example_Brinkman2D::get_permeablity() const
{
  double K = this->example_database["permeability"];
  /*if (this->example_database["read_permeability_from_file"])
    K = -1;
   */
  return K;
}

//double Example_Brinkman2D::get_stab() const
//{
//  double equal_order_stab_weight = this->example_database["equal_order_stab_weight"];
//  return equal_order_stab_weight;
//}
