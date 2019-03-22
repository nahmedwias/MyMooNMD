#include <Example_Brinkman3D.h>

#include <Database.h>
#include <FEDatabase3D.h>
#include <SquareMatrix3D.h>
#include <string.h>
#include <MainUtilities.h>
#include <Brinkman3D.h>



#ifdef _MPI
#include <mpi.h>
#endif


/* examples */


namespace ansatz_lin_const //0
{
#include "Brinkman_3D/AnsatzLinConst.h"
}

namespace ansatz_quad_lin //0
{
#include "Brinkman_3D/AnsatzQuadLin.h"
}

namespace poiseuille //0
{
#include "Brinkman_3D/Poiseuille.h"
}



Example_Brinkman3D::Example_Brinkman3D(
  const ParameterDatabase& user_input_parameter_db) 
 : Example3D(user_input_parameter_db)
{
    int example_code = this->example_database["example"];
    switch( example_code )
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
        
        //            /** some variables to change values in the example */
        //            ansatz_lin_const::DIMENSIONLESS_VISCOSITY = this->get_nu();
       
       // read parameters from local database
       ansatz_lin_const::viscosity = get_viscosity();
       ansatz_lin_const::effective_viscosity = get_effective_viscosity();
       ansatz_lin_const::permeability = get_permeablity();

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
        
        //            /** some variables to change values in the example */
        //            ansatz_lin_const::DIMENSIONLESS_VISCOSITY = this->get_nu();
       
       // read parameters from local database
       ansatz_quad_lin::viscosity = get_viscosity();
       ansatz_quad_lin::effective_viscosity = get_effective_viscosity();
       ansatz_quad_lin::permeability = get_permeablity();

        ansatz_quad_lin::ExampleFile();
        break;
        
            
        case 2:
            /** exact_solution */
            exact_solution.push_back( poiseuille::ExactU1 );
            exact_solution.push_back( poiseuille::ExactU2 );
            exact_solution.push_back( poiseuille::ExactU3 );
            exact_solution.push_back( poiseuille::ExactP );
            
            /* boundary condition */
            boundary_conditions.push_back( poiseuille::BoundCondition );
            boundary_conditions.push_back( poiseuille::BoundCondition );
            boundary_conditions.push_back( poiseuille::BoundCondition );
            boundary_conditions.push_back( BoundConditionNoBoundCondition );
            
            /* boundary values */
            boundary_data.push_back( poiseuille::U1BoundValue );
            boundary_data.push_back( poiseuille::U2BoundValue );
            boundary_data.push_back( poiseuille::U3BoundValue );
            boundary_data.push_back( BoundaryValueHomogenous );
            
            /* coefficients */
            problem_coefficients = poiseuille::LinCoeffs;
            
            //            /** some variables to change values in the example */
            //            ansatz_lin_const::DIMENSIONLESS_VISCOSITY = this->get_nu();
            
            // read parameters from local database
            poiseuille::viscosity = get_viscosity();
            poiseuille::effective_viscosity = get_effective_viscosity();
            poiseuille::permeability = get_permeablity();

            poiseuille::ExampleFile();
            break;

        default:
        ErrThrow("Unknown Brinkman example!");
        
    }
    
    
    //        case 100:
    //            /** exact_solution */
    //            exact_solution.push_back( poiseuille::ExactU1 );
    //            exact_solution.push_back( poiseuille::ExactU2 );
    //            exact_solution.push_back( poiseuille::ExactU3 );
    //            exact_solution.push_back( poiseuille::ExactP );
    //
    //            /** boundary condition */
    //            boundary_conditions.push_back( poiseuille::BoundCondition );
    //            boundary_conditions.push_back( poiseuille::BoundCondition );
    //            boundary_conditions.push_back( poiseuille::BoundCondition );
    //            boundary_conditions.push_back( BoundConditionNoBoundCondition );
    //
    //            /** boundary values */
    //            boundary_data.push_back( poiseuille::U1BoundValue );
    //            boundary_data.push_back( poiseuille::U2BoundValue );
    //            boundary_data.push_back( poiseuille::U3BoundValue );
    //            boundary_data.push_back( BoundaryValueHomogenous );
    //
    //            /** coefficients */
    //            problem_coefficients = poiseuille::LinCoeffs;
    //
    //            //this->example_database["equal_order_stab_weight"] = 500;
    //            ////poiseuille::stab_weight = this->example_database["equal_order_stab_weight"];
    //            //double equal_order_stab_weight = this->example_database["equal_order_stab_weight"];
    //
    //            poiseuille::ExampleFile();
    //
    //            break;
    //        case 1:
    //            /** exact_solution */
    //            exact_solution.push_back( Poiseuille_Hannukainen::ExactU1 );
    //            exact_solution.push_back( Poiseuille_Hannukainen::ExactU2 );
    //            exact_solution.push_back( Poiseuille_Hannukainen::ExactU3 );
    //            exact_solution.push_back( Poiseuille_Hannukainen::ExactP );
    //
    //            /** boundary condition */
    //            boundary_conditions.push_back( Poiseuille_Hannukainen::BoundCondition );
    //            boundary_conditions.push_back( Poiseuille_Hannukainen::BoundCondition );
    //            boundary_conditions.push_back( Poiseuille_Hannukainen::BoundCondition );
    //            boundary_conditions.push_back( BoundConditionNoBoundCondition );
    //
    //            /** boundary values */
    //            boundary_data.push_back( Poiseuille_Hannukainen::U1BoundValue );
    //            boundary_data.push_back( Poiseuille_Hannukainen::U2BoundValue );
    //            boundary_data.push_back( Poiseuille_Hannukainen::U3BoundValue );
    //            boundary_data.push_back( BoundaryValueHomogenous );
    //
    //            /** coefficients */
    //            problem_coefficients = Poiseuille_Hannukainen::LinCoeffs;
    //
    //            Poiseuille_Hannukainen::ExampleFile();
    //            break;
    //        case 2:
    //            /** exact_solution */
    //            exact_solution.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::ExactU1 );
    //            exact_solution.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::ExactU2 );
    //            exact_solution.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::ExactU3 );
    //            exact_solution.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::ExactP );
    //
    //            /** boundary condition */
    //            boundary_conditions.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::BoundCondition );
    //            boundary_conditions.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::BoundCondition );
    //            boundary_conditions.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::BoundCondition );
    //            boundary_conditions.push_back( BoundConditionNoBoundCondition );
    //
    //            /** boundary values */
    //            boundary_data.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::U1BoundValue );
    //            boundary_data.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::U2BoundValue );
    //            boundary_data.push_back( Poiseuille_Hannukainen_with_inscribed_physical_sphere::U3BoundValue );
    //            boundary_data.push_back( BoundaryValueHomogenous );
    //
    //            /** coefficients */
    //            problem_coefficients = Poiseuille_Hannukainen_with_inscribed_physical_sphere::LinCoeffs;
    //
    //            Poiseuille_Hannukainen_with_inscribed_physical_sphere::ExampleFile();
    //            break;
    //        case 3:
    //            /** exact_solution */
    //            exact_solution.push_back( sine_cosine::ExactU1 );
    //            exact_solution.push_back( sine_cosine::ExactU2 );
    //            exact_solution.push_back( sine_cosine::ExactU3 );
    //            exact_solution.push_back( sine_cosine::ExactP );
    //
    //            /** boundary condition */
    //            boundary_conditions.push_back( sine_cosine::BoundCondition );
    //            boundary_conditions.push_back( sine_cosine::BoundCondition );
    //            boundary_conditions.push_back( sine_cosine::BoundCondition );
    //            boundary_conditions.push_back( BoundConditionNoBoundCondition );
    //
    //            /** boundary values */
    //            boundary_data.push_back( sine_cosine::U1BoundValue );
    //            boundary_data.push_back( sine_cosine::U2BoundValue );
    //            boundary_data.push_back( sine_cosine::U3BoundValue );
    //            boundary_data.push_back( BoundaryValueHomogenous );
    //
    //            /** coefficients */
    //            problem_coefficients = sine_cosine::LinCoeffs;
    //
    //            sine_cosine::ExampleFile();
    //            break;
    //        case 4:
    //            /** exact_solution */
    //            exact_solution.push_back( sine2_sine2::ExactU1 );
    //            exact_solution.push_back( sine2_sine2::ExactU2 );
    //            exact_solution.push_back( sine2_sine2::ExactU3 );
    //            exact_solution.push_back( sine2_sine2::ExactP );
    //
    //            /** boundary condition */
    //            boundary_conditions.push_back( sine2_sine2::BoundCondition );
    //            boundary_conditions.push_back( sine2_sine2::BoundCondition );
    //            boundary_conditions.push_back( sine2_sine2::BoundCondition );
    //            boundary_conditions.push_back( BoundConditionNoBoundCondition );
    //
    //            /** boundary values */
    //            boundary_data.push_back( sine2_sine2::U1BoundValue );
    //            boundary_data.push_back( sine2_sine2::U2BoundValue );
    //            boundary_data.push_back( sine2_sine2::U3BoundValue );
    //            boundary_data.push_back( BoundaryValueHomogenous );
    //
    //            /** coefficients */
    //            problem_coefficients = sine2_sine2::LinCoeffs;
    //
    //            sine_cosine::ExampleFile();
    //            break;
    //
    //            //        case 3:
    //            //            /** exact_solution */
    //            //            exact_solution.push_back( driven_cavity::ExactU1 );
    //            //            exact_solution.push_back( driven_cavity::ExactU2 );
    //            //            exact_solution.push_back( driven_cavity::ExactP );
    //            //
    //            //            /** boundary condition */
    //            //            boundary_conditions.push_back( driven_cavity::BoundCondition );
    //            //            boundary_conditions.push_back( driven_cavity::BoundCondition );
    //            //            boundary_conditions.push_back( BoundConditionNoBoundCondition );
    //            //
    //            //            /** boundary values */
    //            //            boundary_data.push_back( driven_cavity::U1BoundValue );
    //            //            boundary_data.push_back( driven_cavity::U2BoundValue );
    //            //            boundary_data.push_back( BoundaryValueHomogenous );
    //            //
    //            //            /** coefficients */
    //            //            problem_coefficients = driven_cavity::LinCoeffs;
    //            //
    //            //            driven_cavity::ExampleFile();
    //            //            break;
    //            //        case 4:
    //            //            /** exact_solution */
    //            //            exact_solution.push_back( flow_around_cylinder::ExactU1 );
    //            //            exact_solution.push_back( flow_around_cylinder::ExactU2 );
    //            //            exact_solution.push_back( flow_around_cylinder::ExactP );
    //            //
    //            //            /** boundary condition */
    //            //            boundary_conditions.push_back( flow_around_cylinder::BoundCondition );
    //            //            boundary_conditions.push_back( flow_around_cylinder::BoundCondition );
    //            //            boundary_conditions.push_back( BoundConditionNoBoundCondition );
    //            //
    //            //            /** boundary values */
    //            //            boundary_data.push_back( flow_around_cylinder::U1BoundValue );
    //            //            boundary_data.push_back( flow_around_cylinder::U2BoundValue );
    //            //            boundary_data.push_back( BoundaryValueHomogenous );
    //            //
    //            //            /** coefficients */
    //            //            problem_coefficients = flow_around_cylinder::LinCoeffs;
    //            //            
    //            //            flow_around_cylinder::ExampleFile();
    //            //            break;
    //            
    //        default:
    //            ErrThrow("Unknown Brinkman example!");
}


void Example_Brinkman3D::do_post_processing(Brinkman3D& brinkman3d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(brinkman3d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_Brinkman3D","No post processing done for the current example.");
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// functions to avoid global parameters (TDatabase) in the example file

double Example_Brinkman3D::get_viscosity() const
{
  double viscosity = this->example_database["viscosity"];
  return viscosity;
}
double Example_Brinkman3D::get_effective_viscosity() const
{
  double effective_viscosity = this->example_database["effective_viscosity"];
  return effective_viscosity;
}
double Example_Brinkman3D::get_permeablity() const
{
  double K = this->example_database["permeability"];
  /*if (this->example_database["read_permeability_from_file"])
    K = -1;
   */
  return K;
}

//double Example_Brinkman3D::get_stab() const
//{
//  double equal_order_stab_weight = this->example_database["equal_order_stab_weight"];
//  return equal_order_stab_weight;
//}
