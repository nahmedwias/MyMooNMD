/**
 * @brief A test program for solving coupled (in 1 direction) Brinkman2D-CD2D problems.
 *
 * This serves as a test for solving Thermohydraulic_Brinkman2D problems. It is intended to
 * perform calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes wrong)
 * the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then this tests
 * shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 *
 * @date 12/09/2018
 * @author Laura Blank
 *
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Chrono.h>
#include <LocalAssembling2D.h>
#include <MainUtilities.h> //for error measuring

#include <Brinkman2D.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <Example_Brinkman2D.h>
#include <BoundaryAssembling2D.h>


#include <CD2D.h>
#include <Example_CD2D.h>
#include <AlgebraicFluxCorrection.h>
#include <Multigrid.h>
#include <algorithm>



//=======================================================================
// compare the computed bulk and boundary errors (apparent for Nitsche approach) 
// in the Brinkman2D object with the given ones in the array
void compareAllErrors(const Brinkman2D& brinkman2d, std::array<double, 8> reference_errors)
{
  const double eps = 1e-11; //2e-9;

  Output::print(setprecision(14),reference_errors[0]);
  Output::print(setprecision(14),reference_errors[1]);
  Output::print(setprecision(14),reference_errors[2]);
  Output::print(setprecision(14),reference_errors[3]);
  Output::print(setprecision(14),reference_errors[4]);
  Output::print(setprecision(14),reference_errors[5]); //boundary error of u
  Output::print(setprecision(14),reference_errors[6]); // boundary error of u.n


  // check the errors
  if( fabs(brinkman2d.getL2VelocityError() - reference_errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 velocity error not correct. ",
        brinkman2d.getL2VelocityError(), " and ",  reference_errors[0]);
  }
  if( fabs(brinkman2d.getL2DivergenceError() - reference_errors[1]) > eps)
  {
    ErrThrow("Program 1: L2 divergence velocity error not correct. ",
        brinkman2d.getL2DivergenceError() - reference_errors[1]);
  }
  if( fabs(brinkman2d.getH1SemiVelocityError() - reference_errors[2]) > eps)
  {
    ErrThrow("Program 1: H1-semi velocity error not correct. ",
        brinkman2d.getH1SemiVelocityError() - reference_errors[2]);
  }
 
  if( fabs(brinkman2d.getL2BoundaryError() - reference_errors[5]) > eps)
  {
    ErrThrow("Program 1: L2 velocity error at the boundary is not correct. ",
        brinkman2d.getL2BoundaryError() - reference_errors[5]);
  }
 if( fabs(brinkman2d.getL2NormNormalComponentError() - reference_errors[6]) > eps)
  {
    ErrThrow("Program 1: L2 error of the normal velocity at the boundayr is not correct. ",
        brinkman2d.getL2NormNormalComponentError() - reference_errors[6]);
  }

 
 if( fabs(brinkman2d.getL2PressureError() - reference_errors[3]) > eps )
  {
    ErrThrow("Program 1: L2 pressure error not correct.",
        brinkman2d.getL2PressureError() - reference_errors[3]);
  }
  if( fabs(brinkman2d.getH1SemiPressureError() - reference_errors[4]) > eps )
  {
    ErrThrow("Program 1: H1-semi pressure error not correct.",
        brinkman2d.getH1SemiPressureError() - reference_errors[4]);
  }
}

//=======================================================================
// compare the computed bulk errors of the CD2D solution with the given ones in
// the array
void compareErrors(const CD2D& cd2d, std::array <double, 4> reference_errors)
{
  const double eps = 1e-11;
  
  // check the errors
  if( fabs(cd2d.get_L2_error() - reference_errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 error not correct. ",
             cd2d.get_L2_error() - reference_errors[0]);
  }
  if( fabs(cd2d.get_H1_semi_error() - reference_errors[1]) > eps )
  {
    ErrThrow("Program 1: H1-semi error not correct. ",
             cd2d.get_H1_semi_error() - reference_errors[1]);
  }
  if( fabs(cd2d.get_SD_error() - reference_errors[2]) > eps )
  {
    ErrThrow("Program 1: sd error not correct.",
             cd2d.get_SD_error() - reference_errors[2]);
  }
  if( fabs(cd2d.get_L_inf_error() - reference_errors[3]) > eps )
  {
    ErrThrow("Program 1: L_inf error not correct.",
             cd2d.get_L_inf_error() - reference_errors[3]);
  }
}


//=======================================================================
// Here the actual computations for Brinkman2D take place.
// This version includes boundary errors.
void check_brinkman2d(TDomain & domain, ParameterDatabase& db, int velocityCode,int pressureCode,
    std::array<double, 8> reference_errors, unsigned int nRefinements)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocityCode;
  TDatabase::ParamDB->PRESSURE_SPACE = pressureCode;

  Output::print("FEDATABASE");
  Output::set_outfile("Test.out");
  TDatabase::WriteParamDB((char*)"FEDATA");

  Brinkman2D brinkman2d(domain, db);
  brinkman2d.assemble();
  brinkman2d.solve();
  brinkman2d.output(nRefinements);

  // compare computed with given errors
  compareAllErrors(brinkman2d, reference_errors); 
}

//=======================================================================
// Here the actual computations for CD2D take place
void check_cd2d(TDomain & domain, ParameterDatabase& db, int element_code, 
                std::array<double, 4> errors)
{
  TDatabase::ParamDB->ANSATZ_ORDER = element_code;
  
  CD2D cd2d(domain, db);
  cd2d.assemble();
  cd2d.solve();
  cd2d.output();
  // compare computed with given errors
  compareErrors(cd2d, errors); // throws upon a difference
}


//=======================================================================
//============ EXAMPLE 10 - Geothermal_Energy_Brinkman2D.h ==============
//=======================================================================

void tests_on_triangles_P1P1_Example10(unsigned int nRefinements, ParameterDatabase& db)
{ 

db["geo_file"].set("/Users/blank/ParMooN/ParMooN/data/mesh/geothermal2d.mesh", false);
db["boundary_file"].set("/Users/blank/ParMooN/ParMooN/data/mesh/geothermal2d.PRM", false);

// default construct a domain object
  TDomain domain(db);

  // Initialization of the default parameters
  TDatabase::SetDefaultParameters(); 

  // refine grid up to the finest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }


db["problem_type"] = 3;

db["verbosity"] = 2;

    std::array<double, 8>  reference_errors_brinkman2d;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 10; // Geothermal_Energy_Brinkman2D.h
  db["permeability"] = 0.0000001;
  db["effective_viscosity"] = 0.;
  db["viscosity"] = 0.001;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.5;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.5;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.01;

  db["corner_stab_weight"] = 0.1;

  db["coefficient_function_type"] = 3;
  db["use_source_sink_function"] = true;

  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  
  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {0., 0.};

  TDatabase::ParamDB->n_nitsche_boundary = 2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  //TDatabase::ParamDB->l_T = -1;
  //--------------------------------------------------------------------------------------------------------------------
  
  Output::print("\nstarting with Brinkman2D on TwoTriangles (non-symmetric Galerkin formulation with non-symmetric GLS stabilization, Grad-Div Stabilization, non-symmetric Nitsche method ), Example 10, P1/P1, with Nitsche and Neumann bcs and with visc_eff 0, visc = 0.001, perm = 0.0000001");

  reference_errors_brinkman2d = {{0.0024036452085024, 0.044524856675523, 0.044953295857573, 8.7256326558545, 24.135888171285, 
                       0.00015015057080671, 6.3264594781094e-05}};
  check_brinkman2d(domain, db, 1,1, reference_errors_brinkman2d, nRefinements);


//===================  CD2D  =====================
  // Re-Initialization of the default parameters
  TDatabase::SetDefaultParameters(); 

  // refine grid up to the finest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }

  std::array<double, 4> reference_errors_cd2d;

db["problem_type"] = 1;

db["verbosity"] = 2;

  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 4; // Geothermal_Energy_CD2D.h

db["space_discretization_type"] = "galerkin";  //supg;  //

//Peclet number eps=1/PE_NR (diffusion coefficient)
db["diffusion_coefficient"]=  0.001;
///TDatabase::ParamDB->DISCTYPE = 1;
TDatabase::ParamDB->RE_NR = 1;
TDatabase::ParamDB->PE_NR = 1000;

  db["coefficient_function_type"] = 4;
  db["use_source_sink_function"] = true;

  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;

//change how the SUPG parameter is computed
TDatabase::ParamDB->SDFEM_TYPE = 0;

// factors for SDFEM (delta_K=h_K*DELTAi)
//   DELTA0 for high Peclet number
//   DELTA1 for low Peclet number
TDatabase::ParamDB->DELTA0 =  0.1;
TDatabase::ParamDB->DELTA1 =  0.0;

// LPS, stabilization parameters
// choose ANSATZ_ORDER 100, 201, 302, 403, 504,... 
TDatabase::ParamDB->LP_FULL_GRADIENT = 1;
TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0.5;
TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;
TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = 1;

TDatabase::ParamDB->LP_STREAMLINE = 0;
TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0.1;
TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = 1;
TDatabase::ParamDB->Axial3D = 0;

  
  //---------------------------------------------------------------------------------


  Output::print("\nstarting with P2");
  reference_errors_cd2d = {{ 1073.5001901761, 125.37410635258, 3.9646773568233, 149.99925787946 }};
  check_cd2d(domain, db, 2, reference_errors_cd2d);

  
  if(db["solver_type"].is("iterative") && db["preconditioner"].is("multigrid"))
  {
    // the higher order elements do not/very slowly converge.
    // @todo find out why multigrid does not allow high order elements for CD2D
    return;
  }

}


// ========================================================================
// =======================================================================
// main program
// =======================================================================
// ========================================================================
int main(int argc, char* argv[])
{
  //  declaration of databases
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  // high order quadrature for computing errors
  // TDatabase::ParamDB->INPUT_QUAD_RULE = 99;

  unsigned int nRefinements = 0;
  Output::setVerbosity(2);

  //---------------------------------------------------------------------
  // merge TDatabase with Problem specific ParameterDatabase db
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(ParameterDatabase::default_output_database(), true);
  db.merge(Example2D::default_example_database(), true);
  
  db.merge(Brinkman2D::get_default_Brinkman2D_parameters(),true);
  db.merge(CD2D::get_default_CD2D_parameters(),true);

  db.merge(TDomain::default_domain_parameters(),true);

  db["output_compute_errors"] = true;
  db["output_write_vtk"] = false;
 
 // db.add("geo_file", "TwoTriangles", "", {});
 // db.add("boundary_file", "Default_UnitSquare", "");

  nRefinements = 0;

//  tests_on_triangles_P1P1_Example10(nRefinements, db);

  

/*// Tests including boundary errors
nRefinements = 7;
tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example8(nRefinements, db);
tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(nRefinements, db);
tests_on_triangles_P2P2_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(nRefinements, db);
*/

return 0;
}
