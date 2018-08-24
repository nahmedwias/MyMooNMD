/**
 * @brief A test program for solving 2D Brinkman problems.
 *
 * This serves as a test for solving Brinkman problems in 2D. It is intended to
 * perform Brinkman calculations with different examples in different setups to test
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
 * @date 25/07/2016
 * @author Laura Blank
 *
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Brinkman2D.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <Example_Brinkman2D.h>
#include "LocalAssembling.h"
#include <BoundaryAssembling2D.h>
#include <MainUtilities.h> //for error measuring
#include <Chrono.h>
#include <algorithm>


//----------------------------------------------------------------------------
// compare the computed errors in the Brinkman2D object with the given ones in
// the array
void compareErrors(const Brinkman2D& brinkman2d, std::array<double, 5> reference_errors)
{
  const double eps = 2e-9;

  Output::print(setprecision(14),reference_errors[0]);
  Output::print(setprecision(14),reference_errors[1]);
  Output::print(setprecision(14),reference_errors[2]);
  Output::print(setprecision(14),reference_errors[3]);
  Output::print(setprecision(14),reference_errors[4]);

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
//----------------------------------------------------------------------------
// compare the computed boundary errors in the Brinkman2D object with the given ones in
// the array
void compareAllErrors(const Brinkman2D& brinkman2d, std::array<double, 8> reference_errors)
{
  const double eps = 2e-9;

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


//----------------------------------------------------------------------
// Here the actual computations take place
void check_brinkman2d(TDomain & domain, ParameterDatabase& db, int velocityCode,int pressureCode,
    std::array<double, 5> reference_errors, unsigned int nRefinements)
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
  compareErrors(brinkman2d, reference_errors); 
}
//----------------------------------------------------------------------
// Here the actual computations take place
// This version includes boundary errors
void check_brinkman2d_New(TDomain & domain, ParameterDatabase& db, int velocityCode,int pressureCode,
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
////======================================================================================================
//void tests_on_quads(unsigned int nRefinements, ParameterDatabase& db)
//{
//    std::array<double, 4> reference_errors;
//    
//    // default construct a domain object
//    TDomain domain(db);
//    
//    // Initialization of the default parameters
//    TDatabase::SetDefaultParameters();
//    
//    // refine grid up to the coarsest level
//    for(unsigned int i = 0; i < nRefinements; i++)
//    {
//        domain.RegRefineAll();
//    }
//    
//    //------------------------------------------------------------------------------------------------------------
//    Output::print("\nStarting with P2/P1 on quads for Poiseuille flow with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
//    reference_errors = {{ 2.0257088643869e-16, 2.1906901043565e-15, 5.9796305144209e-15, 4.8385464788293e-14 }};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    
//    //    //------------------------------------------------------------------------------------------------------------
//    //    db["effective_viscosity"]=0.01;
//    //    db["viscosity"]=10;
//    //    db["permeability"]=2;
//    //
//    //    Output::print("\nStarting with P2/P1 on quads for Poiseuille flow with visc_eff=0.01, visc=10, perm=2");
//    //    reference_errors = {{1.9508695699745e-16,3.6646227316896e-15,1.3277195222811e-16,1.0734187033819e-15}};
//    //
//    //    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    //
//    //    db["effective_viscosity"]=1;
//    //    db["viscosity"]=1;
//    //    db["permeability"]=1;
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
//    reference_errors = {{2.2608321672551e-05, 0.00059223693639448, 9.7288545102352e-06,9.0848530457398e-05}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    db["example"] = 0; // Poiseuille
//    
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 1; // Poiseuille_Hannukainen
//    db["effective_viscosity"]=0.001;
//    db["viscosity"]=1;
//    db["permeability"]=1;
//    
//    Output::print("\nStarting with P2/P1 on quads for Poiseuille flow with visc_eff=0.001, visc=1, perm=2 for example ",db["example"] );
//    reference_errors = {{0.097304080869912,2.9931806337455,0.00323505572483,0.027631651143564}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    
//    db["effective_viscosity"]=1;
//    db["viscosity"]=1;
//    db["permeability"]=1;
//    db["example"] = 0; // Poiseuille
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 2; // Poiseuille_Hannukainen with inscribed sphere
//    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
//    reference_errors = {{0.81900218507518,11.068008588367,81.992399830346, 294.15272705334}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    db["example"] = 0; // Poiseuille
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 3; // sincos
//    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
//    reference_errors = {{0.0040823711607968, 0.10522135481493,0.017649828876871, 0.51271046930702}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    db["example"] = 0; // Poiseuille
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 4; // sincos2
//    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"]);
//    reference_errors = {{0.0466869785743,1.2606661474721, 0.15644878749027,2.5764786606906}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    db["example"] = 0; // Poiseuille
//}
//
////======================================================================================================
//void tests_on_quads_stab(unsigned int nRefinements, ParameterDatabase& db)
//{
//    std::array<double, 4>  reference_errors;
//    
//    // default construct a domain object
//    TDomain domain(db);
//    
//    // Initialization of the default parameters
//    TDatabase::SetDefaultParameters();
//    
//    // refine grid up to the coarsest level
//    for(unsigned int i = 0; i < nRefinements; i++)
//    {
//        domain.RegRefineAll();
//    }
//    
//    //------------------------------------------------------------------------------------------------------------
//    //P1P1
//    Output::print("\nStarting with residual-based equal-order stabilization for P1/P1 (parameter=0.4) on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"]);
//    reference_errors = {{ 0.20649934279357,0.97728617227336, 1.3352827196758, 5.1734925650333}};
//    
//    db["P1P1_stab"] = true;
//    TDatabase::ParamDB->equal_order_stab_weight_PkPk=0.4;
//    
//    check_brinkman2d(domain, db, 1,1, reference_errors, nRefinements);
//    
//    db["P1P1_stab"] = false;
//    
//    //------------------------------------------------------------------------------------------------------------
//    //P2P2
//    Output::print("\nStarting with residual-based equal-order stabilization for P2/P2 (parameter=0.4) on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"]);
//    reference_errors = {{0.11209308296508,1.6714743655252, 0.88411612007029, 10.192320142372}};
//    
//    db["P2P2_stab"] = true;
//    TDatabase::ParamDB->equal_order_stab_weight_PkPk=0.4;
//    
//    check_brinkman2d(domain, db, 2,2, reference_errors, nRefinements);
//    
//    db["P2P2_stab"] = false;
//    
//    //------------------------------------------------------------------------------------------------------------
//    //P2P2
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nStarting with residual-based equal-order stabilization for P2/P2 (parameter=0.4) on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"]);
//    reference_errors = {{0.012909128829698,0.19363961444426,0.10019110454835,1.1661459244268}};
//    
//    db["P2P2_stab"] = true;
//    TDatabase::ParamDB->equal_order_stab_weight_PkPk=0.4;
//    
//    check_brinkman2d(domain, db, 2,2, reference_errors, nRefinements);
//    
//    db["P2P2_stab"] = false;
//    db["example"] = 0; // Poiseuille
//    
//}
//
////======================================================================================================
//void tests_on_quads_Nitsche(unsigned int nRefinements, ParameterDatabase& db)
//{
//    // default construct a domain object
//    TDomain domain(db);
//    
//    // Initialization of the default parameters
//    TDatabase::SetDefaultParameters();
//    
//    std::array<double, 4>  reference_errors;
//    
//    // refine grid up to the coarsest level
//    for(unsigned int i = 0; i < nRefinements; i++)
//    {
//        domain.RegRefineAll();
//    }
//    //--------------------------------------------------------------------------------------------------------------------
//    
//    TDatabase::ParamDB->n_nitsche_boundary=4;
//    TDatabase::ParamDB->nitsche_boundary_id={0, 1, 2, 3};
//    TDatabase::ParamDB->nitsche_penalty={1000, 1000, 1000, 1000};
//    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P1 on quads with visc_eff=1, visc=1, perm= 1");
//    reference_errors = {{0.00071104275913606,0.0086539449600541,0.0137946704876,0.11201103666807}};
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    
//    
//    //--------------------------------------------------------------------------------------------------------------------
//    TDatabase::ParamDB->n_nitsche_boundary=2;
//    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
//    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
//    TDatabase::ParamDB->n_neumann_boundary=2;
//    TDatabase::ParamDB->neumann_boundary_id={1,3};
//    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
//    Output::print("\nstarting with Nitsche (parameter 1000)  on P2/P1 on quads with Neumann and visc_eff=1, visc=1, perm= 1");
//    reference_errors = {{0.0018492680116681, 0.00052482584609959, 0.5, 2.0166230982544e-14}};
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    
//    
//    //--------------------------------------------------------------------------------------------------------------------
//    TDatabase::ParamDB->n_nitsche_boundary=2;
//    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
//    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
//    TDatabase::ParamDB->n_neumann_boundary=2;
//    TDatabase::ParamDB->neumann_boundary_id={1,3};
//    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
//    db["effective_viscosity"]=0.001;
//    db["viscosity"]=10;
//    db["permeability"]=1;
//    Output::print("\nstarting with Nitsche (parameter 1000)  on P2/P1 on quads with Neumann and visc_eff=0.001, visc=10, perm= 1");
//    reference_errors = {{4.8144012875101e-05,0.0010541202262723, 0.5,1.184219258062e-15}};
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    
//    //--------------------------------------------------------------------------------------------------------------------
//    TDatabase::ParamDB->n_nitsche_boundary=2;
//    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
//    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
//    TDatabase::ParamDB->n_neumann_boundary=2;
//    TDatabase::ParamDB->neumann_boundary_id={1,3};
//    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
//    db["effective_viscosity"]=0.001;
//    db["viscosity"]=10;
//    db["permeability"]=1;
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P1 on quads with Neumann and visc_eff=0.001, visc=10, perm= 1");
//    reference_errors = {{ 0.19394824424421,3.8205668689833,0.5,1.0532389294782e-15}};
//    check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
//    
//    db["example"] = 0; // Poiseuille
//    
//    //--------------------------------------------------------------------------------------------------------------------
//    TDatabase::ParamDB->n_nitsche_boundary=2;
//    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
//    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
//    TDatabase::ParamDB->n_neumann_boundary=2;
//    TDatabase::ParamDB->neumann_boundary_id={1,3};
//    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
//    db["effective_viscosity"]=0.001;
//    db["viscosity"]=10;
//    db["permeability"]=1;
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P2 on quads with Neumann and visc_eff=0.001, visc=10, perm= 1");
//    reference_errors = {{0.19394824424421,3.8205668689833,2.8822850639982,87.950352486737}};
//    check_brinkman2d(domain, db, 2,2, reference_errors, nRefinements);
//    
//    db["example"] = 0; // Poiseuille
//    
//    //--------------------------------------------------------------------------------------------------------------------
//    TDatabase::ParamDB->n_nitsche_boundary=2;
//    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
//    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
//    TDatabase::ParamDB->n_neumann_boundary=2;
//    TDatabase::ParamDB->neumann_boundary_id={1,3};
//    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
//    db["effective_viscosity"]=0.001;
//    db["viscosity"]=10;
//    db["permeability"]=1;
//    db["PkPk_stab"]=true;
//    TDatabase::ParamDB->equal_order_stab_weight_PkPk=0.4;
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nstarting with Nitsche (parameter 1000) on P1/P1 on quads with Stab (parameter=0.4), Neumann and visc_eff=0.001, visc=10, perm= 1 ");
//    reference_errors = {{0.3582641143964,3.3275876201918,0.5119324115733,1.2550049173583}};
//    check_brinkman2d(domain, db, 1,1, reference_errors, nRefinements);
//    
//    db["example"] = 0; // Poiseuille
//    db["PkPk_stab"]=false;
//    
//    //--------------------------------------------------------------------------------------------------------------------
//    TDatabase::ParamDB->n_nitsche_boundary=2;
//    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
//    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
//    TDatabase::ParamDB->n_neumann_boundary=2;
//    TDatabase::ParamDB->neumann_boundary_id={1,3};
//    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
//    db["effective_viscosity"]=0.001;
//    db["viscosity"]=10;
//    db["permeability"]=1;
//    db["PkPk_stab"]=true;
//    TDatabase::ParamDB->equal_order_stab_weight_PkPk=0.4;
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P2 on quads with Stab (parameter=0.4), Neumann and visc_eff=0.001, visc=10, perm= 1");
//    reference_errors = {{2.7665723232465,38.383469514259,3.8290546915193,26.920248715798}};
//    check_brinkman2d(domain, db, 2,2, reference_errors, nRefinements);
//    
//    db["example"] = 0; // Poiseuille
//    db["P2P2_stab"]=false;
//}




//################# EXAMPLE 1 - Exponential Flow (Hannukainen & Co) ####################################
//======================================================================================================
void tests_on_triangles_P2P1_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ // default construct a domain object
  TDomain domain(db);
  // Initialization of the default parameters
  TDatabase::SetDefaultParameters(); 
  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
   
  db["viscosity"] = 0.004;
  db["PkPk_stab"] = false;
  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "symmetric Galerkin formulation";
  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};
  db["GradDiv_stab"] = false;

  TDatabase::ParamDB->n_nitsche_boundary = 0; //2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles (symmetric Galerkin formulation), Example 1, P2/P1, with Dirichlet and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{0.0054680957600273, 0.014317962190331, 0.1443805008035, 5.246987770953e-05, 0.00061083280165159}};
  check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);

  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  Output::print("\nstarting with Brinkman2D on TwoTriangles (nonsymmetric Galerkin formulation), Example 1, P2/P1, with Dirichlet and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{0.0054680957600273, 0.014317962190331, 0.1443805008035, 5.246987770953e-05, 0.00061083280165159}};
  check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
}

//======================================================================================================
void tests_on_triangles_P2P1_PenaltyFreeNonSymmetricNitsche_Example1(unsigned int nRefinements, ParameterDatabase& db)
{
  TDomain domain(db);
  TDatabase::SetDefaultParameters();

  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;
  db["PkPk_stab"] = false;
  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["GradDiv_stab"] = false;

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1,3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2; //2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 1, P2/P1, with penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{0.049354569383435, 0.049746071781475 , 0.20428824240268, 9.3291684594421e-05, 0.0010263458082476}};
  check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_Example1(unsigned int nRefinements, ParameterDatabase& db)
{
  TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.01;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.01;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by h_T";
  db["GradDiv_stab"] = false;

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 0; //2;                     
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (symmetric GLS), Dirichlet and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{  1.2424692993215, 0.031419104874062, 16.668762604057, 0.0018335349292131, 0.034019179940349 }}; 
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (nonsymmetric GLS), Dirichlet and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{ 1.2347001456568, 1.3166645102027, 16.750021539147,0.018674085154587, 0.22010507412489}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ 
  TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.01;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.01;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by h_T";
  db["GradDiv_stab"] = false;

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;                     
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (symmetric GLS), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{10.192177714504, 5.8853778774783, 25.146773888043, 0.16699522666312, 3.1493307465155 }};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";

  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (nonsymmetric GLS), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{9.599443219491, 4.1492978030163,  23.384929664058, 0.022130923948982, 0.27177748647831}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab100_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.01;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.01;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 100;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by h_T";


  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1,3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;                     
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{ 10.869461373038, 0.0015555188075453, 24.270540284686, 0.1885807925005, 2.3578866662985 }};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (nonsymmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{10.581773844742, 0.29752211831635, 24.022855667394, 0.10196998348134, 0.59373617246339}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{10.666137970649, 1.2062306735909, 24.166125351829, 0.21575466767375, 3.3261411856147}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (nonsymmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{9.8125894591994, 8.617299671234, 25.051206032843, 0.045806976649214, 0.32576244958451}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}



//############################## EXAMPLE 8 - SinCos_BadiaCodina ###########################################
//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example8(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 8; // Poiseuille_Hannukainen
  db["permeability"] = 1;
  db["effective_viscosity"] = 0.;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  TDatabase::ParamDB->n_neumann_boundary = 0;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 4;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 1, 2, 3};
  TDatabase::ParamDB->nitsche_penalty = {0, 0, 0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 1");
  reference_errors = {{5.967728815198, 27.218645176007, 45.026548460105, 2.5622671365087, 33.995673079637}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (nonsymmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 1");
  reference_errors = {{2.3166504721005, 24.851364206268, 27.233327688498, 0.00084168519671418, 0.010828366364366}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_smallK_Example8(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 8; // Poiseuille_Hannukainen
  db["permeability"] = 0.00001;
  db["effective_viscosity"] = 0.;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  TDatabase::ParamDB->n_neumann_boundary = 0;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 4;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 1, 2, 3};
  TDatabase::ParamDB->nitsche_penalty = {0, 0, 0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  TDatabase::ParamDB->l_T = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{2.3166504721003, 24.851364206268, 27.233327688498, 84.168519671385, 1082.8366364365}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001 ");
  reference_errors = {{5.6528468593088, 33.609960101853, 43.609769730799, 362.86296558258, 2403.2148842124}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{1.3167265706952, 21.848230410711, 28.026777395638, 93.875214873957, 2044.9301910272}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{5.6528468593088, 33.609960101853, 43.609769730799, 362.86296558258, 2403.2148842124}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  //l_T=1
  TDatabase::ParamDB->l_T = 1;
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";


  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{2.3166504721003, 24.851364206268, 27.233327688498, 84.168519671385, 1082.8366364365}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{5.6528468593088, 33.609960101853, 43.609769730799, 362.86296558258, 2403.2148842124}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{1.3167265706952, 21.848230410711, 28.026777395638, 93.875214873957, 2044.9301910272}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{5.6528468593088, 33.609960101853, 43.609769730799, 362.86296558258, 2403.2148842124}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 



  db["permeability"] = 0.01;
  //l_T=-1
  TDatabase::ParamDB->l_T = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01 ");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{1.3167265706952, 21.848230410711, 28.026777395638, 0.093875214873957, 2.0449301910272}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  //l_T=1
  TDatabase::ParamDB->l_T = 1;
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";


  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{2.3166504721003, 24.851364206268, 27.233327688498, 0.084168519671385, 1.0828366364365}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{1.3167265706952, 21.848230410711, 28.026777395638, 0.093875214873957, 2.0449301910272}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 


}

// includes boundary errors
void tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example8(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 8>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 8; // SinCos BadiaCodina
  db["permeability"] = 0.001;
  db["effective_viscosity"] = 0.;
  db["viscosity"] = 1;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  db["corner_stab_weight"] = 1;

  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  
  TDatabase::ParamDB->n_neumann_boundary = 0;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 4;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 1, 2, 3};
  TDatabase::ParamDB->nitsche_penalty = {0, 0, 0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  //TDatabase::ParamDB->l_T = -1;

 Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS) (0.1), Grad-Div stab (0.1), corner stab (1), scaling by L_0 (0.1), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 1, perm = 0.001");
  reference_errors = {{0.0009227556862607, 0.68537870569107, 0.97251316225737, 0.10670734510353, 110.50582765597, 0.0049705167016144, 0.0023616321724773}};
  check_brinkman2d_New(domain, db, 1, 1, reference_errors, nRefinements); 

}

// includes boundary errors
void tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 8>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 10;
  db["effective_viscosity"] = 0.00001;
  db["viscosity"] = 1;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  db["corner_stab_weight"] = 0;

  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  
  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  //TDatabase::ParamDB->l_T = -1;

 Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS) (0.1), Grad-Div stab (0.1), corner stab (1), scaling by L_0 (0.1), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 1, perm = 0.001");
  reference_errors = {{0.26758348461979, 0.27660180763361, 28.548747008085, 3.0184064668482e-5, 0.012638134363589, 2.708449611399, 0.0065121848533241}};
  check_brinkman2d_New(domain, db, 1, 1, reference_errors, nRefinements); 

}

// includes boundary errors
void tests_on_triangles_P2P2_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 8>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 10;
  db["effective_viscosity"] = 0.00001;
  db["viscosity"] = 1;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  db["corner_stab_weight"] = 0;

  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  
  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  //TDatabase::ParamDB->l_T = -1;

 Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS) (0.1), Grad-Div stab (0.1), corner stab (1), scaling by L_0 (0.1), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 1, perm = 0.001");
  reference_errors = {{0.014505696395323, 0.050048815242381, 2.8694217752035, 3.2893923393833e-6, 0.0029729670540557, 0.10649331439528, 0.0015928285252614}};
  check_brinkman2d_New(domain, db, 2, 2, reference_errors, nRefinements); 

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

  unsigned int nRefinements = 2;
  Output::setVerbosity(2);

  //---------------------------------------------------------------------
  // merge TDatabase with Problem specific ParameterDatabase db
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(ParameterDatabase::default_output_database(),true);
  db.merge(Example2D::default_example_database(),true);
  
  db.merge(Brinkman2D::get_default_Brinkman2D_parameters(),true);

  db["example"] = 1;
  db["output_compute_errors"] = true;
  db["output_write_vtk"] = false;
  //db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "TwoTriangles", "", {"UnitSquare","Default_UnitSquare", "TwoTriangles"});
  //db.add("P1P1_stab", (bool) false, "" );
  //db.add("PkPk_stab", (bool) false, "", {true, false} );
  //db.add("equal_order_stab_weight_PkPk", (double) 0., "", (double) -1000, (double) 1000 );
  //db.add("refinement_n_initial_steps", (size_t) 2.0 , "", (size_t) 0, (size_t) 10000);

  //----------------------------------------
  //tests_on_quads(nRefinements, db);

  //----------------------------------------
  // tests_on_quads_stab(nRefinements, db);

  //----------------------------------------
  //tests_on_quads_Nitsche(nRefinements, db);

  tests_on_triangles_P2P1_Example1(nRefinements, db);

  tests_on_triangles_P2P1_PenaltyFreeNonSymmetricNitsche_Example1(nRefinements, db);

  tests_on_triangles_P1P1_GLSStab_Example1(nRefinements, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_Example1(nRefinements, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab100_Example1(nRefinements, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(nRefinements, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example8(nRefinements, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_smallK_Example8(nRefinements, db);

// Tests including boundary errors
nRefinements = 7;
tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example8(nRefinements, db);
tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(nRefinements, db);
tests_on_triangles_P2P2_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(nRefinements, db);

return 0;
}
