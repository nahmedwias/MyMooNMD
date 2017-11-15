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
#include <LocalAssembling2D.h>
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

//----------------------------------------------------------------------
// Here the actual computations take place
void check_brinkman2d(TDomain & domain, ParameterDatabase& db, int velocityCode,int pressureCode,
    std::array<double, 5> reference_errors)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocityCode;
  TDatabase::ParamDB->PRESSURE_SPACE = pressureCode;

  Output::print("FEDATABASE");
  Output::set_outfile("Test.out");
  TDatabase::WriteParamDB("FEDATA");

  Brinkman2D brinkman2d(domain, db);
  brinkman2d.assemble();
  brinkman2d.solve();
  brinkman2d.output();

  // compare computed with given errors
  compareErrors(brinkman2d, reference_errors); 
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
//    check_brinkman2d(domain, db, 2,1, reference_errors);
//    
//    //    //------------------------------------------------------------------------------------------------------------
//    //    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.01;
//    //    TDatabase::ParamDB->VISCOSITY=10;
//    //    TDatabase::ParamDB->PERMEABILITY=2;
//    //
//    //    Output::print("\nStarting with P2/P1 on quads for Poiseuille flow with visc_eff=0.01, visc=10, perm=2");
//    //    reference_errors = {{1.9508695699745e-16,3.6646227316896e-15,1.3277195222811e-16,1.0734187033819e-15}};
//    //
//    //    check_brinkman2d(domain, db, 2,1, reference_errors);
//    //
//    //    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=1;
//    //    TDatabase::ParamDB->VISCOSITY=1;
//    //    TDatabase::ParamDB->PERMEABILITY=1;
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
//    reference_errors = {{2.2608321672551e-05, 0.00059223693639448, 9.7288545102352e-06,9.0848530457398e-05}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors);
//    db["example"] = 0; // Poiseuille
//    
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 1; // Poiseuille_Hannukainen
//    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
//    TDatabase::ParamDB->VISCOSITY=1;
//    TDatabase::ParamDB->PERMEABILITY=1;
//    
//    Output::print("\nStarting with P2/P1 on quads for Poiseuille flow with visc_eff=0.001, visc=1, perm=2 for example ",db["example"] );
//    reference_errors = {{0.097304080869912,2.9931806337455,0.00323505572483,0.027631651143564}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors);
//    
//    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=1;
//    TDatabase::ParamDB->VISCOSITY=1;
//    TDatabase::ParamDB->PERMEABILITY=1;
//    db["example"] = 0; // Poiseuille
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 2; // Poiseuille_Hannukainen with inscribed sphere
//    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
//    reference_errors = {{0.81900218507518,11.068008588367,81.992399830346, 294.15272705334}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors);
//    db["example"] = 0; // Poiseuille
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 3; // sincos
//    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
//    reference_errors = {{0.0040823711607968, 0.10522135481493,0.017649828876871, 0.51271046930702}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors);
//    db["example"] = 0; // Poiseuille
//    
//    //------------------------------------------------------------------------------------------------------------
//    db["example"] = 4; // sincos2
//    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"]);
//    reference_errors = {{0.0466869785743,1.2606661474721, 0.15644878749027,2.5764786606906}};
//    
//    check_brinkman2d(domain, db, 2,1, reference_errors);
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
//    check_brinkman2d(domain, db, 1,1, reference_errors);
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
//    check_brinkman2d(domain, db, 2,2, reference_errors);
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
//    check_brinkman2d(domain, db, 2,2, reference_errors);
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
//    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
//    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P1 on quads with visc_eff=1, visc=1, perm= 1");
//    reference_errors = {{0.00071104275913606,0.0086539449600541,0.0137946704876,0.11201103666807}};
//    check_brinkman2d(domain, db, 2,1, reference_errors);
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
//    check_brinkman2d(domain, db, 2,1, reference_errors);
//    
//    
//    //--------------------------------------------------------------------------------------------------------------------
//    TDatabase::ParamDB->n_nitsche_boundary=2;
//    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
//    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
//    TDatabase::ParamDB->n_neumann_boundary=2;
//    TDatabase::ParamDB->neumann_boundary_id={1,3};
//    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
//    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
//    TDatabase::ParamDB->VISCOSITY=10;
//    TDatabase::ParamDB->PERMEABILITY=1;
//    Output::print("\nstarting with Nitsche (parameter 1000)  on P2/P1 on quads with Neumann and visc_eff=0.001, visc=10, perm= 1");
//    reference_errors = {{4.8144012875101e-05,0.0010541202262723, 0.5,1.184219258062e-15}};
//    check_brinkman2d(domain, db, 2,1, reference_errors);
//    
//    //--------------------------------------------------------------------------------------------------------------------
//    TDatabase::ParamDB->n_nitsche_boundary=2;
//    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
//    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
//    TDatabase::ParamDB->n_neumann_boundary=2;
//    TDatabase::ParamDB->neumann_boundary_id={1,3};
//    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
//    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
//    TDatabase::ParamDB->VISCOSITY=10;
//    TDatabase::ParamDB->PERMEABILITY=1;
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P1 on quads with Neumann and visc_eff=0.001, visc=10, perm= 1");
//    reference_errors = {{ 0.19394824424421,3.8205668689833,0.5,1.0532389294782e-15}};
//    check_brinkman2d(domain, db, 2,1, reference_errors);
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
//    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
//    TDatabase::ParamDB->VISCOSITY=10;
//    TDatabase::ParamDB->PERMEABILITY=1;
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P2 on quads with Neumann and visc_eff=0.001, visc=10, perm= 1");
//    reference_errors = {{0.19394824424421,3.8205668689833,2.8822850639982,87.950352486737}};
//    check_brinkman2d(domain, db, 2,2, reference_errors);
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
//    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
//    TDatabase::ParamDB->VISCOSITY=10;
//    TDatabase::ParamDB->PERMEABILITY=1;
//    db["PkPk_stab"]=true;
//    TDatabase::ParamDB->equal_order_stab_weight_PkPk=0.4;
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nstarting with Nitsche (parameter 1000) on P1/P1 on quads with Stab (parameter=0.4), Neumann and visc_eff=0.001, visc=10, perm= 1 ");
//    reference_errors = {{0.3582641143964,3.3275876201918,0.5119324115733,1.2550049173583}};
//    check_brinkman2d(domain, db, 1,1, reference_errors);
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
//    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
//    TDatabase::ParamDB->VISCOSITY=10;
//    TDatabase::ParamDB->PERMEABILITY=1;
//    db["PkPk_stab"]=true;
//    TDatabase::ParamDB->equal_order_stab_weight_PkPk=0.4;
//    db["example"] = 1; // Poiseuille_Hannukainen
//    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P2 on quads with Stab (parameter=0.4), Neumann and visc_eff=0.001, visc=10, perm= 1");
//    reference_errors = {{2.7665723232465,38.383469514259,3.8290546915193,26.920248715798}};
//    check_brinkman2d(domain, db, 2,2, reference_errors);
//    
//    db["example"] = 0; // Poiseuille
//    db["P2P2_stab"]=false;
//
//    
//
//}
//

//======================================================================================================
void tests_on_triangles_P2P1(unsigned int nRefinements, ParameterDatabase& db)
{   // default construct a domain object
  TDomain domain(db);
  // Initialization of the default parameters
  TDatabase::SetDefaultParameters(); 
  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  TDatabase::ParamDB->PERMEABILITY=1.;
  TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.004;
  TDatabase::ParamDB->VISCOSITY=0.004;
  TDatabase::ParamDB->PkPk_stab=false;
  TDatabase::ParamDB->GradDiv_stab=false;
  db["Galerkin_type"] = "symmetric Galerkin formulation";
  TDatabase::ParamDB->n_neumann_boundary=2;
  TDatabase::ParamDB->neumann_boundary_id={1,3};
  TDatabase::ParamDB->neumann_boundary_value={-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary=0; //2;                     
  TDatabase::ParamDB->nitsche_boundary_id={0,2};
  TDatabase::ParamDB->nitsche_penalty={0, 0};
  TDatabase::ParamDB->s1=-1;
  TDatabase::ParamDB->s2=-1;

  //TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;↲                       
  Output::print("\nstarting with Brinkman2D on TwoTriangles (symmetric Galerkin formulation), P2/P1, with Dirichlet and Neumann bcs and with visc_eff=visc=0.004, perm= 1");
  reference_errors = {{0.0054680957600273, 0.014317962190331, 0.1443805008035, 5.246987770953e-05, 0.00061083280165159}};
  check_brinkman2d(domain, db, 2,1, reference_errors); 

  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  Output::print("\nstarting with Brinkman2D on TwoTriangles (nonsymmetric Galerkin formulation), P2/P1, with Dirichlet and Neumann bcs and with visc_eff=visc=0.004, perm= 1");
  reference_errors = {{0.0054680957600273, 0.014317962190331, 0.1443805008035, 5.246987770953e-05, 0.00061083280165159}};
  check_brinkman2d(domain, db, 2,1, reference_errors); 
}

//======================================================================================================
void tests_on_triangles_P2P1_PenaltyFreeNonSymmetricNitsche(unsigned int nRefinements, ParameterDatabase& db)
{  // default construct a domain object
  TDomain domain(db);

  // Initialization of the default parameters
  TDatabase::SetDefaultParameters();

  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  TDatabase::ParamDB->PERMEABILITY=1.;
  TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.004;
  TDatabase::ParamDB->VISCOSITY=0.004;
  TDatabase::ParamDB->PkPk_stab=false;
  TDatabase::ParamDB->GradDiv_stab=false;
  db["Galerkin_type"] = "symmetric Galerkin formulation";

  TDatabase::ParamDB->n_neumann_boundary=2;
  TDatabase::ParamDB->neumann_boundary_id={1,3};
  TDatabase::ParamDB->neumann_boundary_value={-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary=2; //2;                     
  TDatabase::ParamDB->nitsche_boundary_id={0,2};
  TDatabase::ParamDB->nitsche_penalty={0, 0};
  TDatabase::ParamDB->s1=-1;
  TDatabase::ParamDB->s2=-1;

  //TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;↲                       
  Output::print("\nstarting with Brinkman2D on TwoTriangles, P2/P1, with penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff=visc=0.004, perm= 1");
  reference_errors = {{0.049354569383435, 0.049746071781475 , 0.20428824240268, 9.3291684594421e-05, 0.0010263458082476}};
  check_brinkman2d(domain, db, 2,1, reference_errors); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSstab(unsigned int nRefinements, ParameterDatabase& db)
{  // default construct a domain object
  TDomain domain(db);
  // Initialization of the default parameters
  TDatabase::SetDefaultParameters();
  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  TDatabase::ParamDB->PERMEABILITY = 1.;
  TDatabase::ParamDB->EFFECTIVE_VISCOSITY = 0.004;
  TDatabase::ParamDB->VISCOSITY = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase
  db["PkPk_stab"]=true;
  db["equal_order_stab_weight_PkPk"]=0.01;
  TDatabase::ParamDB->PkPk_stab=true;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk=0.01;
  TDatabase::ParamDB->GradDiv_stab=false;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1,3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 0; //2;                     
  TDatabase::ParamDB->nitsche_boundary_id = {0,2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;↲                       
  Output::print("\nstarting with Brinkman2D on TwoTriangles, P1/P1-Stab (symmetric GLS), Dirichlet and Neumann bcs and with visc_eff=visc=0.004, perm= 1");
  reference_errors = {{  1.2424692993215, 0.031419104874062, 16.668762604057, 0.0018335349292131, 0.034019179940349 }}; 
  check_brinkman2d(domain, db, 1, 1, reference_errors); 

  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  Output::print("\nstarting with Brinkman2D on TwoTriangles, P1/P1-Stab (nonsymmetric GLS), Dirichlet and Neumann bcs and with visc_eff=visc=0.004, perm= 1");
  reference_errors = {{ 1.2347001456568, 1.3166645102027, 16.750021539147,0.018674085154587, 0.22010507412489}};
  check_brinkman2d(domain, db, 1, 1, reference_errors); 
}


//======================================================================================================
void tests_on_triangles_P1P1_GLSstab_PenaltyFreeNonSymmetricNitsche(unsigned int nRefinements, ParameterDatabase& db)
{  // default construct a domain object
  TDomain domain(db);
  // Initialization of the default parameters
  TDatabase::SetDefaultParameters();
  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  //--------------------------------------------------------------------------------------------------------------------↲
  db["example"] = 1; // Poiseuille_Hannukainen
  TDatabase::ParamDB->PERMEABILITY = 1.;
  TDatabase::ParamDB->EFFECTIVE_VISCOSITY = 0.004;
  TDatabase::ParamDB->VISCOSITY = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["PkPk_stab"]=true;
  db["equal_order_stab_weight_PkPk"]=0.01;
  TDatabase::ParamDB->PkPk_stab=true;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk=0.01;
  TDatabase::ParamDB->GradDiv_stab=false;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1,3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;                     
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;↲                       
  Output::print("\nstarting with Brinkman2D on TwoTriangles, P1/P1-Stab (symmetric GLS), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff=visc=0.004, perm= 1");
  reference_errors = {{10.192177714504, 5.8853778774783, 25.146773888043, 0.16699522666312, 3.1493307465155 }};
  check_brinkman2d(domain, db, 1, 1, reference_errors); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  Output::print("\nstarting with Brinkman2D on TwoTriangles, P1/P1-Stab (nonsymmetric GLS), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff=visc=0.004, perm= 1");
  reference_errors = {{9.599443219491, 4.1492978030163,  23.384929664058, 0.022130923948982, 0.27177748647831}};
  check_brinkman2d(domain, db, 1, 1, reference_errors); 



}

//========================================================================
// =======================================================================
// main program
// =======================================================================
//========================================================================
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
  db.add("equal_order_stab_weight_PkPk", (double) 0., "", (double) -1000, (double) 1000 );
  db.add("refinement_n_initial_steps", (size_t) 2.0 , "", (size_t) 0, (size_t) 10000);

  //----------------------------------------
  //tests_on_quads(nRefinements, db);

  //----------------------------------------
  // tests_on_quads_stab(nRefinements, db);

  //----------------------------------------
  //tests_on_quads_Nitsche(nRefinements, db);

  //----------------------------------------
  tests_on_triangles_P2P1(nRefinements, db);

  //----------------------------------------
  tests_on_triangles_P2P1_PenaltyFreeNonSymmetricNitsche(nRefinements,db);

  //----------------------------------------
  tests_on_triangles_P1P1_GLSstab(nRefinements, db);

  //-----------------------------------------
  tests_on_triangles_P1P1_GLSstab_PenaltyFreeNonSymmetricNitsche(nRefinements, db);

  return 0;
}
