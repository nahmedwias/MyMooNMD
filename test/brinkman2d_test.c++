/**
 * @brief A test program for solving Brinkman problems.
 *
 * This serves as a test for solving Brinkman problems. It is intended to
 * perform Brinkman calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
 * So far only one such test is implemented.
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
void compareErrors(const Brinkman2D& brinkman2d, std::array<double, 4> reference_errors)
{
    const double eps = 2e-9;
    
    
    Output::print(setprecision(14),reference_errors[0]);
    Output::print(setprecision(14),reference_errors[1]);
    Output::print(setprecision(14),reference_errors[2]);
    Output::print(setprecision(14),reference_errors[3]);
    //    std::array<double, int(4)> current_errors;
    //    current_errors = brinkman2d.get_errors();
    
    
    // check the errors
    if( fabs(brinkman2d.getL2VelocityError() - reference_errors[0]) > eps )
    {
        ErrThrow("Program 1: L2 velocity error not correct. ",
                 brinkman2d.getL2VelocityError() - reference_errors[0]);
    }
    if( fabs(brinkman2d.getH1SemiVelocityError() - reference_errors[1]) > eps)
    {
        ErrThrow("Program 1: H1-semi velocity error not correct. ",
                 brinkman2d.getH1SemiVelocityError() - reference_errors[1]);
    }
    
    if( fabs(brinkman2d.getL2PressureError() - reference_errors[2]) > eps )
    {
        ErrThrow("Program 1: L2 pressure error not correct.",
                 brinkman2d.getL2PressureError() - reference_errors[2]);
    }
    if( fabs(brinkman2d.getH1SemiPressureError() - reference_errors[3]) > eps )
    {
        ErrThrow("Program 1: H1-semi pressure error not correct.",
                 brinkman2d.getH1SemiPressureError() - reference_errors[3]);
        
    }
}

//----------------------------------------------------------------------
// Here the actual computations take place
void check_brinkman2d(TDomain & domain, ParameterDatabase& db, int velocityCode,int pressureCode,
                      std::array<double, 4> reference_errors)
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
    compareErrors(brinkman2d, reference_errors); // throws upon a difference
}

//======================================================================================================
void tests_on_quads(unsigned int nRefinements, ParameterDatabase& db)
{
    std::array<double, 4> reference_errors;
    
    // default construct a domain object
    TDomain domain(db);
    
    // Initialization of the default parameters
    TDatabase::SetDefaultParameters();
    
    // refine grid up to the coarsest level
    for(unsigned int i = 0; i < nRefinements; i++)
    {
        domain.RegRefineAll();
    }
    
    //------------------------------------------------------------------------------------------------------------
    Output::print("\nStarting with P2/P1 on quads for Poiseuille flow with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
    reference_errors = {{ 2.0257088643869e-16, 2.1906901043565e-15, 5.9796305144209e-15, 4.8385464788293e-14 }};
    
    check_brinkman2d(domain, db, 2,1, reference_errors);
    
    //    //------------------------------------------------------------------------------------------------------------
    //    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.01;
    //    TDatabase::ParamDB->VISCOSITY=10;
    //    TDatabase::ParamDB->PERMEABILITY=2;
    //
    //    Output::print("\nStarting with P2/P1 on quads for Poiseuille flow with visc_eff=0.01, visc=10, perm=2");
    //    reference_errors = {{1.9508695699745e-16,3.6646227316896e-15,1.3277195222811e-16,1.0734187033819e-15}};
    //
    //    check_brinkman2d(domain, db, 2,1, reference_errors);
    //
    //    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=1;
    //    TDatabase::ParamDB->VISCOSITY=1;
    //    TDatabase::ParamDB->PERMEABILITY=1;
    
    //------------------------------------------------------------------------------------------------------------
    db["example"] = 1; // Poiseuille_Hannukainen
    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
    reference_errors = {{2.2608321672551e-05, 0.00059223693639448, 9.7288545102352e-06,9.0848530457398e-05}};
    
    check_brinkman2d(domain, db, 2,1, reference_errors);
    db["example"] = 0; // Poiseuille
    
    
    //------------------------------------------------------------------------------------------------------------
    db["example"] = 1; // Poiseuille_Hannukainen
    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
    TDatabase::ParamDB->VISCOSITY=1;
    TDatabase::ParamDB->PERMEABILITY=1;
    
    Output::print("\nStarting with P2/P1 on quads for Poiseuille flow with visc_eff=0.001, visc=1, perm=2 for example ",db["example"] );
    reference_errors = {{0.097304080869912,2.9931806337455,0.00323505572483,0.027631651143564}};
    
    check_brinkman2d(domain, db, 2,1, reference_errors);
    
    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=1;
    TDatabase::ParamDB->VISCOSITY=1;
    TDatabase::ParamDB->PERMEABILITY=1;
    db["example"] = 0; // Poiseuille
    
    //------------------------------------------------------------------------------------------------------------
    db["example"] = 2; // Poiseuille_Hannukainen with inscribed sphere
    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
    reference_errors = {{0.094501334268564,1.2200819222258,8.9575878583122,32.442781399324}};
    
    check_brinkman2d(domain, db, 2,1, reference_errors);
    db["example"] = 0; // Poiseuille
    
    //------------------------------------------------------------------------------------------------------------
    db["example"] = 3; // sincos
    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"] );
    reference_errors = {{0.0040823711607968, 0.10522135481493,0.017649828876871, 0.51271046930702}};
    
    check_brinkman2d(domain, db, 2,1, reference_errors);
    db["example"] = 0; // Poiseuille
    
    //------------------------------------------------------------------------------------------------------------
    db["example"] = 4; // sincos2
    Output::print("\nStarting with P2/P1 on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"]);
    reference_errors = {{1.4228987330444,10.513375915835, 12.734621500535, 85.258423292316}};
    
    check_brinkman2d(domain, db, 2,1, reference_errors);
    db["example"] = 0; // Poiseuille
}

//======================================================================================================
void tests_on_quads_stab(unsigned int nRefinements, ParameterDatabase& db)
{
    std::array<double, 4>  reference_errors;
    
    // default construct a domain object
    TDomain domain(db);
    
    // Initialization of the default parameters
    TDatabase::SetDefaultParameters();
    
    // refine grid up to the coarsest level
    for(unsigned int i = 0; i < nRefinements; i++)
    {
        domain.RegRefineAll();
    }
    
    //------------------------------------------------------------------------------------------------------------
    //P1P1
    Output::print("\nStarting with residual-based equal-order stabilization for P1/P1 (parameter=0.4) on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"]);
    reference_errors = {{ 0.20649934279357,0.97728617227336, 1.3352827196758, 5.1734925650333}};
    
    db["P1P1_stab"] = true;
    TDatabase::ParamDB->equal_order_stab_weight_P1P1=0.4;
    
    check_brinkman2d(domain, db, 1,1, reference_errors);
    
    db["P1P1_stab"] = false;
    
    //------------------------------------------------------------------------------------------------------------
    //P2P2
    Output::print("\nStarting with residual-based equal-order stabilization for P2/P2 (parameter=0.4) on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"]);
    reference_errors = {{0.11209308296508,1.6714743655252, 0.88411612007029, 10.192320142372}};
    
    db["P2P2_stab"] = true;
    TDatabase::ParamDB->equal_order_stab_weight_P2P2=0.4;
    
    check_brinkman2d(domain, db, 2,2, reference_errors);
    
    db["P2P2_stab"] = false;
    
    //------------------------------------------------------------------------------------------------------------
    //P2P2
    db["example"] = 1; // Poiseuille_Hannukainen
    Output::print("\nStarting with residual-based equal-order stabilization for P2/P2 (parameter=0.4) on quads with visc_eff=1, visc=1, perm=1 for example ",db["example"]);
    reference_errors = {{0.012909128829698,0.19363961444426,0.10019110454835,1.1661459244268}};
    
    db["P2P2_stab"] = true;
    TDatabase::ParamDB->equal_order_stab_weight_P2P2=0.4;
    
    check_brinkman2d(domain, db, 2,2, reference_errors);
    
    db["P2P2_stab"] = false;
    db["example"] = 0; // Poiseuille
    
}

//======================================================================================================
void tests_on_quads_Nitsche(unsigned int nRefinements, ParameterDatabase& db)
{
    // default construct a domain object
    TDomain domain(db);
    
    // Initialization of the default parameters
    TDatabase::SetDefaultParameters();
    
    std::array<double, 4>  reference_errors;
    
    // refine grid up to the coarsest level
    for(unsigned int i = 0; i < nRefinements; i++)
    {
        domain.RegRefineAll();
    }
    //--------------------------------------------------------------------------------------------------------------------
    
    TDatabase::ParamDB->n_nitsche_boundary=4;
    TDatabase::ParamDB->nitsche_boundary_id={0, 1, 2, 3};
    TDatabase::ParamDB->nitsche_penalty={1000, 1000, 1000, 1000};
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P1 on quads with visc_eff=1, visc=1, perm= 1");
    reference_errors = {{0.00071104275913606,0.0086539449600541,0.0137946704876,0.11201103666807}};
    check_brinkman2d(domain, db, 2,1, reference_errors);
    
    
    //--------------------------------------------------------------------------------------------------------------------
    TDatabase::ParamDB->n_nitsche_boundary=2;
    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
    TDatabase::ParamDB->n_neumann_boundary=2;
    TDatabase::ParamDB->neumann_boundary_id={1,3};
    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
    Output::print("\nstarting with Nitsche (parameter 1000)  on P2/P1 on quads with Neumann and visc_eff=1, visc=1, perm= 1");
    reference_errors = {{0.0018492680116681, 0.00052482584609959, 0.5, 2.0166230982544e-14}};
    check_brinkman2d(domain, db, 2,1, reference_errors);
    
    
    //--------------------------------------------------------------------------------------------------------------------
    TDatabase::ParamDB->n_nitsche_boundary=2;
    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
    TDatabase::ParamDB->n_neumann_boundary=2;
    TDatabase::ParamDB->neumann_boundary_id={1,3};
    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
    TDatabase::ParamDB->VISCOSITY=10;
    TDatabase::ParamDB->PERMEABILITY=1;
    Output::print("\nstarting with Nitsche (parameter 1000)  on P2/P1 on quads with Neumann and visc_eff=0.001, visc=10, perm= 1");
    reference_errors = {{4.8144012875101e-05,0.0010541202262723, 0.5,1.184219258062e-15}};
    check_brinkman2d(domain, db, 2,1, reference_errors);
    
    //--------------------------------------------------------------------------------------------------------------------
    TDatabase::ParamDB->n_nitsche_boundary=2;
    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
    TDatabase::ParamDB->n_neumann_boundary=2;
    TDatabase::ParamDB->neumann_boundary_id={1,3};
    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
    TDatabase::ParamDB->VISCOSITY=10;
    TDatabase::ParamDB->PERMEABILITY=1;
    db["example"] = 1; // Poiseuille_Hannukainen
    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P1 on quads with Neumann and visc_eff=0.001, visc=10, perm= 1");
    reference_errors = {{ 0.89129631559785,4.821553845982,0.5,7.4513009195699e-16}};
    check_brinkman2d(domain, db, 2,1, reference_errors);
    
    db["example"] = 0; // Poiseuille
    
    //--------------------------------------------------------------------------------------------------------------------
    TDatabase::ParamDB->n_nitsche_boundary=2;
    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
    TDatabase::ParamDB->n_neumann_boundary=2;
    TDatabase::ParamDB->neumann_boundary_id={1,3};
    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
    TDatabase::ParamDB->VISCOSITY=10;
    TDatabase::ParamDB->PERMEABILITY=1;
    db["example"] = 1; // Poiseuille_Hannukainen
    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P2 on quads with Neumann and visc_eff=0.001, visc=10, perm= 1");
    reference_errors = {{ 0.89129631559785,4.821553845982, 0.50408661422354,2.342959540112}};
    check_brinkman2d(domain, db, 2,2, reference_errors);
    
    db["example"] = 0; // Poiseuille
    
    //--------------------------------------------------------------------------------------------------------------------
    TDatabase::ParamDB->n_nitsche_boundary=2;
    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
    TDatabase::ParamDB->n_neumann_boundary=2;
    TDatabase::ParamDB->neumann_boundary_id={1,3};
    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
    TDatabase::ParamDB->VISCOSITY=10;
    TDatabase::ParamDB->PERMEABILITY=1;
    db["P1P1_stab"]=true;
    TDatabase::ParamDB->equal_order_stab_weight_P1P1=0.4;
    db["example"] = 1; // Poiseuille_Hannukainen
    Output::print("\nstarting with Nitsche (parameter 1000) on P1/P1 on quads with Stab (parameter=0.4), Neumann and visc_eff=0.001, visc=10, perm= 1 ");
    reference_errors = {{0.90295941738137,2.1622048531632, 0.50378342568305,0.55512036706005}};
    check_brinkman2d(domain, db, 1,1, reference_errors);
    
    db["example"] = 0; // Poiseuille
    db["P1P1_stab"]=false;
    
    //--------------------------------------------------------------------------------------------------------------------
    TDatabase::ParamDB->n_nitsche_boundary=2;
    TDatabase::ParamDB->nitsche_boundary_id={0, 2};
    TDatabase::ParamDB->nitsche_penalty={1000, 1000};
    TDatabase::ParamDB->n_neumann_boundary=2;
    TDatabase::ParamDB->neumann_boundary_id={1,3};
    TDatabase::ParamDB->neumann_boundary_value={-1, 0};
    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=0.001;
    TDatabase::ParamDB->VISCOSITY=10;
    TDatabase::ParamDB->PERMEABILITY=1;
    db["P2P2_stab"]=true;
    TDatabase::ParamDB->equal_order_stab_weight_P2P2=0.4;
    db["example"] = 1; // Poiseuille_Hannukainen
    Output::print("\nstarting with Nitsche (parameter 1000) on P2/P2 on quads with Stab (parameter=0.4), Neumann and visc_eff=0.001, visc=10, perm= 1");
    reference_errors = {{0.95894542662398, 10.718540461577, 0.52113216529381, 3.0875564733177}};
    check_brinkman2d(domain, db, 2,2, reference_errors);
    
    db["example"] = 0; // Poiseuille
    db["P2P2_stab"]=false;

    

}


//======================================================================================================
void tests_on_triangles(unsigned int nRefinements, ParameterDatabase& db)
{//Not yet Nitsche (normals)
    // default construct a domain object
    TDomain domain(db);
    
    // Initialization of the default parameters
    TDatabase::SetDefaultParameters();
    
    // refine grid up to the coarsest level
    for(unsigned int i = 0; i < nRefinements; i++)
    {
        domain.RegRefineAll();
    }
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
    db.merge(ParameterDatabase::default_output_database());
    db.merge(Example2D::default_example_database());
    db["example"] = 0; // known Poiseuille solution
    db["output_compute_errors"] = true;
    db["output_write_vtk"] = true;
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare","Default_UnitSquare", "TwoTriangles"});
    db.add("P1P1_stab", (bool) false, "" );
    db.add("P2P2_stab", (bool) false, "" );
    db.add("refinement_n_initial_steps", (size_t) 2.0 , "", (size_t) 0, (size_t) 10000);
    
    //----------------------------------------
    tests_on_quads(nRefinements, db);
    
    //----------------------------------------
    tests_on_quads_stab(nRefinements, db);
    
    
    //----------------------------------------
    tests_on_quads_Nitsche(nRefinements, db);
    
    
    return 0;
}
