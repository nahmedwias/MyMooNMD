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

// compare the computed errors in the Brinkman2D object with the given ones in
// the array
void compareErrors(const Brinkman2D& brinkman2d, std::array<double, 4> errors)
{
    const double eps = 2e-9;
    
    // check the errors
    if( fabs(brinkman2d.getL2VelocityError() - errors[0]) > eps )
    {
        ErrThrow("Program 1: L2 velocity error not correct. ",
                 brinkman2d.getL2VelocityError() - errors[0]);
    }
    if( fabs(brinkman2d.getH1SemiVelocityError() - errors[2]) > eps )
    {
        ErrThrow("Program 1: H1-semi velocity error not correct. ",
                 brinkman2d.getH1SemiVelocityError() - errors[2]);
    }
    if( fabs(brinkman2d.getL2PressureError() - errors[3]) > eps )
    {
        ErrThrow("Program 1: L2 pressure error not correct.",
                 brinkman2d.getL2PressureError() - errors[3]);
    }
    if( fabs(brinkman2d.getH1SemiPressureError() - errors[4]) > 2*eps )
    {
        ErrThrow("Program 1: H1-semi pressure error not correct.",
                 brinkman2d.getH1SemiPressureError() - errors[4]);
    }
}

// Here the actual computations take place
void check_brinkman2d(TDomain & domain, ParameterDatabase& db, int velocityCode,int pressureCode,
                   std::array<double, 4> errors)
{
    TDatabase::ParamDB->VELOCITY_SPACE = velocityCode;
    // automatically choose pressure space to get inf-sup stable pair
    TDatabase::ParamDB->PRESSURE_SPACE = pressureCode;
    
    Brinkman2D brinkman2d(domain, db);
    brinkman2d.assemble();
    brinkman2d.solve();
    brinkman2d.output();
    // compare computed with given errors
    compareErrors(brinkman2d, errors); // throws upon a difference
}

void tests_on_quads(unsigned int nRefinements, ParameterDatabase& db)
{
    // default construct a domain object
    TDomain domain(db);
    
    
    // refine grid up to the coarsest level
    for(unsigned int i = 0; i < nRefinements; i++)
    {
        domain.RegRefineAll();
    }
    
    std::array<double, 4> errors;
    
    Output::print("\nstarting with P2/P1 on quads");
    errors = {{ 2.0257088643869e-16, 2.1906901043565e-15, 5.9796305144209e-15, 4.8385464788293e-14 }};
    check_brinkman2d(domain, db, 2,1, errors);
 
    
    // TODO
//    Output::print("\nstarting with residual-based equal-order stabilization on P1/P1 on quads");
//    errors = {{ 0.36989913525384, 8.7083491818683,
//        0.1118419830706, 2.6029572935706 }};
//    check_brinkman2d(domain, db, 1012, errors);
//    
//    Output::print("\nstarting with residual-based equal-order stabilization on P2/P2 on quads");
//    errors = {{ 0.36989913525384, 8.7083491818683,
//        0.1118419830706, 2.6029572935706 }};
//    check_brinkman2d(domain, db, 1012, errors);

}

void tests_on_quadsNitsche(unsigned int nRefinements, ParameterDatabase& db)
{
    // default construct a domain object
    TDomain domain(db);
    
//    db["n_nitsche_boundary"]=4;
//    db["nitsche_boundary_id"]= (0, 1, 2, 3);
//    db["nitsche_penalty"]= (1, 1, 1, 1);
//    
    // refine grid up to the coarsest level
    for(unsigned int i = 0; i < nRefinements; i++)
    {
        domain.RegRefineAll();
    }
    
    std::array<double, 4> errors;
    
    Output::print("\nstarting with Nitsche on P2/P1 on quads");
    errors = {{ 0.27161360354054,1.6902055317218, 3.6216269331927, 8.4556818185136}};
    check_brinkman2d(domain, db, 2,1, errors);
    
    
    // TODO
    //    Output::print("\nstarting with Nitsche on quads");
    //    errors = {{ 0.36989913525384, 8.7083491818683,
    //        0.1118419830706, 2.6029572935706 }};
    //    check_brinkman2d(domain, db, 2, errors);
    //
    //    Output::print("\nstarting with residual-based equal-order stabilization on P1/P1 on quads");
    //    errors = {{ 0.36989913525384, 8.7083491818683,
    //        0.1118419830706, 2.6029572935706 }};
    //    check_brinkman2d(domain, db, 1012, errors);
    //
    //    Output::print("\nstarting with residual-based equal-order stabilization on P2/P2 on quads");
    //    errors = {{ 0.36989913525384, 8.7083491818683,
    //        0.1118419830706, 2.6029572935706 }};
    //    check_brinkman2d(domain, db, 1012, errors);
    
}



void tests_on_triangles(unsigned int nRefinements, ParameterDatabase& db)
{//Not yet Nitsche (normals)
    // default construct a domain object
    TDomain domain(db);
    
    // refine grid up to the coarsest level
    for(unsigned int i = 0; i < nRefinements; i++)
    {
        domain.RegRefineAll();
    }
    
            // TODO
//    std::array<double, 4> errors;
    

//    Output::print("\nstarting with P2/P1 on triangles");
//    errors = {{ 2.1136884064519, 23.120239110875,
//        0.30277518654981, 4.4428829381584 }};
//    check_brinkman2d(domain, db, 2, errors);
//    

//    Output::print("\nstarting with P1Mini on triangles");
//    errors = {{ 1.7301620785317, 23.120239110875,
//        0.32075488021636, 4.4428829381584 }};
//    check_brinkman2d(domain, db, 101, errors);
//    cout << "P1MINI works only on triangles" << endl;
//    
//    Output::print("\nstarting with residual-based equal-rder stabilization on P1/P1 on quads");
//    errors = {{ 0.36989913525384, 8.7083491818683,
//        0.1118419830706, 2.6029572935706 }};
//    check_brinkman2d(domain, db, 1012, errors);
//    
//    Output::print("\nstarting with residual-based equal-order stabilization on P2/P2 on quads");
//    errors = {{ 0.36989913525384, 8.7083491818683,
//        0.1118419830706, 2.6029572935706 }};
//    check_brinkman2d(domain, db, 1012, errors);
}



// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
    //  declaration of databases
    TDatabase Database;
    TFEDatabase2D FEDatabase;
    
    // velocity space code for Raviart-Thomas (RT) and
    // Brezzi-Douglas-Marini(BDM) elements:
    // 1000    RT_0
    // 1001    RT_1
    // 1002    RT_2
    // 1003    RT_3
    // 1011    BDM_1
    // 1012    BDM_2
    // 1013    BDM_3
    
    // 101 P1Mini
    //cout << "P1MINI works only on triangles" << endl;
    
    TDatabase::ParamDB->VELOCITY_SPACE = 2;
    // automatically choose pressure space to get inf-sup stable pair
    TDatabase::ParamDB->PRESSURE_SPACE = 1;
    // high order quadrature for computing errors
    TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
    
    TDatabase::ParamDB->VISCOSITY=1;
    TDatabase::ParamDB->EFFECTIVE_VISCOSITY=1;
    TDatabase::ParamDB->PERMEABILITY=1;
    
    

    
    unsigned int nRefinements = 2;
    
    Output::setVerbosity(2);
    
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Solver<>::default_solver_database());
    db.merge(ParameterDatabase::default_output_database());
    db.merge(Example2D::default_example_database());


    
    db["example"] = 0; // known Poiseuille Hannukainen solution
    db["residual_tolerance"] = 1.0e-13;
    
    db["output_compute_errors"] = true;
    
    Output::print("\n\n ----------- direct solver -----------\n");
    db["solver_type"] = "direct";
    
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    
    tests_on_quads(nRefinements, db);
//    db["geo_file"] = "TwoTriangles";
//  tests_on_triangles(nRefinements, db);
    
    
    db["solver_type"] = "iterative";
    db["max_n_iterations"] = 10000;
    
    Output::print("\n\n --------- fgmres+lsc solver ---------\n");
    db["preconditioner"] = "least_squares_commutator";
    
//    db["geo_file"] = "UnitSquare";
//    tests_on_quads(nRefinements, db);
//    db["geo_file"] = "TwoTriangles";
//    tests_on_triangles(nRefinements, db);
    
    Output::print("\n\n -------- fgmres+simple solver -------\n");
    db["preconditioner"] = "semi_implicit_method_for_pressure_linked_equations";
    
//    db.add("n_nitsche_boundary",4);
//    db.add("nitsche_boundary_id", 0 1 2 3);
//           db.add("nitsche_penalty", 1 1 1 1);
//    tests_on_quadsNitsche(nRefinements, db);
    

    


   
    return 0;
}
