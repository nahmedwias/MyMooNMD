/**
 * @brief A test program for the solving of Darcy2D problems.
 *
 * This serves as a test for the solving of Darcy2D problems. It is intended to
 * perform Darcy2D calculations with different examples in different setups to 
 * test a wide variety of ParMooN core functionality.
 * So far only one such test is implemented.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes 
 * wrong) the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then 
 * this tests shows you how many other program parts are affected by your 
 * changes. If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Darcy2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_CD2D.h>

#include <MainUtilities.h> //for error measuring

// compare the computed errors in the Darcy2D object with the given ones in 
// the array
void compareErrors(const Darcy2D& darcy2d, std::array<double, 5> errors)
{
  const double eps = 1e-9;
  // check the errors
  if( fabs(darcy2d.getL2VelocityError() - errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 velocity error not correct. ",
             darcy2d.getL2VelocityError() - errors[0]);
  }
  if( fabs(darcy2d.getL2DivergenceError() - errors[1]) > eps )
  {
    ErrThrow("Program 1: L2 velocity divergence error not correct. ",
             darcy2d.getL2DivergenceError() - errors[1]);
  }
  if( fabs(darcy2d.getH1SemiVelocityError() - errors[2]) > eps )
  {
    ErrThrow("Program 1: H1-semi velocity error not correct. ",
             darcy2d.getH1SemiVelocityError() - errors[2]);
  }
  if( fabs(darcy2d.getL2PressureError() - errors[3]) > eps )
  {
    ErrThrow("Program 1: L2 pressure error not correct.",
             darcy2d.getL2PressureError() - errors[3]);
  }
  if( fabs(darcy2d.getH1SemiPressureError() - errors[4]) > eps )
  {
    ErrThrow("Program 1: H1-semi pressure error not correct.",
             darcy2d.getH1SemiPressureError() - errors[4]);
  }
}

// Here the actual computations take place
void check(TDomain & domain, int velocityCode, std::array<double, 5> errors)
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db["problem_type"] = 0; // problem type is not needed
  db["example"] = 0;
  db.add("solver_type", (size_t)2, "");
  
  TDatabase::ParamDB->VELOCITY_SPACE = velocityCode;
  // automatically choose pressure space to get inf-sup stable pair
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  
  Darcy2D darcy2d(domain, db);
  darcy2d.assemble();
  darcy2d.solve();
  darcy2d.output();
  // compare computed with given errors
  compareErrors(darcy2d, errors); // throws upon a difference
}

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of databases
  TDatabase Database;
  TFEDatabase2D FEDatabase;

  // Set Database values (this is what is usually done by the input-file)
  TDatabase::ParamDB->PROBLEM_TYPE = 0; // problem type is not needed
  TDatabase::ParamDB->EXAMPLE = 0; // known sin-cos solution
  // velocity space code for Raviart-Thomas (RT) and 
  // Brezzi-Douglas-Marini(BDM) elements:
  // 1000    RT_0
  // 1001    RT_1
  // 1002    RT_2
  // 1003    RT_3
  // 1011    BDM_1
  // 1012    BDM_2
  // 1013    BDM_3
  TDatabase::ParamDB->VELOCITY_SPACE = 1000;
  // automatically choose pressure space to get inf-sup stable pair
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->MEASURE_ERRORS = 1; // compute errors
  // high order quadrature for computing errors 
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  TDatabase::ParamDB->WRITE_VTK = 0; // do not write vtk files
  TDatabase::ParamDB->SIGMA_PERM = 1.; // permeability
  TDatabase::ParamDB->SC_VERBOSE = 0; // supress solver output
  TDatabase::ParamDB->SOLVER_TYPE = 2; // use direct solver

  unsigned int nRefinements = 4;
  
  // cells are squares
  {
    // default construct a domain object
    TDomain domain;
    // the domain is initialised with default description and default
    // initial mesh
    domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");

    // refine grid up to the coarsest level
    for(unsigned int i = 0; i < nRefinements; i++)
    {
      domain.RegRefineAll();
    }

    std::array<double, 5> errors;
    
    Output::print("\nstarting with RT0 on quads");
    errors = {{ 0.50616706358332, 6.2926762302213, 28.267854806247, 
                0.079978980936799, 4.4428829381584 }};
    check(domain, 1000, errors);
    
    Output::print("\nstarting with RT1 on quads");
    errors = {{ 0.025524493858475, 0.32016058303363, 3.1653722663713, 
                0.0040560305244671, 0.50312981240831 }};
    check(domain, 1001, errors);
    
    Output::print("\nstarting with RT2 on quads");
    errors = {{ 0.00084661911007719, 0.010629840538212, 0.17550885275941, 
                0.00013465340078883, 0.027920142668719 }};
    check(domain, 1002, errors);
    
    Output::print("\nstarting with RT3 on quads");
    errors = {{ 2.094182550122e-05, 0.00026303115293237, 0.0063544526017813, 
                3.3317619753011e-06, 0.0010111221715242 }};
    check(domain, 1003, errors);
    
    Output::print("\nstarting with BDM1 on quads");
    errors = {{ 0.12631404854406, 6.2926762302213, 7.0233297819314, 
                0.080703320512897, 4.4428829381584 }};
    check(domain, 1011, errors);
    
    Output::print("\nstarting with BDM2 on quads");
    errors = {{ 0.0046193106648853, 0.59661934959593, 0.54355308167768, 
                0.00755750437127, 0.708252543173 }};
    check(domain, 1012, errors);
    
    Output::print("\nstarting with BDM3 on quads");
    errors = {{ 0.00025013087012341, 0.037688202895706, 0.036176275134976, 
                0.00047733398111123, 0.068099895638466 }};
    check(domain, 1013, errors);
  }
  
  // cells are triangles
  {
    // default construct a domain object
    TDomain domain;
    // the domain is initialised with default description and default
    // initial mesh
    domain.Init((char*)"Default_UnitSquare", (char*)"TwoTriangles");

    // refine grid up to the coarsest level
    for(unsigned int i = 0; i < nRefinements; i++)
    {
      domain.RegRefineAll();
    }

    std::array<double, 5> errors;

    Output::print("\nstarting with RT0 on triangles");
    errors = {{ 0.50408358816775, 5.1429079193159, 28.151331748344, 
                0.0652645752817, 4.4428829381584 }};
    check(domain, 1000, errors);
    
    Output::print("\nstarting with RT1 on triangles");
    errors = {{ 0.028196037109588, 0.390978706619, 3.7142101823095, 
                0.0049554126577833, 0.63339791851833 }};
    check(domain, 1001, errors);
    
    Output::print("\nstarting with RT2 on triangles");
    errors = {{ 0.0011486412161111, 0.017485986235863, 0.2702544676034, 
                2.7194995909334e-05, 0.055705709939193 }};
    check(domain, 1002, errors);
    
    Output::print("\nstarting with RT3 on triangles");
    errors = {{ 4.5103729086427e-05, 0.0010136026623928, 0.012803427384749, 
                1.2004518891142e-05, 0.0035296757619195 }};
    check(domain, 1003, errors);
    
    Output::print("\nstarting with BDM1 on triangles");
    errors = {{ 0.099862267417049, 5.1429079193159, 7.4906379811947,
                0.0657343839084, 4.4428829381584 }};
    check(domain, 1011, errors);
    
    Output::print("\nstarting with BDM2 on triangles");
    errors = {{ 0.0035151373346507, 0.390978706619, 0.58694587503024, 
                0.004954555083239, 0.63428676374772 }};
    check(domain, 1012, errors);
    
    Output::print("\nstarting with BDM3 on triangles");
    errors = {{ 0.00017651588561095, 0.017485986236131, 0.029975320141965, 
                2.6130552218382e-05, 0.055748396623119 }};
    check(domain, 1013, errors);
  }
  return 0;
}
