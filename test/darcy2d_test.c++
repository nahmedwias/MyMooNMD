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
  TDatabase::ParamDB->VELOCITY_SPACE = velocityCode;
  // automatically choose pressure space to get inf-sup stable pair
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  
  Darcy2D darcy2d(domain);
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
    for(int i = 0; i < nRefinements; i++)
    {
      domain.RegRefineAll();
    }

    std::array<double, 5> errors;
    
    Output::print("\nstarting with RT0 on quads");
    errors = { 0.50665741516721, 6.2926762302213, 28.267995185207, 
               0.080296618540803, 4.4428829381584 };
    check(domain, 1000, errors);
    
    Output::print("\nstarting with RT1 on quads");
    errors = { 0.027836487104237, 0.32016058303363, 3.1656276117037, 
               0.0047912634003738, 0.50325154768943 };
    check(domain, 1001, errors);
    
    Output::print("\nstarting with RT2 on quads");
    errors = { 0.00084993853637184, 0.010629840538212, 0.17542990320099, 
               0.0001346764830322, 0.027920315010037 };
    check(domain, 1002, errors);
    
    Output::print("\nstarting with RT3 on quads");
    errors = { 2.1078861709082e-05, 0.00026303115293237, 0.0063440847932246, 
               3.3567074487812e-06, 0.0010111343061848 };
    check(domain, 1003, errors);
    
    Output::print("\nstarting with BDM1 on quads");
    errors = { 0.12677477212599, 6.2926762302213, 7.024091423718, 
               0.080589930926966, 4.4428829381584 };
    check(domain, 1011, errors);
    
    Output::print("\nstarting with BDM2 on quads");
    errors = { 0.0046500390802298, 0.59661934959593, 0.54462179448715, 
               0.0075575013964879, 0.70825261815867 };
    check(domain, 1012, errors);
    
    Output::print("\nstarting with BDM3 on quads");
    errors = { 0.00025076541593067, 0.037688202895721, 0.036241187043479, 
               0.00047733415522407, 0.068099895944846 };
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
    for(int i = 0; i < nRefinements; i++)
    {
      domain.RegRefineAll();
    }

    std::array<double, 5> errors;

    Output::print("\nstarting with RT0 on triangles");
    errors = {0.50446546401145, 5.1429079193159, 28.151331748344, 
              0.065512804062722, 4.4428829381584 };
    check(domain, 1000, errors);
    
    Output::print("\nstarting with RT1 on triangles");
    errors = {0.030316560421862, 0.390978706619, 3.7129997012541, 
              0.0055749146617156, 0.63333061754082 };
    check(domain, 1001, errors);
    
    Output::print("\nstarting with RT2 on triangles");
    errors = {0.0011515922414381, 0.017485986235774, 0.27018639952401, 
              2.7335455463106e-05, 0.055706379020305 };
    check(domain, 1002, errors);
    
    Output::print("\nstarting with RT3 on triangles");
    errors = {4.5272848786446e-05, 0.001013602663064, 0.012806644670806, 
              1.2011482390322e-05, 0.0035296782886598 };
    check(domain, 1003, errors);
    
    Output::print("\nstarting with BDM1 on triangles");
    errors = { 0.10074771842405, 5.1429079193159, 7.4896219862889, 
               0.065660844395943, 4.4428829381584 };
    check(domain, 1011, errors);
    
    Output::print("\nstarting with BDM2 on triangles");
    errors = { 0.0035229371407354, 0.390978706619, 0.58724813677087, 
               0.0049545500597403, 0.63428187104069 };
    check(domain, 1012, errors);
    
    Output::print("\nstarting with BDM3 on triangles");
    errors = { 0.00017672344384691, 0.017485986235533, 0.029985240418542, 
               2.6133409945217e-05, 0.055748288399196 };
    check(domain, 1013, errors);
  }
  return 0;
}
