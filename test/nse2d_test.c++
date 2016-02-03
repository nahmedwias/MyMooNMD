/**
 * @brief A test program for the solving of NSE2D problems.
 *
 * This serves as a test for the solving of NSE2D problems. It is intended to
 * perform NSE2D calculations with different examples in different setups to test
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
 * @author Naveed, Ulrich, Clemens
 *
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <NSE2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_NSE2D.h>

#include <MainUtilities.h> //for error measuring

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  /** Program 1
   *  This program tests direct solve with galerkin discretization
   * direct solver and multigrid 
   */
  {

    //  declaration of databases
    TDatabase Database;
    TFEDatabase2D FEDatabase;

    // default construct a domain object
    TDomain domain;

    TDatabase::ParamDB->PROBLEM_TYPE = 5; //NSE Problem
    TDatabase::ParamDB->EXAMPLE = 2; 
    TDatabase::ParamDB->UNIFORM_STEPS = 3;
    TDatabase::ParamDB->RE_NR=1;
    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LEVELS =1;
    TDatabase::ParamDB->SOLVER_TYPE = 2;
    TDatabase::ParamDB->LAPLACETYPE = 0;
    TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE= 1e-10;
    TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE=50;
    TDatabase::ParamDB->MEASURE_ERRORS = 1;
    
    TDatabase::ParamDB->VELOCITY_SPACE = 2;
    // if pressure space is -4711, the depending on 
    // velcoity space, the pressure space is auto chosen in
    // the NSE2D class or set to 1 for triangles
    TDatabase::ParamDB->PRESSURE_SPACE = -4711; 
    
    // possibly parameters in the database
    Database.CheckParameterConsistencyNSE();
    // the domain is initialised with default description and default
    // initial mesh
    domain.Init((char*)"Default_UnitSquare", (char*)"TwoTriangles");
    // refine grid up to the coarsest level
    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    {
      domain.RegRefineAll();
    }

    //=========================================================================
    // creat an object 
    NSE2D nse2d(domain);
    // assemble all 
    nse2d.assemble();
    // check stopping criterion
    nse2d.stopIt(0);

    if(  (nse2d.getImpulsResidual() - 20.57064552  > 1e-6)
      || (nse2d.getMassResidual() -   0.5821843959 > 1e-6)
      || (nse2d.getFullResidual() -   20.57888227  > 1e-6)
    )
    {
      ErrThrow("Residual at iteration 0 is not correct");
    }
    
    for(unsigned int k=1;; k++)
    {
      Output::print<1>("nonlinear step " , setw(3), k-1, "\t",
        nse2d.getResiduals());
      nse2d.solve();
      
      if(TDatabase::ParamDB->PROBLEM_TYPE == 3)
        break;
      // checking the first nonlinear iteration
      if( (k==2) && 
         ( (nse2d.getImpulsResidual() -   0.3965057248    >1e-6)        
         ||  (nse2d.getMassResidual() -   2.194728449e-16 > 1e-6)
         ||  (nse2d.getFullResidual() -   0.3965057248    >1e-6) ) )
      {
        ErrThrow("Residual at iteration 1 is not correct");
      }
      
      nse2d.assemble_nonlinear_term();;
      if(nse2d.stopIt(k))
        break;
    }
    
    double err[4];
    TAuxParam2D NSEaux_error;
    MultiIndex2D NSAllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D *velocity_space = &nse2d.get_velocity_space();
    const TFESpace2D *pressure_space = &nse2d.get_pressure_space();
    
    const TFEFunction2D *u1 = nse2d.get_velocity_component(0);
    const TFEFunction2D *u2 = nse2d.get_velocity_component(1);
    
    u1->GetErrors(nse2d.get_example().get_exact(0), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &velocity_space, err);
    // errors in second velocity component
    u2->GetErrors(nse2d.get_example().get_exact(1), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &velocity_space, err + 2);
    // check the L2-error of the velcoity
    if( fabs( sqrt(err[0]*err[0] + err[2]*err[2]) - 0.000610487) > 1e-6 )
    {
      ErrThrow("Program 1: L2 norm of velocity is not correct.");
    }
    // check the H1-error of the velcoity
    if( fabs( sqrt(err[1]*err[1] + err[3]*err[3]) - 0.0389713) > 1e-6 )
    {
      ErrThrow("Program 1: H1-semi norm of velocity not correct.");
    }
    
    // errors in pressure
    const TFEFunction2D &p = nse2d.get_pressure();
    p.GetErrors(nse2d.get_example().get_exact(2), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &pressure_space, err);
    // check the L2-error of the pressure
    if( fabs(err[0] - 0.0107332) > 1e-6)
    {
      ErrThrow("Program 1: L2 norm of pressure is not correct.");
    }
    // check the H1-error of the pressure
    if(fabs( err[1] - 0.512621) > 1e-6 )
    {
      ErrThrow("Program 1: H1-norm of pressure is not correct.");
    } 
    Output::print<1>("Elements: P2/P1: test passed for the Direct Solver: ");
  } // end program 1
  
  //////////////////////////////////////////////////////////////////////////////////
  {

    //  declaration of databases
    TDatabase Database;
    TFEDatabase2D FEDatabase;

    // default construct a domain object
    TDomain domain;

    TDatabase::ParamDB->PROBLEM_TYPE = 5; //NSE Problem
    TDatabase::ParamDB->EXAMPLE = 2; 
    TDatabase::ParamDB->UNIFORM_STEPS = 4;
    TDatabase::ParamDB->RE_NR=1;
    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LEVELS =1;
    TDatabase::ParamDB->SOLVER_TYPE = 2;
    TDatabase::ParamDB->LAPLACETYPE = 0;
    TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE= 1e-10;
    TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE=50;
    TDatabase::ParamDB->MEASURE_ERRORS = 1;
    
    TDatabase::ParamDB->VELOCITY_SPACE = 12;
    // Q_2/P_1^disc
    TDatabase::ParamDB->PRESSURE_SPACE = -4711; 
    
    // possibly parameters in the database
    Database.CheckParameterConsistencyNSE();

    // the domain is initialised with default description and default
    // initial mesh
    domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");
    // refine grid up to the coarsest level
    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    {
      domain.RegRefineAll();
    }

    //=========================================================================
    // creat an object 
    NSE2D nse2d(domain);
    // assemble all 
    nse2d.assemble();
    // check stopping criterion
    nse2d.stopIt(0);
    
    OutPut("NOTE: There is something wrong in the nse.getImpulsResidual()"
        " One have to correct that and update the test accordingly");
    
    if(  (nse2d.getImpulsResidual() - 29.42572519 > 1e-6)
      || (nse2d.getMassResidual() - 1.183051948   > 1e-6)
      || (nse2d.getFullResidual() - 29.4494977   > 1e-6)
    )
    { 
      ErrThrow("Residual at iteration 0 is not correct");
    }
    
    for(unsigned int k=1;; k++)
    {
      Output::print<1>("nonlinear step " , setw(3), k-1, "\t",
        nse2d.getResiduals());
      nse2d.solve();
      
      if(TDatabase::ParamDB->PROBLEM_TYPE == 3)
        break;
      // checking the first nonlinear iteration
      if( (k==2) && ( (nse2d.getImpulsResidual() - 0.1984203854 >1e-6)        
         ||  (nse2d.getMassResidual() - 6.936837905e-16     > 1e-6)
         ||  (nse2d.getFullResidual() - 0.1984203854        >1e-6) ) )
      {
        ErrThrow("Residual at iteration 1 is not correct");
      }
      
      nse2d.assemble_nonlinear_term();;
      if(nse2d.stopIt(k))
        break;      
    }
    nse2d.output();
    double err[4];
    TAuxParam2D NSEaux_error;
    MultiIndex2D NSAllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D *velocity_space = &nse2d.get_velocity_space();
    const TFESpace2D *pressure_space = &nse2d.get_pressure_space();
    
    const TFEFunction2D *u1 = nse2d.get_velocity_component(0);
    const TFEFunction2D *u2 = nse2d.get_velocity_component(1);
    
    u1->GetErrors(nse2d.get_example().get_exact(0), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &velocity_space, err);
    // errors in second velocity component
    u2->GetErrors(nse2d.get_example().get_exact(1), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &velocity_space, err + 2);
    // check the L2-error of the velcoity
    if( fabs( sqrt(err[0]*err[0] + err[2]*err[2]) - 6.37755e-05) > 1e-6 )
    {
      ErrThrow("Program 1: L2 norm of velocity is not correct.");
    }
    // check the H1-error of the velcoity
    if( fabs( sqrt(err[1]*err[1] + err[3]*err[3]) - 0.00661116) > 1e-6 )
    {
      ErrThrow("Program 1: H1-semi norm of velocity not correct.");
    }
    
    // errors in pressure
    const TFEFunction2D &p = nse2d.get_pressure();
    p.GetErrors(nse2d.get_example().get_exact(2), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &pressure_space, err);
    // check the L2-error of the pressure
    if( fabs(err[0] - 0.00190367) > 1e-6)
    {
      ErrThrow("Program 1: L2 norm of pressure is not correct.");
    }
    // check the H1-error of the pressure
    if(fabs( err[1] - 0.17806) > 1e-6 )
    {
      ErrThrow("Program 1: H1-norm of pressure is not correct.");
    } 
    Output::print<1>("Elements: Q_2/P_1^disc: test passed for the Direct Solver: ");
  } // end program 1
  //===================================================================================
  /** @brief Multigrid Test: for P2/P1 elements**/
  //===================================================================================
  
  return 0;
}
