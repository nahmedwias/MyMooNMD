/**
 * @brief A test program for the solving of NSE2D problems using Multigrid.
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
 * @todo FIXME This test seems to test a direct solver, although it is named a
 * multigrid test. I suppose the reason is that it was not correctly transferred
 * to using the new Database. Fix this!
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
  ErrThrow("This test is named multigrid but actually tests a direct solver!"
      "Please set up a new (NSE2D) multigrid test in a better framework.")
  // test 1;
  //===================================================================================
  /** @brief Multigrid Test: for P2/P1 elements**/
  //===================================================================================
  {

    //  declaration of databases
    TDatabase Database;
    TFEDatabase2D FEDatabase;

    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(ParameterDatabase::default_nonlinit_database());
    db["problem_type"].set<size_t>(5);
    db["example"] = 2;
    
    db.add("refinement_n_initial_steps", (size_t) 4, "");
    db.add("multigrid_n_levels", (size_t) 3, "");

    db["nonlinloop_maxit"] = 100;
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_slowfactor"] = 1.;

    // default construct a domain object
    TDomain domain(db);

    TDatabase::ParamDB->RE_NR=1;
    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->NSTYPE = 4;
    
    TDatabase::ParamDB->LAPLACETYPE = 0;
    
    TDatabase::ParamDB->VELOCITY_SPACE = 2;
    // if pressure space is -4711, the depending on 
    // velocity space, the pressure space is auto chosen in
    // the NSE2D class or set to 1 for triangles
    TDatabase::ParamDB->PRESSURE_SPACE = -4711; 
    
    
    TDatabase::ParamDB->SOLVER_TYPE = 1;
    TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE= 10;
    TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE = 1e-11;
    TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE = 0.0;
    TDatabase::ParamDB->SC_GMRES_RESTART= 20;
    TDatabase::ParamDB->SC_MG_CYCLE_SADDLE= 0;
    TDatabase::ParamDB->SC_PRE_SMOOTH_SADDLE= 2;
    TDatabase::ParamDB->SC_POST_SMOOTH_SADDLE= 2;
    TDatabase::ParamDB->SC_COARSE_MAXIT_SADDLE= 100;
    TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SADDLE= 0.1;
    TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE= 0;
    TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE= 0;
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE= 0.8;
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE= 0.8;
    TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE= 5;
    TDatabase::ParamDB->SC_SOLVER_SADDLE= 16;
    TDatabase::ParamDB->SC_SMOOTHER_SADDLE= 4;
    TDatabase::ParamDB->SC_MG_TYPE_SADDLE= 0;
    TDatabase::ParamDB->SC_COARSE_SMOOTHER_SADDLE = 17;
    
    // possibly parameters in the database
    Database.CheckParameterConsistencyNSE();
    // the domain is initialised with default description and default
    // initial mesh
    domain.Init((char*)"Default_UnitSquare", (char*)"TwoTriangles");
    // refine grid up to the coarsest level
    size_t n_ref = domain.get_n_initial_refinement_steps();
    for(size_t i=0; i < n_ref; i++)
    {
      domain.RegRefineAll();
    }

    //=========================================================================
    // creat an object 
    NSE2D nse2d(domain, db);
    // assemble all 
    nse2d.assemble();
    // check stopping criterion
    nse2d.stopIt(0);

    if(  (fabs(nse2d.getImpulsResidual() - 27.96300423) > 1e-6)
      || (fabs(nse2d.getMassResidual()   - 0.442647637)   > 1e-6)
      || (fabs(nse2d.getFullResidual()   - 27.96650752)   > 1e-6) )
    {
      ErrThrow("Residual at iteration 0 is not correct");
    }
    
    // nonlinear iterations     
    for(int k=0; ; k++)
    {
      nse2d.solve();
      Output::print<1>("nonlinear step " , setw(3), k, "\t",
        setprecision(10),nse2d.getResiduals());
      
      if( (k==1) 
        && ((fabs(nse2d.getImpulsResidual()  - 0.2040084814) > 1e-6)
        || (fabs(nse2d.getMassResidual()  - 4.095420856e-14)   > 1e-6)
        || (fabs(nse2d.getFullResidual()     - 0.2040084814)   > 1e-6)) )
     {
       ErrThrow("Residual at iteration ", k, " is not correct");
     }
      if(db["problem_type"].is(3))
        break;
      
      nse2d.assemble_nonlinear_term();
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
    if( fabs( sqrt(err[0]*err[0] + err[2]*err[2]) - 7.5635e-05) > 1e-6 )
    {
      ErrThrow("Program 1: L2 norm of velocity is not correct.");
    }
    // check the H1-error of the velcoity
    if( fabs( sqrt(err[1]*err[1] + err[3]*err[3]) - 0.00973175) > 1e-6 )
    {
      ErrThrow("Program 1: H1-semi norm of velocity not correct.");
    }
    
    // errors in pressure
    const TFEFunction2D &p = nse2d.get_pressure();
    p.GetErrors(nse2d.get_example().get_exact(2), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &pressure_space, err);
    // check the L2-error of the pressure
    if( fabs(err[0] - 0.00194342) > 1e-6)
    {
      ErrThrow("Program 1: L2 norm of pressure is not correct.");
    }
    // check the H1-error of the pressure
    if(fabs( err[1] - 0.228822) > 1e-6 )
    {
      ErrThrow("Program 1: H1-norm of pressure is not correct.");
    } 
    Output::print<1>("Elements: P2/P1: test passed for multigrid: ");
  } // end program 1
  
  //===================================================================================
  /** @brief Multigrid Test: for Q2/P1^disc elements**/
  //===================================================================================
  {
    //  declaration of databases
    TDatabase Database;
    TFEDatabase2D FEDatabase;

    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(ParameterDatabase::default_nonlinit_database());
    db["problem_type"].set<size_t>(5);
    db["example"] = 2;
    
    db["nonlinloop_maxit"] = 100;
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_slowfactor"] = 1.;

    db.add("refinement_n_initial_steps", (size_t) 4, "");
    db.add("multigrid_n_levels", (size_t) 2, "");

    // default construct a domain object
    TDomain domain(db);

    TDatabase::ParamDB->RE_NR=1;
    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE = 0;

    TDatabase::ParamDB->VELOCITY_SPACE = 12;
    // Q_2/P_1^disc
    TDatabase::ParamDB->PRESSURE_SPACE = -4711; 
    
    
    TDatabase::ParamDB->SOLVER_TYPE = 1;
    TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE= 10;
    TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE = 1e-11;
    TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE = 0.0;
    TDatabase::ParamDB->SC_GMRES_RESTART= 20;
    TDatabase::ParamDB->SC_MG_CYCLE_SADDLE= 0;
    TDatabase::ParamDB->SC_PRE_SMOOTH_SADDLE= 2;
    TDatabase::ParamDB->SC_POST_SMOOTH_SADDLE= 2;
    TDatabase::ParamDB->SC_COARSE_MAXIT_SADDLE= 100;
    TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SADDLE= 0.1;
    TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE= 0;
    TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE= 0;
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE= 0.8;
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE= 0.8;
    TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE= 5;
    TDatabase::ParamDB->SC_SOLVER_SADDLE= 16;
    TDatabase::ParamDB->SC_SMOOTHER_SADDLE= 2;
    TDatabase::ParamDB->SC_MG_TYPE_SADDLE= 0;
    TDatabase::ParamDB->SC_COARSE_SMOOTHER_SADDLE = 17;
    
    // possibly parameters in the database
    Database.CheckParameterConsistencyNSE();

    // the domain is initialised with default description and default
    // initial mesh
    domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");
    // refine grid
    size_t n_ref = domain.get_n_initial_refinement_steps();
    for(size_t i=0; i< n_ref; i++)
    {
      domain.RegRefineAll();
    }

    //=========================================================================
    // creat an object 
    NSE2D nse2d(domain, db);
    // assemble all 
    nse2d.assemble();
    // check stopping criterion
    nse2d.stopIt(0);
    
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
      
      if(db["problem_type"].is(3))
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
    Output::print<1>("Elements: Q_2/P_1^disc: test passed for the Multigrid: ");
  } // end program 1
  
  return 0;
}
