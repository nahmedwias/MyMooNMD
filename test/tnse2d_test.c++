/**
 * @brief A test program to test Time Stokes/Navier Stokes program
 * 
 * Residuals at the first nonlinear iteration are compared with the 
 * Volker main program (current MooNMD code)
 * Also the errors are compared with that code
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Time_NSE2D.h>
#include <TimeDiscRout.h>
#include <MainUtilities.h>
#include <list>





void testBE()
{
  
}

void testCN(Time_NSE2D &tnse, int m)
{
  cout<<"testCN: " << m<< "  " ;
  double errors[6];
  errors[0]=errors[1]=errors[2]=errors[3]=errors[4]=errors[5]=0.;  
  TAuxParam2D aux;
  MultiIndex2D NSAllDerivatives[3] = {D00, D10, D01};
  const TFESpace2D *velocity_space = &tnse.get_velocity_space();
  const TFESpace2D *pressure_space = &tnse.get_pressure_space();
  
  TFEFunction2D * u1 = tnse.get_velocity_component(0);
  TFEFunction2D * u2 = tnse.get_velocity_component(1);
  TFEFunction2D p = tnse.get_pressure();
  
    
  u1->GetErrors(tnse.get_example().get_exact(0), 3, NSAllDerivatives, 2, L2H1Errors,nullptr,
                &aux,1, &velocity_space,errors);   
  
  u2->GetErrors(tnse.get_example().get_exact(1), 3, NSAllDerivatives, 2, L2H1Errors,nullptr,
                  &aux,1, &velocity_space,errors+2);
  p.GetErrors(tnse.get_example().get_exact(2), 3, NSAllDerivatives, 2, L2H1Errors, nullptr, 
               &aux, 1, &pressure_space, errors+4);
  
  double eps = 1E-6, err=0;
  if(m==1)
  {
    err = sqrt(errors[0]*errors[0] + errors[2]*errors[2]);
    if( fabs(err - 0.00110029) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of velocity is not correct. ", 
               err);
    err = sqrt(errors[1]*errors[1] + errors[3]*errors[3]);
    if( fabs(err -  0.014616) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of velocity is not correct.");
    if( fabs(errors[4] - 0.0198907) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of pressure is not correct.");
    if( fabs(errors[5] - 0.0974169) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of pressure is not correct.");
  }
  else if(m==2)
  {
    err = sqrt(errors[0]*errors[0] + errors[2]*errors[2]);
    if( fabs(err - 0.00222401) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of velocity is not correct. ", 
               err);
    err = sqrt(errors[1]*errors[1] + errors[3]*errors[3]);
    if( fabs(err - 0.0290912) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of velocity is not correct.");
    
    if( fabs(errors[4] - 0.0208086) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of pressure is not correct.");
    if( fabs(errors[5] - 0.153391) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of pressure is not correct.");
  }
  else if(m==18)
  {
    err = sqrt(errors[0]*errors[0] + errors[2]*errors[2]);
    if( fabs(err - 0.0175156) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of velocity is not correct. ", 
               err);
    err = sqrt(errors[1]*errors[1] + errors[3]*errors[3]);
    if( fabs(err - 0.228106) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of velocity is not correct.");
    
    if( fabs(errors[4] - 0.0446204) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of pressure is not correct.");
    if( fabs(errors[5] - 1.109628) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of pressure is not correct.");
  }
  else if(m==19)
  {
    err = sqrt(errors[0]*errors[0] + errors[2]*errors[2]);
    if( fabs(err - 1.818933e-02) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of velocity is not correct. ", 
               err);
    err = sqrt(errors[1]*errors[1] + errors[3]*errors[3]);
    if( fabs(err - 2.368664e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of velocity is not correct.");

    if( fabs(errors[4] -   4.584916e-02) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of pressure is not correct.", errors[4]);
    if( fabs(errors[5] - 1.152722e+00) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of pressure is not correct.");
  }
  else if(m==20)
  {
    err = sqrt(errors[0]*errors[0] + errors[2]*errors[2]);
    if( fabs(err - 1.881760e-02) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of velocity is not correct. ", 
               err);
    err = sqrt(errors[1]*errors[1] + errors[3]*errors[3]);
    if( fabs(err - 2.450352e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of velocity is not correct.");

    if( fabs(errors[4] - 4.699491e-02) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm of pressure is not correct.", errors[4]);
    if( fabs(errors[5] - 1.192964e+00) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm of pressure is not correct.");
  }
}

void time_integration(int td, Time_NSE2D& tnse)
{
  TDatabase::TimeDB->TIME_DISC = td;
  
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  
  tnse.assemble_initial_time();
  
  int step=0;
  int image=0;
  
  while(TDatabase::TimeDB->CURRENTTIME < 
    TDatabase::TimeDB->ENDTIME-1e-10)
  {
    step ++;
    TDatabase::TimeDB->INTERNAL_STARTTIME 
       = TDatabase::TimeDB->CURRENTTIME;
    SetTimeDiscParameters(1);
    
    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    TDatabase::TimeDB->CURRENTTIME += tau;
    
    Output::print<1>("\nCURRENT TIME: ", 
                     TDatabase::TimeDB->CURRENTTIME);
    // prepare the right hand side vector
    // only needed once per time step
    tnse.assemble_rhs();
    // assemble the nonlinear matrices
    tnse.assemble_nonlinear_term();
    // prepare the matrices for defect computations
    // and solvers
    tnse.assemble_system();
    // nonlinear iteration
    for(unsigned int k=0;; k++)
    {
      if(tnse.stopIte(k))
        break;
      tnse.solve();
      // assemble the nonlinear matrices 
      tnse.assemble_nonlinear_term();
      // prepare the matrices for next nonlinear iteration
      tnse.assemble_system();         
    }
    // post processing: error computations
    // and solutions for visualization
    // tnse.output(step,image);
    testCN(tnse, step);
  }
  
}

int main(int argc, char* argv[])
{
  
  // test with Crank Nicolson euler 
  {
    TDatabase Database;
    TFEDatabase2D FEDatabase;
    
    TDatabase::ParamDB->MEASURE_ERRORS=1;
    TDatabase::ParamDB->EXAMPLE =101;
    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->RE_NR = 1;
    TDatabase::ParamDB->VELOCITY_SPACE=12;
    TDatabase::ParamDB->PRESSURE_SPACE  =-4711;
    // FLOW_PROBLEM_TYPE: 3 Stokes
    // FLOW_PROBLEM_TYPE: 5 NSE
    TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 5;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE= 1;
    TDatabase::ParamDB->NSE_NONLINEAR_FORM= 0;
    TDatabase::ParamDB->SOLVER_TYPE = 2;
    
    TDatabase::TimeDB->TIME_DISC= 2;
    TDatabase::TimeDB->STARTTIME=0;
    TDatabase::TimeDB->ENDTIME=1;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.05;
    TDatabase::ParamDB->UNIFORM_STEPS = 1;
    TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE = 0;
    TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE = 1e-10;
    TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 2;
    //  declaration of databases
    TDomain domain;
    SetTimeDiscParameters(0);
    // some parameters
       
    // the domain is initialised with default description and default
    // initial mesh
    domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");
    for(int i=0; i< TDatabase::ParamDB->UNIFORM_STEPS; ++i)
    domain.RegRefineAll();
    TDatabase::ParamDB->SOLVER_TYPE=2;
    
    Time_NSE2D tnse(domain);
    // direct SOLVER_TYPE = 2
    time_integration(2,tnse);
  }
  
  
  Output::print<1>("TEST SUCCESFULL: ");
  return 0;
}
