/**
 * @brief A test program to test Time_CD2D program
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Time_CD2D.h>
#include <TimeDiscRout.h>
#include <MainUtilities.h>
#include <list>


void testBE()
{
  
}

void testCN(Time_CD2D &tcd, int m)
{
  double errors[5];
  errors[0]=errors[1]=errors[2]=errors[3]=0.;
  TAuxParam2D aux;
  MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
  const TFEFunction2D& function = tcd.get_function();
  const TFESpace2D* space = function.GetFESpace2D();
  
  function.GetErrors(tcd.get_example().get_exact(0), 3, AllDerivatives, 4,
                       SDFEMErrors, tcd.get_example().get_coeffs(), &aux, 1,
                       &space, errors);
  double eps = 1E-6;
  if(m==0)
  {
    if( fabs(errors[0] - 3.360188e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct. ", errors[0]);
    if( fabs(errors[1] - 2.521492e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==1)
  {
    if( fabs(errors[0] - 1.541462e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 2.647214e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==2)
  {
    if( fabs(errors[0] - 2.207736e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 2.782276e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==3)
  {
    if( fabs(errors[0] - 2.116156e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 2.924973e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  } 
  else if(m==18)
  {
    if( fabs(errors[0] - 4.577653e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.192082e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==19)
  {
    if( fabs(errors[0] - 4.812381e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.509557e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==20)
  {
    if( fabs(errors[0] - 5.059094e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.843309e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
}

void time_integration(int td, Time_CD2D& tcd)
{
  TDatabase::TimeDB->TIME_DISC = td;
  
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  
  tcd.assemble_initial_time();
  
  int step=0;
  testCN(tcd, step);
  
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
    tcd.assemble();
    tcd.solve();
    
    tcd.descale_stiffness(tau, TDatabase::TimeDB->THETA1);

    testCN(tcd, step);
  }
  
}

int main(int argc, char* argv[])
{
  
  // test with Crank Nicolson euler 
  {
    TDatabase Database;
    TFEDatabase2D FEDatabase;
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db["example"] = 0;
    db.add("reynolds_number",1,"");

    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->ANSATZ_ORDER=1;
    
    TDatabase::TimeDB->STARTTIME=0;
    TDatabase::TimeDB->ENDTIME=1;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.05;
    //  declaration of databases
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    TDomain domain(db);
    SetTimeDiscParameters(0);
    // some parameters
       
    for(int i=0; i< 5; ++i)
    domain.RegRefineAll();
    
    Time_CD2D tcd(domain, db);
    time_integration(2,tcd);
  }
  
  { 
    TDatabase Database;
    TFEDatabase2D FEDatabase;
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(ParameterDatabase::default_output_database());
    db["example"] = 0;
    db.add("reynolds_number",1,"");

    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->ANSATZ_ORDER=1;
    
    TDatabase::TimeDB->STARTTIME=0;
    TDatabase::TimeDB->ENDTIME=1;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.05;
    
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    TDomain domain(db);
    SetTimeDiscParameters(0);
    
    for(int i=0; i< 5; ++i)
      domain.RegRefineAll();
    TDatabase::ParamDB->SC_GMRES_RESTART= 20;
    
    Time_CD2D tcd(domain, db);
    time_integration(2,tcd);
  }
  Output::print<1>("TEST SUCCESFULL: ");
  return 0;
}
