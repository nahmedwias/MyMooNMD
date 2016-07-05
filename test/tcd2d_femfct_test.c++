/**
 * @brief A test program for a time dependent CD2D problem
 * algorithm "linear Crank-Nicolson FEM-FCT" as implemented
 * in the namespace AlgebraicFluxCorrection.
 * Uses LeVeque's rotating bodies example.
 *
 * @date 2016/01/13
 * @author Clemens Bartsch, Naveed Ahmed
 */
#include <Time_CD2D.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Time_CD2D.h>
#include <TimeDiscRout.h>
#include <MainUtilities.h>


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
  double eps1 = 1E-6;
  double eps2 = 1E-5;
  if(m==0)
  {
    if( fabs(errors[0] - 0.0705721) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct. ", errors[0]);
    if( fabs(errors[1] - 6.91068) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");

  }
  else if(m==1)
  {
    if( fabs(errors[0] - 0.0707957) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.85906) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==2)
  {
    if( fabs(errors[0] - 0.0711781) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.81221) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==3)
  {
    if( fabs(errors[0] - 0.0713873) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.76706) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==18)
  {
    if( fabs(errors[0] - 0.0791431) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.17222) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==19)
  {
    if( fabs(errors[0] - 0.0798896) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.13685) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==20)
  {
    if( fabs(errors[0] - 0.0799563) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.10199) > eps2 )
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

 //test crank-nicolson linear fem-fct scheme
int main(int argc, char* argv[])
{
  {
    TDatabase Database;
    TFEDatabase2D FEDatabase;
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db["example"] = 3;

    TDatabase::ParamDB->RE_NR = 1e-20;

    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->ANSATZ_ORDER=1;
    TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION = 2;
    TDatabase::ParamDB->FEM_FCT_PRELIMITING = 0;

    TDatabase::ParamDB->SOLVER_TYPE=2;

    TDatabase::TimeDB->STARTTIME=0;
    TDatabase::TimeDB->ENDTIME=0.02;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.001;

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
}



