/**
 * @brief A test program to test Time_CD2D program
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Time_CD2D.h>
#include <TimeDiscRout.h>
#include <TimeDiscretizations.h>
#include <MainUtilities.h>
#include <AuxParam2D.h>
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
  tcd.output();
  
  const TFESpace2D* space = function.GetFESpace2D();
  
  function.GetErrors(tcd.get_example().get_exact(0), 3, AllDerivatives, 4,
                       SDFEMErrors, tcd.get_example().get_coeffs(), &aux, 1,
                       &space, errors);
  if(m==1){cout << errors[0] << "  "<<
  errors[1]<< "  " <<endl;}
  double eps = 1E-6;
  if(m==0)
  {cout << "M: " << m << endl;
    if( fabs(errors[0] - 3.360188e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct. ", errors[0]);
    if( fabs(errors[1] - 2.521492e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==1)
  {cout << "M: " << m << endl;
    if( fabs(errors[0] - 1.541462e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 2.647214e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==2)
  {cout << "M: " << m << endl;
    if( fabs(errors[0] - 2.207736e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 2.782276e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==3)
  {cout << "M: " << m << endl;
    if( fabs(errors[0] - 2.116156e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 2.924973e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  } 
  else if(m==18)
  {cout << "M: " << m << endl;
    if( fabs(errors[0] - 4.577653e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.192082e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==19)
  {cout << "M: " << m << endl;
    if( fabs(errors[0] - 4.812381e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.509557e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==20)
  {cout << "M: " << m << endl;
    if( fabs(errors[0] - 5.059094e-03) > eps )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( fabs(errors[1] - 6.843309e-01) > eps )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
}

void time_integration(int td, Time_CD2D& tcd, TimeDiscretization& tss)
{
  TDatabase::TimeDB->TIME_DISC = td;
    
  tcd.assemble_initial_time();
  
  tss.current_step_=0;
  testCN(tcd, tss.current_step_);

  while(tss.current_time_ < tss.get_end_time()-1e-10)
  {
    tss.current_step_++;
    TDatabase::TimeDB->INTERNAL_STARTTIME 
       = tss.current_time_;
    tss.set_time_disc_parameters();
    SetTimeDiscParameters(1);
    
    tss.current_time_ += tss.get_step_length();;
    TDatabase::TimeDB->CURRENTTIME += tss.get_step_length();
    
    Output::print<1>("\nCURRENT TIME: ", tss.current_time_);
    tcd.assemble();
    tcd.solve();
    testCN(tcd, tss.current_step_);
  }
  
}

int main(int argc, char* argv[])
{
  
  // test with Crank Nicolson euler 
  {
    TDatabase Database;
    TFEDatabase2D FEDatabase;
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database());
    db["example"] = 0;
    db["reynolds_number"] = 1;

    db["space_discretization_type"] = "galerkin";
    db["time_discretization"] = "crank_nicolson";
    db["time_step_length"] = 0.05;
    TDatabase::ParamDB->ANSATZ_ORDER=1;
    
    db["time_end"]=1.;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.05;
    //  declaration of databases
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    TDomain domain(db);
    
    SetTimeDiscParameters(0);
    // some parameters
       
    for(int i=0; i< 5; ++i)
    domain.RegRefineAll();
    
    db.add("solver_type", "direct", "", {"direct", "petsc"});
    TDatabase::TimeDB->CURRENTTIME = db["time_start"];
    Time_CD2D tcd(domain, db);
    TimeDiscretization& tss = tcd.get_time_stepping_scheme();
    tss.current_step_ = 0;
    tss.current_time_ = db["time_start"];
    
    time_integration(2, tcd, tss);
    
    // I don't know what has changed but with setting the default parameters in
    // the old database here, it does not work.
    Database.SetDefaultParameters();
    TDatabase::ParamDB->ANSATZ_ORDER=1;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.05;
    
    Output::print("\n\nTesting PETSc\n");
    db["solver_type"] = "petsc";
    TDatabase::TimeDB->CURRENTTIME = db["time_start"];
    tss.current_time_ = db["time_start"];
    Time_CD2D tcd_petsc(domain, db);
    //TimeDiscretization& tss = tcd.get_time_stepping_scheme();
    tss.current_step_ = 0;
    tss.set_time_disc_parameters();
    time_integration(2, tcd_petsc, tss);
  }
  return 0;
}
