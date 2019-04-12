/**
 * @brief A test program for a time dependent CD2D problem
 * algorithm "linear Crank-Nicolson FEM-FCT" as implemented
 * in the namespace AlgebraicFluxCorrection.
 * Uses LeVeque's rotating bodies example.
 *
 * @date 2016/01/13
 * @author Clemens Bartsch, Naveed Ahmed
 */
#include "all_defines_external_libraries.h"
#include <AlgebraicFluxCorrection.h>
#include "TimeConvectionDiffusion.h"
#include <Database.h>
#include <FEDatabase2D.h>
#include <TimeDiscRout.h>
#include <MainUtilities.h>
#include <ConvDiff.h>
#include <AuxParam2D.h>
#include <TimeDiscretizations.h>

void testCN(TimeConvectionDiffusion<2> &tcd, int m)
{
  double errors[5];
  errors[0]=errors[1]=errors[2]=errors[3]=0.;
  TAuxParam2D aux;
  MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
  const TFEFunction2D& function = tcd.get_function();
  const TFESpace2D* space = function.GetFESpace2D().get();

  function.GetErrors(tcd.get_example().get_exact(0), 3, AllDerivatives, 4,
                     conv_diff_l2_h1_linf_error<2>,
                     tcd.get_example().get_coeffs(), &aux, 1, &space, errors);
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

void time_integration(int td, TimeConvectionDiffusion<2>& tcd, TimeDiscretization& tss)
{
  TDatabase::TimeDB->TIME_DISC = td;
  TDatabase::TimeDB->CURRENTTIME = tcd.get_db()["time_start"];

  tcd.assemble_initial_time();

  int step=0;
  testCN(tcd, step);

  double end_time = tcd.get_db()["time_end"];
  while(TDatabase::TimeDB->CURRENTTIME < end_time -1e-10)
  {
    step ++;
    TDatabase::TimeDB->INTERNAL_STARTTIME
       = TDatabase::TimeDB->CURRENTTIME;
    tss.set_time_disc_parameters();
    SetTimeDiscParameters(1);

    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    TDatabase::TimeDB->CURRENTTIME += tau;

    Output::print<1>("\nCURRENT TIME: ",
         TDatabase::TimeDB->CURRENTTIME);
    tcd.assemble();
    tcd.solve();


    testCN(tcd, step);
  }

}

 //test crank-nicolson linear fem-fct scheme
int main(int, char**)
{
  {
    TDatabase Database;
    TFEDatabase2D FEDatabase;
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database(),true);
    db["example"] = 3;
    db["reynolds_number"] = 1e-20;

    db.merge(AlgebraicFluxCorrection::default_afc_database(), true);
    db["algebraic_flux_correction"].set("fem-fct-cn");

    db["space_discretization_type"] = "galerkin";
    db["time_discretization"] = "crank_nicolson";
    db["time_start"]=0.;
    db["time_end"]=0.02;
    db["time_step_length"] = 0.001;
    TDatabase::ParamDB->ANSATZ_ORDER=1;

    db["time_end"] = 0.02;
    TDatabase::TimeDB->TIMESTEPLENGTH = 0.001;

    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    db.add("refinement_n_initial_steps", (size_t) 5,"");
    TDomain domain(db);
    SetTimeDiscParameters(0);
    // refine grid
    domain.refine_and_get_hierarchy_of_collections(db);

    db.add("solver_type", "direct", "", {"direct", "petsc"});
    TimeConvectionDiffusion<2> tcd(domain, db);
    
    TimeDiscretization& tss = tcd.get_time_stepping_scheme();
    tss.current_step_ = 0;
    tss.set_time_disc_parameters();
    
    time_integration(2,tcd, tss);
    
#ifdef PARMOON_WITH_PETSC
    db["solver_type"] = "petsc";
    TimeConvectionDiffusion<2> tcd_petsc(domain, db);
    time_integration(2, tcd_petsc, tss);
#endif // PARMOON_WITH_PETSC
  }
}



