#include <VOF_TwoPhase2D.h>


/** @brief constructor*/
VOF_TwoPhase2D::VOF_TwoPhase2D(const TDomain& domain,
                               const ParameterDatabase& param_db_tnse,
                               const ParameterDatabase& param_db_tcd
                               )
: example_tnse2d_(param_db_tnse),
  tnse2d_(domain,param_db_tnse,example_tnse2d_),
  example_tcd2d_(param_db_tcd),
  phaseconvection2d_(domain,param_db_tcd,example_tcd2d_)
{
  if (param_db_tnse["example"].value_as_string()
      != param_db_tcd["example"].value_as_string())
  {
    ErrThrow("The example number of TNSE and TCD should be the same!");
  }
  else
  {
    example_number_ = param_db_tnse["example"];
  }

  /* Properties of liquid phase, eg. rhol=1000 and mul=1.e-3 */
  rhol_ = param_db_tnse["fluid_density"];
  mul_  = param_db_tnse["fluid_dynamic_viscosity"];

  if (example_number_ == 40 || example_number_ == 50 )
  {
    rhog_ = TDatabase::ParamDB->P7; // rho of gas = rho2=rhomin
    mug_  = TDatabase::ParamDB->P8; // mu of gas  = mu2 =mumin
//    double density_ratio   = rhol_/rhog_;
//    double viscosity_ratio = mul_/mug_;
  }


  cout << "Constructor of VOF TWO PHASE 2D IS OK!" << endl;

}

