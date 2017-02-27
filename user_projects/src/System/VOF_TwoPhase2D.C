#include <VOF_TwoPhase2D.h>


/** @brief constructor*/
VOF_TwoPhase2D::VOF_TwoPhase2D(const TDomain& domain,
                               const ParameterDatabase& param_db_tnse,
                               const ParameterDatabase& param_db_tcd
                               )
: example_tnse2d_(param_db_tnse), example_tcd2d_(param_db_tcd),
  tnse2d_(domain,param_db_tnse,example_tnse2d_),
  phaseconvection2d_(domain,param_db_tcd,example_tcd2d_)
{
  cout << "Constructor of VOF TWO PHASE 2D IS OK!" << endl;

}

