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
  /* First, check and set consistent parameters  */
  manage_example_parameters();

  /* Properties of liquid phase, eg. rhol=1000 and mul=1.e-3 */
  rhol_ = param_db_tnse["fluid_density"];
  mul_  = param_db_tnse["fluid_dynamic_viscosity"];



  cout << "Constructor of VOF TWO PHASE 2D IS OK!" << endl;
}

/** Check and set input parameters */
void VOF_TwoPhase2D::manage_example_parameters()
{
  // Check example number consistency
  if (example_tnse2d_.get_database()["example"].value_as_string()
        != example_tcd2d_.get_database()["example"].value_as_string())
  {
      ErrThrow("The example number of TNSE and TCD should be the same!");
  }
  else
    example_number_ = example_tnse2d_.get_database()["example"];

  switch(example_number_)
  {
    case 40: case 50:
      rhog_ = TDatabase::ParamDB->P7; // rho of gas = rho2=rhomin
      mug_  = TDatabase::ParamDB->P8; // mu of gas  = mu2 =mumin
//    double density_ratio   = rhol_/rhog_;
//    double viscosity_ratio = mul_/mug_;
      break;

    default:
      break;
  }
}

/** Update the vectors mu and rho with the phase fraction vector */
void VOF_TwoPhase2D::update_field_vectors()
{
  this->rho_vector_ = this->phaseconvection2d_.get_solution();
  this->mu_vector_ = this->phaseconvection2d_.get_solution();
  BlockVector unity = this->rho_vector_; unity = 1;
  this->rho_vector_.scale(this->rhol_-this->rhog_);
  this->mu_vector_.scale(this->mul_-this->mug_);
  this->rho_vector_.add_scaled(unity,this->rhog_);
  this->mu_vector_.add_scaled(unity,this->mug_);
}

/** Write output of vectors in a file */
void VOF_TwoPhase2D::output_vectors(std::string filename_phi,
                                    std::string filename_rho,
                                    std::string filename_mu)
{
  this->phaseconvection2d_.get_solution().write(filename_phi);
  this->rho_vector_.write(filename_rho);
  this->mu_vector_.write(filename_mu);
}





