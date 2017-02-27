#include <VOF_TwoPhase2D.h>


/** @brief constructor*/
VOF_TwoPhase2D::VOF_TwoPhase2D(const TDomain& domain,
                               const ParameterDatabase& param_db_tnse,
                               const ParameterDatabase& param_db_tcd
                               )
: example_tnse2d_(param_db_tnse),
  tnse2d_(domain,param_db_tnse,example_tnse2d_),
  example_tcd2d_(param_db_tcd),
  phaseconvection2d_(domain,param_db_tcd,example_tcd2d_),
  tnse_variable_fluid_(param_db_tnse["dimensional_nse"]),
  solve_convection_(param_db_tcd["solve_cd"]),
  nse2cd_coupling_(param_db_tcd["coupling_nse_cd"]),
  cd2nse_coupling_(param_db_tnse["coupling_cd_nse"]),
  rho_vector_(phaseconvection2d_.get_solution()),
  mu_vector_(phaseconvection2d_.get_solution())
{
  /* First, check and set consistent parameters  */
  manage_example_parameters();

  /* Properties of liquid phase, eg. rhol=1000 and mul=1.e-3 */
  rhol_ = param_db_tnse["fluid_density"];
  mul_  = param_db_tnse["fluid_dynamic_viscosity"];

  /* Initialize vectors to constant values of one phase, e.g. liquid  */
  this->rho_vector_ = rhol_;
  this->mu_vector_ = mul_;


}

/** Check and set input parameters */
void VOF_TwoPhase2D::manage_example_parameters()
{
  /********************************************************************
     * Check example number consistency!
     ********************************************************************/
  if (example_tnse2d_.get_database()["example"].value_as_string()
        != example_tcd2d_.get_database()["example"].value_as_string())
  {
      ErrThrow("The example number of TNSE and TCD should be the same!");
  }
  else
    example_number_ = example_tnse2d_.get_database()["example"];

  /********************************************************************
   * MANAGE PARAMETERS FOR BENCHMARK PROBLEMS!
   ********************************************************************/
  switch(example_number_)
  {
    case 10:
      this->nse2cd_coupling_  = false;
      this->cd2nse_coupling_  = false;
      break;
    case 20: case 21: case 22:
      this->tnse_variable_fluid_ = true;
      this->nse2cd_coupling_     = true;
      this->cd2nse_coupling_     = false;
      break;
    case 30: case 31: case 32:
//      this->tnse_variable_fluid_ = true;
      this->nse2cd_coupling_  = false;
//      this->cd2nse_coupling_  = true;
//      this->solve_convection_ = true;
      break;
    case 40: case 41:
//      this->tnse_variable_fluid_ = true;
//      this->nse2cd_coupling_  = true;
//      this->cd2nse_coupling_  = true;
//      this->solve_convection_ = true;
      break;
  }


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

/* Print some info, mostly useful after the constructor */
void VOF_TwoPhase2D::output_initial_info()
{
  Output::print<3>("The velocity space is ", TDatabase::ParamDB->VELOCITY_SPACE);
  Output::print<3>("The pressure space is ", TDatabase::ParamDB->PRESSURE_SPACE);
  Output::print<3>("The ansatz   space is ", TDatabase::ParamDB->ANSATZ_ORDER);
  Output::print<3>("Convection_example number     ",
                   this->example_tcd2d_.get_database()["example"]);
  Output::print<3>("Time NSE Example number       ",
                   this->example_tnse2d_.get_database()["example"]);
  Output::print<3>("Fluid properties are non constant : "
      + std::to_string(this->tnse_variable_fluid_));
  Output::print<3>("TCD is solved            : " + std::to_string(this->solve_convection_));
  Output::print<3>("Coupling NSE >> CD is    : " + std::to_string(this->nse2cd_coupling_));
  Output::print<3>("Coupling CD >> NSE is    : " + std::to_string(this->cd2nse_coupling_));
  Output::print<3>("The density   of liquid is ", this->rhol_);
  Output::print<3>("The viscosity of liquid is ", this->mul_);
  Output::print<3>("The density   of gas    is ", this->rhog_);
  Output::print<3>("The viscosity of gas    is ", this->mug_);
}

