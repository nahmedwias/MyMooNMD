#include <VOF_TwoPhase3D.h>

#ifdef _MPI
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif


/** @brief constructor*/
VOF_TwoPhase3D::VOF_TwoPhase3D(const std::list< TCollection* > grid_collections,
                               const ParameterDatabase& param_db_tnse,
                               const ParameterDatabase& param_db_tcd
#ifdef _MPI
                               , int maxSubDomainPerDof
#endif
                               )
: example_tnse3d_(param_db_tnse),
  tnse3d_(grid_collections,param_db_tnse,example_tnse3d_
#ifdef _MPI
          , maxSubDomainPerDof
#endif
          ),
  example_tcd3d_(param_db_tcd),
  phaseconvection3d_(grid_collections,param_db_tcd,example_tcd3d_
#ifdef _MPI
                     , maxSubDomainPerDof
#endif
                    )
//  /* Properties of liquid phase, eg. rhol=1000 and mul=1.e-3 */
//  rhol_(param_db_tnse["fluid_density"]),
//  mul_(param_db_tnse["fluid_dynamic_viscosity"]),
//  /* List of booleans */
//  tnse_variable_fluid_(param_db_tnse["dimensional_nse"]),
//  solve_convection_(param_db_tcd["solve_cd"]),
//  nse2cd_coupling_(param_db_tcd["coupling_nse_cd"]),
//  cd2nse_coupling_(param_db_tnse["coupling_cd_nse"]),
//  /* Copy the block vector structure */
//  rho_vector_(phaseconvection2d_.get_solution()),
//  mu_vector_(phaseconvection2d_.get_solution()),
//  unity_vector_(phaseconvection2d_.get_solution()),
//  /* Construct FeFunctions */
//  rho_fefunction_(&this->phaseconvection2d_.get_space(),(char*)"r",(char*)"r",
//                  rho_vector_.block(0),rho_vector_.length(0)),
//  mu_fefunction_(&this->phaseconvection2d_.get_space(),(char*)"m",(char*)"m",
//                  mu_vector_.block(0),mu_vector_.length(0))
{
//  /* Initialize unity vector  */
//  this->unity_vector_ = 1;
//  this->tnse2d_.set_rho_mu_fefunct(&this->rho_fefunction_,
//                                   &this->mu_fefunction_,
//                                   &this->phaseconvection2d_.get_function());
///* at the end of the constructor, rho and mu are constant equal to phase fraction,
// * whatever the examples and booleans are. The call to "update_field_vectors" will
// * calculate their value, depending on the example number
// */
}



///** Check and set input parameters */
//void VOF_TwoPhase2D::manage_example_parameters()
//{
//  /********************************************************************
//     * Check example number consistency!
//     ********************************************************************/
//  if (example_tnse2d_.get_database()["example"].value_as_string()
//        != example_tcd2d_.get_database()["example"].value_as_string())
//  {
//      ErrThrow("The example number of TNSE and TCD should be the same!");
//  }
//  else
//    example_number_ = example_tnse2d_.get_database()["example"];
//
//  /********************************************************************
//   * MANAGE PARAMETERS FOR BENCHMARK PROBLEMS!
//   ********************************************************************/
//  switch(example_number_)
//  {
//    case 10:
//      this->nse2cd_coupling_  = false;
//      this->cd2nse_coupling_  = false;
//      Output::info<5>("Example " + std::to_string(example_number_) +
//                      ":the couplings must be deactivated, but you are free "
//                      "to activate 'tnse_variable_properties' and 'solve tcd'.");
//      break;
//    case 20: case 21: case 22:
//      this->tnse_variable_fluid_ = true;
//      this->nse2cd_coupling_     = true;
//      this->cd2nse_coupling_     = false;
//      this->solve_convection_    = true;
//      Output::info<5>("Example " + std::to_string(example_number_) +
//                      ":all booleans must be true but cd2nse_coupling is false.");
//      break;
//    case 30: case 31: case 32:
////      this->tnse_variable_fluid_ = true;
//      this->nse2cd_coupling_  = false;
////      this->cd2nse_coupling_  = true;
////      this->solve_convection_ = true;
//      Output::info<5>("Example " + std::to_string(example_number_) +
//                      ":coupling nse2cd must be false, but the other booleans"
//                      " are free to be changed.");
//      break;
//    case 40: case 41: case 42: case 43: case 50:
////      this->tnse_variable_fluid_ = true;
////      this->nse2cd_coupling_  = true;
////      this->cd2nse_coupling_  = true;
////      this->solve_convection_ = true;
//      break;
//  }
//
//  if (this->tnse_variable_fluid_) // using variable rho and mu requires
//  {                               // NSTYPE 3 or 4, LAPLACETYPE 1
//    if (TDatabase::ParamDB->LAPLACETYPE != 1
//        || (TDatabase::ParamDB->NSTYPE != 4
//        &&  TDatabase::ParamDB->NSTYPE != 3))
//    {
//      ErrThrow("In order to assemble TNSE2D with variable fluid properties,"
//          "LAPLACTYPE must be 1 and NSTYPE 3 or 4. If the problem contains"
//          " slip conditions, you must use NSTYPE 4.");
//    }
//  }
//
//  // Set the boolean of TNSE2D taken from VOF
//  this->tnse2d_.set_bool_variable_properties(this->tnse_variable_fluid_);
//}
//
///** Update the vectors mu and rho with the phase fraction vector,
// * depending on the example number */
//void VOF_TwoPhase2D::update_field_vectors()
//{
//  // note : the following two lines are redundant in the initial
//  // time...because already equal in the constructor...
//  // but necessary to recalculate it before updating..
//  // it may be not nice to calculate the property fields by
//  // using only the values at nodes with blockvector...any better way?
//  this->rho_vector_ = this->phaseconvection2d_.get_solution();
//  this->mu_vector_ = this->phaseconvection2d_.get_solution();
//  //At this point, rhol and mul = phase fraction
//  switch(example_number_)
//  {
//    case 10: case 20: case 3:
//    case 21: case 22:
//      this->rho_vector_ = this->rhol_;
//      this->mu_vector_ = this->mul_;
//      // note that in case of standard tnse2d, these values
//      // won't be used or read in the tnse2d object
//      break;
//    case 30: case 31:
//    case 32:
//      // here rho=1 and viscosity = phase field
//      this->rho_vector_ = this->unity_vector_;
//      break;
//    case 41:
//      // rayleigh taylor: visco = 1 and rho=phase field
//      this->mu_vector_ = this->unity_vector_;
//      break;
//    case 43:
//      // drop pressure for CSF test
//      this->rhog_ = TDatabase::ParamDB->P7; // rho of gas = rho2=rhomin
//      this->rho_vector_.scale(this->rhol_-this->rhog_);
//      this->rho_vector_.add_scaled(this->unity_vector_,this->rhog_);
//      this->mu_vector_ = this->unity_vector_;
//      break;
//    case 40: case 42: case 50:
//      // dambreak and 2-phase flows
//      this->rhog_ = TDatabase::ParamDB->P7; // rho of gas = rho2=rhomin
//      this->mug_  = TDatabase::ParamDB->P8; // mu of gas  = mu2 =mumin
////    double density_ratio   = rhol_/rhog_;
////    double viscosity_ratio = mul_/mug_;
//      this->rho_vector_.scale(this->rhol_-this->rhog_);
//      this->rho_vector_.add_scaled(this->unity_vector_,this->rhog_);
//      this->mu_vector_.scale(this->mul_-this->mug_);
//      this->mu_vector_.add_scaled(this->unity_vector_,this->mug_);
//      break;
//    default:
//      ErrThrow("Unknown example number");
//      break;
//  }
//  // transfer fe functions to tnse2d after updates
//  this->tnse2d_.set_rho_mu_fefunct(&this->rho_fefunction_,
//                                   &this->mu_fefunction_,
//                                   &this->phaseconvection2d_.get_function());
//}
//
///** Write output of vectors in a file */
//void VOF_TwoPhase2D::output_vectors(std::string filename_phi,
//                                    std::string filename_rho,
//                                    std::string filename_mu)
//{
//  this->phaseconvection2d_.get_solution().write(filename_phi);
//  this->rho_vector_.write(filename_rho);
//  this->mu_vector_.write(filename_mu);
//}
//
///* Print some info, mostly useful after the constructor */
//void VOF_TwoPhase2D::output_initial_info()
//{
//  Output::print<3>("The velocity space is ", TDatabase::ParamDB->VELOCITY_SPACE);
//  Output::print<3>("The pressure space is ", TDatabase::ParamDB->PRESSURE_SPACE);
//  Output::print<3>("The ansatz   space is ", TDatabase::ParamDB->ANSATZ_ORDER);
//  Output::print<3>("Convection_example number     ",
//                   this->example_tcd2d_.get_database()["example"]);
//  Output::print<3>("Time NSE Example number       ",
//                   this->example_tnse2d_.get_database()["example"]);
//  Output::print<3>("Fluid properties are non constant : "
//      + std::to_string(this->tnse_variable_fluid_));
//  Output::print<3>("TCD is solved            : " + std::to_string(this->solve_convection_));
//  Output::print<3>("Coupling NSE >> CD is    : " + std::to_string(this->nse2cd_coupling_));
//  Output::print<3>("Coupling CD >> NSE is    : " + std::to_string(this->cd2nse_coupling_));
//  Output::print<3>("The density   of liquid is ", this->rhol_);
//  Output::print<3>("The viscosity of liquid is ", this->mul_);
//  Output::print<3>("The density   of gas    is ", this->rhog_);
//  Output::print<3>("The viscosity of gas    is ", this->mug_);
//}
//
