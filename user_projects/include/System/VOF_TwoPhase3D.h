/** ************************************************************************
 * 
 * @name         VOF_TwoPhase3D
 * @brief        store everything needed to solve a two-phase problem
 *               modelled with Volume Of Fluid
 *               
 * @author       Najib Alia
 * @History      25.04.2017
**************************************************************************/

#ifndef __VOF_TwoPhase3D__
#define __VOF_TwoPhase3D__

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <FESpace3D.h>
#include <FEFunction3D.h>
#include <Time_NSE3D.h>
#include <Time_CD3D.h>
#include <Example_TimeNSE3D.h>
#include <Example_TimeCD3D.h>

#include <MainUtilities.h>
#include <ParameterDatabase.h>
#include <vector>
#include <deque>
#include <utility>
#include <array>

class VOF_TwoPhase3D
{
  public:
    Example_TimeNSE3D example_tnse3d_;
    Time_NSE3D tnse3d_;
    Example_TimeCD3D example_tcd3d_;
    Time_CD3D phaseconvection3d_;

    /* example number of vof=example of tnse=tcd */
    int example_number_;
    /* rhol = constant density of liquid phase
     * default value is 1 */
    double rhol_ = 1;
    /* mul = constant dyn. visco of liquid phase
     * default value is 1 */
    double mul_  = 1;
    /* rhog = constant density of gas phase
     * default value is 0 */
    double rhog_ = 0;
    /* mug = constant dyn. visco of gas phase
     * default value is 0 */
    double mug_  = 0;

    /* Boolean parameters which activates and
     * controls different features of this class
     */
    bool tnse_variable_fluid_ = false;
    bool solve_convection_    = false;
    bool nse2cd_coupling_     = false;
    bool cd2nse_coupling_     = false;

    /* Vector equal to property fields at the nodes */
    BlockVector rho_vector_;
    BlockVector mu_vector_;
    BlockVector unity_vector_; // same structure and equals to 1

    /* FEFunction equal to property fields, construction
     * based on the above BlockVectors
     */
    TFEFunction3D rho_fefunction_;
    TFEFunction3D mu_fefunction_;
//
  public:
    /** @brief constructor*/
    VOF_TwoPhase3D(const std::list< TCollection* > grid_collections,
                   const ParameterDatabase& param_db_tnse,
                   const ParameterDatabase& param_db_tcd
#ifdef _MPI
                   , int maxSubDomainPerDof
#endif
                   );

    /*************************************************************/
   /**
    * Special member functions mostly deleted
    * ...needs to be optimized
    */
   //! Delete copy constructor.
    VOF_TwoPhase3D(const VOF_TwoPhase3D&) = delete;

   //! Delete move constructor.
    VOF_TwoPhase3D(VOF_TwoPhase3D&&) = delete;

   //! Delete copy assignment operator.
    VOF_TwoPhase3D& operator=(const VOF_TwoPhase3D&) = delete;

   //! Delete move assignment operator.
    VOF_TwoPhase3D& operator=(VOF_TwoPhase3D&&) = delete;

   //! Default destructor. Most likely causes memory leaks.
   ~VOF_TwoPhase3D() = default;

   /*************************************************************/

   /* Check that the input parameters are consistent,
    * and correct/throw errors if not the case
    */
   void manage_example_parameters();

   /* Update the BlockVectors rho and mu with the phase fraction
    * vector, via the equation: rho = rhol.phi + rhog.(1-phi)
    * idem for mu
    */
   void update_field_vectors();

   /* Write the vectors in a file for output */
   void output_vectors(std::string filename_phi,
                       std::string filename_rho,
                       std::string filename_mu);

   /* Print some info, mostly useful after the constructor */
   void output_initial_info();

};

#endif // __VOF_TwoPhase3D__


