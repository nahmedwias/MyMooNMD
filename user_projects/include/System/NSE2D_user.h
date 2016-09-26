/** ***************************************************************************
 *
 * @name   NSE2D
 * @brief  store everything needed to solve a Navier-Stokes problem
 *
 *         Store matrix, right hand side, FE spaces, FE functions and the 
 *         solution vector of a Stokes problem. This wraps up everything which 
 *         is necessary to solve a Stokes problem in 2D.
 *
 * @author     Ulrich Wilbrandt
 * @date       06.09.13
 *
 ******************************************************************************/

#ifndef NSE2D_H_
#define NSE2D_H_

#include <FEVectFunct2D.h>
#include <BlockFEMatrix.h>
#include <BlockVector_user.h>
#include <Residuals.h>
#include <ParameterDatabase.h>
#include <Solver.h>
#include <Example_NSE2D.h>
#include <PostProcessing2D.h>
#include <MainUtilities.h> // FixedSizeQueue
#include <PostProcessing2D.h>
#include <utility>
#include <array>


class NSE2D
{
  public:

    
    /** @brief assemble matrix, 
     * 
     * This assembles everything which is not related to the nonlinear term.
     * I.e. it assembles a Stokes matrix.
     */
    void assemble(TFEFunction2D* rho_field=nullptr,
                  TFEFunction2D* mu_field=nullptr);
    
    /** @brief assemble nonlinear term
     * 
     * The matrix blocks to which the nonlinear term contributes are reset to 
     * zero and then completely reassembled, including the linear and nonlinear
     * terms. If this->assemble() has been called before, the matrix is now set
     * up correctly. 
     */
    void assemble_nonlinear_term(TFEFunction2D* rho_field=nullptr,
                                 TFEFunction2D* mu_field=nullptr);

    /* @brief initialize phase fraction, rho_field and mu_field FEfunctions
     * ( = same value in all the domain)
     */
    void update_multiphase(BlockVector new_phase_fraction);



    /** @brief Set the solution vector of a NSE2D object
    * @details can be used to prescribe a non-zero initial solution
    * or to adapt solution between iterations (case of multiphase flow)
    * ATTENTION: it is dangerous, no check if new_solution has a different block structure
    * or length...!!!
    */
    void set_solution(BlockVector new_solution)
    { this->get_solution() = new_solution; }

    /** @brief Update solution vector of a NSE2D object
    * @details used exclusively for multiphase flow (weight = phase fraction,
    * or level set function or any marker function)
    */
    void update_solution(BlockVector weight_vector);

//    TFEFunction2D & get_rho_field()
//    { return this->systems.front().rho_field; }
//
//    TFEFunction2D & get_mu_field()
//    { return this->systems.front().mu_field; }
//
//    TFEFunction2D & get_phase_fraction()
//    { return this->systems.front().phase_fraction; }


};



#endif /* NSE2D_H_ */
