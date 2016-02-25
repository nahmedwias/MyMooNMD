#ifndef __PRROBUSTTIMENSE__
#define __PRROBUSTTIMENSE__

#include <Time_NSE2D.h>
#include <BlockFEMatrix.h>

class PrRobustTime_NSE2D : Time_NSE2D
{
  protected:
    struct SystemPerGrid
    {
      /** @brief Finite Element space for the projection*/
      TFESpace2D projection_space;
      
      /** @brief projection matrix for reconstruciton*/      
      BlockFEMatrix ProjectionMatrix;
     
      /** @brief Mass Matrix for vector valued space*/
      BlockFEMatrix MassMatrix;
      
      /** @brief right hand side
       * this right hand side vector depends on the 
       * vector space
       */
      BlockVector rhsXh;
      /** @brief Modifed_Mass matrices
       * these blocks comes from the multiplication
       * of the MassMatrix of RT/BDM and the ProjectionMatrix
       * matrix: 
       * (P0*M*P0^T P0*M*P1^T )
       * (P1*M*P0^T P1*M*P1^T )
       */
      BlockFEMatrix Modifed_Mass;
      /** @brief solution vector*/
      BlockVector solutionVV;
      /** @brief Finite element function*/
      TFEFunction2D fefctVV;
      /** @brief constructor */
      SystemPerGrid(const Example_NSE2D& example, TCollection& coll, 
                    const TFESpace2D& velocity_space, 
                    const TFESpace2D& pressure_space, int order);
    };
    /** @brief a complete System on each grid 
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<SystemPerGrid> Systems;
    
    /** @brief right hand side vector from previous time step (on finest mesh)*/
      BlockVector old_rhs_modified;
  public:
    PrRobustTime_NSE2D(const TDomain& domain, Example_NSE2D& _example, 
                       int reference_id = -4711);
  
    /** @brief
     * Assemble the projection matrix only:
     * This will be used later for the post-processing 
     * of the nonlinear term:
     */
    bool assembleProjMat();
    
    /** @brief Assemble all the matrices and rhs before the time iterations
     * 
     * This includes the assembling of: Stiff_matrix, ProjectionMatrix, MassMatrix
     * for RT/BDM elements, and the Modifed_Mass matrix which comes from the 
     * multiplication of the projection matrix and MassMatrix of RT/BDM
     * results in 2 by 2 square blocks
     */
    void assemble_initial_time();
    // will be replaced latter
    void assemble_initial_timeNS();
    
    // assemble the projection matrix
    // FIXME: temporary function which is for the checking of previous working 
    // code: delete after the complete check
    void assembleProjectionMatrix();
    /** @brief 
     * 1. assembling the right hand side for special reconstruciton of the 
     * test function: Pressure robust method
     * 2. scaling of the B-blocks due to time stepping
     * this function will prepare the right hand side during the time 
     * discretization
     */
    void assemble_rhs();
    void assemble_rhsNS();
    /** @brief assemble nonlinear term
     * 
     * The matrix blocks to which the nonlinear term contributes are reset to 
     * zero and then completely reassembled, including the linear and nonlinear
     * terms. If this->assemble() has been called before, the matrix is now set
     * up correctly. 
     */
    void assemble_nonlinearNS();
    /** @brief Assemble the system matrix
     * This function will prepare the system which will be 
     * used for solvers
     */
    void assemble_system_matrix();
    
    /** @brief check if one of the stopping criteria is fulfilled
     * 
     * either converged, maximun number of iterations reached, or slow 
     * convergence
     * 
     * @param it_counter current iterate
     */
    bool stopIte(unsigned int it_counter);
    
    /** @brief solve the system*/
    void solve();
    /** @brief descale matrices
     * This function will reset all A-blocks due to the addition
     * of mass matrices and scaling factor
     * during the function call PrRobustTime_NSE2D::assemble_system_matrix()
     */
    void descaleMatrices();
    // post-processing 
    void output(int m, int& image);
    
    
    // getters and setters
    const TFESpace2D & get_projection_space() const 
    { 
      return this->Systems.front().projection_space; 
    }
    
};
#endif