#ifndef __PRROBUSTNSE2D__
#define __PRROBUSTNSE2D__

#include <NSE2D.h>
#include <BlockFEMatrix.h>

class PrRobustNSE2D : NSE2D
{
  protected:  
    struct SystemPerGrid
    {
      /** @brief Finite Element space for the projection*/
      TFESpace2D projection_space;      
      
      /** @brief projection matrix for reconstruciton*/      
      BlockFEMatrix ProjectionMatrix;
      
      /** @brief right hand side for reconstruciton*/
      BlockVector rhsXh;
      
      /** @brief constructor */
      SystemPerGrid(const Example_NSE2D& example, TCollection& coll, 
                    unsigned int order, const TFESpace2D& testSpace, 
                    const TFESpace2D& presSpace);
    };
    /** @brief a complete System on each grid 
     * 
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<SystemPerGrid> Systems;
    
    /**
       * @brief store the square root of the residual from previous iteration
       */
      double oldResidual;
      /** @brief store the initial residual so that the nonlinear iteration can 
       *         be stopped as soon as a desired reduction is achieved
       */
      double initial_residual;
    
  public:
    // 
    PrRobustNSE2D(const TDomain& domain, const Example_NSE2D& _example, 
                  unsigned int reference_id = -4711);
    
    /** @brief assemvle 
     * assemble the system matrix and right hand side of the NSE2D.
     * In addition: 
     * Inside the definition of this function if the reconstruciton
     * is used then additionally the right hand side is assembled 
     * using vector value spaces. In the next, the projection matrix 
     * is assembled via the function "assembleProjectionMatrix". 
     * 
     * The system right hand side is then obtained my multiplyig 
     * the "ProjectionMatrix' with right hand side vector "rhsXh".
     * 
     * obtained by usual assembling with the matrix 
     * obtained from the RT0 projection
     */
    void assembleMatrixRhs();
    
    // assemble the projection matrix
    void assembleProjectionMatrix();
    
    /** @brief stopping creterion and residuals*/
    bool stopIteration(unsigned int it);
    /** @brief assembling of nonlinear matrices
     * no right hand side and pressure blocks
     */
    void assembleNL();
    
    /** @brief 
     * solve the resulting system
     */
    void solve();
    /** @brief 
     * post processing
     */
    void output(int i);
    
    // getters and setters
    const TFESpace2D & get_projection_space() const 
    { 
      return this->Systems.front().projection_space; 
    }
    
    const TFESpace2D & get_velocity_space() const
    { return this->NSE2D::systems.front().velocity_space; }
    
    const TFESpace2D & get_pressure_space() const
    { return this->NSE2D::systems.front().pressure_space; }
    
    const TFEVectFunct2D & get_velocity() const
    { return this->NSE2D::systems.front().u; } 
    
    TFEFunction2D *get_velocity_component(int i)
    { return (i==0) ? this->NSE2D::systems.front().u.GetComponent(0)
                    : this->NSE2D::systems.front().u.GetComponent(1); }
                    
    TFEFunction2D &get_pressure()
    {return this->NSE2D::systems.front().p; }    
};

#endif