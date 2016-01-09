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
      /** @brief matrix used for the reconstruciton*/
      BlockFEMatrix projMat;    
      
      BlockVector rhsXh;
      /** */
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
    
  public:
    // 
    PrRobustNSE2D(const TDomain& domain, const Example_NSE2D& _example, 
                  unsigned int reference_id = -4711);
    
    /** @brief update the right hand side
     * this function will multiply the right hand side
     * obtained by usual assembling with the matrix 
     * obtained from the RT0 projection
     */
    void assembleMatrixRhs();
    
    /** @brief 
     * solve the resulting system
     */
    void solve();
    
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
    
  private:
    /** @brief assemble the right hand side for vector space
     * This assemble's only the right hand side 
     * which uses the vector-space
     */
    void assemblerhs();
};

#endif