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
      
      BlockFEMatrix Modifed_Mass;
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
  
    
    void assemble_initial_time();
    
    // assemble the projection matrix
    void assembleProjectionMatrix();
    
    void assemble_rhs();
    
    void AssembleRHS_Only();
    
    void assemble_system_matrix();
    
    void solve();
    
    void output(int m);
    
    
    // getters and setters
    const TFESpace2D & get_projection_space() const 
    { 
      return this->Systems.front().projection_space; 
    }
    
};
#endif