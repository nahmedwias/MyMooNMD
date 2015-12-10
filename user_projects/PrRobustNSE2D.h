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
    BlockFEMatrix matrix;    
    /** */
    /** @brief constructor */
    SystemPerGrid(const Example_NSE2D& example, TCollection& coll, 
                    unsigned int order);
  };
  
  std::deque<SystemPerGrid> Systems;
  
public:
    PrRobustNSE2D(const TDomain& domain, const Example_NSE2D& _example, 
                  unsigned int reference_id = -4711);
};

#endif