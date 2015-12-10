/** ************************************************************************
 * 
 * @class BlockFEMatrix
 * @brief store the information 
 * @author Naveed, Ulrich, Clemens 
 * @date   04.12.2015  
 * @History Added methods
 * 
 * 
* *************************************************************************/

#ifndef __BLOCKFEMATRIX__
#define __BLOCKFEMATRIX__


#include <FEMatrix.h>
#include <BlockMatrix.h>


class BlockFEMatrix : public BlockMatrix
{
  public:
    /** @brief a constructor
     * 
     * @param  testSpace: the 
     * @param  ansatzSpace: the space to project on
     */
    BlockFEMatrix(const TFESpace2D* testSpace, const TFESpace2D* ansatzSpace);

  private:
};

#endif