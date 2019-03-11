#ifndef BLOCKFEMATRIX_PR_H
#define BLOCKFEMATRIX_PR_H

#include <BlockFEMatrix.h>
/**
 * @todo write docs
 */
class BlockFEMatrix_PR : public BlockFEMatrix
{
public:
    /**
     * @todo write docs
     */
    BlockFEMatrix_PR(const TFESpace2D& velocity);
    
    static BlockFEMatrix ProJMat(const TFESpace2D& velocity); 

};

#endif // BLOCKFEMATRIX_PR_H
