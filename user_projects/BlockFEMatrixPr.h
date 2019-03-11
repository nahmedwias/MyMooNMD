#ifndef BLOCKFEMATRIXPR_H
#define BLOCKFEMATRIXPR_H

#include "../include/Matrix/BlockFEMatrix.h"

/**
 * @todo write docs
 */
class BlockFEMatrixPr : public BlockFEMatrix
{
public:
    /** Default constructor */
#ifdef __2D__
   explicit BlockFEMatrixPr(std::vector< const TFESpace2D* > spaces_rows, 
                    std::vector< const TFESpace2D* > spaces_cols);
#endif

  static BlockFEMatrixPr Projection_NSE2D(const TFESpace2D& velocity, 
                   const TFESpace2D& projection, const TFESpace2D& pressure);
};

#endif // BLOCKFEMATRIXPR_H
