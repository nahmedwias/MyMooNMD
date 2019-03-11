#include "BlockFEMatrix_PR.h"

BlockFEMatrix_PR::BlockFEMatrix_PR(const TFESpace2D& velocity)
:BlockFEMatrix({&velocity, &velocity})
{
  
}

BlockFEMatrix_PR BlockFEMatrix::ProJMat(const TFESpace2D& velocity)
{
  BlockFEMatrix my_matrix({&velocity, &velocity});
  
  my_matrix.replace_blocks(FEMatrix(&velocity, &velocity), {{0,0}, {1, 1}}, 
                           {false, false});

  return my_matrix;
}
