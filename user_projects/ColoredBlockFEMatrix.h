/** ************************************************************************
*
* @class      ColoredBlockFEMatrix
* @brief      extends ColoredBlockMatrix by handling active degrees of freedom
*
*             A BlockMatrix of subtype BlockFEMatrix stores FEMatrices
*             instead of simple algebraic TMatrices. Thus it can access
*             information on finite element spaces and active degrees of
*             freedoms and exploits this algoritmically. The majority of
*             BlockMatrices in ParMooN are in fact BlockFEMatrices.
*
*             \todo This class is still in its infancy.
*
*
* @author     Naveed Ahmed, Clemens Bartsch, Ulrich Wilbrandt
* @date       2015/12/08
*
*
*
****************************************************************************/

#ifndef USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_
#define USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_

#include <ColoredBlockMatrix.h>

class ColoredBlockFEMatrix : public ColoredBlockMatrix
{

};



#endif /* USER_PROJECTS_COLOREDBLOCKFEMATRIX_H_ */
