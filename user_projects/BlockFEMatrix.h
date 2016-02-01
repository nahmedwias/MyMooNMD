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
     * result a matrix used for the pressure robust methods
     * @param  testSpace: the 
     * @param  ansatzSpace: the space to project on
     */
    BlockFEMatrix(const TFESpace2D* testSpace, const TFESpace2D* ansatzSpace, 
                  const TFESpace2D* pressureSpace);

    /** @brief compute y = y + a * Mx 
     * M is this matrix
     * 
     * add the matrix-vector product "M*x", scaled by "a", to y: only active
     * 
     * 
     * @param x the vector multiplied by M
     * @param y targeted vector: result will be stored in this vector
     * @param factor scaling factor (optional) default 1.
     */
    double addScaleActive(const double *x, double *y, double factor = 1.0);
    
    /** @brief*/
    const TFESpace2D * get_space_of_block(unsigned int b, bool test) const
    {
      return test ? TestSpace : AnsatzSpace;
    }
  private:
    const TFESpace2D *TestSpace;
    const TFESpace2D *AnsatzSpace;
    const TFESpace2D *PressureSpace;
};

#endif