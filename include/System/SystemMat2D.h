/** ************************************************************************ 
*
* @class     SystemMat2D
* @brief     stores the information of a 2D system matrix
* 
* This is a base class with very little functionality. Derived classes for 
* differnet problem types exist and should be used rather than using this class
* directly.
* 
* @author    Ulrich Wilbrandt,
* @date      17.03.15
 ************************************************************************  */


#ifndef __SYSTEMMAT2D__
#define __SYSTEMMAT2D__

#include <SquareMatrix2D.h>
#include <vector>

/**class for 2D  NSE system matrix */
class SystemMat2D
{
  protected:

    /** @brief vector of finite element fespaces 
     * 
     * Each space appears exactly once in this vector. Typically this vector has
     * size 1 (scalar) or two (Saddle-point problems).
     */
    std::vector<TFESpace2D*> fe_spaces;
    
    /** @brief vector of the square matrices */
    std::vector<TSquareMatrix2D*> sq_matrices;
  
    /** @brief vector of the rectangular matrices */
    std::vector<TMatrix2D*> rect_matrices;
    
    /** @brief method for resudual calculation */
    DefectProc *defect;
    
  public:
    /** constructor */
     SystemMat2D(unsigned int n_spaces, unsigned int n_sq_matrices, 
                 unsigned int n_rect_matrices)
      : fe_spaces(n_spaces, NULL), sq_matrices(n_sq_matrices, NULL),
        rect_matrices(n_rect_matrices, NULL), defect(NULL)
     {};
};

#endif // __SYSTEMMAT2D__