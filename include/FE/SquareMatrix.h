#ifndef __SQUAREMATRIX__
#define __SQUAREMATRIX__

#include <FEMatrix.h>

/** @brief deprecated, use the FEMatrix class instead */
class TSquareMatrix : public FEMatrix
{
  protected:
    /** generate the matrix, called from derived classes */
    TSquareMatrix(const TFESpace1D * space);
    TSquareMatrix(const TFESpace2D * space);
    #ifdef __3D__
    TSquareMatrix(const TFESpace3D * space);
    #endif // 3D
  public:

    /** destructor: free Entries array */
    ~TSquareMatrix() = default;
    
    TSquareMatrix(const TSquareMatrix & m) = default;
};

#endif // __SQUAREMATRIX__
