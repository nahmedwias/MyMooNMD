#include <SquareMatrix.h>

TSquareMatrix::TSquareMatrix(const TFESpace1D * space)
 : FEMatrix(space)
{
  
}
TSquareMatrix::TSquareMatrix(const TFESpace2D * space)
 : FEMatrix(space)
{
  
}
#ifdef __3D__
TSquareMatrix::TSquareMatrix(const TFESpace3D * space)
 : FEMatrix(space)
{
  
}
#endif // 3D

