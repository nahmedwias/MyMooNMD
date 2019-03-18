#include <Matrix2D.h>

TMatrix2D::TMatrix2D(const TFESpace2D * testspace,
                     const TFESpace2D * ansatzspace)
 : FEMatrix(testspace, ansatzspace)
{
  
}
