#include <Matrix3D.h>

TMatrix3D::TMatrix3D(const TFESpace3D * testspace,
                     const TFESpace3D * ansatzspace)
 : FEMatrix(testspace, ansatzspace)
{
  
}
