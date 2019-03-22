#ifndef __REFTETRABARYDESC__
#define __REFTETRABARYDESC__

#include <RefDesc.h>

/** refinement descriptor for barycentric refinement of a tetrahedron */
class TRefTetraBaryDesc : public TRefDesc
{
  public:
    explicit TRefTetraBaryDesc(TShapeDesc *shape);
};

#endif // __REFTETRABARYDESC__
