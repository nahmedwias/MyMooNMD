// =======================================================================
// @(#)RefTriBaryDesc.h        1.1 10/30/98
//
// Class:       TRefTriBaryDesc
// Purpose:     refinement descriptor for regular refinement of a triangle
//
// Author:      Volker Behns  17.07.97
//
// =======================================================================

#ifndef __REFTRIBARYDESC__
#define __REFTRIBARYDESC__

#include <RefDesc.h>

/** refinement descriptor for regular refinement of a triangle */
class TRefTriBaryDesc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a regular refinement of a triangle */
    TRefTriBaryDesc(TShapeDesc *shape);
};

#endif
