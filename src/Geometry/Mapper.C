// =======================================================================
// @(#)Mapper.C        1.3 11/15/99
// 
// Class:       TMapper
// Purpose:     mapper for geometric objects
//
// Author:      Volker Behns  30.07.97
//
// =======================================================================

#include <Mapper.h>

#ifndef NULL
#define NULL 0
#endif

// MapTriReg0: map tri, reg, (0,0)
static const int DatMapTriReg0RefVerts[] = {0,2,1,5,4,3};
static const int DatMapTriReg0RefEdges[] = {5,4,3,2,1,0,7,6,8};
static const int DatMapTriReg0RefFaces[] = {0,2,1,3};
static const int DatMapTriReg0OrigVerts[] = {0,2,1};
static const int DatMapTriReg0OrigEdges[] = {2,1,0};

// MapTriReg1: map tri, reg, (0,1)
static const int DatMapTriReg1RefVerts[] = {1,0,2,3,5,4};
static const int DatMapTriReg1RefEdges[] = {1,0,5,4,3,2,8,7,6};
static const int DatMapTriReg1RefFaces[] = {1,0,2,3};
static const int DatMapTriReg1OrigVerts[] = {1,0,2};
static const int DatMapTriReg1OrigEdges[] = {0,2,1};

// MapTriReg2: map tri, reg, (0,2)
static const int DatMapTriReg2RefVerts[] = {2,1,0,4,3,5};
static const int DatMapTriReg2RefEdges[] = {3,2,1,0,5,4,6,8,7};
static const int DatMapTriReg2RefFaces[] = {2,1,0,3};
static const int DatMapTriReg2OrigVerts[] = {2,1,0};
static const int DatMapTriReg2OrigEdges[] = {1,0,2};

// MapQuadReg0: map quad, reg, (0,0)
static const int DatMapQuadReg0RefVerts[] = {0,3,2,1,7,6,5,4,8};
static const int DatMapQuadReg0RefEdges[] = {7,6,5,4,3,2,1,0,11,10,9,8};
static const int DatMapQuadReg0RefFaces[] = {0,3,2,1};
static const int DatMapQuadReg0OrigVerts[] = {0,3,2,1};
static const int DatMapQuadReg0OrigEdges[] = {3,2,1,0};

// MapQuadReg1: map quad, reg, (0,1)
static const int DatMapQuadReg1RefVerts[] = {1,0,3,2,4,7,6,5,8};
static const int DatMapQuadReg1RefEdges[] = {1,0,7,6,5,4,3,2,8,11,10,9};
static const int DatMapQuadReg1RefFaces[] = {1,0,3,2};
static const int DatMapQuadReg1OrigVerts[] = {1,0,3,2};
static const int DatMapQuadReg1OrigEdges[] = {0,3,2,1};

// MapQuadReg2: map quad, reg, (0,2)
static const int DatMapQuadReg2RefVerts[] = {2,1,0,3,5,4,7,6,8};
static const int DatMapQuadReg2RefEdges[] = {3,2,1,0,7,6,5,4,9,8,11,10};
static const int DatMapQuadReg2RefFaces[] = {2,1,0,3};
static const int DatMapQuadReg2OrigVerts[] = {2,1,0,3};
static const int DatMapQuadReg2OrigEdges[] = {1,0,3,2};

// MapQuadReg3: map quad, reg, (0,3)
static const int DatMapQuadReg3RefVerts[] = {3,2,1,0,6,5,4,7,8};
static const int DatMapQuadReg3RefEdges[] = {5,4,3,2,1,0,7,6,10,9,8,11};
static const int DatMapQuadReg3RefFaces[] = {3,2,1,0};
static const int DatMapQuadReg3OrigVerts[] = {3,2,1,0};
static const int DatMapQuadReg3OrigEdges[] = {2,1,0,3};

//Constructor
TMapper::TMapper(Mapper which)
{
  switch (which)
  {
  //
  // maps for triangles
  //
    case MapTriReg0: // map tri, reg, (0,0)
         MapRefVerts = (const int *) DatMapTriReg0RefVerts;
         MapRefEdges = (const int *) DatMapTriReg0RefEdges;
         MapRefFaces = (const int *) DatMapTriReg0RefFaces;

         MapOrigVerts = (const int *) DatMapTriReg0OrigVerts;
         MapOrigEdges = (const int *) DatMapTriReg0OrigEdges;

         break;

    case MapTriReg1: // map tri, reg, (0,1)
         MapRefVerts = (const int *) DatMapTriReg1RefVerts;
         MapRefEdges = (const int *) DatMapTriReg1RefEdges;
         MapRefFaces = (const int *) DatMapTriReg1RefFaces;

         MapOrigVerts = (const int *) DatMapTriReg1OrigVerts;
         MapOrigEdges = (const int *) DatMapTriReg1OrigEdges;

         break;

    case MapTriReg2: // map tri, reg, (0,2)
         MapRefVerts = (const int *) DatMapTriReg2RefVerts;
         MapRefEdges = (const int *) DatMapTriReg2RefEdges;
         MapRefFaces = (const int *) DatMapTriReg2RefFaces;

         MapOrigVerts = (const int *) DatMapTriReg2OrigVerts;
         MapOrigEdges = (const int *) DatMapTriReg2OrigEdges;

         break;

    case MapQuadReg0: // map quad, reg, (0,0)
         MapRefVerts = (const int *) DatMapQuadReg0RefVerts;
         MapRefEdges = (const int *) DatMapQuadReg0RefEdges;
         MapRefFaces = (const int *) DatMapQuadReg0RefFaces;

         MapOrigVerts = (const int *) DatMapQuadReg0OrigVerts;
         MapOrigEdges = (const int *) DatMapQuadReg0OrigEdges;

         break;

    case MapQuadReg1: // map quad, reg, (0,1)
         MapRefVerts = (const int *) DatMapQuadReg1RefVerts;
         MapRefEdges = (const int *) DatMapQuadReg1RefEdges;
         MapRefFaces = (const int *) DatMapQuadReg1RefFaces;

         MapOrigVerts = (const int *) DatMapQuadReg1OrigVerts;
         MapOrigEdges = (const int *) DatMapQuadReg1OrigEdges;

         break;

    case MapQuadReg2: // map quad, reg, (0,2)
         MapRefVerts = (const int *) DatMapQuadReg2RefVerts;
         MapRefEdges = (const int *) DatMapQuadReg2RefEdges;
         MapRefFaces = (const int *) DatMapQuadReg2RefFaces;

         MapOrigVerts = (const int *) DatMapQuadReg2OrigVerts;
         MapOrigEdges = (const int *) DatMapQuadReg2OrigEdges;

         break;

    case MapQuadReg3: // map quad, reg, (0,3)
         MapRefVerts = (const int *) DatMapQuadReg3RefVerts;
         MapRefEdges = (const int *) DatMapQuadReg3RefEdges;
         MapRefFaces = (const int *) DatMapQuadReg3RefFaces;

         MapOrigVerts = (const int *) DatMapQuadReg3OrigVerts;
         MapOrigEdges = (const int *) DatMapQuadReg3OrigEdges;

         break;


    default:
         MapRefVerts = NULL;
         MapRefEdges = NULL;
         MapRefFaces = NULL;

         MapOrigVerts = NULL;
         MapOrigEdges = NULL;
  }
}
