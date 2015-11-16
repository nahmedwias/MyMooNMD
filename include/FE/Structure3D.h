// =======================================================================
// @(#)Structure3D.h        1.3 09/14/99
// 
// Class:       TStructure3D
//
// Purpose:     build and store a matrix Structure3D
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation
//
//              start of reimplementation 26.08.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __STRUCTURE3D__
#define __STRUCTURE3D__

#include <FESpace2D.h>
#include <FESpace3D.h>
#include <Structure.h>

class TStructure3D : public TStructure
{
  protected:

  public:
    /** generate the matrix Structure3D, both space with 3D collection */
    TStructure3D(const TFESpace3D *testspace, const TFESpace3D *ansatzspace)
  : TStructure(testspace, ansatzspace){}

};

#endif
