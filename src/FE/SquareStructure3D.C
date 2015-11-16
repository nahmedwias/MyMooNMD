// =======================================================================
// @(#)SquareStructure3D.C        1.6 09/17/99
// 
// Class:       TStructure
//
// Purpose:     build and store a structure for a square matrix in 3d
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
//              02.09.1999 adding connections over mortar edges
//                         use __ADD_LINK__ (Volker Behns)
//
// =======================================================================

#include <DefineParams.h>
#include <FEDatabase3D.h>
#include <SquareStructure3D.h>
#include <HNDesc.h>
#include <HangingNode.h>
#include <string.h>
#include <Database.h>


