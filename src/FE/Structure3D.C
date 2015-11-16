// =======================================================================
// @(#)Structure3D.C        1.10 11/24/99
// 
// Class:       TStructure3D
//
// Purpose:     build and store a matrix structure
//
// Author:      Gunar Matthies
//
// History:     24.11.97 start implementation (Gunar Matthies)
//
//              04.08.1998 start reimplementation (Gunar Matthies)
//              02.09.1999 adding connections over mortar edges
//                         use __ADD_LINK__ (Volker Behns)
//
// =======================================================================

#include <FEDatabase3D.h>
#include <Structure3D.h>
#include <string.h>
#include <Database.h>

