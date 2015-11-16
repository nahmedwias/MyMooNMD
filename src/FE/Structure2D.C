// =======================================================================
// @(#)Structure2D.C        1.10 11/24/99
// 
// Class:       TStructure
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
#ifdef _MPI
#  include "mpi.h"
#endif

#include <DefineParams.h>

#include <FEDatabase2D.h>
#include <MooNMD_Io.h>
#include <Structure.h>
#include <string.h>
#include <stdlib.h>

#ifdef __MORTAR__
#ifdef __ADD_LINK__
  #include <Database.h>
  #include <It_Mortar.h>
  #include <MortarBaseJoint.h>
  #include <stdlib.h>
#endif // __ADD_LINK__
#endif // __MORTAR__
#include <Database.h>
