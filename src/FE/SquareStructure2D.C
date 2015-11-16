// =======================================================================
// @(#)SquareStructure2D.C        1.6 09/17/99
//
// Class:       TStructure
//
// Purpose:     build and store a structure for a square matrix in 2d
//
// Author:      Gunar Matthies
//
// History:     06.08.1998 start implementation
//
//              02.09.1999 adding connections over mortar edges
//                         use __ADD_LINK__ (Volker Behns)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <DefineParams.h>
#include <FEDatabase2D.h>
#include <SquareStructure2D.h>
#include <Database.h>
#include <string.h>

#include <Joint.h>

#ifdef __MORTAR__
#ifdef __ADD_LINK__
#include <Database.h>
#include <It_Mortar.h>
#include <MortarJoint.h>
#include <stdlib.h>
#endif                                            // __ADD_LINK__
#endif                                            // __MORTAR__

#include <Database.h>

