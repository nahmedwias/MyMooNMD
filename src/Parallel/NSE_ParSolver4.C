// =======================================================================
// @(#)NSE_ParSolver4.C
//
// Class:      TNSE_ParSolver4
// Purpose:    Class containing all info needed for communication between subdomains
//             for square matrices
//
// Author:     Sashikumaar Ganesan (19.10.10)
//
// History:    Start of implementation 19.10.10 (Sashikumaar Ganesan)
//
// =======================================================================

#ifdef _MPI
#  include "mpi.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <NSE_ParSolver4.h>
#include <Database.h>
#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>
#include <ParVectorNSE3D.h>


    /** constructor */
TNSE_ParSolver4::TNSE_ParSolver4(TParFECommunicator3D *velo_FEComm, 
        TSquareStructure3D *squareStructure, TParFECommunicator3D *p_FECom, 
        TStructure3D * structureBT, TStructure3D *structureB)
            : TNSE_ParSolver(velo_FEComm, squareStructure, p_FECom,
                            structureBT, structureB, 4)
 {

 }

#endif