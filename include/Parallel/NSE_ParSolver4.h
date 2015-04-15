// =======================================================================
// @(#)NSE_ParSolver4.h
//
// Class:      TNSE_ParSolver4
// Purpose:    Class containing all info needed for  parallel solvers
//
// Author:     Sashikumaar Ganesan (19.10.10)
//
// History:    Start of implementation 19.10.10 (Sashikumaar Ganesan)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"


#include <NSE_ParSolver.h>

#include <SquareStructure3D.h>
#include <SquareMatrix3D.h>
#include <Matrix3D.h>
#include <ParFECommunicator3D.h>


#ifndef __PARSSOLVERNSE3D4__
#define __PARSSOLVERNSE3D4__

class TNSE_ParSolver4 : public TNSE_ParSolver
{
  protected:


  public:
    /** constructor */
TNSE_ParSolver4(TParFECommunicator3D *velo_FEComm, TSquareStructure3D *squareStructure, TParFECommunicator3D *p_FECom, 
                TStructure3D * structureBT, TStructure3D *structureB);

};
#endif
#endif
