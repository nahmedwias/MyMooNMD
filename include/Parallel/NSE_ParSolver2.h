// =======================================================================
// @(#)NSE_ParSolver2.h
//
// Class:      TNSE_ParSolver2
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


#ifndef __PARSSOLVERNSE3D2__
#define __PARSSOLVERNSE3D2__

class TNSE_ParSolver2 : public TNSE_ParSolver
{
  protected:

  private:

   void InitMumps();

   //void InitMumps_2();

   void InitSuperLU();
    
   void InitHips();

   void InitpARMS();

   void AssembleDistMatrix_A(TSquareMatrix3D *MatA);

   
//    void AssembleMatrix(TSquareMatrix3D *MatA, TMatrix3D *MatB1T, 
//                        TMatrix3D *MatB2T, TMatrix3D *MatB3T, TMatrix3D  *MatB1, 
//                        TMatrix3D *MatB2, TMatrix3D * MatB3, double *Rhs);
   
   void AssembleDistMatrix(TSquareMatrix3D *MatA, TMatrix3D *MatB1T, TMatrix3D *MatB2T,
               TMatrix3D *MatB3T, TMatrix3D  *MatB1, TMatrix3D *MatB2, TMatrix3D * MatB3,
               double *Rhs);


  public:
    /** constructor */
   TNSE_ParSolver2(TParFECommunicator3D *velo_FEComm, TSquareStructure3D *squareStructure,
                   TParFECommunicator3D *p_FECom, TStructure3D * structureBT,
                   TStructure3D *structureB);


    void Solve(TSquareMatrix3D *MatA, TMatrix3D *MatB1T, TMatrix3D *MatB2T,
               TMatrix3D *MatB3T, TMatrix3D  *MatB1, TMatrix3D *MatB2, TMatrix3D * MatB3,
               TParVectorNSE3D  *rhs, TParVectorNSE3D  *sol);

    void Solve(TSquareMatrix3D *MatA, TParVectorNSE3D  *rhs, TParVectorNSE3D  *sol); 
      

};
#endif
#endif

