// =======================================================================
// @(#)NSE_ParSolver.h
//
// Class:      TNSE_ParSolver 
// Purpose:    Super class for all NSE_ParSolver* 
//
// Author:     Sashikumaar Ganesan (19.10.10)
//
// History:    Start of implementation 19.10.10 (Sashikumaar Ganesan)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"

#ifdef __3D__
#include <SquareStructure3D.h>
#include <SquareMatrix3D.h>
#include <Matrix3D.h>
#include <ParFECommunicator3D.h>

#endif

#include <ParSolver.h>

#ifndef __NSEPARSOLVER__
#define __NSEPARSOLVER__

class TNSE_ParSolver : public TParSolver
{
  protected:

   /** global no. dof of velo */
   int Global_N_U;

   /** global no. dof of pressure */
   int Global_N_P;

#ifdef __3D__
   /** matrix strcuture */
   TSquareStructure3D *SquareStructure;

   /** matrix strcuture */
   TStructure3D *Structure;

   /** matrix strcuture */
   TStructure3D *StructureT;

   /** ParCommunicator of the corresponding Fespace*/
   TParFECommunicator3D *Velo_FEComm, *P_FEComm;

#endif

   /** No. Dist dof of velocity*/
   int N_distU;

   /** No. Dist dof of velocity*/
   int N_distActiveU;

   /** Starting of Global number of the Dist dof of velocity*/
   int Global_beginU;

   /** Rowptr for the dist rectangle A matrix (N_DistDof x N) block of (N x N) matrix  
    i.e., N = sum(N_DistDof) in all process */
   int *DistA_RowPtr;

   /** Global column index of dist mat A */
   int *DistA_KCol;

   /** No. entries for A mat(ONE component) in the rectangle system matrix  */
   int DistA_N_Entries;

   /** Rowptr for the dist rectangle BT matrix (N_DistDof x N) block of (N x N) matrix  
    i.e., N = sum(N_DistDof) in all process */
   int *DistBT_RowPtr;

   /** Global column index of dist mat BT */
   int *DistBT_KCol;

   /** No. entries for BT mat in the rectangle system matrix  */
   int DistBT_N_Entries;

   /** No. Dist dof of pressure*/
   int N_distP;

   /** Starting of Global number of the Dist dof of pressure*/
   int Global_beginP;

// //    int *distu_coldisp, *distp_coldisp;
// //    int *distu_range, *distp_range;

   /** Rowptr for the dist rectangle BT matrix (N_DistDof x N) block of (N x N) matrix  
    i.e., N = sum(N_DistDof) in all process */
   int *DistB_RowPtr;

   /** Global column index of dist mat BT */
   int *DistB_KCol;

   /** No. entries for B mat in the rectangle system matrix  */
   int DistB_N_Entries;

   /**  DistDof of Loc Dof velocity */
   int *DistDofofLocDof_U;

   /**  DistDof of Loc Dof pressure */
   int *DistDofofLocDof_P;

   /** OwnDof for Pressure Constant */
   int Pressure_Const_DistDof;

  /**  Aii position in the system matrix */
  int *Aii_pos, *Aii_Deptpos;

  int MaxN_NzInDepRow_All_A, MaxN_NzInDepRow_All_BT, MaxN_NzInDepRow_All_B;

  void GetDistArray(double *RHS, int nu, int np, double *tmp_rhs);
   
  void GetLocArray(double *tmp_rhs, int nu, int np, double *RHS);  
     
  
  public:
    /** constructor */
    TNSE_ParSolver(TParFECommunicator3D *velo_FEComm, TSquareStructure3D *squareStructure,
        TParFECommunicator3D *p_FECom, TStructure3D * structureBT, TStructure3D *structureB,
        int type);

    /**methods */
    void SetUpDistMat();

    /** destructor */
    ~TNSE_ParSolver();
};
#endif
#endif
