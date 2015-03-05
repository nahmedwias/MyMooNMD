// =======================================================================
// @(#)Scalar_ParSolver.h
//
// Class:      TScalar_ParSolver 
// Purpose:    Super class for all Scalar_ParSolver
//
// Author:     Sashikumaar Ganesan (22.12.10)
//
// History:    Start of implementation 22.12.10 (Sashikumaar Ganesan)
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"

#ifdef __3D__
#include <SquareStructure3D.h>
#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>

#endif

#include <ParSolver.h>

#ifndef __SCALARPARSOLVER__
#define __SCALARPARSOLVER__

class TScalar_ParSolver : public TParSolver
{
  protected:

   /** global no. dof of velo */
   int Global_DOF;

   /** No. rhs */
   int N_Rhs;
   
#ifdef __3D__
   /** matrix strcuture */
   TSquareStructure3D *SquareStructure;

   /** ParCommunicator of the corresponding Fespace*/
   TParFECommunicator3D *FEComm;

#endif

   /** No. Dist dof */
   int N_distDof;

   /** No. Dist dof */
   int N_distActiveDof;

   /** Starting of Global number of the Dist dof*/
   int Global_beginDof;

   /** Rowptr for the dist rectangle A matrix (N_DistDof x N) block of (N x N) matrix  
    i.e., N = sum(N_DistDof) in all process */
   int *DistA_RowPtr;

   /** Global column index of dist mat A */
   int *DistA_KCol;

   /** No. entries for A mat(ONE component) in the rectangle system matrix  */
   int DistA_N_Entries;

   /**  DistDof of Loc Dof velocity */
   int *DistDofofLocDof;

  int MaxN_NzInDepRow_All_A;

//   void GetDistArray(double *RHS, int nu, int np, double *tmp_rhs);
   
//   void GetLocArray(double *tmp_rhs, int nu, int np, double *RHS);  

  private:

   void AssembleDistMatrix(TSquareMatrix3D *Mat);
   
   void SetUpIterDistMat();
   
   void InitMumps();

   void AssembleIterDistMatrix(TSquareMatrix3D *Mat);
   
  public:
    /** constructor */
    TScalar_ParSolver(TParFECommunicator3D *fEComm, TSquareStructure3D *squareStructure,
                      int type, TParVector3D *parSol, TParVector3D *parRhs);

    /**methods */
    void SetUpDistMat();
    
    void Solve(TSquareMatrix3D *Mat, bool FACTORIZE);    
    
    
    
//     int GetIndex(int *Array, int Length, int pos);

    /** destructor */
    ~TScalar_ParSolver();
};
#endif
#endif
