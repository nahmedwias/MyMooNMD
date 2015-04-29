// =======================================================================
// @(#)ParDirectSolver.h
//
// Class:      TParDirectSolver 
// Purpose:    Class for interfacing ParMooN with MumpsSolver
//
// Author:     Sashikumaar Ganesan & Abdus Shamim (27.04.15)
//
// History:    Start of implementation 27.04.15 (Sashikumaar Ganesan & Abdus Shamim)
//
// =======================================================================
#ifdef _MPI
 #  include "mpi.h"
 
#include <SquareStructure3D.h>
#include <SquareMatrix3D.h>
#include <ParFECommunicator3D.h>

#include <MumpsSolver.h>

#ifndef __DIRECTSOLVER_MUMPS__
#define __DIRECTSOLVER_MUMPS__

class TParDirectSolver
{
protected:
  TParFECommunicator3D *ParComm;
  MPI_Comm Comm;
  
  TMumpsSolver *Mumps;
  
  //Mumps Solver Parameters
  int N_Master, NDof;
  
  //row ptr and col ptr of system matrix
  int *RowPtr, *KCol;
  
  // number of non zeroes in system matrix
  int N_Nz;
  
  // number of global degrees of freedom over all ranks
  int N_Eqns;
  
  // number of rhs, OwnRhs(only master dofs)
  int N_rhs;
  double *OwnRhs;
  
  double *GlobalRhs;
  int GlobalRhsSize;
  
  //(I,J) for MUMPS input
  int *I_rn, *J_cn;
  
  //SqMatrix
  TSquareMatrix3D *Mat;
  
  //MatLoc is the system matrix entry values(only master rows)
  double *MatLoc;
  
public:
  TParDirectSolver(TParFECommunicator3D *parcomm,TSquareMatrix3D *mat);
  
  ~TParDirectSolver();
  
  void AssembleLocMatrix();
  
  void InitMumps();
  
  void GetRhs(double *Rhs);
  
  void UpdateSol(double *Sol);
  
  void Solve(double *Sol, double *Rhs, bool Factorize);
  
};
 
#endif
#endif
