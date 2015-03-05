// =======================================================================
// @(#)ParSolver.h
//
// Class:      ParSolver 
// Purpose:    Super class for all Parallel Solver classes 
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

#include <ParFECommunicator3D.h>
#endif

#include <MumpsSolver.h>
#include <ParVector3D.h>
// #include <ParmsSolver.h>
// #include <HipsSolver.h>

#ifndef __PARSOLVER__
#define __PARSOLVER__

class TParSolver
{
  protected:

   /** MPI_Comm for which the fespace communications are needed */
   MPI_Comm Comm;

   /** total number No.of equations in the system matrix ( same in all processors )*/
   int N_Eqns;

    /** rhs */
    double *RHS;

     /** rhs */
    double *TCLS_SOL;   
    
   /** total number of equations in the distributed system matrix*/
   int N_DistMatEqns;

   /** total number of non-zero entries distributed system matrix*/
   int N_DistMatEntries;

   /** row index of  distributed system matrix MooNMD based*/
   int *DistMat_RowPtr;
   
   /** row index of  distributed system matrix solvers based */
   int *DistMat_Irn;

   /** column index of distributed system matrix solvers based  */
   int  *DistMat_Jcn;

   /** distributed system matrix values of length N_DistMatEntries */
    double *DistMat_Entries;

   /** parallel solution vector */    
    TParVector3D *ParSol;
 
   /** parallel Rhs vector */    
    TParVector3D *ParRhs;    
    
//    /** No. of Master Dofs which is needed in orther processors in iterative solver in each neib*/
//     int *TCLS_N_MastDofToNeib;
//     
//    /** No. of Master Dofs which is needed in orther processors in iterative solver in each neib*/
//     int TCLS_N_MastDofToNeib_All;
//     
//    /** list of Master Dofs which is needed in orther processors in iterative solver in each neib*/
//     int *TCLS_MastDofToNeib;

   /** external solver MUMPS */
  TMumpsSolver *MUMPS_Solver;
  
   /** external solver SuperLU_Solver */
//   TSuperLUSolver *SuperLU_Solver;
  
   /** external solver pARMS */
//   TParmsSolver *Parms_Solver;

   /** external solver Hips */
//   THipsSolver *Hips_Solver;


#ifdef __3D__
   /** matrix strcuture */
   TSquareStructure3D *SquareStructure;

#endif

  public:
    
   void  Sort(int *Array, int length);

   int GetIndex(int *Array, int Length, int pos);

};
#endif
#endif
