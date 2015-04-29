/** ************************************************************************ 
*
* @class     TSystemMatScalar3D
* @brief     stores the information of a 3D scalar system matrix 
* @author    Sashikumaar Ganesan 
* @date      23.01.15
* @History    
 ************************************************************************  */


#ifndef __SYSTEMMATSCALAR3D__
#define __SYSTEMMATSCALAR3D__

#include <SquareMatrix3D.h>
#include <ItMethod.h>

#ifdef _MPI
//#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>

#include <ParDirectSolver.h>
#endif

/** class for 3D scalar system matrix */
class TSystemMatScalar3D
{
  protected:
    /** own fespace and parallel FE Communicator */
    #ifdef _MPI
    TParFEMapper3D **ParMapper;
    TParFECommunicator3D **ParComm;
    MPI_Comm Comm;
    TParDirectSolver *DS;
    #endif
    
    /** Number of multigrid levels */
    int N_Levels;
    
    /** starting level for the solver, e.g. Start_Level = N_Levels-1 for direct solver */
    int Start_Level;    
    
    /** fespace */
    TFESpace3D **FeSpaces, *fesp[1], *ferhs[1];
    
    /** Boundary condition and Boundary Value */
    BoundCondFunct3D *BoundaryConditions[1]; 
    
    /** Boundary condition and Boundary Value */
    BoundValueFunct3D *BoundaryValues[1];    
    
    /** N_DOF at fine level */
    int N_DOF;
    
    /** Discretization type */
    int Disctype;
    
    /** SOLVER type */
    int SOLVER;
       
    /** number of matrices in the system matrix*/
    int N_Matrices;

    /** sqstructure of the system matrix */
    TSquareStructure3D **sqstructure;

    /** A is the stiffness/system mat for stationary problem   */
    TSquareMatrix3D **sqmatrixA, *SQMATRICES[1];
    TSquareMatrix **sqmatrices;
        
    /** Discrete form for the equation */
    TDiscreteForm3D *DiscreteFormARhs;
    
    /** rhs for assemble */
    double *RHSs[1];
    
   /** variables for multigrid */
    double Parameters[2], N_aux, *Itmethod_sol, *Itmethod_rhs;
    TMultiGrid3D *MG;
    TMGLevel3D *MGLevel;
    TItMethod *Itmethod, *prec;
   
  public:
    /** Constructor*/
  TSystemMatScalar3D(int N_levels, TFESpace3D **fespaces, int disctype, int solver);
    
    /** destrcutor */
    ~TSystemMatScalar3D();
    
    /** Initilize the discrete forms and the matrices */
    void Init(CoeffFct3D *BilinearCoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *BoundValue);
 
    /** assemble the system matrix */
    void Assemble(TAuxParam3D *aux, double **sol, double **rhs);

    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
    
#ifdef _MPI
    TParFECommunicator3D *Get_ParComm(int level){
      return ParComm[level];
    }
#endif
};

#endif
