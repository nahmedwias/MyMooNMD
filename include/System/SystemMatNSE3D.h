/** ************************************************************************ 
*
* @class     TSystemMatNSE3D
* @brief     stores the information of a 3D NSE system matrix 
* @author    Sashikumaar Ganesan, 
* @date      27.01.15
* @History    
 ************************************************************************  */


#ifndef __SYSTEMMATNSE3D__
#define __SYSTEMMATNSE3D__

#include <SquareMatrix3D.h>
#include <FEVectFunct3D.h>
#include <NSE_MultiGrid.h>
#include <ItMethod.h>

/**class for 3D  NSE system matrix */
class TSystemMatNSE3D
{
  protected:

    /** Number of multigrid levels */
    int N_Levels;
    
    /** starting level for the solver, e.g. Start_Level = N_Levels-1 for direct solver */
    int Start_Level;  
    
    /** DOFs of velocity and pressure spaces */    
    int N_TotalDOF, N_U, N_P, N_Active, N_DirichletDof;
    
    /** velocity and pressure fespace */
    TFESpace3D **U_Space, **P_Space, *FeSpaces[4];
    
    /** velo FE function */
    TFEVectFunct3D **Velocity;
    
    /** Fe functions of NSE */
    TFEFunction3D **Pressure;    
    
    /** Discretization type */
    int Disctype;
       
    /** NSE type */
    int NSEType; 
           
    /** Bilinear coefficient   */
    CoeffFct3D *LinCoeffs[1];    
    
    /** NSEaux is used to pass additional fe functions (eg. mesh velocity) that is nedded for assembling */
    TAuxParam3D *NSEaux, *NSEaux_error;
    
    /** method for matrix vector mult */
    MatVecProc *MatVect;  
    
    /** method for resudual calculation */
    DefectProc *Defect;   
    
    /** Solver type */
    int SOLVER;
       
    /** number of matrices in the system matrix*/
    int N_Matrices;

    /** sqstructureA of the system matrix */
    TSquareStructure3D **sqstructureA;

    /** structure of the system matrix */
    TStructure3D **structureB, **structureBT;
    
    /** A is the stiffness/system mat for NSE velocity component   */
    TSquareMatrix3D **SqmatrixA11, **SqmatrixA12, **SqmatrixA13;
    TSquareMatrix3D **SqmatrixA21, **SqmatrixA22, **SqmatrixA23;
    TSquareMatrix3D **SqmatrixA31, **SqmatrixA32, **SqmatrixA33, *SQMATRICES[9];
    TSquareMatrix **sqmatrices;
  
    /** B is the  system mat for NSE pressure component   */
    TMatrix3D **MatrixB1, **MatrixB2, **MatrixB3, **MatrixB1T, **MatrixB2T, **MatrixB3T, *MATRICES[6];
    TMatrix **matrices;
    
    /** Boundary conditon */
    BoundCondFunct3D *BoundaryConditions[3];
  
     /** Boundary values */   
    BoundValueFunct3D *BoundaryValues[3];
        
    /** Discrete form for the equation */
    TDiscreteForm3D *DiscreteFormARhs, *DiscreteFormNL;

    /** variables for multigrid */
    int N_aux;
    double Parameters[2], *Itmethod_sol, *Itmethod_rhs;
    TNSE_MultiGrid *MG;
    TNSE_MGLevel *MGLevel;
    TItMethod *Itmethod, *prec;
    
  public:
    /** constructor */
     TSystemMatNSE3D(int N_levels, TFESpace3D **velocity_fespace, TFESpace3D **presssure_fespace, TFEVectFunct3D **velocity, 
                     TFEFunction3D **pressure, int disctype, int nsetype, int solver);

    /** destrcutor */
//     ~TSystemMatNSE3D();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct3D *lincoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *U1BoundValue, 
              BoundValueFunct3D *U2BoundValue, BoundValueFunct3D *U3BoundValue);
    
    /** assemble the system matrix */
    void Assemble(double **sol, double **rhs);
    
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear(double **sol, double **rhs);
    
    /** get the resudual of the NSE system */
    void GetResidual(double *sol, double *rhs, double *res);
    
    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
  
    /** measure the error in the NSE */
    void MeasureErrors(DoubleFunct3D *ExactU1, DoubleFunct3D *ExactU2, DoubleFunct3D *ExactU3, DoubleFunct3D *ExactP,
                       double *u_error, double *p_error);
};

#endif
