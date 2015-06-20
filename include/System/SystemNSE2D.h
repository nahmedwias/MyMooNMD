/** ************************************************************************ 
*
* @class     TSystemNSE2D
* @brief     stores the information of a 2D NSE system matrix 
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History   Added methods (Sashi, 23.08.14)
 ************************************************************************  */


#ifndef __SYSTEMNSE2D__
#define __SYSTEMNSE2D__

#include <SquareMatrix2D.h>

/**class for 2D  NSE system matrix */
class TSystemNSE2D
{
  protected:

    /** DOFs of velocity and pressure spaces */    
    int N_U, N_P, N_Active, N_DirichletDof;
    
    /** velocity fespace */
    TFESpace2D *FeSpaces[5];
    
    /** Fe functions of NSE */
    TFEFunction2D *FeFct[5];    
    
    /** Discretization type */
    int Disctype;
       
    /** NSE type */
    int NSEType; 
           
    /** Bilinear coefficient   */
    CoeffFct2D *LinCoeffs[1];    
    
    /** NSEaux is used to pass additional fe functions (eg. mesh velocity) that is nedded for assembling */
    TAuxParam2D *NSEaux, *NSEaux_error;
    
    /** method for resudual calculation */
    DefectProc *Defect;   
    
    /** Solver type */
    int Solver;
       
    /** number of matrices in the system matrix*/
    int N_Matrices;

    /** sqstructureA of the system matrix */
    TSquareStructure2D *sqstructureA;

    /** structure of the system matrix */
    TStructure2D *structureB, *structureBT;
    
    /** A is the stiffness/system mat for NSE velocity component   */
    TSquareMatrix2D *SqmatrixA11, *SqmatrixA12, *SqmatrixA21, *SqmatrixA22, *SQMATRICES[9];
    TSquareMatrix **sqmatrices;
  
    /** B is the  system mat for NSE pressure component   */
    TMatrix2D *MatrixB1, *MatrixB2, *MatrixB1T, *MatrixB2T, *MATRICES[8];
    TMatrix **matrices;
    
    /** Boundary conditon */
    BoundCondFunct2D *BoundaryConditions[2];
  
     /** Boundary values */   
    BoundValueFunct2D *BoundaryValues[2];
        
    /** Discrete form for the equation */
    TDiscreteForm2D *DiscreteFormARhs, *DiscreteFormNL;

    
  public:
    /** constructor */
     TSystemNSE2D(TFESpace2D *velocity_fespace, TFESpace2D *presssure_fespace, TFEVectFunct2D *Velocity, 
                      TFEFunction2D *p, int disctype, int nsetype, int solver);

    /** destrcutor */
    ~TSystemNSE2D();

    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue,
              TAuxParam2D *aux);
    
    
    /** assemble the system matrix */
    void Assemble(double *sol, double *rhs);
    
    /** assemble the nonlinear part of the NSE system */
    void AssembleNonLinear(double *sol, double *rhs);
    
    /** get the resudual of the NSE system */
    void GetResidual(double *sol, double *rhs, double *res);
    
    /** solve the system matrix */
    void  Solve(double *sol, double *rhs);
  
    /** measure the error in the NSE */
    void MeasureErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
                                    double *u_error, double *p_error);
};

#endif
