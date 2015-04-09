/** ************************************************************************ 
*
* @class     TSystemMatTimeScalar2D_ALE
* @brief     stores the information of a timedependent part of a 2D scalar system matrix 
* @author    Sashikumaar Ganesan
* @date      08.08.14
* @History 
 ************************************************************************  */


#ifndef __SYSTEMMATTIMESCALAR2D_ALE__
#define __SYSTEMMATTIMESCALAR2D_ALE__

#include <SquareMatrix2D.h>
#include <SystemMatTimeScalar2D.h>

/**class for 2D scalar system matrix */
class TSystemMatTimeScalar2D_ALE : public TSystemMatTimeScalar2D
{
  protected:
    /** No. of Grid DOFs */
    int N_GridDOFs, N_GridActive, *GridKCol, *GridRowPtr;
    
    /** Grid Posistions */
    double *MeshVelo, *gridpos, *gridpos_old, *gridpos_ref, *griddisp, *GridRhs, *Entries[4];   
    
    /** grid fespace */
    TFESpace2D *GridFESpace;    
    
    /** Fgrid BC */    
    BoundValueFunct2D *GridBoundValue[1];
     
    /** grid pos vector */
    TFEVectFunct2D *GridPos, *RefGridPos;
     
    /** Fe functions of NSE */
    TFEFunction2D *MeshVeloFct[2];    
    
    /** Discrete form for moving mesh */
    TDiscreteForm2D *DiscreteFormMARhs, *DiscreteFormGrid;  
    
    /** Old M mass matrix */
    TSquareMatrix2D *sqmatrixM_old;  
    
    /** marices for the moving grid **/
    TSquareMatrix2D *SqmatrixG11, *SqmatrixG12, *SqmatrixG21, *SqmatrixG22, *SQMATRICES_GRID[4];
    
    /** structure for the moving grid **/
    TSquareStructure2D *SquareStructureG;
    
    /** aux for mesh */
    TAuxParam2D *Aux_ALE, *Meshaux;

    /** Grid bounadry conditions **/ 
    BoundCondFunct2D *GridBoundaryConditions[1];
    BoundValueFunct2D *GridBoundValues[1];
      
    /** method for Modify Mesh Coords */
    ModifyMeshCoords *ModifyCoord;
    
     /** method for Modify Mesh Coords */
    ModifyBoundCoords *ModifyBoudary;   
    
    /** */
    bool SolveLinearElastic, CONSERVATIVEALE;
    
  public:
    /** constructor */
    TSystemMatTimeScalar2D_ALE(TFESpace2D *fespace, int disctype, int solver, TFESpace2D *gridFESpace, TFEVectFunct2D *MeshVelocity, 
                                                       bool conservativeale);

//     /** destrcutor */
//     ~TSystemMatTimeScalar2D();

    /** methods */
    void Init(CoeffFct2D *BilinearCoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue, 
              CoeffFct2D *GridBilinearCoeffs, BoundCondFunct2D *GridBoundCond, BoundValueFunct2D *gridBoundValue,
              TAuxParam2D *aux);

    void AddMeshModifyFunction(ModifyMeshCoords *modifyCoord)
    {  ModifyCoord = modifyCoord; SolveLinearElastic = FALSE; }
    
    void AddBoundModifyFunction(ModifyBoundCoords *modifyboudary)
    {  ModifyBoudary = modifyboudary; SolveLinearElastic = TRUE; }  
    
    
//     /** return the stiffness matric */
//     TSquareMatrix2D *GetAMatrix()
//     { return sqmatrixA; }

    /** store M mat for next time step */
    void StoreMmat();

    /** move mesh  */   
    void MoveMesh(int N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
                                             double * Iso_refX, double Currtime);
    
    /** move mesh  */
    void MoveMesh(double Currtime);
    
    /** Get Mesh Velo */
    void GetMeshVelo(int N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
                     double * Iso_refX,  double Currtime, double tau);
    
    /** Get Mesh Velo */
    void GetMeshVelo(double Currtime, double tau);
    
    /** assemble the Mesh mat and rhs */ 
    void AssembleMeshMat();
    
    /** assemble the Mass mat and rhs */
    void AssembleMRhs(double *sol, double *rhs);   
    
    /** assemble the Mass, stiff mat and rhs */
    void AssembleMARhs(double *sol, double *rhs); 

    /** M = M + (tau*TDatabase::TimeDB->THETA1)*A */ 
    /** B = (tau*TDatabase::TimeDB->THETA1)*rhs +(tau*TDatabase::TimeDB->THETA2)*oldrhs + [ M - (tau*TDatabase::TimeDB->THETA2)A]*oldsol */  
    void AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol);
 
     /** solve the system matrix */
    void Solve(double *sol, double *rhs);  
    
//     /** return the residual of the system for the given sol*/
//     double GetResidual(double *sol);
//     
};

#endif
