/** ************************************************************************ 
*
* @class     TSystemMatTNSE2D_ALE
* @brief     stores the information of a 2D TNSE system matrix 
* @author    Sashikumaar Ganesan,Birupaksha Pal 
* @date      03.09.14
* @History   Added methods (Sashi, 03.09.14)
 ************************************************************************  */

#ifndef __SYSTEMMATTNSE2D_ALE__
#define __SYSTEMMATTNSE2D_ALE__

#include <SystemMatTNSE2D.h>

/**class for 2D  TNSE system matrix */
class TSystemMatTNSE2D_ALE : public TSystemMatTNSE2D
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
    
//     /** M - mass/system mat for TNSE velocity component   */
//     TSquareMatrix2D *OldSqmatrixM11, *OldSqmatrixM12, *OldSqmatrixM21, *OldSqmatrixM22;

       
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
     TSystemMatTNSE2D_ALE(TFESpace2D *velocity_fespace, 
                          TFESpace2D *presssure_fespace,
                          TFEVectFunct2D *Velocity,
                          TFEFunction2D *p, 
#ifdef __PRIVATE__
                          TFESpace2D *Projection_space,
#endif
                          TFESpace2D *gridFESpace,
                          TFEVectFunct2D *MeshVelocity, bool conservativeale);
// 
//     /** destrcutor */
//     ~TSystemMatTNSE2D_ALE();
// 
    /** methods */
    /** Initilize the discrete forms and the matrices */    
    void Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, 
              BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue,
              CoeffFct2D *GridBilinearCoeffs, BoundCondFunct2D *GridBoundCond,
              BoundValueFunct2D *gridBoundValue, TAuxParam2D *aux,
              TAuxParam2D *nseaux_error);
    
    
    void AddMeshModifyFunction(ModifyMeshCoords *modifyCoord)
    {  ModifyCoord = modifyCoord; SolveLinearElastic = FALSE; }
    
    void AddBoundModifyFunction(ModifyBoundCoords *modifyboudary)
    {  ModifyBoudary = modifyboudary; SolveLinearElastic = TRUE; }  
 
    /** assemble the M, A and rhs */
    void Assemble(double *sol, double *rhs);
  
//     /** store M mat for next time step */
//     void StoreMmat();
    
//     /** move mesh  */   
//     void MoveMesh(int N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
//                                              double * Iso_refX, double Currtime);
//     
//     /** move mesh  */
//     void MoveMesh(double Currtime);
    
    /** Get Mesh Velo */
    void GetMeshVeloAndMove(int N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
                     double * Iso_refX,  double Currtime, double tau);
    
    /** Get Mesh Velo */
    void GetMeshVeloAndMove(double Currtime, double tau);
    
    /** assemble the Mesh mat and rhs */ 
    void AssembleMeshMat();
    
    /** scale B matices and assemble rhs based on the \theta scheme  */
    void AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol);

    /** scale B matices and assemble rhs based on the \theta scheme  */
    void AssembleSystMatNonLinear();
     
    /** restoring the mass matrix */
    void RestoreMassMat();
        
    /** assemble the nonlinear part of the NSE system */
    void AssembleANonLinear(double *sol, double *rhs);

    /** solve the system matrix */
    void  Solve(double *sol);
    
    /** get the resudual of the NSE system */
    void GetTNSEResidual(double *sol, double *res);   
 
//     /** measure the error in the NSE */
//     void MeasureTNSEErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP, double *AllErrors);
//     
//     void GridMovingScheme(TFESpace2D *GridFESpace);
};

#endif
