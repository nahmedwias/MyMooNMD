/** ************************************************************************ 
* @brief     source file for TSystemMatTimeScalar2D_ALE
* @author    Sashikumaar Ganesan
* @date      19.02.15
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemMatTimeScalar2D_ALE.h>
#include <SystemMatTimeScalar2D.h>
#include <FEVectFunct2D.h>
#include <TimeConvDiff2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
#include <LinAlg.h>
#include <FreeSurface2D.h>

TSystemMatTimeScalar2D_ALE::TSystemMatTimeScalar2D_ALE(TFESpace2D *fespace, int disctype, int solver, TFESpace2D *gridFESpace, TFEVectFunct2D *MeshVelocity, 
                                                       bool conservativeale) : TSystemMatTimeScalar2D(fespace,  disctype, solver)
{
  char WString[] = "w";  
  TFESpace2D *fesp[1];
    
  
  /** old M mass matrix */
  sqmatrixM_old = new TSquareMatrix2D(sqstructure);  
  
  GridFESpace = gridFESpace;
  N_GridDOFs = gridFESpace->GetN_DegreesOfFreedom();
  N_GridActive = gridFESpace->GetActiveBound();
  
  // grid 
  SquareStructureG= new TSquareStructure2D(GridFESpace); 
  SquareStructureG->Sort();
   
  // for mesh
  SqmatrixG11 = new TSquareMatrix2D(SquareStructureG); // G11
  SqmatrixG12 = new TSquareMatrix2D(SquareStructureG); // G12
  SqmatrixG21 = new TSquareMatrix2D(SquareStructureG); // G21
  SqmatrixG22 = new TSquareMatrix2D(SquareStructureG); // G22
  
  SQMATRICES_GRID[0] = SqmatrixG11;
  SQMATRICES_GRID[1] = SqmatrixG12;
  SQMATRICES_GRID[2] = SqmatrixG21;
  SQMATRICES_GRID[3] = SqmatrixG22;

  
  Entries[0] = SqmatrixG11->GetEntries();
  Entries[1] = SqmatrixG12->GetEntries();
  Entries[2] = SqmatrixG21->GetEntries();
  Entries[3] = SqmatrixG22->GetEntries();

  GridKCol = SquareStructureG->GetKCol();
  GridRowPtr = SquareStructureG->GetRowPtr();
     
  fesp[0] = GridFESpace;
  Meshaux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);  
    
  MeshVeloFct[0] = MeshVelocity->GetComponent(0);
  MeshVeloFct[1] = MeshVelocity->GetComponent(1);
  MeshVelo =  MeshVelocity->GetValues();
  
  gridpos = new double[2*N_GridDOFs];
  gridpos_old = new double[2*N_GridDOFs];   
  gridpos_ref = new double[2*N_GridDOFs];   
  griddisp = new double[2*N_GridDOFs];   
  
  
   memset(gridpos, 0, 2*N_GridDOFs*SizeOfDouble);
   GridPos = new TFEVectFunct2D(GridFESpace, WString, WString, gridpos, N_GridDOFs, 2);  
   GridPos->GridToData();
   
   GridRhs = new double[2*N_GridDOFs];
   
   RefGridPos = new TFEVectFunct2D(GridFESpace, WString, WString, gridpos_ref, N_GridDOFs, 2);     
 
   memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);
   memcpy(gridpos_ref, gridpos, 2*N_GridDOFs*SizeOfDouble); 
   
   Aux_ALE = NULL;
   SolveLinearElastic = TRUE;
   CONSERVATIVEALE = conservativeale;
} // TSystemMatTimeScalar2D_ALE


void TSystemMatTimeScalar2D_ALE::Init(CoeffFct2D *BilinearCoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue, 
                                      CoeffFct2D *GridBilinearCoeffs, BoundCondFunct2D *GridBoundCond, BoundValueFunct2D *gridBoundValue)
{
  BoundaryConditions[0] =  BoundCond;
  BoundaryValues[0] = BoundValue;
  
  GridBoundValue[0] = gridBoundValue;
  GridBoundaryConditions[0] = GridBoundCond;
  
  TDiscreteForm2D *DiscreteFormMRhs_Galerkin;
  TDiscreteForm2D *DiscreteFormARhs_Galerkin; 
  TDiscreteForm2D *DiscreteFormMRhs_SUPG;
  TDiscreteForm2D *DiscreteFormARhs_SUPG;
  
  TDiscreteForm2D *DiscreteFormMARhs_Galerkin;
  TDiscreteForm2D *DiscreteFormMatrixMARhs_SUPG;
     
   InitializeDiscreteForms_ScalarMoving(DiscreteFormMARhs_Galerkin, DiscreteFormGrid, DiscreteFormMatrixMARhs_SUPG,
                                        BilinearCoeffs, GridBilinearCoeffs);
   
  
  InitializeDiscreteFormsScalar(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormMRhs_SUPG,
                                  DiscreteFormARhs_SUPG, BilinearCoeffs);
  
    switch(Disctype)
     {
      case GALERKIN:
      case LOCAL_PROJECTION:
           DiscreteFormARhs = DiscreteFormARhs_Galerkin;
           DiscreteFormMRhs = DiscreteFormMRhs_Galerkin;
           DiscreteFormMARhs = DiscreteFormMARhs_Galerkin;
      break;
      
      case SUPG:
           DiscreteFormARhs = DiscreteFormARhs_SUPG;
           DiscreteFormMRhs = DiscreteFormMRhs_SUPG;
           DiscreteFormMARhs = DiscreteFormMatrixMARhs_SUPG;          
      break;
      
      default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
     }  
       
} // Init


void TSystemMatTimeScalar2D_ALE::StoreMmat()
 { 
   
   sqmatrixM_old->Reset(); 
   MatAdd(sqmatrixM_old, sqmatrixM, 1.); 
   
 } // TSystemMatTimeScalar2D_ALE::StoreMmat()

 
void TSystemMatTimeScalar2D_ALE::MoveMesh(int N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
                                             double * Iso_refX, double Currtime)
{
 int i, N_GridBDDOFs;
 
  GridPos->GridToData();   
  memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);  
 
  // modyfy the boundary 
  RefGridPos->DataToGrid();  
  ModifyBoudary(N_MovVert, MovBoundVert, Free_Joint, Iso_refX, Currtime);    
  
  // data with updated BD values
  GridPos->GridToData();  

  N_GridBDDOFs = N_GridDOFs - N_GridActive;
  
  memset(GridRhs, 0, 2*N_GridDOFs*SizeOfDouble);     
  memcpy(GridRhs+N_GridActive, gridpos+N_GridActive, N_GridBDDOFs*SizeOfDouble);     //rhs1  
  memcpy(GridRhs+(N_GridDOFs+N_GridActive), gridpos+(N_GridDOFs+N_GridActive), N_GridBDDOFs*SizeOfDouble);   //rhs2       
    
  Daxpy(N_GridBDDOFs, -1., gridpos_old+N_GridActive, GridRhs+N_GridActive);
  Daxpy(N_GridBDDOFs, -1., gridpos_old+(N_GridDOFs+N_GridActive), GridRhs+(N_GridDOFs+N_GridActive));    
     
  memcpy(griddisp, GridRhs, 2*N_GridDOFs*SizeOfDouble);   
    
  SolveGridEquation(Entries, griddisp, GridRhs, GridKCol, GridRowPtr, N_GridDOFs);
           
  memcpy(gridpos, gridpos_old, 2*N_GridDOFs*SizeOfDouble);
  Daxpy(2*N_GridDOFs, 1., griddisp, gridpos);
  GridPos->DataToGrid();     

} // MoveMesh
 
void TSystemMatTimeScalar2D_ALE::MoveMesh(double Currtime)
{
 int i;
    
   // domain at the given time  
   for(i=0;i<N_GridDOFs;i++)
     ModifyCoord(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], Currtime);   
  
   // move the mesh
   GridPos->DataToGrid(); 
} //TSystemMatTimeScalar2D_ALE::MoveMesh(


void TSystemMatTimeScalar2D_ALE::GetMeshVelo(int N_MovVert, TVertex **MovBoundVert, TIsoBoundEdge **Free_Joint, 
                                             double * Iso_refX,  double Currtime, double tau)
{
 int i, N_GridBDDOFs;
 
  GridPos->GridToData();   
  memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);  
 
  // modyfy the boundary 
  ModifyBoudary(N_MovVert, MovBoundVert, Free_Joint, Iso_refX, Currtime);    
  
  // data with updated BD values
  GridPos->GridToData();  

  N_GridBDDOFs = N_GridDOFs - N_GridActive;
  
  memset(GridRhs, 0, 2*N_GridDOFs*SizeOfDouble);     
  memcpy(GridRhs+N_GridActive, gridpos+N_GridActive, N_GridBDDOFs*SizeOfDouble);     //rhs1  
  memcpy(GridRhs+(N_GridDOFs+N_GridActive), gridpos+(N_GridDOFs+N_GridActive), N_GridBDDOFs*SizeOfDouble);   //rhs2       
    
  Daxpy(N_GridBDDOFs, -1., gridpos_old+N_GridActive, GridRhs+N_GridActive);
  Daxpy(N_GridBDDOFs, -1., gridpos_old+(N_GridDOFs+N_GridActive), GridRhs+(N_GridDOFs+N_GridActive));    
     
  memcpy(MeshVelo, GridRhs, 2*N_GridDOFs*SizeOfDouble);   
    
  SolveGridEquation(Entries, MeshVelo, GridRhs, GridKCol, GridRowPtr, N_GridDOFs);
         
  Dscal(2*N_GridDOFs, -1./tau, MeshVelo); // - sign due to -w\cdot\nabla C in the equation    
  
  memcpy(gridpos, gridpos_old, 2*N_GridDOFs*SizeOfDouble);
  //   Daxpy(2*N_GridDOFs, 1., MeshVelo, gridpos);
  // move the mesh back to original pos, as we calculate only the mesh velocity
  GridPos->DataToGrid();   

} // TSystemMatTimeScalar2D_ALE::GetMeshVelo


void TSystemMatTimeScalar2D_ALE::GetMeshVelo(double Currtime, double tau)
{
 int i;
  
  GridPos->GridToData();   
  memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);  

   //  move velo in current time  
   for(i=0;i<N_GridDOFs;i++)
     ModifyCoord(gridpos_ref[i], gridpos_ref[i+N_GridDOFs], gridpos[i], gridpos[i+N_GridDOFs], Currtime);   

   //compute mesh velocity
   memcpy(MeshVelo, gridpos, 2*N_GridDOFs*SizeOfDouble);     
   Daxpy(2*N_GridDOFs, -1., gridpos_old, MeshVelo);        
   Dscal(2*N_GridDOFs, -1./tau, MeshVelo); // - sign du*/ //e to -w\cdot\nabla C in the equation   
   memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble); 
} //TSystemMatTimeScalar2D_ALE::GetMeshVelo(
  
  
  
void TSystemMatTimeScalar2D_ALE::AssembleMeshMat()
{
  
 TFESpace2D *fesp[1];

   fesp[0] = GridFESpace; 
   SQMATRICES_GRID[0]->Reset();
   SQMATRICES_GRID[1]->Reset();
   SQMATRICES_GRID[2]->Reset();
   SQMATRICES_GRID[3]->Reset();    

     Assemble2D(1, fesp,
             4, SQMATRICES_GRID,
             0, NULL,
             0, NULL, NULL,
             DiscreteFormGrid,
             GridBoundaryConditions,
             GridBoundValue,
             Meshaux);
     
  // for Dirichlet rows in off-diagonal matrices
  memset(Entries[1] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);
  memset(Entries[2] + GridRowPtr[N_GridActive], 0, (GridRowPtr[N_GridDOFs] - GridRowPtr[N_GridActive])*SizeOfDouble);     
  
}

void TSystemMatTimeScalar2D_ALE::AssembleMRhs(double *sol, double *rhs)
{
  int N_DOF, N_Active, N_SquareMatrices;
  double *RHSs[1];

  TFESpace2D *fesp[2], *ferhs[1];
  TFEFunction2D  *fefct[4];
  TAuxParam2D *aux;
  
    N_DOF = FeSpace->GetN_DegreesOfFreedom();
    N_Active =  FeSpace->GetActiveBound();
    
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);
  
    fesp[0] = FeSpace;
    fesp[1] = GridFESpace;    
    ferhs[0] = FeSpace;
    
    // initialize matrices
    SQMATRICES[0] = sqmatrixM;
    SQMATRICES[0]->Reset();   

    N_SquareMatrices =1;
    
   if(Disctype == SUPG)
    {
     N_SquareMatrices = 2;  
     SQMATRICES[1] = sqmatrixS;
     SQMATRICES[1]->Reset();  
    }      

    
   aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); 

    // assemble
    Assemble2D(2, fesp,
               N_SquareMatrices, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteFormMRhs,
               BoundaryConditions,
               BoundaryValues,
               aux);
     
    delete aux;
    
   // copy Dirichlet values from rhs into sol
   memcpy(sol+N_Active, rhs+N_Active, (N_DOF - N_Active)*SizeOfDouble);  

   if(Disctype==SUPG)
    {     
     MatAdd(sqmatrixM, sqmatrixS, 1.);
    }
     
} // TSystemMatScalar2D::AssembleMRhs 

void TSystemMatTimeScalar2D_ALE::AssembleMARhs(double *sol, double *rhs)
{
  int N_DOF, N_Active, N_SquareMatrices;
  double *RHSs[1];

  TFESpace2D *fesp[2], *ferhs[1];
  TFEFunction2D  *fefct[4];
    
    // assemble the mass mat and rhs 
    N_DOF = FeSpace->GetN_DegreesOfFreedom();
    N_Active =  FeSpace->GetActiveBound();
    
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);
  
    fesp[0] = FeSpace;
    fesp[1] = GridFESpace;    
    ferhs[0] = FeSpace;
    
    // initialize matrices
    SQMATRICES[0] = sqmatrixA;
    SQMATRICES[0]->Reset();   
    SQMATRICES[1] = sqmatrixM;
    SQMATRICES[1]->Reset();   
    
    N_SquareMatrices =2;
    
   if(Disctype == SDFEM)
    {
     N_SquareMatrices = 3;  
     SQMATRICES[2] = sqmatrixS;
     SQMATRICES[2]->Reset();  
    }      

    fefct[0] = MeshVeloFct[0];
    fefct[1] = MeshVeloFct[1];
    fefct[2] = MeshVeloFct[0]; //for calculating divergence of w
    fefct[3] = MeshVeloFct[1]; //for calculating divergence of w    

    if(Aux_ALE==NULL)
     { 
      // defined in TimeConvDiff2D.h
      Aux_ALE =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
                            TimeCDParamsVeloFieldN_Fct,
                            TimeCDParamsVeloFieldN_ParamFct,
                            TimeCDParamsVeloFieldN_FEValues_ALE,
                            fesp+1, fefct,
                            TimeCDParamsVeloFieldFct_ALE,
                            TimeCDParamsVeloFieldFEFctIndex_ALE,
                            TimeCDParamsVeloFieldFEMultiIndex_ALE,
                            TimeCDParamsVeloFieldN_Params_ALE,
                            TimeCDParamsVeloFieldBeginParam);

     }

    // assemble
    Assemble2D(2, fesp,
               N_SquareMatrices, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteFormMARhs,
               BoundaryConditions,
               BoundaryValues,
               Aux_ALE);
     
      // copy Dirichlet values from rhs into sol
      memcpy(sol+N_Active, rhs+N_Active, (N_DOF - N_Active)*SizeOfDouble);  
      
      
    if(Disctype==SUPG)
     {     
      MatAdd(sqmatrixM, sqmatrixS, 1.);
     }
     
//                  memset(defect, 0, N_DOF*SizeOfDouble);   
//         MatVectActive(sqmatrixA, sol, defect);
//         OutPut("MatA-Sol " << Ddot(N_DOF, defect, defect ) << endl);  
//        OutPut("MatA-Sol " << Ddot(N_DOF, sol, sol ) << endl); 
//     
//     cout << " AssembleMARhs " << endl;
//   exit(0);   
  
} // TSystemMatScalar2D::AssembleMARhs 

void TSystemMatTimeScalar2D_ALE::AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol)
{
    int N_DOF, N_Active, N_SquareMatrices;
    double tau;
    
    N_DOF = FeSpace->GetN_DegreesOfFreedom();
    N_Active =  FeSpace->GetActiveBound();   
    
    tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  
    memset(B, 0, N_DOF*SizeOfDouble);      
    
    // old rhs multiplied with current subtime step and theta3 on B
    Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3,  oldrhs, B);    
    
    // add rhs from current sub time step to rhs array B
    Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4,  rhs, B);    
    
    // defect = M * oldsol
    memset(defect, 0, N_DOF*SizeOfDouble);      

    if(CONSERVATIVEALE)        
     {
      MatAdd(sqmatrixM_old, sqmatrixA, - tau*TDatabase::TimeDB->THETA2);
      gamma = 0.;

      MatVectActive(sqmatrixM_old, oldsol, defect);
      
      // store the mass matrix for next time step (only for linear problem)
      this->StoreMmat();
     }
    else
     {
      MatAdd(sqmatrixM, sqmatrixA, - tau*TDatabase::TimeDB->THETA2);
      gamma = -tau*TDatabase::TimeDB->THETA2;  
     
      MatVectActive(sqmatrixM, oldsol, defect);  
     }

     // B:= B + defec    
     Daxpy(N_Active, 1, defect, B);
     
//    OutPut("B  " << Ddot(N_DOF, oldsol, oldsol) << endl); 
//       OutPut("B  " << Ddot(N_DOF, B, B ) << endl);  
//       exit(0);
     // set Dirichlet values
     memcpy(B+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
     memcpy(sol+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
     
     //assemble the system matrix
     MatAdd(sqmatrixM, sqmatrixA, -gamma + tau*TDatabase::TimeDB->THETA1);
     
//                OutPut("B  " << Ddot(N_DOF, B, B ) << endl);
} // AssembleSystMat
 
void TSystemMatTimeScalar2D_ALE::Solve(double *sol, double *rhs)
{
  
    switch(SOLVER)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:
        cout << "GMG solver not yet implemented " <<endl;
      break;

      case DIRECT:
        DirectSolver(sqmatrixM, B, sol);
      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
}


// double TSystemMatTimeScalar2D_ALE::GetResidual(double *sol)
// {
//   int N_DOF = FeSpace->GetN_DegreesOfFreedom(); 
//   double residual_scalar;
//   
// //   if(SystMatAssembled)
//    {
//     memset(defect, 0, N_DOF*SizeOfDouble);         
//     ScalarDefect(sqmatrixM, sol, B, defect, residual_scalar);   
//    }
// //   else
//    {
// //     OutPut("Assemble the System Matrix before calculating the GetResidual" << endl);
// //     exit(4711);;   
//    }
//    
//    return residual_scalar;
//       
// }
// 
// 
