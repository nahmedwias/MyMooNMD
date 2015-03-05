/** ************************************************************************ 
* @brief     source file for TSystemMatTimeScalar2D
* @author    Sashikumaar Ganesan
* @date      08.08.14
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemMatTimeScalar2D.h>
#include <SystemMatScalar2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>

TSystemMatTimeScalar2D::TSystemMatTimeScalar2D(TFESpace2D *fespace, int disctype, int solver): TSystemMatScalar2D(fespace,  disctype, solver)
{
  /** need it for solver */
  sqmatrices = (TSquareMatrix **)SQMATRICES;
  
  /** M mass matrix */
  sqmatrixM = new TSquareMatrix2D(sqstructure);  
  N_Matrices++;
  
  /** working rhs, used in AssembleSystMat() */
  int N_DOF = FeSpace->GetN_DegreesOfFreedom(); 
  B = new double[N_DOF];
  defect = new double[N_DOF];
  
  gamma =0.;
  
  /** time-consistent part of the SUPG matrix */
  if(Disctype==SDFEM || Disctype==SUPG)
   {
    sqmatrixS = new TSquareMatrix2D(sqstructure); 
    N_Matrices++;
    
    sqmatrixK = new TSquareMatrix2D(sqstructure); 
    N_Matrices++; 
   }
   
  SystMatAssembled  = FALSE;
} // constructor


TSystemMatTimeScalar2D::~TSystemMatTimeScalar2D()
{
  delete sqstructure;
  delete sqmatrixA; 
  
  delete sqmatrixM;
  if(Disctype==SDFEM || Disctype==SUPG)
   {
     delete sqmatrixS;
     delete sqmatrixK;    
   }
}


void TSystemMatTimeScalar2D::Init(CoeffFct2D *BilinearCoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue)
{
  BoundaryConditions[0] =  BoundCond;
  BoundaryValues[0] = BoundValue;
    
  TDiscreteForm2D *DiscreteFormMRhs_Galerkin;
  TDiscreteForm2D *DiscreteFormARhs_Galerkin; 
  TDiscreteForm2D *DiscreteFormMRhs_SUPG;
  TDiscreteForm2D *DiscreteFormARhs_SUPG;

  
  InitializeDiscreteFormsScalar(DiscreteFormMRhs_Galerkin, DiscreteFormARhs_Galerkin, DiscreteFormMRhs_SUPG,
                                  DiscreteFormARhs_SUPG, BilinearCoeffs);
  
    switch(Disctype)
     {
      case GALERKIN:
      case LOCAL_PROJECTION:
           DiscreteFormARhs = DiscreteFormARhs_Galerkin;
           DiscreteFormMRhs = DiscreteFormMRhs_Galerkin;
      break;
      
      case SUPG:
           DiscreteFormARhs = DiscreteFormARhs_SUPG;
           DiscreteFormMRhs = DiscreteFormMRhs_SUPG;
      break;
      
      default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
     }  
       
  
} // Init


void TSystemMatTimeScalar2D::AssembleMRhs(TAuxParam2D *aux, double *sol, double *rhs)
{
  int N_DOF, N_Active, N_SquareMatrices;
  double *RHSs[1];
 

  TFESpace2D *fesp[1], *ferhs[1];
//   TSquareMatrix2D *SQMATRICES[2];

    
    N_DOF = FeSpace->GetN_DegreesOfFreedom();
    N_Active =  FeSpace->GetActiveBound();
//     N_DirichletDof = N_DOF - N_Active;
    
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);
  
    fesp[0] = FeSpace;
    ferhs[0] = FeSpace;
    
    // initialize matrices
    SQMATRICES[0] = sqmatrixM;
    SQMATRICES[0]->Reset();    
    N_SquareMatrices =1;
    
   if(Disctype == SDFEM)
    {
     N_SquareMatrices = 2;  
     SQMATRICES[1] = sqmatrixS;
     SQMATRICES[1]->Reset();  
    }      
    
    if(aux==NULL)
     { aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); }
    
    // assemble
    Assemble2D(1, fesp,
               N_SquareMatrices, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteFormMRhs,
               BoundaryConditions,
               BoundaryValues,
               aux);
     
      // copy Dirichlet values from rhs into sol
      memcpy(sol+N_Active, rhs+N_Active, (N_DOF - N_Active)*SizeOfDouble);  
      
} // TSystemMatScalar2D::AssembleMRhs 



void TSystemMatTimeScalar2D::AssembleARhs(TAuxParam2D *aux, double *sol, double *rhs)
{
  int N_DOF, N_Active, N_SquareMatrices;
  double *RHSs[1];
 

  TFESpace2D *fesp[1], *ferhs[1];
//   TSquareMatrix2D *SQMATRICES[2];

//     BoundCondFunct2D *BoundaryConditions[1] = { BoundaryCond };
//     BoundValueFunct2D *BoundaryValues[1] = { BoundaryValue };
    
    N_DOF = FeSpace->GetN_DegreesOfFreedom();
    N_Active =  FeSpace->GetActiveBound();
    
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);
  
    fesp[0] = FeSpace;
    ferhs[0] = FeSpace;
    
    // initialize matrices
    SQMATRICES[0] = sqmatrixA;
    SQMATRICES[0]->Reset();    
    N_SquareMatrices =1;
    
   if(Disctype == SDFEM)
    {
     N_SquareMatrices = 2;  
     SQMATRICES[1] = sqmatrixK;
     SQMATRICES[1]->Reset();  
    }      
    
    if(aux==NULL)
     { aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); }
    
    // assemble
    Assemble2D(1, fesp,
               N_SquareMatrices, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteFormARhs,
               BoundaryConditions,
               BoundaryValues,
               aux);
     
      // copy Dirichlet values from rhs into sol
      memcpy(sol+N_Active, rhs+N_Active, (N_DOF - N_Active)*SizeOfDouble);  
      
} // TSystemMatScalar2D::AssembleARhs 

void TSystemMatTimeScalar2D::AssembleSystMat(double *oldrhs, double *oldsol, double *rhs, double *sol)
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
    
//        cout << "AssembleSystMat  gamma " << gamma << endl; 
    // M = M + (- tau*TDatabase::TimeDB->THETA2)
    MatAdd(sqmatrixM, sqmatrixA, - tau*TDatabase::TimeDB->THETA2);

    if(Disctype==SUPG)
     {     
      MatAdd(sqmatrixM, sqmatrixS, 1.);
     }
     // set current factor of steady state matrix
     gamma = -tau*TDatabase::TimeDB->THETA2;  
    
     // defect = M * oldsol
     memset(defect, 0, N_DOF*SizeOfDouble);  
     MatVectActive(sqmatrixM, oldsol, defect); 
     
     // B:= B + defec    
     Daxpy(N_Active, 1, defect, B);
    
     // set Dirichlet values
     memcpy(B+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
     memcpy(sol+N_Active, rhs+N_Active, (N_DOF-N_Active)*SizeOfDouble);
     
     //assemble the system matrix
     MatAdd(sqmatrixM, sqmatrixA, -gamma + tau*TDatabase::TimeDB->THETA1);
     gamma = tau*TDatabase::TimeDB->THETA1;    
   
     SystMatAssembled  = TRUE;     
     
//      cout << "AssembleSystMat  gamma " << gamma << endl;  
     
} // AssembleSystMat

void TSystemMatTimeScalar2D::RestoreMassMat()
{

//   cout << "RestoreMassMat  gamma " << gamma << endl;
  if(SystMatAssembled)
   {
     // restore the mass matrix
     MatAdd(sqmatrixM, sqmatrixA, -gamma);
     gamma = 0.;
     
     if(Disctype==SUPG)
     {     
      MatAdd(sqmatrixM, sqmatrixS, -1.);
     }
  
     SystMatAssembled  = FALSE;  
   }
  else
  {
    cout << "System is not assembled to restore " <<endl;
  }
   
//   cout << "RestoreMassMat" << endl;
//   exit(0);
  
}

void TSystemMatTimeScalar2D::Solve(double *sol, double *rhs)
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


double TSystemMatTimeScalar2D::GetResidual(double *sol)
{
  int N_DOF = FeSpace->GetN_DegreesOfFreedom(); 
  double residual_scalar;
  
//   if(SystMatAssembled)
   {
    memset(defect, 0, N_DOF*SizeOfDouble);         
    ScalarDefect(sqmatrixM, sol, B, defect, residual_scalar);   
   }
//   else
   {
//     OutPut("Assemble the System Matrix before calculating the GetResidual" << endl);
//     exit(4711);;   
   }
   
   return residual_scalar;
      
}


