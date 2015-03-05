/** ************************************************************************ 
* @brief     source file for TSystemMatTNSE2D
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemMatNSE2D.h>
#include <SystemMatTNSE2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
// #include <TNSE2D_ParamRout.h>
#include <MainUtilities.h>
#include <Upwind.h>

#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>

#ifdef __PRIVATE__  
 #include <VMS.h>
#endif

TSystemMatTNSE2D::TSystemMatTNSE2D(TFESpace2D *velocity_fespace, TFESpace2D *presssure_fespace, TFEVectFunct2D *Velocity, 
                                   TFEFunction2D *p, int disctype, int nsetype, int solver
#ifdef __PRIVATE__  
                                   ,TFESpace2D *Projection_space
#endif    
                                   ) : TSystemMatNSE2D(velocity_fespace, presssure_fespace, Velocity, p, disctype, nsetype, solver)
{
  B = new double[2*N_U+N_P];
  defect = new double[2*N_U+N_P];

  gamma =0.;  
  
    // allocate the mass matrices in addition
    switch(NSEType)
     {
      case 1:
        SqmatrixM11 = new TSquareMatrix2D(sqstructureA);


      break;

      case 2:
        SqmatrixM11 = new TSquareMatrix2D(sqstructureA);
      break;

      case 3:
        SqmatrixM11 = new TSquareMatrix2D(sqstructureA);
        SqmatrixM12 = new TSquareMatrix2D(sqstructureA);
        SqmatrixM21 = new TSquareMatrix2D(sqstructureA);
        SqmatrixM22 = new TSquareMatrix2D(sqstructureA);
      break;

      case 4:
        SqmatrixM11 = new TSquareMatrix2D(sqstructureA);
        SqmatrixM12 = new TSquareMatrix2D(sqstructureA);
        SqmatrixM21 = new TSquareMatrix2D(sqstructureA);
        SqmatrixM22 = new TSquareMatrix2D(sqstructureA);
      break;
      
      default:
            OutPut("Unknown NSETYPE, it must be 1 to 4" << endl);
            exit(4711);;      
      
     }  
 
#ifdef __PRIVATE__ 
   if (Disctype == VMS_PROJECTION)
   { 
     if(NSEType==1 || NSEType==2)
     {
      Error("NSETYPE should be 3 or 4 for VMS_PROJECTION !!!" << endl);
      exit(-1);          
     }
    FeSpaces[2] = velocity_fespace; //  to be included the convolution space if needed
    FeSpaces[3] = Projection_space; 
     
    sqstructureL = new TSquareStructure2D(Projection_space);
    sqstructureL->Sort();
    structure_tilde_G = new TStructure2D(velocity_fespace, Projection_space);
    structure_G = new TStructure2D(Projection_space, velocity_fespace);    

    MatricesL = new TSquareMatrix2D(sqstructureL);
    Matrices_tilde_G11 = new TMatrix2D(structure_tilde_G);
    Matrices_tilde_G22 = new TMatrix2D(structure_tilde_G);    
    Matrices_G11 = new TMatrix2D(structure_G);
    Matrices_G22 = new TMatrix2D(structure_G);
   }
#endif

 NSE_Rhsaux = NULL;
 SystMatAssembled  = FALSE;
 olderror_l_2_l_2u = 0.;
}

TSystemMatTNSE2D::~TSystemMatTNSE2D()
{
//     delete NSEaux; 
//        
//     if(NSEaux_error)
//       delete NSEaux_error;

}


void TSystemMatTNSE2D::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue,
                            TAuxParam2D *aux, TAuxParam2D *nseaux_error)
{
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormColetti;
  TDiscreteForm2D *DiscreteFormGL00Convolution;
  TDiscreteForm2D *DiscreteFormGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormSmagorinsky;
  TDiscreteForm2D *DiscreteFormVMSProjection;

  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLColetti;
  TDiscreteForm2D *DiscreteFormNLGL00Convolution;
  TDiscreteForm2D *DiscreteFormNLGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky;
  TDiscreteForm2D *DiscreteFormNLVMSProjection;

  TDiscreteForm2D *DiscreteFormRHS;
  TDiscreteForm2D *DiscreteFormRHSColetti;
  TDiscreteForm2D *DiscreteFormRHSSmagorinskyExpl;
  TDiscreteForm2D *DiscreteFormMatrixGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormGL00AuxProblemRHS;
  TDiscreteForm2D *DiscreteFormRHSLESModel;
  TDiscreteForm2D *DiscreteFormRHSAuxProblemU;
  TDiscreteForm2D *DiscreteFormMatrixAuxProblemU;
  
  // save the boundary condition
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  

  // save the boundary values  
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
 
  // save the nse bilinear coefficient   
  LinCoeffs[0] = lincoeffs;
  
  //default, i.e., velocity for nonlinear term
  NSEaux = aux;
  
  NSE_Rhsaux = new TAuxParam2D(1, 0, 0, 0, FeSpaces, NULL, NULL, NULL, NULL, 0, NULL);
  
  // aux for calculating the error
  NSEaux_error = nseaux_error;
  
  // set the Discreteforms
  InitializeDiscreteForms(DiscreteFormGalerkin,DiscreteFormUpwind,
              DiscreteFormSmagorinsky,DiscreteFormColetti,
              DiscreteFormGL00Convolution,DiscreteFormGL00AuxProblem,
              DiscreteFormVMSProjection,
              DiscreteFormNLGalerkin,
              DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky,
              DiscreteFormNLColetti,DiscreteFormNLGL00Convolution,
              DiscreteFormNLGL00AuxProblem,
              DiscreteFormNLVMSProjection,
              DiscreteFormRHS,
              DiscreteFormRHSColetti,
              DiscreteFormRHSLESModel,
              DiscreteFormMatrixGL00AuxProblem,
              DiscreteFormGL00AuxProblemRHS,
              DiscreteFormRHSSmagorinskyExpl,
              DiscreteFormMatrixAuxProblemU,
              DiscreteFormRHSAuxProblemU,
              LinCoeffs[0], NSEType);

    // find discrete form
    switch(Disctype)
       {
          case GALERKIN:
            DiscreteFormARhs = DiscreteFormGalerkin;
            DiscreteFormNL = DiscreteFormNLGalerkin;
            DiscreteFormRhs = DiscreteFormRHS;
          break;

          case UPWIND:
            DiscreteFormARhs = DiscreteFormUpwind;
            DiscreteFormNL = DiscreteFormNLUpwind;    
            DiscreteFormRhs = DiscreteFormRHS;
            break;

//           case SMAGORINSKY:
//             DiscreteFormARhs = DiscreteFormSmagorinsky;
//             DiscreteFormNL = DiscreteFormNLSmagorinsky;              
//             break;

#ifdef __PRIVATE__ 
          case VMS_PROJECTION:
    
            DiscreteFormARhs = DiscreteFormVMSProjection;
            DiscreteFormNL = DiscreteFormNLVMSProjection;
            DiscreteFormRhs = DiscreteFormRHS;    
            break;
#endif
	    
          default:
            Error("Unknown DISCTYPE" << Disctype << endl);
            exit(-1);
        } 
     
     // set the discrete form for the Stokes equation
      if (TDatabase::ParamDB->STOKES_PROBLEM)
       {
        DiscreteFormARhs = DiscreteFormUpwind;     
        DiscreteFormNL = NULL;
       }

} // TSystemMatTNSE2D::Init

 
 
/* Assemble M, A and rhs */ 
void TSystemMatTNSE2D::Assemble(double *sol, double *rhs)
{
  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  
  double *RHSs[3];

  TFESpace2D *fesprhs[3];
  
  N_Rhs = 2;
  N_FESpaces = 2;   
      
     // initialize matrices
     switch(NSEType)
      {
        case 1:
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixM11;
         MATRICES[0] = MatrixB1;
         MATRICES[1] = MatrixB2;

         SQMATRICES[0]->Reset();
         SQMATRICES[1]->Reset();
         MATRICES[0]->Reset();
         MATRICES[1]->Reset();

         N_SquareMatrices = 2;
         N_RectMatrices = 2;

         break;

        case 2:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[0] = SqmatrixM11;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 2;
          N_RectMatrices = 4;
        break;

        case 3:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          SQMATRICES[4] = SqmatrixM11;
          SQMATRICES[5] = SqmatrixM22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 6;
          N_RectMatrices = 2;
  
#ifdef __PRIVATE__  
        if(Disctype == VMS_PROJECTION)
        {
          N_SquareMatrices = 7;
          SQMATRICES[6] =  MatricesL;
          SQMATRICES[6]->Reset();

          N_RectMatrices = 6;
          MATRICES[2] = Matrices_tilde_G11;
          MATRICES[3] = Matrices_tilde_G22;
          MATRICES[4] = Matrices_G11;
          MATRICES[5] = Matrices_G22;
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();

          N_FESpaces = 4;
        }
#endif    
  
        break;

        case 4:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          SQMATRICES[4] = SqmatrixM11;
          SQMATRICES[5] = SqmatrixM22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 6;
          N_RectMatrices = 4;

#ifdef __PRIVATE__  
        if(Disctype == VMS_PROJECTION)
        {
          N_SquareMatrices = 7;
          SQMATRICES[6] =  MatricesL;
          SQMATRICES[6]->Reset();

          N_RectMatrices = 8;
          MATRICES[4] = Matrices_tilde_G11;
          MATRICES[5] = Matrices_tilde_G22;
          MATRICES[6] = Matrices_G11;
          MATRICES[7] = Matrices_G22;
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();
          MATRICES[6]->Reset();
          MATRICES[7]->Reset();

          N_FESpaces = 4;
        }
#endif        
          break;
      } //  switch(NSEType)
     
      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;
      RHSs[2] = rhs + 2*N_U;
      memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
     
      fesprhs[0] = FeSpaces[0];
      fesprhs[1] = FeSpaces[0];
      fesprhs[2] = FeSpaces[1];
  
      // assemble
      Assemble2D(N_FESpaces, FeSpaces,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormARhs,
        BoundaryConditions,
        BoundaryValues,
        NSEaux);

       
      if( (Disctype==UPWIND) && (!TDatabase::ParamDB->STOKES_PROBLEM) )
       {
        switch(NSEType)
         {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << endl;
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            UpwindForNavierStokes(LinCoeffs[0],SQMATRICES[3], FeFct[0], FeFct[1]);
           break;
         }                        // endswitch
       }                          // endif     
            
      // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {
        if(NSEType <4)
         {
          OutPut("For slip with friction bc NSTYPE 4 is ");
          OutPut("necessary !!!!! " << endl);
          exit(4711);
         }
 
        N_FESpaces = 1;
        N_SquareMatrices = 8;
        N_RectMatrices = 2;
        N_Rhs = 2;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;
        SQMATRICES[2] = SqmatrixA12;
        SQMATRICES[3] = SqmatrixA21;
        SQMATRICES[4] = SqmatrixM11;
        SQMATRICES[5] = SqmatrixM22;
        SQMATRICES[6] = SqmatrixM12;
        SQMATRICES[7] = SqmatrixM21;

        MATRICES[0] = MatrixB1T;
        MATRICES[1] = MatrixB2T;


        Assemble2DSlipBC(N_FESpaces, FeSpaces,
                         N_SquareMatrices, SQMATRICES,
                         N_RectMatrices, MATRICES,
                         N_Rhs, RHSs, fesprhs,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         NSEaux,
                         FeFct[0], FeFct[1]);

      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=      
    
     // set rhs for Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
     
#ifdef __PRIVATE__   
      // update matrices
      if (Disctype == VMS_PROJECTION)
        {
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixA12;
         SQMATRICES[2] = SqmatrixA21;
         SQMATRICES[3] = SqmatrixA22;
         SQMATRICES[6] =  MatricesL;
         MATRICES[2] = Matrices_tilde_G11;
         MATRICES[3] = Matrices_tilde_G22;
         MATRICES[4] = Matrices_G11;
         MATRICES[5] = Matrices_G22;

         VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[3]->GetN_DegreesOfFreedom(),
                                     SQMATRICES, MATRICES);
        }
#endif  
     
} // TSystemMatTNSE2D::Assemble(T

 
void TSystemMatTNSE2D::AssembleRhs(double *sol, double *rhs)
{
  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  
  double *RHSs[3];
  
  TFESpace2D *fesprhs[3];
  
  
      N_Rhs = 2;
      N_FESpaces = 1;
      N_SquareMatrices = 0;
      N_RectMatrices = 0;
  
      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;

      memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
     
      fesprhs[0] = FeSpaces[0];
      fesprhs[1] = FeSpaces[0];
      fesprhs[2] = FeSpaces[1];  
  
      // assemble
      Assemble2D(N_FESpaces, FeSpaces,
        N_SquareMatrices, NULL,
        N_RectMatrices, NULL,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormRhs,
        BoundaryConditions,
        BoundaryValues,
        NSE_Rhsaux);      
}



void TSystemMatTNSE2D::AssembleSystMat(double scale, double *oldrhs, double *rhs, double *sol)
{
 double tau, val = TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT;
  
  tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
  memset(B, 0, (2*N_U+N_P)*SizeOfDouble);    
     
  // old rhs multiplied with current subtime step and theta3 on B
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs, B);
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, oldrhs+N_U, B+N_U);   
    
  // add rhs from current sub time step to rhs array B
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs, B);
  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs+N_U, B+N_U);   
 
  // scale the pressure matrices, not in nonlinear step 
  if (scale != 1.0)
  {
   switch(NSEType)
    {
     case 1:
     case 3:
        Dscal(MatrixB1->GetN_Entries(), scale, MatrixB1->GetEntries());
        Dscal(MatrixB2->GetN_Entries(), scale, MatrixB2->GetEntries());
     break;

    case 2:
    case 4:     
         Dscal(MatrixB1T->GetN_Entries(), scale, MatrixB1T->GetEntries());
         Dscal(MatrixB2T->GetN_Entries(), scale, MatrixB2T->GetEntries());
	 
      // scale divergence constraint
      if(val>0) 
       {
        Dscal(MatrixB1->GetN_Entries(), val*scale, MatrixB1->GetEntries());
        Dscal(MatrixB2->GetN_Entries(), val*scale, MatrixB2->GetEntries());
       }
      break;
    } // switch(NSETyp
  } //  if (scale != 1.0)
    
   // Also currently : M := M + gamma A
   // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A 
   // defect = M * sol
   // B:= B + defect 
   memset(defect, 0, (2*N_U+N_P)*SizeOfDouble);  
   switch(NSEType)
    {
     case 1:
     case 2:
       MatAdd(SqmatrixM11, SqmatrixA11, -tau*TDatabase::TimeDB->THETA2);          
       gamma = - tau*TDatabase::TimeDB->THETA2;
   
       MatVectActive(SqmatrixM11, sol, defect);
       MatVectActive(SqmatrixM11, sol+N_U, defect+N_U);
       Daxpy(N_Active, 1, defect, B);
       Daxpy(N_Active, 1, defect+N_U, B+N_U);
 
       // assembling of system matrix       
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma + tau*TDatabase::TimeDB->THETA1);   
       gamma = tau*TDatabase::TimeDB->THETA1;
     break;

     case 3:
     case 4:
       MatAdd(SqmatrixM11, SqmatrixA11, - tau*TDatabase::TimeDB->THETA2);
       MatAdd(SqmatrixM12, SqmatrixA12, - tau*TDatabase::TimeDB->THETA2);
       MatAdd(SqmatrixM21, SqmatrixA21, - tau*TDatabase::TimeDB->THETA2);
       MatAdd(SqmatrixM22, SqmatrixA22, - tau*TDatabase::TimeDB->THETA2);       
       gamma = - tau*TDatabase::TimeDB->THETA2;

       MatVectActive(SqmatrixM11, sol, defect);
       Daxpy(N_Active, 1, defect, B);
       MatVectActive(SqmatrixM12, sol+N_U, defect);
       Daxpy(N_Active, 1, defect, B);
       MatVectActive(SqmatrixM21, sol, defect+N_U);
       Daxpy(N_Active, 1, defect+N_U, B+N_U);
       MatVectActive(SqmatrixM22, sol+N_U, defect+N_U);
       Daxpy(N_Active, 1, defect+N_U, B+N_U);

       //assembling system matrix
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma + tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM12, SqmatrixA12, -gamma + tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM21, SqmatrixA21, -gamma + tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM22, SqmatrixA22, -gamma + tau*TDatabase::TimeDB->THETA1);       
       gamma = tau*TDatabase::TimeDB->THETA1;     
     break;
    } 
  
   // set rhs for Dirichlet nodes
   memcpy(B+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
   memcpy(B+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
//               cout<< "B " <<Ddot((2*N_U+N_P), B, B)<<endl; 
   SystMatAssembled  = TRUE;
} // AssembleSystMat


/* assemble only LHS, not rhs */
void TSystemMatTNSE2D::AssembleSystMatNonLinear()
{
 double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
     
  if(SystMatAssembled)
   {
    cout << "Restore System mat before calling AssembleSystMat" <<endl;
   }

   switch(NSEType)
    {
     case 1:
     case 2: 
       // assembling of system matrix       
       MatAdd(SqmatrixM11, SqmatrixA11,   tau*TDatabase::TimeDB->THETA1);   
       gamma = tau*TDatabase::TimeDB->THETA1;
     break;

     case 3:
     case 4:       
       //assembling system matrix
       MatAdd(SqmatrixM11, SqmatrixA11, tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM12, SqmatrixA12,  tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM21, SqmatrixA21,  tau*TDatabase::TimeDB->THETA1);
       MatAdd(SqmatrixM22, SqmatrixA22,  tau*TDatabase::TimeDB->THETA1);       
       gamma = tau*TDatabase::TimeDB->THETA1;     
     break;
    } 
  
   SystMatAssembled  = TRUE;
} // AssembleSystMatNonLinear




void TSystemMatTNSE2D::RestoreMassMat()
{

//   cout << "RestoreMassMat  gamma " << gamma << endl;
  if(SystMatAssembled)
   {
    // restore the mass matrix
    switch(NSEType)
     {
      case 1:
      case 2:
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma);          
       gamma = 0.;
      break;

     case 3:
     case 4:
       //assembling system matrix
       MatAdd(SqmatrixM11, SqmatrixA11, -gamma);
       MatAdd(SqmatrixM12, SqmatrixA12, -gamma);
       MatAdd(SqmatrixM21, SqmatrixA21, -gamma);
       MatAdd(SqmatrixM22, SqmatrixA22, -gamma);       
       gamma = 0.;     
     break;
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


void TSystemMatTNSE2D::AssembleANonLinear(double *sol, double *rhs)
{
 int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces, last_sq;

     N_RectMatrices = 0;          
     N_Rhs = 0;
     N_FESpaces = 1;
 
     // set the nonliner matrices
      switch(TDatabase::ParamDB->NSTYPE)
       {
        case 1:
        case 2:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[0]->Reset();

          N_SquareMatrices = 1;
        break;

        case 3:
        case 4:
          if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
           {
            SQMATRICES[0] = SqmatrixA11;
            SQMATRICES[1] = SqmatrixA22;
            SQMATRICES[0]->Reset();
            SQMATRICES[1]->Reset();

            N_SquareMatrices = 2;
            last_sq = 1;

#ifdef __PRIVATE__   
            if (Disctype == VMS_PROJECTION)
              {
               SQMATRICES[0] = SqmatrixA11;
               SQMATRICES[1] = SqmatrixA12;
               SQMATRICES[2] = SqmatrixA21;
               SQMATRICES[3] = SqmatrixA22;
               SQMATRICES[0]->Reset();
               SQMATRICES[1]->Reset();
               SQMATRICES[2]->Reset();
               SQMATRICES[3]->Reset();

               N_SquareMatrices = 4;
               last_sq = 3;

               N_RectMatrices = 2;
               MATRICES[0] = Matrices_tilde_G11;
               MATRICES[1] = Matrices_tilde_G22;
               MATRICES[0]->Reset();
               MATRICES[1]->Reset();
       
               N_FESpaces = 4;
              }  
#endif

           }
          else
           {
            // Newton method
            cout<< "Newton method not tested " <<endl;
            exit(0);
           }

         break;
        } // switch(TDatabase::ParamDB->NSTYPE)
         
    
      // assemble the nonlinear part of NSE
      Assemble2D(N_FESpaces, FeSpaces,
                 N_SquareMatrices, SQMATRICES,
                 N_RectMatrices, MATRICES,
                 N_Rhs, NULL, NULL,
                 DiscreteFormNL,
                 BoundaryConditions,
                 BoundaryValues,
                 NSEaux);    

       // apply upwind disc
      if( (Disctype==UPWIND) && (!TDatabase::ParamDB->STOKES_PROBLEM) )
       {
        switch(NSEType)
         {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << endl;
            UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
            UpwindForNavierStokes(LinCoeffs[0],SQMATRICES[last_sq], FeFct[0], FeFct[1]);
          break;
         }                        // endswitch
       }                          // endif     
       
       
      // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      { 
        N_FESpaces = 1;
        N_SquareMatrices = 4;
        N_RectMatrices = 0;
        N_Rhs = 0;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA12;
        SQMATRICES[2] = SqmatrixA21;
        SQMATRICES[3] = SqmatrixA22;

        Assemble2DSlipBC(N_FESpaces, FeSpaces,
                         N_SquareMatrices, SQMATRICES,
                         N_RectMatrices, NULL,
                         N_Rhs, NULL, NULL,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         NSEaux,
                         FeFct[0], FeFct[1]);

      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=         

     // set rhs for Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble);       
     
#ifdef __PRIVATE__   
      // update matrices
      if (Disctype == VMS_PROJECTION)
        {
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixA12;
         SQMATRICES[2] = SqmatrixA21;
         SQMATRICES[3] = SqmatrixA22;
         SQMATRICES[6] =  MatricesL;
         MATRICES[2] = Matrices_tilde_G11;
         MATRICES[3] = Matrices_tilde_G22;
         MATRICES[4] = Matrices_G11;
         MATRICES[5] = Matrices_G22;

         VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[3]->GetN_DegreesOfFreedom(),
                                     SQMATRICES, MATRICES);
        }
#endif       
     
} //TSystemMatTNSE2D::AssembleNonLinear(

 
void TSystemMatTNSE2D::Solve(double *sol)
{  
  if(!SystMatAssembled)
  {
    cout << "System Matrix is not assembled to solve " <<endl;
    exit(0);
  }
  
    switch(Solver)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:
        cout << "GMG solver not yet implemented " <<endl;
      break;

      case DIRECT:
        switch(NSEType)
         {
          case 1:
            DirectSolver(SqmatrixM11, MatrixB1,  MatrixB2, B, sol);
          break;

          case 2:
             DirectSolver(SqmatrixA11, MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, B, sol);
          break;

          case 3:
           cout << "Direct solver not yet implemented for NSTYPE 3 " <<endl;
          break;

          case 4:
             DirectSolver(SqmatrixM11, SqmatrixM12, SqmatrixM21, SqmatrixM22, 
                          MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, B, sol); 
          break;
      } //  switch(NSEType) 

      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    

}



void TSystemMatTNSE2D::GetTNSEResidual(double *sol, double *res)
{
  
  if(!SystMatAssembled)
   {
    cout << "System Matrix is not assembled to calculate residual " <<endl;
    exit(0);
   }
    
    SQMATRICES[0] = SqmatrixM11;
    SQMATRICES[1] = SqmatrixM12;
    SQMATRICES[2] = SqmatrixM21;
    SQMATRICES[3] = SqmatrixM22;
       
    MATRICES[0] = MatrixB1;
    MATRICES[1] = MatrixB2;
    MATRICES[2] = MatrixB1T;
    MATRICES[3] = MatrixB2T;
       
   Defect(sqmatrices, matrices, sol, B, res); 
   
} // TSystemMatTNSE2D::GetResidual


void TSystemMatTNSE2D::MeasureTNSEErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
                                    double *AllErrors)
{
  MultiIndex2D TimeNSAllDerivatives[3] = { D00, D10, D01 };

  double errors[4],  u_error[4];    
    
     // errors in first velocity component
     FeFct[0]->GetErrors(ExactU1, 3, TimeNSAllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces, errors);
      u_error[0] = errors[0];
      u_error[1] = errors[1];
      
     // errors in second velocity component
     FeFct[1]->GetErrors(ExactU2, 3, TimeNSAllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces, errors);
     u_error[2] = errors[0];
     u_error[3] = errors[1];      
      
      // errors in pressure
     FeFct[2]->GetErrors(ExactP, 3, TimeNSAllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces+1, errors);     

     
     // calculate all errors
     AllErrors[0] = sqrt(u_error[0]*u_error[0]+u_error[2]*u_error[2]);
     AllErrors[1] = sqrt(u_error[1]*u_error[1]+u_error[3]*u_error[3]);
     AllErrors[2] = errors[0];
     AllErrors[3] = errors[1];    
     
      // error in L^infty(0,t,L^2)
      if(AllErrors[0] > AllErrors[5])
       {
        AllErrors[5]  = AllErrors[0];
        AllErrors[4]  =  TDatabase::TimeDB->CURRENTTIME;
      }

      
      // error in L^2(0,t,L^2)    
      AllErrors[6] += (u_error[0]*u_error[0] + u_error[2]*u_error[2] +olderror_l_2_l_2u)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;      
      olderror_l_2_l_2u = u_error[0]*u_error[0] + u_error[2]*u_error[2];
     
}
    
