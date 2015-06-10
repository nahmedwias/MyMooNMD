/** ************************************************************************ 
* @brief     source file for TSystemMatNSE2D
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History 
 ************************************************************************  */

#ifdef __2D__

#include <Database.h>
#include <SystemMatNSE2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <NSE2D_ParamRout.h>
#include <MainUtilities.h>
#include <Upwind.h>

#include <stdlib.h>
#include <string.h>


TSystemMatNSE2D::TSystemMatNSE2D(TFESpace2D *velocity_fespace, TFESpace2D *presssure_fespace, TFEVectFunct2D *Velocity, 
                                 TFEFunction2D *p, int disctype, int nsetype, int solver)
{
  //store the FEspaces and fefunct
  FeSpaces[0] = velocity_fespace;
  FeSpaces[1] = presssure_fespace;
  
  N_U = velocity_fespace->GetN_DegreesOfFreedom();
  N_P = presssure_fespace->GetN_DegreesOfFreedom();
    
  N_Active =  velocity_fespace->GetActiveBound();
  N_DirichletDof = N_U - N_Active;  
  
  FeFct[0] = Velocity->GetComponent(0);
  FeFct[1] = Velocity->GetComponent(1); 
  FeFct[2] = p;
   
  //set the discretization type
  Disctype = disctype;
  
  // NSE type
  NSEType = nsetype;
  
  //set the solver type
  Solver = solver;
  
  // build matrices
  // first build matrix structure
  sqstructureA = new TSquareStructure2D(FeSpaces[0]);
  sqstructureA->Sort();  // sort column numbers: numbers are in increasing order
  
  structureB = new TStructure2D(FeSpaces[1], FeSpaces[0]);
  structureBT = new TStructure2D(FeSpaces[0], FeSpaces[1]);
    
    switch(NSEType)
     {
      case 1:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        Defect = Defect_NSE1;
      break;

      case 2:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);
        MatrixB1T = new TMatrix2D(structureBT);
        MatrixB2T = new TMatrix2D(structureBT);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        Defect = Defect_NSE2;
      break;

      case 3:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA22 = new TSquareMatrix2D(sqstructureA);
        Defect = Defect_NSE3;
      break;

      case 4:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);
        MatrixB1T = new TMatrix2D(structureBT);
        MatrixB2T = new TMatrix2D(structureBT);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA22 = new TSquareMatrix2D(sqstructureA);
        Defect = Defect_NSE4;
      break;
      
      default:
            OutPut("Unknown NSETYPE, it must be 1 to 4" << endl);
            exit(4711);;      
      
     }  
 
   // matrices for methods
   sqmatrices = (TSquareMatrix **)SQMATRICES;
   matrices = (TMatrix **)MATRICES;
   
   NSEaux_error = NULL;
   NSEaux = NULL;
}

TSystemMatNSE2D::~TSystemMatNSE2D()
{
    delete NSEaux; 
       
    if(NSEaux_error)
      delete NSEaux_error;

}


void TSystemMatNSE2D::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue, TAuxParam2D *aux)
{
  TDiscreteForm2D *DiscreteFormGalerkin, *DiscreteFormSDFEM, *DiscreteFormUpwind, *DiscreteFormSmagorinsky;
  TDiscreteForm2D *DiscreteFormVMSProjection, *DiscreteFormNLGalerkin, *DiscreteFormNLSDFEM, *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky, *DiscreteFormNLVMSProjection, *DiscreteFormPressSep, *DiscreteFormAuxProbPressSep;
  TDiscreteForm2D *DiscreteFormNSRFBRhs;
    
  // save the boundary condition
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  

  // save the boundary values  
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
 
  // save the nse bilinear coefficient   
  LinCoeffs[0] = lincoeffs;
  
  //default, i.e., velocity for nonlinear term
  if(aux==NULL)
   {    
    NSEaux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo,
                           NSN_FEValuesVelo,
                           FeSpaces, FeFct,
                           NSFctVelo,
                           NSFEFctIndexVelo, NSFEMultiIndexVelo,
                           NSN_ParamsVelo, NSBeginParamVelo);       
   }
  else
  {
   NSEaux = aux;
  }
  
  // aux for calculating the error
   if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
       NSEaux_error =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo,
                             NSN_ParamFctVelo,
                             NSN_FEValuesVelo,
                             FeSpaces, FeFct,
                             NSFctVelo,
                             NSFEFctIndexVelo, NSFEMultiIndexVelo,
                             NSN_ParamsVelo, NSBeginParamVelo);      
       
    }
  
  
  // set the Discreteforms
  InitializeDiscreteForms(
    DiscreteFormGalerkin, DiscreteFormSDFEM,
    DiscreteFormUpwind, DiscreteFormSmagorinsky,
    DiscreteFormVMSProjection,
    DiscreteFormNLGalerkin, DiscreteFormNLSDFEM,
    DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky,
    DiscreteFormNLVMSProjection,
    DiscreteFormPressSep,
    DiscreteFormAuxProbPressSep,
    DiscreteFormNSRFBRhs,
    LinCoeffs[0], NSEType);

    // find discrete form
    switch(Disctype)
       {
          case GALERKIN:
            DiscreteFormARhs = DiscreteFormGalerkin;
            DiscreteFormNL = DiscreteFormNLGalerkin;
          break;

          case SDFEM:
            DiscreteFormARhs = DiscreteFormSDFEM;
            DiscreteFormNL = DiscreteFormNLSDFEM; 
          break;

          case UPWIND:
            DiscreteFormARhs = DiscreteFormUpwind;
            DiscreteFormNL = DiscreteFormNLUpwind;    
            break;

          case SMAGORINSKY:
            DiscreteFormARhs = DiscreteFormSmagorinsky;
            DiscreteFormNL = DiscreteFormNLSmagorinsky;              
            break;

          case VMS_PROJECTION:
	      DiscreteFormARhs = DiscreteFormVMSProjection;
	      if (TDatabase::ParamDB->NSTYPE != 1)
	       {
                OutPut("VMS only for NSTYPE 1 implemented !!!"<<endl);
		exit(4711);
	       }

            DiscreteFormNL = DiscreteFormNLVMSProjection;
            break;

          default:
            Error("Unknown DISCTYPE" << endl);
            exit(-1);
        } 
     
     // set the discrete form for the Stokes equation
      if (TDatabase::ParamDB->PROBLEM_TYPE==STOKES)
       {
        DiscreteFormARhs = DiscreteFormUpwind;     
        DiscreteFormNL = NULL;
       }
} // TSystemMatNSE2D::Init

 
void TSystemMatNSE2D::Assemble(double *sol, double *rhs)
{
  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  
  double *RHSs[3];

  TFESpace2D *fesprhs[3];

     // initialize matrices
     switch(NSEType)
      {
        case 1:
          SQMATRICES[0] = SqmatrixA11;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 2;
        break;

        case 2:
          SQMATRICES[0] = SqmatrixA11;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 4;
        break;

        case 3:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 4;
          N_RectMatrices = 2;
        break;

        case 4:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 4;
          N_RectMatrices = 4;

          break;
      } //  switch(NSEType)
      
      N_Rhs = 2;
      N_FESpaces = 2;   
     
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
 
     
      if( (Disctype==UPWIND) && (!TDatabase::ParamDB->PROBLEM_TYPE == 3) )
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
        N_SquareMatrices = 4;
        N_RectMatrices = 2;
        N_Rhs = 2;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;
        SQMATRICES[2] = SqmatrixA12;
        SQMATRICES[3] = SqmatrixA21;

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
      
//     cout << "Test Assemble " << endl; 
} // TSystemMatNSE2D::Assemble(T

void TSystemMatNSE2D::AssembleNonLinear(double *sol, double *rhs)
{
 int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces, last_sq;

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
           }
          else
           {
            // Newton method
            cout<< "Newton method not tested " <<endl;
            exit(0);
           }

         break;
        } // switch(TDatabase::ParamDB->NSTYPE)
            
      N_RectMatrices = 0;          
      N_Rhs = 0;
      N_FESpaces = 1;
    
      // assemble the nonlinear part of NSE
      Assemble2D(N_FESpaces, FeSpaces,
                 N_SquareMatrices, SQMATRICES,
                 N_RectMatrices, NULL,
                 N_Rhs, NULL, NULL,
                 DiscreteFormNL,
                 BoundaryConditions,
                 BoundaryValues,
                 NSEaux);    

       // apply upwind disc
      if( (Disctype==UPWIND) && (!TDatabase::ParamDB->PROBLEM_TYPE == 3) )
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
      
      
      
} //TSystemMatNSE2D::AssembleNonLinear(


void TSystemMatNSE2D::GetResidual(double *sol, double *rhs, double *res)
{
  
   Defect(sqmatrices, matrices, sol, rhs, res); 
   
} // TSystemMatNSE2D::GetResidual

void TSystemMatNSE2D::Solve(double *sol, double *rhs)
{
  
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
            DirectSolver(SqmatrixA11, MatrixB1,  MatrixB2, rhs, sol);
          break;

          case 2:
             DirectSolver(SqmatrixA11, MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, rhs, sol);
          break;

          case 3:
           cout << "Solver not included for NSTYPE 3 in this version" <<endl;
            cout << "try NSTYPE 4 " <<endl;   
	    exit(0);
          break;

          case 4:
             DirectSolver(SqmatrixA11, SqmatrixA12, SqmatrixA21, SqmatrixA22, 
                          MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, rhs, sol); 
          break;
      } //  switch(NSEType) 

      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
  
}

void TSystemMatNSE2D::MeasureErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
                                    double *u_error, double *p_error)
{
  double errors[4];
  
     // errors in first velocity component
     FeFct[0]->GetErrors(ExactU1, 3, NSAllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces, errors);
      u_error[0] = errors[0];
      u_error[1] = errors[1];
      
     // errors in second velocity component
     FeFct[1]->GetErrors(ExactU2, 3, NSAllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces, errors);
     u_error[2] = errors[0];
     u_error[3] = errors[1];      
      
      // errors in pressure
     FeFct[2]->GetErrors(ExactP, 3, NSAllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces+1, errors);     
     p_error[0] = errors[0];
     p_error[1] = errors[1];        
}
    
    
#endif // #ifdef __3D__
