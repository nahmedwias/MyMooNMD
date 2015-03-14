/** ************************************************************************ 
* @brief     source file for TSystemMatNSE2D
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemMatNSE2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <Upwind.h>

#include <stdlib.h>
#include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>

TSystemMatNSE2D::TSystemMatNSE2D(TFEVectFunct2D *Velocity, TFEFunction2D *p,
                                 int disctype, int nsetype, int solver)
{
  //store the FEspaces and fefunct
  this->FeSpaces[0] = Velocity->GetFESpace2D(); // velocity space
  this->FeSpaces[1] = p->GetFESpace2D(); // pressure space
  
  this->N_U = FeSpaces[0]->GetN_DegreesOfFreedom();
  this->N_P = FeSpaces[1]->GetN_DegreesOfFreedom();
    
  this->N_Active =  FeSpaces[0]->GetActiveBound();
  this->N_DirichletDof = N_U - N_Active;  
  
  this->FeFct[0] = Velocity->GetComponent(0);
  this->FeFct[1] = Velocity->GetComponent(1); 
  this->FeFct[2] = p;
   
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
  sqstructureC = NULL; // only for NSTYPE 14
  
  structureB = new TStructure2D(FeSpaces[1], FeSpaces[0]);
  structureBT = new TStructure2D(FeSpaces[0], FeSpaces[1]);
    
  switch(NSEType)
  {
    case 1:
      MatrixB1 = new TMatrix2D(structureB);
      MatrixB2 = new TMatrix2D(structureB);
      MatrixB1T = NULL;
      MatrixB2T = NULL;

      SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
      SqmatrixA12 = NULL;
      SqmatrixA21 = NULL;
      SqmatrixA22 = NULL;
      sqmatrices[0] = SqmatrixA11;
      matrices[0] = MatrixB1;
      matrices[1] = MatrixB2;
      Defect = Defect_NSE1;
    break;

    case 2:
      MatrixB1 = new TMatrix2D(structureB);
      MatrixB2 = new TMatrix2D(structureB);
      MatrixB1T = new TMatrix2D(structureBT);
      MatrixB2T = new TMatrix2D(structureBT);

      SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
      SqmatrixA12 = NULL;
      SqmatrixA21 = NULL;
      SqmatrixA22 = NULL;
      sqmatrices[0] = SqmatrixA11;
      matrices[0] = MatrixB1;
      matrices[1] = MatrixB2;
      matrices[2] = MatrixB1T;
      matrices[3] = MatrixB2T;
      Defect = Defect_NSE2;
    break;

    case 3:
      MatrixB1 = new TMatrix2D(structureB);
      MatrixB2 = new TMatrix2D(structureB);
      MatrixB1T = NULL;
      MatrixB2T = NULL;

      SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
      SqmatrixA12 = new TSquareMatrix2D(sqstructureA);
      SqmatrixA21 = new TSquareMatrix2D(sqstructureA);
      SqmatrixA22 = new TSquareMatrix2D(sqstructureA);
      sqmatrices[0] = SqmatrixA11;
      sqmatrices[1] = SqmatrixA12;
      sqmatrices[2] = SqmatrixA21;
      sqmatrices[3] = SqmatrixA22;
      matrices[0] = MatrixB1;
      matrices[1] = MatrixB2;
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
      sqmatrices[0] = SqmatrixA11;
      sqmatrices[1] = SqmatrixA12;
      sqmatrices[2] = SqmatrixA21;
      sqmatrices[3] = SqmatrixA22;
      matrices[0] = MatrixB1;
      matrices[1] = MatrixB2;
      matrices[2] = MatrixB1T;
      matrices[3] = MatrixB2T;
      Defect = Defect_NSE4;
    break;
    
    case 14:
      MatrixB1 = new TMatrix2D(structureB);
      MatrixB2 = new TMatrix2D(structureB);
      MatrixB1T = new TMatrix2D(structureBT);
      MatrixB2T = new TMatrix2D(structureBT);

      SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
      SqmatrixA12 = new TSquareMatrix2D(sqstructureA);
      SqmatrixA21 = new TSquareMatrix2D(sqstructureA);
      SqmatrixA22 = new TSquareMatrix2D(sqstructureA);
      sqmatrices[0] = SqmatrixA11;
      sqmatrices[1] = SqmatrixA12;
      sqmatrices[2] = SqmatrixA21;
      sqmatrices[3] = SqmatrixA22;
      //sqmatrices[4] = new TSquareMatrix2D(this->N_P); // empty matrix
      sqstructureC = new TSquareStructure2D(FeSpaces[1]);
      sqmatrices[4] = new TSquareMatrix2D(sqstructureC);
      matrices[0] = MatrixB1;
      matrices[1] = MatrixB2;
      matrices[2] = MatrixB1T;
      matrices[3] = MatrixB2T;
      // a method to compute the defect for NSTYPE 14 does not exist in 2D,
      // therefore we take the one for NSTYPE 4, which excludes C
      OutPut("WARNING: NSTYPE 14 is not fully supported, take NSTYPE 4\n");
      Defect = Defect_NSE4;
      break;
    default:
      OutPut("Unknown NSETYPE, it must be 1 to 4" << endl);
      exit(4711);
  }  
 
  NSEaux_error = NULL;
  NSEaux = NULL;
}

TSystemMatNSE2D::~TSystemMatNSE2D()
{
  delete NSEaux; 
  
  if(NSEaux_error)
    delete NSEaux_error;
  delete SqmatrixA11;
  delete SqmatrixA12;
  delete SqmatrixA21;
  delete SqmatrixA22;
  delete MatrixB1;
  delete MatrixB2;
  delete MatrixB1T;
  delete MatrixB2T;
  delete sqstructureA;
  delete structureB;
  delete structureBT;
}


void TSystemMatNSE2D::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond,
                           BoundValueFunct2D *U1BoundValue,
                           BoundValueFunct2D *U2BoundValue)
{
  // save the boundary condition
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  

  // save the boundary values  
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
 
  // save the nse bilinear coefficient   
  LinCoeffs[0] = lincoeffs;
} // TSystemMatNSE2D::Init

 
void TSystemMatNSE2D::Assemble(LocalAssembling2D& la, double *sol, double *rhs)
{
  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
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
    case 14:
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
      if(NSEType == 14)
        OutPut("WARNING: NSTYPE 14 is not fully supported, take NSTYPE 4\n");
      break;
  } //  switch(NSEType)
  
  N_Rhs = 2;
  N_FESpaces = 2;
  
  double *RHSs[3] = {rhs, rhs+N_U, rhs+2*N_U };
  memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
  
  TFESpace2D *fesprhs[3] = { FeSpaces[0], FeSpaces[0], FeSpaces[1]};
  
  // assemble
//   Assemble2D(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES, 
//              N_RectMatrices, MATRICES, N_Rhs, RHSs, fesprhs, DiscreteFormARhs, 
//              BoundaryConditions, BoundaryValues, NSEaux);
  Assemble2D(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES, 
             N_RectMatrices, MATRICES, N_Rhs, RHSs, fesprhs, 
             BoundaryConditions, BoundaryValues, la);
  
  if( (Disctype==UPWIND) && (!TDatabase::ParamDB->PROBLEM_TYPE==3) )
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
      case 14:
        // do upwinding with two matrices
        cout << "UPWINDING DONE : level " << endl;
        UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
        UpwindForNavierStokes(LinCoeffs[0],SQMATRICES[3], FeFct[0], FeFct[1]);
        break;
    } // endswitch
  } // endif
  
  // slip with boundary condition
  if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
  {
    if(NSEType <4)
    {
      OutPut("For slip with friction bc NSTYPE 4 is ");
      OutPut("necessary !!!!! " << endl);
      exit(4711);
    }
    ErrMsg("Assemble2DSlipBC does not work");
    exit(1);
 
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

    //Assemble2DSlipBC(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES,
    //                 N_RectMatrices, MATRICES, N_Rhs, RHSs, fesprhs,
    //                 NULL, BoundaryConditions, BoundaryValues, NSEaux,
    //                 FeFct[0], FeFct[1]);

  } // (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=      
  
  // set rhs for Dirichlet nodes
  memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
  memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
} // TSystemMatNSE2D::Assemble(...)

void TSystemMatNSE2D::AssembleNonLinear(LocalAssembling2D& la, double *sol,
                                        double *rhs)
{
  int N_SquareMatrices, last_sq;
  
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
    case 14:
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
      if(TDatabase::ParamDB->NSTYPE == 14)
        OutPut("WARNING: NSTYPE 14 is not fully supported, take NSTYPE 4\n");
      break;
  } // switch(TDatabase::ParamDB->NSTYPE)
  
  int N_RectMatrices = 0;
  int N_Rhs = 0;
  int N_FESpaces = 1;
  
  // assemble the nonlinear part of NSE
  //Assemble2D(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES,
  //           N_RectMatrices, NULL, N_Rhs, NULL, NULL, DiscreteFormNL,
  //           BoundaryConditions, BoundaryValues, NSEaux);
  Assemble2D(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES,
             N_RectMatrices, NULL, N_Rhs, NULL, NULL,
             BoundaryConditions, BoundaryValues, la);
  
  // apply upwind disc
  if( (Disctype==UPWIND) && (!TDatabase::ParamDB->PROBLEM_TYPE==3) )
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
      case 14:
        // do upwinding with two matrices
        cout << "UPWINDING DONE : level " << endl;
        UpwindForNavierStokes(LinCoeffs[0], SQMATRICES[0], FeFct[0], FeFct[1]);
        UpwindForNavierStokes(LinCoeffs[0],SQMATRICES[last_sq], FeFct[0], FeFct[1]);
        break;
    } // endswitch
  } // endif     
  
  
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

    ErrMsg("Assemble2DSlipBC does not work");
    exit(1);
    //Assemble2DSlipBC(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES,
    //                 N_RectMatrices, NULL, N_Rhs, NULL, NULL, NULL,
    //                 BoundaryConditions, BoundaryValues, NSEaux, 
    //                 FeFct[0], FeFct[1]);
  }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1
   
  // set rhs for Dirichlet nodes
  memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
  memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
} //TSystemMatNSE2D::AssembleNonLinear(


void TSystemMatNSE2D::GetResidual(double *sol, double *rhs, double *res)
{
  this->Defect(sqmatrices, matrices, sol, rhs, res);
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
          DirectSolver(SqmatrixA11, MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2,
                       rhs, sol);
          break;

        case 3:
          ErrMsg("Solver not included for NSTYPE 3 in this version. "
                 << "try NSTYPE 4");
          exit(1);
          break;

        case 4:
          DirectSolver(SqmatrixA11, SqmatrixA12, SqmatrixA21, SqmatrixA22,
                       MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, rhs, sol); 
          break;
        case 14:
          OutPut("WARNING: NSTYPE 14 is not fully supported, take NSTYPE 4\n");
          DirectSolver(SqmatrixA11, SqmatrixA12, SqmatrixA21, SqmatrixA22,
                       MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, rhs, sol); 
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
  TAuxParam2D NSEaux_error;
  MultiIndex2D NSAllDerivatives[3] = { D00, D10, D01 };
  
  // errors in first velocity component
  FeFct[0]->GetErrors(ExactU1, 3, NSAllDerivatives, 2, L2H1Errors, NULL,
                      &NSEaux_error, 1, FeSpaces, errors);
  u_error[0] = errors[0];
  u_error[1] = errors[1];
  
  // errors in second velocity component
  FeFct[1]->GetErrors(ExactU2, 3, NSAllDerivatives, 2, L2H1Errors, NULL,
                      &NSEaux_error, 1, FeSpaces, errors);
  u_error[2] = errors[0];
  u_error[3] = errors[1];      
  
  // errors in pressure
  FeFct[2]->GetErrors(ExactP, 3, NSAllDerivatives, 2, L2H1Errors, NULL,
                      &NSEaux_error, 1, FeSpaces+1, errors);     
  p_error[0] = errors[0];
  p_error[1] = errors[1];        
}
    
