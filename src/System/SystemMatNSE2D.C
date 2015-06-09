/** ************************************************************************ 
* @brief     source file for TSystemMatNSE2D
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History 
* ************************************************************************  */
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

TSystemMatNSE2D::TSystemMatNSE2D(TFEVectFunct2D *velocity,
                                 TFEFunction2D *pressure)
    : SystemMat2D(2, 5, 4) // 2 spaces, 5 square, 4 rectangular matrices
                           // the last two are overwritten depending on NSTYPE
{
  //store the FEspaces and fefunct
  this->fe_spaces[0] = velocity->GetFESpace2D(); // velocity space
  this->fe_spaces[1] = pressure->GetFESpace2D(); // pressure space
  // save the boundary condition
  this->BoundaryConditions[0] = fe_spaces[0]->GetBoundCondition();
  this->BoundaryConditions[1] = fe_spaces[0]->GetBoundCondition();
  
  // build matrices
  // first build matrix structure
  TSquareStructure2D *sqstructureA = new TSquareStructure2D(fe_spaces[0]);
  sqstructureA->Sort();  // sort column numbers: numbers are in increasing order
      
  TStructure2D *structureB = new TStructure2D(fe_spaces[1], fe_spaces[0]);
  TStructure2D *structureBT = new TStructure2D(fe_spaces[0], fe_spaces[1]);
  
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      /*
       * ( A  0  B1^T )
       * ( 0  A  B2^T )
       * ( B1 B2 0    )
       * 
       * B1^T and B2^T are not explicitly stored.
       */
      this->sq_matrices = { new TSquareMatrix2D(sqstructureA) };
      this->rect_matrices = { new TMatrix2D(structureB),
                              new TMatrix2D(structureB) };
      this->defect = Defect_NSE1;
      break;
      
    case 2:
      /*
       * ( A  0  B1T )
       * ( 0  A  B2T )
       * ( B1 B2 0    )
       * 
       * B1T and B2T are explicitly stored.
       */
      this->sq_matrices = { new TSquareMatrix2D(sqstructureA)};
      this->rect_matrices = { new TMatrix2D(structureB),
                              new TMatrix2D(structureB),
                              new TMatrix2D(structureBT),
                              new TMatrix2D(structureBT) };
      this->defect = Defect_NSE2;
      break;
      
    case 3:
      /*
       * ( A11 A12 B1^T )
       * ( A21 A22 B2^T )
       * ( B1  B2  0    )
       * 
       * B1^T and B2^T are not explicitly stored.
       */
      this->sq_matrices = { new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(sqstructureA) };
      this->rect_matrices = { new TMatrix2D(structureB),
                              new TMatrix2D(structureB) };
      this->defect = Defect_NSE3;
      break;
      
    case 4:
      /*
       * ( A11 A12 B1T )
       * ( A21 A22 B2T )
       * ( B1  B2  0    )
       * 
       * B1T and B2T are explicitly stored.
       */
      this->sq_matrices = { new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(sqstructureA) };
      this->rect_matrices = { new TMatrix2D(structureB),
                              new TMatrix2D(structureB),
                              new TMatrix2D(structureBT),
                              new TMatrix2D(structureBT)};
      this->defect = Defect_NSE4;
      break;
      
    case 14:
      /*
       * ( A11 A12 B1^T )
       * ( A21 A22 B2^T )
       * ( B1  B2  C    )
       * 
       * B1^T and B2^T are not explicitly stored.
       */
      this->sq_matrices = { new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(sqstructureA),
                            new TSquareMatrix2D(new TSquareStructure2D(
                                                fe_spaces[1])) };
      this->rect_matrices = { new TMatrix2D(structureB),
                              new TMatrix2D(structureB),
                              new TMatrix2D(structureBT),
                              new TMatrix2D(structureBT) };
      // a method to compute the defect for NSTYPE 14 does not exist in 2D,
      // therefore we take the one for NSTYPE 4, which excludes C
      OutPut("WARNING: NSTYPE 14 is not fully supported, take NSTYPE 4\n");
      this->defect = Defect_NSE4;
      break;
    default:
      OutPut("Unknown NSETYPE, it must be 1 to 4" << endl);
      exit(4711);
    }
  }

TSystemMatNSE2D::~TSystemMatNSE2D()
{
  delete sq_matrices[0]->GetStructure(); // delete structure of all A matrices
  if(TDatabase::ParamDB->NSTYPE == 14)
    delete sq_matrices.back()->GetStructure(); // delete structure of matrix C
  delete rect_matrices[0]->GetStructure(); // delete structure of matrix B
  if(TDatabase::ParamDB->NSTYPE != 1 && TDatabase::ParamDB->NSTYPE != 3)
    delete rect_matrices.back()->GetStructure(); //delete structure of matrix BT
  // delete matrices
  for(auto mat : sq_matrices)
    delete mat;
  for(auto mat : rect_matrices)
    delete mat;
}

void TSystemMatNSE2D::Init(BoundValueFunct2D *U1BoundValue,
                           BoundValueFunct2D *U2BoundValue)
{
  // save the boundary values  
  this->BoundaryValues[0] = U1BoundValue;
  this->BoundaryValues[1] = U2BoundValue;
} // TSystemMatNSE2D::Init

void TSystemMatNSE2D::Assemble(LocalAssembling2D& la, double *sol, double *rhs)
{
  // reset the matrices to zero
  for(auto mat : sq_matrices)
    mat->Reset();
  for(auto mat : rect_matrices)
    mat->Reset();
  
  int N_Rhs = 2;
  int N_FESpaces = 2;
  int N_U = this->fe_spaces[0]->GetN_DegreesOfFreedom();
  int N_P = this->fe_spaces[1]->GetN_DegreesOfFreedom();
  int N_Active = this->fe_spaces[0]->GetN_ActiveDegrees();
  int N_DirichletDof = N_U - N_Active;
  
  double *RHSs[3] = {rhs, rhs + N_U, rhs + 2 * N_U};
  memset(rhs, 0, (2 * N_U + N_P) * SizeOfDouble);
  
  TFESpace2D *fesprhs[3] = {fe_spaces[0], fe_spaces[0], fe_spaces[1]};
  
  // assemble
  Assemble2D(N_FESpaces, &fe_spaces[0], sq_matrices.size(), &sq_matrices[0],
             rect_matrices.size(), &rect_matrices[0], N_Rhs, RHSs, fesprhs,
             BoundaryConditions, BoundaryValues, la);
  
  if((TDatabase::ParamDB->DISCTYPE == UPWIND) 
     && (!TDatabase::ParamDB->PROBLEM_TYPE == 3))
  {
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        // do upwinding with one matrix
        UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                              la.get_fe_function(0), la.get_fe_function(1));
        cout << "UPWINDING DONE : level " << endl;
        break;
      case 3:
      case 4:
      case 14:
        // do upwinding with two matrices
        cout << "UPWINDING DONE : level " << endl;
        UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                              la.get_fe_function(0), la.get_fe_function(1));
        UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[3],
                              la.get_fe_function(0), la.get_fe_function(1));
        break;
    } // endswitch
  } // endif
  
  // slip with boundary condition
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
  {
    if(TDatabase::ParamDB->NSTYPE < 4)
    {
      OutPut("For slip with friction bc NSTYPE 4 is ");
      OutPut("necessary !!!!! " << endl);
      exit(4711);
    }
    ErrMsg("Assemble2DSlipBC does not work");
    exit(1);
    
    /*
    N_FESpaces = 1;
    int N_SquareMatrices = 4;
    int N_RectMatrices = 2;
    N_Rhs = 2;
    
    TSquareMatrix2D *SQMATRICES[4];
    SQMATRICES[0] = sq_matrices[0];
    SQMATRICES[1] = sq_matrices[3];
    SQMATRICES[2] = sq_matrices[1];
    SQMATRICES[3] = sq_matrices[2];
    
    TMatrix2D *MATRICES[2];
    MATRICES[0] = rect_matrices[2];
    MATRICES[1] = rect_matrices[3];
    TAuxParam2D NSEaux
    
    Assemble2DSlipBC(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES,
    N_RectMatrices, MATRICES, N_Rhs, RHSs, fesprhs,
    NULL, BoundaryConditions, BoundaryValues, &NSEaux,
    la.get_fe_function(0), la.get_fe_function(1));
    */
  } // (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=      
  
  // set rhs for Dirichlet nodes
  memcpy(sol + N_Active, rhs + N_Active, N_DirichletDof * SizeOfDouble);
  memcpy(sol + N_U + N_Active, rhs + N_U + N_Active,
         N_DirichletDof * SizeOfDouble);
} // TSystemMatNSE2D::Assemble(...)

void TSystemMatNSE2D::AssembleNonLinear(LocalAssembling2D& la, double *sol,
                                        double *rhs)
{
  int N_U = this->fe_spaces[0]->GetN_DegreesOfFreedom();
  int N_Active = this->fe_spaces[0]->GetN_ActiveDegrees();
  int N_DirichletDof = N_U - N_Active;
  
  int N_SquareMatrices, last_sq;
  TSquareMatrix2D *SQMATRICES[2];
  
  // set the nonliner matrices
  // the matrix blocks to which the nonlinear term contributes are reset to zero
  // and then completely reassembled, including the linear and nonlinear terms.
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
    case 2:
      SQMATRICES[0] = sq_matrices[0];
      SQMATRICES[0]->Reset();
      
      N_SquareMatrices = 1;
      break;
      
    case 3:
    case 4:
    case 14:
      if(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE == 0)
      {
        SQMATRICES[0] = sq_matrices[0];
        SQMATRICES[1] = sq_matrices[3];
        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        
        N_SquareMatrices = 2;
        last_sq = 1;
      }
      else
      {
        // Newton method
        cout << "Newton method not tested " << endl;
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
  Assemble2D(N_FESpaces, &fe_spaces[0], N_SquareMatrices, SQMATRICES,
             N_RectMatrices, NULL, N_Rhs, NULL, NULL, BoundaryConditions,
             BoundaryValues, la);
  
  // apply upwind disc
  if((TDatabase::ParamDB->DISCTYPE == UPWIND) && (!TDatabase::ParamDB
      ->PROBLEM_TYPE
                                                  == 3))
  {
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        // do upwinding with one matrix
        UpwindForNavierStokes(la.GetCoeffFct(), SQMATRICES[0],
                              la.get_fe_function(0), la.get_fe_function(1));
        cout << "UPWINDING DONE : level " << endl;
        break;
        
      case 3:
      case 4:
      case 14:
        // do upwinding with two matrices
        cout << "UPWINDING DONE : level " << endl;
        UpwindForNavierStokes(la.GetCoeffFct(), SQMATRICES[0],
                              la.get_fe_function(0), la.get_fe_function(1));
        UpwindForNavierStokes(la.GetCoeffFct(), SQMATRICES[last_sq],
                              la.get_fe_function(0), la.get_fe_function(1));
        break;
    } // endswitch
  } // endif
  
  // slip with boundary condition
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
  {
    N_FESpaces = 1;
    N_SquareMatrices = 4;
    N_RectMatrices = 0;
    N_Rhs = 0;
    
    SQMATRICES[0] = sq_matrices[0];
    SQMATRICES[1] = sq_matrices[1];
    SQMATRICES[2] = sq_matrices[2];
    SQMATRICES[3] = sq_matrices[3];
    
    ErrMsg("Assemble2DSlipBC does not work");
    exit(1);
    // TAuxParam2D NSEaux;
    //Assemble2DSlipBC(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES,
    //                 N_RectMatrices, NULL, N_Rhs, NULL, NULL, NULL,
    //                 BoundaryConditions, BoundaryValues, &NSEaux, 
    //                 la.get_fe_function(0), la.get_fe_function(1));
  }    // (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1
  
  // set rhs for Dirichlet nodes
  memcpy(sol + N_Active, rhs + N_Active, N_DirichletDof * SizeOfDouble);
  memcpy(sol + N_U + N_Active, rhs + N_U + N_Active,
         N_DirichletDof * SizeOfDouble);
} //TSystemMatNSE2D::AssembleNonLinear(

void TSystemMatNSE2D::GetResidual(double *sol, double *rhs, double *res)
{
  this->defect((TSquareMatrix**) &sq_matrices[0], (TMatrix**) &rect_matrices[0],
               sol, rhs, res);
} // TSystemMatNSE2D::GetResidual

void TSystemMatNSE2D::Solve(double *sol, double *rhs)
{
  switch(TDatabase::ParamDB->SOLVER_TYPE)
  {
    case AMG_SOLVE:
      cout << "AMG_SOLVE not yet implemented " << endl;
      break;
      
    case GMG:
      cout << "GMG solver not yet implemented " << endl;
      break;
      
    case DIRECT:
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          DirectSolver(sq_matrices[0], rect_matrices[0], rect_matrices[1], rhs,
                       sol);
          break;
          
        case 2:
          DirectSolver(sq_matrices[0], rect_matrices[2], rect_matrices[3],
                       rect_matrices[0], rect_matrices[1], rhs, sol);
          break;
          
        case 3:
          ErrMsg("Solver not included for NSTYPE 3 in this version. " 
                 << "try NSTYPE 4");
          exit(1);
          break;
          
        case 4:
          DirectSolver(sq_matrices[0], sq_matrices[1], sq_matrices[2],
                       sq_matrices[3], rect_matrices[2], rect_matrices[3],
                       rect_matrices[0], rect_matrices[1], rhs, sol);
          break;
        case 14:
          OutPut("WARNING: NSTYPE 14 is not fully supported, take NSTYPE 4\n");
          DirectSolver(sq_matrices[0], sq_matrices[1], sq_matrices[2],
                       sq_matrices[3], rect_matrices[2], rect_matrices[3],
                       rect_matrices[0], rect_matrices[1], rhs, sol);
      } //  switch(TDatabase::ParamDB->NSTYPE)
      break;
    default:
      OutPut("Unknown Solver\n");
      exit(4711);
  }
}

void TSystemMatNSE2D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->SystemMat2D::sq_matrices[0]->GetN_Rows();
  // reset y
  memset(y, 0.0, n_total_rows*SizeOfDouble);
  this->apply_scaled_add(x, y, factor);
}

void TSystemMatNSE2D::apply_scaled_add(const double *x, double *y, 
                                       double factor) const
{
  if(factor == 0.0)
    // nothing needs to be done
    return;
  
  // number of velocity degrees of freedom
  unsigned int n_v = this->SystemMat2D::sq_matrices[0]->GetN_Rows();
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      this->SystemMat2D::sq_matrices[0]->multiply(  x,       y,       factor);
      this->SystemMat2D::rect_matrices[0]->transpose_multiply(
                                                    x+2*n_v, y,       factor);
      
      this->SystemMat2D::sq_matrices[0]->multiply(  x+n_v,   y+n_v,   factor);
      this->SystemMat2D::rect_matrices[1]->transpose_multiply(
                                                    x+2*n_v, y+n_v,   factor);
      
      this->SystemMat2D::rect_matrices[0]->multiply(x,       y+2*n_v, factor);
      this->SystemMat2D::rect_matrices[1]->multiply(x+n_v,   y+2*n_v, factor);
      break;
    case 2:
      this->SystemMat2D::sq_matrices[0]->multiply(  x,       y,       factor);
      this->SystemMat2D::rect_matrices[2]->multiply(x+2*n_v, y,       factor);
      
      this->SystemMat2D::sq_matrices[0]->multiply(  x+n_v,   y+n_v,   factor);
      this->SystemMat2D::rect_matrices[3]->multiply(x+2*n_v, y+n_v,   factor);
      
      this->SystemMat2D::rect_matrices[0]->multiply(x,       y+2*n_v, factor);
      this->SystemMat2D::rect_matrices[1]->multiply(x+n_v,   y+2*n_v, factor);
      break;
    case 3:
      ErrMsg("TSystemMatNSE2D::apply_scaled_add not yet implemented");
      throw("TSystemMatNSE2D::apply_scaled_add not yet implemented");
      this->SystemMat2D::sq_matrices[0]->multiply(  x,       y,       factor);
      this->SystemMat2D::sq_matrices[1]->multiply(  x+n_v,   y,       factor);
      this->SystemMat2D::rect_matrices[0]->transpose_multiply(
                                                    x+2*n_v, y,       factor);
      
      this->SystemMat2D::sq_matrices[2]->multiply(  x,       y+n_v,   factor);
      this->SystemMat2D::sq_matrices[3]->multiply(  x+n_v,   y+n_v,   factor);
      this->SystemMat2D::rect_matrices[1]->transpose_multiply(
                                                    x+2*n_v, y+n_v,   factor);
      
      this->SystemMat2D::rect_matrices[0]->multiply(x,       y+2*n_v, factor);
      this->SystemMat2D::rect_matrices[1]->multiply(x+n_v,   y+2*n_v, factor);
      break;
    case 4:
    case 14:
      this->SystemMat2D::sq_matrices[0]->multiply(  x,       y,       factor);
      this->SystemMat2D::sq_matrices[1]->multiply(  x+n_v,   y,       factor);
      this->SystemMat2D::rect_matrices[2]->multiply(x+2*n_v, y,       factor);
      
      this->SystemMat2D::sq_matrices[2]->multiply(  x,       y+n_v,   factor);
      this->SystemMat2D::sq_matrices[3]->multiply(  x+n_v,   y+n_v,   factor);
      this->SystemMat2D::rect_matrices[3]->multiply(x+2*n_v, y+n_v,   factor);
      
      this->SystemMat2D::rect_matrices[0]->multiply(x,       y+2*n_v, factor);
      this->SystemMat2D::rect_matrices[1]->multiply(x+n_v,   y+2*n_v, factor);
      if(TDatabase::ParamDB->NSTYPE == 14)
      this->SystemMat2D::sq_matrices[4]->multiply(  x+2*n_v, y+2*n_v, factor);
      break;
    default:
      ErrMsg("Unknown NSTYPE, it must be 1 to 4, or 14");
      throw("unknown NSTYPE");
  }
}

unsigned int TSystemMatNSE2D::n_rows() const
{
  return 3;
}

unsigned int TSystemMatNSE2D::n_cols() const
{
  return 3;
}

unsigned int TSystemMatNSE2D::n_total_rows() const
{
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      return 2 * this->sq_matrices[0]->GetN_Rows()
          + this->rect_matrices[0]->GetN_Columns();
      break;
    case 2:
      return 2 * this->sq_matrices[0]->GetN_Rows()
          + this->rect_matrices[2]->GetN_Rows();
      break;
    case 3:
      return this->sq_matrices[0]->GetN_Rows()
          + this->sq_matrices[2]->GetN_Rows()
          + this->rect_matrices[0]->GetN_Columns();
      break;
    case 4:
    case 14:
      return this->sq_matrices[0]->GetN_Rows()
          + this->sq_matrices[2]->GetN_Rows()
          + this->rect_matrices[2]->GetN_Rows();
      break;
    default:
      ErrMsg("Unknown NSTYPE, it must be 1 to 4, or 14");
      throw("unknown NSTYPE");
  }
}

unsigned int TSystemMatNSE2D::n_total_cols() const
{
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
    case 2:
      return 2 * this->sq_matrices[0]->GetN_Columns()
          + this->rect_matrices[0]->GetN_Columns();
      break;
    case 3:
    case 4:
    case 14:
      return this->sq_matrices[0]->GetN_Columns()
          + this->sq_matrices[2]->GetN_Columns()
          + this->rect_matrices[0]->GetN_Columns();
      break;
    default:
      ErrMsg("Unknown NSTYPE, it must be 1 to 4, or 14");
      throw("unknown NSTYPE");
  }
}

