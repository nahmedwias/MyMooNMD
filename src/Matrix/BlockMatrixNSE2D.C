/** ************************************************************************ 
* @brief     source file for BlockMatrixNSE2D
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History 
* ************************************************************************  */

#include <Database.h>
#include <BlockMatrixNSE2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <Upwind.h>

#include <stdlib.h>
#include <string.h>

/** ************************************************************************ */
BlockMatrixNSE2D::BlockMatrixNSE2D(const TFESpace2D& velocity,
                                   const TFESpace2D& pressure,
                                   BoundValueFunct2D * const * const BoundValue)
    : BlockMatrix(Problem_type::NavierStokes, 2),
      boundary_values({{BoundValue[0], BoundValue[1], BoundValue[2]}})
{
  // build matrices
  // first build matrix structure
  TSquareStructure2D *sqstructureA = new TSquareStructure2D(&velocity);
  sqstructureA->Sort(); // sort column numbers: numbers are in increasing order
      
  TStructure2D *structureB = new TStructure2D(&pressure, &velocity);
  
  // number of pressure degrees of freedom
  unsigned int n_p = pressure.GetN_DegreesOfFreedom();
  
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
    {
      /*
       * ( A  0  B1^T )
       * ( 0  A  B2^T )
       * ( B1 B2 0    )
       * 
       * B1^T and B2^T are not explicitly stored.
       */
      unsigned int n_v = velocity.GetN_DegreesOfFreedom();
      
      this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(n_v));
      //this->BlockMatrix::blocks[2] stays a nullptr
      this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(n_v));
      this->BlockMatrix::blocks[4] = this->BlockMatrix::blocks[0];
      //this->BlockMatrix::blocks[5] stays a nullptr
      this->BlockMatrix::blocks[6].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[7].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[8].reset(new TSquareMatrix2D(n_p));
      break;
    }
    case 2:
    {
      /*
       * ( A  0  B1T )
       * ( 0  A  B2T )
       * ( B1 B2 0    )
       * 
       * B1T and B2T are explicitly stored.
       */
      unsigned int n_v = velocity.GetN_DegreesOfFreedom();
      TStructure2D *structureBT = new TStructure2D(&velocity, &pressure);
      
      this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(n_v));
      this->BlockMatrix::blocks[2].reset(new TMatrix2D(structureBT));
      this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(n_v));
      this->BlockMatrix::blocks[4] = this->BlockMatrix::blocks[0];
      this->BlockMatrix::blocks[5].reset(new TMatrix2D(structureBT));
      this->BlockMatrix::blocks[6].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[7].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[8].reset(new TSquareMatrix2D(n_p));
      break;
    }
    case 3:
    {
      /*
       * ( A11 A12 B1^T )
       * ( A21 A22 B2^T )
       * ( B1  B2  0    )
       * 
       * B1^T and B2^T are not explicitly stored.
       */
      this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(sqstructureA));
      //this->BlockMatrix::blocks[2] stays a nullptr
      this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[4].reset(new TSquareMatrix2D(sqstructureA));
      //this->BlockMatrix::blocks[5] stays a nullptr
      this->BlockMatrix::blocks[6].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[7].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[8].reset(new TSquareMatrix2D(n_p));
      break;
    }
    case 4:
    {
      /*
       * ( A11 A12 B1T )
       * ( A21 A22 B2T )
       * ( B1  B2  0   )
       * 
       * B1T and B2T are explicitly stored.
       */
      TStructure2D *structureBT = new TStructure2D(&velocity, &pressure);
      
      this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[2].reset(new TMatrix2D(structureBT));
      this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[4].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[5].reset(new TMatrix2D(structureBT));
      this->BlockMatrix::blocks[6].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[7].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[8].reset(new TSquareMatrix2D(n_p));
      break;
    }
    case 14:
    {
      /*
       * ( A11 A12 B1^T )
       * ( A21 A22 B2^T )
       * ( B1  B2  C    )
       * 
       * B1^T and B2^T are explicitly stored.
       */
      TStructure2D *structureBT = new TStructure2D(&velocity, &pressure);
      TSquareStructure2D *sqstructureC = new TSquareStructure2D(&pressure);
      // sort column numbers: numbers are in increasing order
      sqstructureC->Sort();
  
      this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[2].reset(new TMatrix2D(structureBT));
      this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[4].reset(new TSquareMatrix2D(sqstructureA));
      this->BlockMatrix::blocks[5].reset(new TMatrix2D(structureBT));
      this->BlockMatrix::blocks[6].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[7].reset(new TMatrix2D(structureB));
      this->BlockMatrix::blocks[8].reset(new TSquareMatrix2D(sqstructureC));
      break;
    }
    default:
      ErrThrow("Unknown NSETYPE, it must be 1, 2, 3, 4, or 14");
      break;
  }
}

/** ************************************************************************ */
BlockMatrixNSE2D::~BlockMatrixNSE2D()
{
  // delete structure of all A matrices
  delete this->BlockMatrix::blocks[0]->GetStructure(); 
  // delete structure of all BT matrices
  delete this->BlockMatrix::blocks[2]->GetStructure();
  // delete structure of all B matrices
  delete this->BlockMatrix::blocks[6]->GetStructure();
  // delete structure of matrix C
  delete this->BlockMatrix::blocks[8]->GetStructure();
}

/** ************************************************************************ */
void BlockMatrixNSE2D::Assemble(LocalAssembling2D& la, BlockVector& rhs)
{
  // this really needs to be rewritten, especially the Assemble2D function!
  int N_Rhs = 3;
  int N_FESpaces = 2;
  const TFESpace2D * v_space = this->get_velocity_space();
  const TFESpace2D * p_space = this->get_pressure_space();
  
  rhs.reset();
  double *RHSs[3] = {rhs.block(0), rhs.block(1), rhs.block(2)};
  
  // what follows is done only to call Assemble2D correctly
  const TFESpace2D *fespmat[2] = {v_space, p_space};
  const TFESpace2D *fesprhs[3] = {v_space, v_space, p_space};
  
  unsigned int n_sq_mat = 4;
  if(TDatabase::ParamDB->NSTYPE == 1 || TDatabase::ParamDB->NSTYPE == 2)
    n_sq_mat = 2;
  if(TDatabase::ParamDB->NSTYPE == 14)
  {
    OutPut("WARNING: NSTYPE 14 needs a special local assembling, the C-block "
         << "is now ignored\n");
    //n_sq_mat = 5;
  }
  TSquareMatrix2D * sq_matrices[5] = {this->get_A_block(0), nullptr, nullptr, 
                                      nullptr, nullptr};
  if(TDatabase::ParamDB->NSTYPE == 1 || TDatabase::ParamDB->NSTYPE == 2)
    sq_matrices[1] = this->get_A_block(3); // no more square blocks used
  else
  {
    sq_matrices[1] = this->get_A_block(1);
    sq_matrices[2] = this->get_A_block(2);
    sq_matrices[3] = this->get_A_block(3);
    sq_matrices[4] = this->get_C_block();
  }
  
  unsigned int n_rect_mat;
  TMatrix2D * rect_matrices[4] = {nullptr, nullptr, nullptr, nullptr};
  if(TDatabase::ParamDB->NSTYPE == 1 || TDatabase::ParamDB->NSTYPE == 3)
  {
    n_rect_mat = 2;
    rect_matrices[0] = this->get_B_block(0);
    rect_matrices[1] = this->get_B_block(1);
  }
  else
  {
    n_rect_mat = 4;
    rect_matrices[0] = this->get_B_block(0);
    rect_matrices[1] = this->get_B_block(1);
    rect_matrices[2] = this->get_BT_block(0);
    rect_matrices[3] = this->get_BT_block(1);
  }
  
  BoundCondFunct2D * boundary_conditions[3] = {
    v_space->GetBoundCondition(), v_space->GetBoundCondition(), 
    p_space->GetBoundCondition() };
  
  std::array<BoundValueFunct2D*, 3> non_const_bound_values;
  non_const_bound_values[0] = this->boundary_values[0];
  non_const_bound_values[1] = this->boundary_values[1];
  non_const_bound_values[2] = this->boundary_values[2];
    
  // assemble
  Assemble2D(N_FESpaces, fespmat, n_sq_mat, sq_matrices,
             n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
             boundary_conditions, non_const_bound_values.data(), la);
  
  if((TDatabase::ParamDB->DISCTYPE == UPWIND) 
     && !(TDatabase::ParamDB->PROBLEM_TYPE == 3))
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
  
  // remove ones in non-active rows in non-diagonal blocks
  this->get_A_block(1)->reset_non_active();
  this->get_A_block(2)->reset_non_active();
  this->get_BT_block(0)->reset_non_active();
  this->get_BT_block(1)->reset_non_active();
} // BlockMatrixNSE2D::Assemble(...)

/** ************************************************************************ */
void BlockMatrixNSE2D::AssembleNonLinear(LocalAssembling2D& la)
{
  const TFESpace2D * v_space = this->get_velocity_space();
  unsigned int n_sq_mat;
  TSquareMatrix2D *sq_mat[2] = {this->get_A_block(0), this->get_A_block(3)};
  // reset the matrices, because all terms are assembled again (including the 
  // linear ones)
  sq_mat[0]->reset();
  sq_mat[1]->reset();
  
  if(TDatabase::ParamDB->NSTYPE == 1 || TDatabase::ParamDB->NSTYPE == 2)
    n_sq_mat = 1;
  else
    n_sq_mat = 2;
  
  unsigned int n_rect_mat = 0;
  unsigned int n_rhs = 0;
  unsigned int n_fe_spaces = 1;
  BoundCondFunct2D * boundary_conditions[1] = { v_space->GetBoundCondition() };
  
  std::array<BoundValueFunct2D*, 3> non_const_bound_values;
  non_const_bound_values[0] = this->boundary_values[0];
  non_const_bound_values[1] = this->boundary_values[1];
  non_const_bound_values[2] = this->boundary_values[2];
  
  // assemble the nonlinear part of NSE
  Assemble2D(n_fe_spaces, &v_space, n_sq_mat, sq_mat, n_rect_mat, nullptr,
             n_rhs, nullptr, nullptr, boundary_conditions, 
             non_const_bound_values.data(), la);
  
  // apply upwind disc
  if((TDatabase::ParamDB->DISCTYPE == UPWIND) 
    && !(TDatabase::ParamDB->PROBLEM_TYPE == 3))
  {
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        // do upwinding with one matrix
        UpwindForNavierStokes(la.GetCoeffFct(), sq_mat[0],
                              la.get_fe_function(0), la.get_fe_function(1));
        cout << "UPWINDING DONE : level " << endl;
        break;
        
      case 3:
      case 4:
      case 14:
        // do upwinding with two matrices
        cout << "UPWINDING DONE : level " << endl;
        UpwindForNavierStokes(la.GetCoeffFct(), sq_mat[0],
                              la.get_fe_function(0), la.get_fe_function(1));
        UpwindForNavierStokes(la.GetCoeffFct(), sq_mat[1],
                              la.get_fe_function(0), la.get_fe_function(1));
        break;
    } // endswitch
  } // endif
  
  // slip with boundary condition
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
  {
    ErrThrow("Assemble2DSlipBC does not work");
    //n_fe_spaces = 1;
    //n_sq_mat = 4;
    //n_rect_mat = 0;
    //n_rhs = 0;
    //TSquareMatrix2D *sq_mat[4];
    //sq_mat[0] = this->get_A_block(0);
    //sq_mat[1] = this->get_A_block(1);
    //sq_mat[2] = this->get_A_block(2);
    //sq_mat[3] = this->get_A_block(3);
    //
    //TAuxParam2D NSEaux;
    //Assemble2DSlipBC(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES,
    //                 N_RectMatrices, NULL, N_Rhs, NULL, NULL, NULL,
    //                 BoundaryConditions, BoundaryValues, &NSEaux, 
    //                 la.get_fe_function(0), la.get_fe_function(1));
  }    // (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1
  
  if(this->combined_matrix)
  {
    // remove possibly defined combined matrix, because it is no longer up to
    // date, after nonlinear terms are assembled.
    delete this->combined_matrix->GetStructure();
    this->combined_matrix.reset();
  }
} //BlockMatrixNSE2D::AssembleNonLinear(

/** ************************************************************************ */
void BlockMatrixNSE2D::Solve(double *sol, double *rhs)
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
          DirectSolver(this->get_A_block(0), this->get_B_block(0),
                       this->get_B_block(1), rhs, sol);
          break;
          
        case 2:
          DirectSolver(this->get_A_block(0), this->get_BT_block(0), 
                       this->get_BT_block(1), this->get_B_block(0), 
                       this->get_B_block(1), rhs, sol);
          break;
          
        case 3:
          ErrThrow("Solver not included for NSTYPE 3 in this version. "
                   + "try NSTYPE 4");
          break;
          
        case 4:
          DirectSolver(this->get_A_block(0), this->get_A_block(1), 
                       this->get_A_block(2), this->get_A_block(3),
                       this->get_BT_block(0),  this->get_BT_block(1),
                       this->get_B_block(0), this->get_B_block(1), rhs, sol);
          break;
        case 14:
          OutPut("WARNING: NSTYPE 14 is not fully supported, take NSTYPE 4\n");
          DirectSolver(this->get_A_block(0), this->get_A_block(1), 
                       this->get_A_block(2), this->get_A_block(3),
                       this->get_BT_block(0), this->get_BT_block(1),
                       this->get_B_block(0), this->get_B_block(1), rhs, sol);
      } //  switch(TDatabase::ParamDB->NSTYPE)
      break;
    default:
      ErrThrow("Unknown Solver");
      break;
  }
}

/** ************************************************************************ */
void BlockMatrixNSE2D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->BlockMatrix::n_total_rows();
  // reset y
  memset(y, 0.0, n_total_rows*SizeOfDouble);
  this->apply_scaled_add(x, y, factor);
}

/** ************************************************************************ */
void BlockMatrixNSE2D::apply_scaled_add(const double *x, double *y, 
                                       double factor) const
{
  if(factor == 0.0)
    // nothing needs to be done
    return;
  
  // number of velocity degrees of freedom
  unsigned int n_v = this->get_velocity_space()->GetN_DegreesOfFreedom();
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      this->BlockMatrix::blocks[0]->multiply(x,       y,       factor);
      this->BlockMatrix::blocks[6]->transpose_multiply(
                                             x+2*n_v, y,       factor);
      
      this->BlockMatrix::blocks[0]->multiply(x+n_v,   y+n_v,   factor);
      this->BlockMatrix::blocks[7]->transpose_multiply(
                                             x+2*n_v, y+n_v,   factor);
      
      this->BlockMatrix::blocks[6]->multiply(x,       y+2*n_v, factor);
      this->BlockMatrix::blocks[7]->multiply(x+n_v,   y+2*n_v, factor);
      break;
    case 2:
      this->BlockMatrix::blocks[0]->multiply(  x,     y,       factor);
      this->BlockMatrix::blocks[2]->multiply(x+2*n_v, y,       factor);
      
      this->BlockMatrix::blocks[0]->multiply(  x+n_v, y+n_v,   factor);
      this->BlockMatrix::blocks[5]->multiply(x+2*n_v, y+n_v,   factor);
      
      this->BlockMatrix::blocks[6]->multiply(x,       y+2*n_v, factor);
      this->BlockMatrix::blocks[7]->multiply(x+n_v,   y+2*n_v, factor);
      break;
    case 3:
      this->BlockMatrix::blocks[0]->multiply(x,       y,       factor);
      this->BlockMatrix::blocks[1]->multiply(x+n_v,   y,       factor);
      this->BlockMatrix::blocks[6]->transpose_multiply(
                                             x+2*n_v, y,       factor);
      
      this->BlockMatrix::blocks[2]->multiply(x,       y+n_v,   factor);
      this->BlockMatrix::blocks[3]->multiply(x+n_v,   y+n_v,   factor);
      this->BlockMatrix::blocks[7]->transpose_multiply(
                                             x+2*n_v, y+n_v,   factor);
      
      this->BlockMatrix::blocks[6]->multiply(x,       y+2*n_v, factor);
      this->BlockMatrix::blocks[7]->multiply(x+n_v,   y+2*n_v, factor);
      break;
    case 4:
    case 14:
      this->BlockMatrix::blocks[0]->multiply(  x,     y,       factor);
      this->BlockMatrix::blocks[1]->multiply(  x+n_v, y,       factor);
      this->BlockMatrix::blocks[2]->multiply(x+2*n_v, y,       factor);
      
      this->BlockMatrix::blocks[3]->multiply(  x,     y+n_v,   factor);
      this->BlockMatrix::blocks[4]->multiply(  x+n_v, y+n_v,   factor);
      this->BlockMatrix::blocks[5]->multiply(x+2*n_v, y+n_v,   factor);
      
      this->BlockMatrix::blocks[6]->multiply(x,       y+2*n_v, factor);
      this->BlockMatrix::blocks[7]->multiply(x+n_v,   y+2*n_v, factor);
      if(TDatabase::ParamDB->NSTYPE == 14)
      this->BlockMatrix::blocks[8]->multiply(x+2*n_v, y+2*n_v, factor);
      break;
    default:
      ErrThrow("Unknown NSETYPE, it must be 1, 2, 3, 4, or 14");
      break;
  }
}

/** ************************************************************************ */
const TFESpace2D* BlockMatrixNSE2D::get_space_of_block(unsigned int b,
                                                       bool test) const
{
  unsigned int space_index = 0;
  if(test)
    space_index = this->block_pattern->test_space_index_of_block(b);
  else
    space_index = this->block_pattern->ansatz_space_index_of_block(b);
  if(space_index == 0)
    return this->get_velocity_space();
  else if(space_index == 1)
    return this->get_pressure_space();
  else
  {
    ErrThrow("could not find space for block " + std::to_string(b));
    throw; // only to avoid compiler warnings
  }
}

/** ************************************************************************ */
TSquareMatrix2D * BlockMatrixNSE2D::get_A_block(unsigned int i)
{
  if(i < 2)
    return (TSquareMatrix2D*)this->BlockMatrix::blocks[i].get();
  else if(i < 4)
    return (TSquareMatrix2D*)this->BlockMatrix::blocks[i+1].get();
  else
    ErrThrow("there are only four A-blocks! " + std::to_string(i));
  throw; // only to avoid a compiler warning
}

/** ************************************************************************ */
const TSquareMatrix2D * BlockMatrixNSE2D::get_A_block(unsigned int i) const 
{
  if(i < 2)
    return (TSquareMatrix2D*)this->BlockMatrix::blocks[i].get();
  else if(i < 4)
    return (TSquareMatrix2D*)this->BlockMatrix::blocks[i+1].get();
  else
    ErrThrow("there are only four A-blocks! " + std::to_string(i));
  throw; // only to avoid a compiler warning
}

/** ************************************************************************ */
TMatrix2D * BlockMatrixNSE2D::get_BT_block(unsigned int i)
{
  if(i < 2)
    return (TMatrix2D*)this->BlockMatrix::blocks[2+3*i].get();
  else
    ErrThrow("There are only two BT-blocks! " + std::to_string(i));
  throw; // only to avoid a compiler warning
}

/** ************************************************************************ */
const TMatrix2D * BlockMatrixNSE2D::get_BT_block(unsigned int i) const 
{
  if(i < 2)
    return (TMatrix2D*)this->BlockMatrix::blocks[2+3*i].get();
  else
    ErrThrow("There are only two BT-blocks! " + std::to_string(i));
  throw; // only to avoid a compiler warning
}

/** ************************************************************************ */
TMatrix2D * BlockMatrixNSE2D::get_B_block(unsigned int i)
{
  if(i < 2)
    return (TMatrix2D*)this->BlockMatrix::blocks[i+6].get();
  else
    ErrThrow("There are only two B-blocks! " + std::to_string(i));
  throw; // only to avoid a compiler warning
}

/** ************************************************************************ */
const TMatrix2D * BlockMatrixNSE2D::get_B_block(unsigned int i) const 
{
  if(i < 2)
    return (TMatrix2D*)this->BlockMatrix::blocks[i+6].get();
  else
    ErrThrow("There are only two B-blocks! " + std::to_string(i));
  throw; // only to avoid a compiler warning
}

