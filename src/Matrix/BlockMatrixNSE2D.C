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
                                   BoundValueFunct2D * const * const BoundValue,
                                   bool massMatrix)
    : BlockMatrix(Problem_type::NavierStokes, 2, massMatrix),
      boundary_values({{BoundValue[0], BoundValue[1], BoundValue[2]}})
{ 
  // build matrices
  this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(&velocity));
  // get that matrix A11, which will be used to construct the others if NSTYPE>2
  const TSquareMatrix2D & A11 = 
    (const TSquareMatrix2D &) *this->BlockMatrix::blocks[0];
  
  this->BlockMatrix::blocks[6].reset(new TMatrix2D(&pressure, &velocity));
  // get that matrix B1, which will be used to construct the other
  const TMatrix2D & B1 = (const TMatrix2D &) *this->BlockMatrix::blocks[6];
  this->BlockMatrix::blocks[7].reset(new TMatrix2D(B1));
  
  
  // number of pressure degrees of freedom
  unsigned int n_p = pressure.GetN_DegreesOfFreedom();
  // number of velocity degrees of freedom
  unsigned int n_v = velocity.GetN_DegreesOfFreedom();
  if(!massMatrix)
  {
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
        std::shared_ptr<TStructure> empty_block_velocity(new TStructure(n_v));
        std::shared_ptr<TStructure> empty_block_pressure(new TStructure(n_p));
        
        this->BlockMatrix::blocks[1].reset(new TMatrix(empty_block_velocity));
        //this->BlockMatrix::blocks[2] stays a nullptr
        this->BlockMatrix::blocks[3].reset(new TMatrix(empty_block_velocity));
        this->BlockMatrix::blocks[4] = this->BlockMatrix::blocks[0];
        //this->BlockMatrix::blocks[5] stays a nullptr
        this->BlockMatrix::blocks[8].reset(new TMatrix(empty_block_pressure));
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
        std::shared_ptr<TStructure> empty_block_velocity(new TStructure(n_v));
        std::shared_ptr<TStructure> empty_block_pressure(new TStructure(n_p));
        
        this->BlockMatrix::blocks[1].reset(new TMatrix(empty_block_velocity));
        this->BlockMatrix::blocks[2].reset(new TMatrix2D(&velocity, &pressure));
        // get that matrix B1T, which will be used to construct the other
        const TMatrix2D & B1T = (const TMatrix2D &) *this->BlockMatrix::blocks[2];
        this->BlockMatrix::blocks[3].reset(new TMatrix(empty_block_velocity));
        this->BlockMatrix::blocks[4] = this->BlockMatrix::blocks[0];
        this->BlockMatrix::blocks[5].reset(new TMatrix2D(B1T));
        this->BlockMatrix::blocks[8].reset(new TMatrix(empty_block_pressure));
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
        std::shared_ptr<TStructure> empty_block_pressure(new TStructure(n_p));
        
        this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(A11));
        //this->BlockMatrix::blocks[2] stays a nullptr
        this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(A11));
        this->BlockMatrix::blocks[4].reset(new TSquareMatrix2D(A11));
        //this->BlockMatrix::blocks[5] stays a nullptr
        this->BlockMatrix::blocks[8].reset(new TMatrix(empty_block_pressure));
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
        std::shared_ptr<TStructure> structureBT(new TStructure(&velocity,
                                                               &pressure));
        std::shared_ptr<TStructure> empty_block_pressure(new TStructure(n_p));
        
        this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(A11));
        this->BlockMatrix::blocks[2].reset(new TMatrix2D(&velocity, &pressure));
        // get that matrix B1T, which will be used to construct the other
        const TMatrix2D & B1T = (const TMatrix2D &) *this->BlockMatrix::blocks[2];
        this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(A11));
        this->BlockMatrix::blocks[4].reset(new TSquareMatrix2D(A11));
        this->BlockMatrix::blocks[5].reset(new TMatrix2D(B1T));
        this->BlockMatrix::blocks[8].reset(new TMatrix(empty_block_pressure));
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
        
        this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(A11));
        this->BlockMatrix::blocks[2].reset(new TMatrix2D(&velocity, &pressure));
        // get that matrix B1T, which will be used to construct the other
        const TMatrix2D & B1T = (const TMatrix2D &) *this->BlockMatrix::blocks[2];
        this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(A11));
        this->BlockMatrix::blocks[4].reset(new TSquareMatrix2D(A11));
        this->BlockMatrix::blocks[5].reset(new TMatrix2D(B1T));
        this->BlockMatrix::blocks[8].reset(new TSquareMatrix2D(&pressure));
        break;
      }
      default:
        ErrThrow("Unknown NSETYPE, it must be 1, 2, 3, 4, or 14");
        break;
    }
  }
  else
  {
    // Here we need the off diagonal A blocks and all B's with
    // no entries. But the method which is used in the BlockMatrix
    // class needs with zero entries so that's why they are 
    // created with zero entries here.
    this->BlockMatrix::blocks[2].reset(new TMatrix2D(&velocity, &pressure));
    this->BlockMatrix::blocks[5].reset(new TMatrix2D(&velocity, &pressure));
    
    this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(A11));
    this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(A11));
    this->BlockMatrix::blocks[4].reset(new TSquareMatrix2D(A11));
    
    this->BlockMatrix::blocks[8].reset(new TSquareMatrix2D(&pressure));
  }
}

/** ************************************************************************ */
BlockMatrixNSE2D::~BlockMatrixNSE2D()
{
}

/** ************************************************************************ */
void BlockMatrixNSE2D::Assemble(LocalAssembling2D& la, BlockVector& rhs, 
                                BlockMatrixNSE2D* mass)
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
  
  unsigned int n_sq_mat = 6;
  if(TDatabase::ParamDB->NSTYPE == 1 || TDatabase::ParamDB->NSTYPE == 2)
    n_sq_mat = 2;
  if(TDatabase::ParamDB->NSTYPE == 14)
  {
    Output::print("WARNING: NSTYPE 14 needs a special local assembling, the ",
                  "C-block is now ignored");
    //n_sq_mat = 5;
  }
  TSquareMatrix2D * sq_matrices[5] = {this->get_A_block(0), nullptr, nullptr, 
                                      nullptr, nullptr};
  if(TDatabase::ParamDB->NSTYPE == 1 || TDatabase::ParamDB->NSTYPE == 2)
    sq_matrices[1] = mass->get_A_block(0); // mass matrix
  else
  {
    sq_matrices[1] = this->get_A_block(1);
    sq_matrices[2] = this->get_A_block(2);
    sq_matrices[3] = this->get_A_block(3);
    // mass matrices 
    sq_matrices[4] = mass->get_A_block(0);
    sq_matrices[5] = mass->get_A_block(3);
    if(TDatabase::ParamDB->NSTYPE == 14)
    {
      n_sq_mat += 1;
      sq_matrices[6] = this->get_C_block();
    }
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
        Output::print<3>("UPWINDING DONE : level ");
        break;
      case 3:
      case 4:
      case 14:
        // do upwinding with two matrices
        Output::print<3>("UPWINDING DONE : level ");
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
      ErrThrow("For slip with friction bc NSTYPE 4 is necessary !!!!! ");
    }
    ErrThrow("Assemble2DSlipBC does not work");
    
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
  this->get_A_block(1)->resetNonActive();
  this->get_A_block(2)->resetNonActive();
  this->get_BT_block(0)->resetNonActive();
  this->get_BT_block(1)->resetNonActive();
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
        Output::print<3>("UPWINDING DONE : level ");
        break;
        
      case 3:
      case 4:
      case 14:
        // do upwinding with two matrices
        Output::print<3>("UPWINDING DONE : level ");
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
    this->combined_matrix.reset();
  }
} //BlockMatrixNSE2D::AssembleNonLinear(

/** ************************************************************************ */
void BlockMatrixNSE2D::Solve(double *sol, double *rhs)
{
  switch(TDatabase::ParamDB->SOLVER_TYPE)
  {
    case AMG_SOLVE:
      Output::print("AMG_SOLVE not yet implemented");
      break;
      
    case GMG:
      Output::print("GMG solver not yet implemented");
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
          ErrThrow("Solver not included for NSTYPE 3 in this version. Try ",
                   "NSTYPE 4");
          break;
          
        case 4:
          DirectSolver(this->get_A_block(0), this->get_A_block(1), 
                       this->get_A_block(2), this->get_A_block(3),
                       this->get_BT_block(0),  this->get_BT_block(1),
                       this->get_B_block(0), this->get_B_block(1), rhs, sol);
          break;
        case 14:
          Output::print("WARNING: NSTYPE 14 is not fully supported, take ",
                        "NSTYPE 4");
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
void BlockMatrixNSE2D::scaleActive(const double factor)
{
  if(factor == 1.)
    return;
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
    case 2:
      this->get_A_block(0)->scaleActive(factor);
      break;
    case 3:
    case 4:
    case 14:
      this->get_A_block(0)->scaleActive(factor);
      this->get_A_block(1)->scaleActive(factor);
      this->get_A_block(2)->scaleActive(factor);
      this->get_A_block(3)->scaleActive(factor);      
      break;
  }

}

/** ************************************************************************ */
void BlockMatrixNSE2D::applyScaledAddActive(const double* x, double* y, double factor) const
{
  if(factor = 0.)
    return;
  unsigned int n_v = this->get_velocity_space()->GetN_DegreesOfFreedom();
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
    case 2:
      this->get_A_block(0)->multiplyActive(x,     y, factor);
      this->get_A_block(0)->multiplyActive(x+n_v, y+n_v, factor);
      break;    
    case 3:
    case 4:
    case 14:
      this->get_A_block(0)->multiplyActive(x,     y,     factor);
      this->get_A_block(1)->multiplyActive(x+n_v, y,     factor);
      this->get_A_block(2)->multiplyActive(x,     y+n_v,  factor);
      this->get_A_block(3)->multiplyActive(x+n_v, y+n_v, factor);
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
    ErrThrow("could not find space for block ", b);
  }
}

/** ************************************************************************ */
TSquareMatrix2D * BlockMatrixNSE2D::get_A_block(unsigned int i)
{
  if(i == 0)
    return (TSquareMatrix2D*)this->BlockMatrix::blocks[0].get();
  else if(i == 1)
  {
    if(TDatabase::ParamDB->NSTYPE > 2)
      return (TSquareMatrix2D*)this->BlockMatrix::blocks[1].get();
    else
      ErrThrow("block A12 does not exist for this NSTYPE");
  }
  else if(i == 2)
  {
    if(TDatabase::ParamDB->NSTYPE > 2)
      return (TSquareMatrix2D*)this->BlockMatrix::blocks[3].get();
    else
      ErrThrow("block A21 does not exist for this NSTYPE");
  }
  else if(i == 3)
    return (TSquareMatrix2D*)this->BlockMatrix::blocks[4].get();
  else
    ErrThrow("there are only four A-blocks! ", i);
}

/** ************************************************************************ */
const TSquareMatrix2D * BlockMatrixNSE2D::get_A_block(unsigned int i) const 
{
  if(i == 0)
    return (TSquareMatrix2D*)this->BlockMatrix::blocks[0].get();
  else if(i == 1)
  {
    if(TDatabase::ParamDB->NSTYPE > 2)
      return (TSquareMatrix2D*)this->BlockMatrix::blocks[1].get();
    else
      ErrThrow("block A12 does not exist for this NSTYPE");
  }
  else if(i == 2)
  {
    if(TDatabase::ParamDB->NSTYPE > 2)
      return (TSquareMatrix2D*)this->BlockMatrix::blocks[3].get();
    else
      ErrThrow("block A21 does not exist for this NSTYPE");
  }
  else if(i == 3)
    return (TSquareMatrix2D*)this->BlockMatrix::blocks[4].get();
  else
    ErrThrow("there are only four A-blocks! ", i);
}

/** ************************************************************************ */
TSquareMatrix2D * BlockMatrixNSE2D::get_C_block()
{
  if(TDatabase::ParamDB->NSTYPE != 14)
  {
    ErrThrow("the C block only exists in case of NSTYPE 14");
  }
  return (TSquareMatrix2D*)this->BlockMatrix::blocks.at(8).get();
}

/** ************************************************************************ */
const TSquareMatrix2D * BlockMatrixNSE2D::get_C_block() const
{
  if(TDatabase::ParamDB->NSTYPE != 14)
  {
    ErrThrow("the C block only exists in case of NSTYPE 14");
  }
  return (TSquareMatrix2D*)this->BlockMatrix::blocks.at(8).get();
}

/** ************************************************************************ */
TMatrix2D * BlockMatrixNSE2D::get_BT_block(unsigned int i)
{
  if(i < 2)
    return (TMatrix2D*)this->BlockMatrix::blocks[2+3*i].get();
  else
    ErrThrow("There are only two BT-blocks! ", i);
}

/** ************************************************************************ */
const TMatrix2D * BlockMatrixNSE2D::get_BT_block(unsigned int i) const 
{
  if(TDatabase::ParamDB->NSTYPE == 1 || TDatabase::ParamDB->NSTYPE == 3)
  {
    ErrThrow("For NSTYPE 1 and 3 the BT blocks are not explicitly stored");
  }
  if(i < 2)
    return (TMatrix2D*)this->BlockMatrix::blocks[2+3*i].get();
  else
    ErrThrow("There are only two BT-blocks! ", i);
}

/** ************************************************************************ */
TMatrix2D * BlockMatrixNSE2D::get_B_block(unsigned int i)
{
  if(i < 2)
    return (TMatrix2D*)this->BlockMatrix::blocks[i+6].get();
  else
    ErrThrow("There are only two B-blocks! ", i);
}

/** ************************************************************************ */
const TMatrix2D * BlockMatrixNSE2D::get_B_block(unsigned int i) const 
{
  if(i < 2)
    return (TMatrix2D*)this->BlockMatrix::blocks[i+6].get();
  else
    ErrThrow("There are only two B-blocks! ", i);
}

