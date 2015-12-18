/** ************************************************************************ 
* @brief     source file for BlockMatrixDarcy2D
* @author    Ulrich Wilbrandt,
* @date      15.03.15
 ************************************************************************  */
#include <Database.h>
#include <BlockMatrixDarcy2D.h>
#include <Darcy2DMixed.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>

/** ************************************************************************ */
BlockMatrixDarcy2D::BlockMatrixDarcy2D(const TFESpace2D& velocity, 
                                       const TFESpace2D& pressure,
                                       BoundValueFunct2D*const*const BoundValue)
 : BlockFEMatrix(Problem_type::Darcy, 2),
   boundary_values({BoundValue[0], BoundValue[1]})
{
  // ( A  B1' )   ( 0 2 )
  // ( B2 C   )   ( 3 1 )
  // create the velocity-velocity coupling matrix 
  this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(&velocity));
  // create the pressure-pressure coupling matrix
  this->BlockMatrix::blocks[1].reset(new TMatrix2D(&velocity, &pressure));
  // create velocity-pressure and pressure-velocity coupling matrices
  this->BlockMatrix::blocks[2].reset(new TMatrix2D(&pressure, &velocity));
  this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(&pressure));
}

/** ************************************************************************ */
BlockMatrixDarcy2D::~BlockMatrixDarcy2D()
{
}

/** ************************************************************************ */
void BlockMatrixDarcy2D::Assemble(LocalAssembling2D& la, BlockVector& rhs)
{
  const TFESpace2D * v_space = this->get_velocity_space();
  const TFESpace2D * p_space = this->get_pressure_space();
  
  rhs.reset();
  double *rhs_blocks[2] = { rhs.block(0), rhs.block(1) };
 
  // reset matrices to zero
  this->BlockMatrix::reset();
  
  const TFESpace2D * pointers_to_spaces[2] = { v_space, p_space };
  // the following cast works, because these blocks were constructed as 
  // TSquareMatrix2D in the constructor of this class 
  TSquareMatrix2D * sq_matrices[2] = {this->get_A_block(), this->get_C_block()};
  TMatrix2D * rect_matrices[2] = { this->get_BT_block(), this->get_B_block() };
  BoundCondFunct2D * boundary_conditions[2] = { 
    v_space->GetBoundCondition(), p_space->GetBoundCondition()};
  // assemble
  Assemble2D_VectFE(2, pointers_to_spaces, 2, sq_matrices, 2, 
                    rect_matrices, 2, rhs_blocks, pointers_to_spaces, la, 
                    boundary_conditions, this->boundary_values.data());
} // void BlockMatrixDarcy2D::Assemble

/** ************************************************************************ */
void BlockMatrixDarcy2D::Solve(BlockVector& sol, BlockVector& rhs)
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
   {
     // the following cast works, because these blocks were constructed as 
     // TSquareMatrix2D in the constructor of this class 
     DirectSolver(this->get_A_block(), this->get_C_block(), 
                  this->get_BT_block(), this->get_B_block(), rhs.get_entries(),
                  sol.get_entries());
   }
   break;
 
   default:
     ErrThrow("Unknown Solver\n");
  }
}

/** ************************************************************************ */
void BlockMatrixDarcy2D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->BlockMatrix::n_total_rows();
  // reset y
  memset(y, 0.0, n_total_rows * SizeOfDouble);
  this->apply_scaled_add(x, y, factor);
}

/** ************************************************************************ */
void BlockMatrixDarcy2D::apply_scaled_add(const double *x, double *y,
                                  double factor) const
{
  // number of velocity degrees of freedom
  unsigned int n_v = this->get_A_block()->GetN_Rows();
  
  this->block(0)->multiply(  x,     y,     factor);
  this->block(1)->multiply(x+n_v,   y,     factor);
  
  this->block(2)->multiply(x,       y+n_v, factor);
  this->block(3)->multiply(  x+n_v, y+n_v, factor);
}

/** ************************************************************************ */


