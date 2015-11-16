/** ************************************************************************ 
* @brief     source file for BlockMatrixDarcy2D
* @author    Ulrich Wilbrandt,
* @date      15.03.15
 ************************************************************************  */
#include <Database.h>
#include <BlockMatrixDarcy2D.h>
#include <Darcy2DMixed.h>
#include <SquareStructure2D.h>
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
 : BlockMatrix(Problem_type::Darcy, 2), 
   boundary_values({BoundValue[0], BoundValue[1]})
{
  // first build matrix structures
  // velocity-velocty coupling
  TSquareStructure2D *sqstructureA = new TSquareStructure2D(&velocity);
  sqstructureA->Sort();  // sort column numbers: numbers are in increasing order
  // pressure-pressure coupling
  TSquareStructure2D *sqstructureC = new TSquareStructure2D(&pressure);
  sqstructureC->Sort();  // sort column numbers: numbers are in increasing order
  // velocity-pressure and pressure-velocity coupling
  TStructure *structureB = new TStructure(&pressure, &velocity);
  TStructure *structureBT = new TStructure(&velocity, &pressure);
  
  // ( A  B1' )   ( 0 2 )
  // ( B2 C   )   ( 3 1 )
  // create the velocity-velocity coupling matrix 
  this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(sqstructureA));
  // create the pressure-pressure coupling matrix
  this->BlockMatrix::blocks[1].reset(new TMatrix2D(structureBT));
  // create velocity-pressure and pressure-velocity coupling matrices
  this->BlockMatrix::blocks[2].reset(new TMatrix2D(structureB));
  this->BlockMatrix::blocks[3].reset(new TSquareMatrix2D(sqstructureC));
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
     cout << "AMG_SOLVE not yet implemented " <<endl;
   break;

   case GMG:
     cout << "GMG solver not yet implemented " <<endl;
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
     OutPut("Unknown Solver\n");
     throw("Unknown Solver");
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
const TFESpace2D* BlockMatrixDarcy2D::get_space_of_block(unsigned int b,
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
  }
}

/** ************************************************************************ */


