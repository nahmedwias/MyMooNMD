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

BlockMatrixDarcy2D::BlockMatrixDarcy2D(TFESpace2D *velocity, TFESpace2D* pressure,
                                     BoundValueFunct2D **BoundValue)
 : BlockMatrix2D(2,2,2)
{
  //store the FEspace
  this->BlockMatrix2D::fe_spaces[0] = velocity;
  this->BlockMatrix2D::fe_spaces[1] = pressure;
  
  // build matrices
  // first build matrix structures
  // velocity-velocty coupling
  TSquareStructure2D *sqstructureA = new TSquareStructure2D(fe_spaces[0]);
  sqstructureA->Sort();  // sort column numbers: numbers are in increasing order
  // pressure-pressure coupling
  TSquareStructure2D *sqstructureC = new TSquareStructure2D(fe_spaces[1]);
  sqstructureC->Sort();  // sort column numbers: numbers are in increasing order
  // velocity-pressure and pressure-velocity coupling
  TStructure2D *structureB = new TStructure2D(fe_spaces[1], fe_spaces[0]);
  TStructure2D *structureBT = new TStructure2D(fe_spaces[0], fe_spaces[1]);
  
  // ( A  B1' )   ( 0 2 )
  // ( B2 C   )   ( 3 1 )
  /** create the velocity-velocity coupling matrix */
  this->BlockMatrix2D::sq_matrices[0] = new TSquareMatrix2D(sqstructureA);
  /** create the pressure-pressure coupling matrix */
  this->BlockMatrix2D::sq_matrices[1] = new TSquareMatrix2D(sqstructureC);
  /** create velocity-pressure and pressure-velocity coupling matrices */
  this->BlockMatrix2D::rect_matrices[0] = new TMatrix2D(structureBT);
  this->BlockMatrix2D::rect_matrices[1] = new TMatrix2D(structureB);
  
  // store the boundary conditions and boundary data
  this->BoundaryConditions[1] = fe_spaces[1]->GetBoundCondition();
  this->BoundaryConditions[0] = fe_spaces[0]->GetBoundCondition();
  this->BoundaryValues[0] = BoundValue[0];
  this->BoundaryValues[1] = BoundValue[1];
  
  this->BlockMatrix2D::defect = NULL; // to be implemented
}

BlockMatrixDarcy2D::~BlockMatrixDarcy2D()
{
  delete this->BlockMatrix2D::sq_matrices[0]->GetStructure();
  delete this->BlockMatrix2D::sq_matrices[1]->GetStructure();
  delete this->BlockMatrix2D::rect_matrices[0]->GetStructure();
  delete this->BlockMatrix2D::rect_matrices[1]->GetStructure();
  delete this->BlockMatrix2D::sq_matrices[0];
  delete this->BlockMatrix2D::sq_matrices[1];
  delete this->BlockMatrix2D::rect_matrices[0];
  delete this->BlockMatrix2D::rect_matrices[1];
}

void BlockMatrixDarcy2D::Assemble(LocalAssembling2D& la, double *sol,
                                 double *rhs)
{
  int N_U = this->BlockMatrix2D::fe_spaces[0]->GetN_DegreesOfFreedom();
  int N_U_Active = this->BlockMatrix2D::fe_spaces[0]->GetActiveBound();
  int N_P = this->BlockMatrix2D::fe_spaces[1]->GetN_DegreesOfFreedom();
  int N_DOF = N_U + N_P;
  
  memset(rhs, 0, N_DOF*SizeOfDouble);
  double *rhs_blocks[2] = { rhs, rhs+N_U };
 
  // reset matrices to zero
  this->BlockMatrix2D::sq_matrices[0]->Reset();
  this->BlockMatrix2D::sq_matrices[1]->Reset();
  this->BlockMatrix2D::rect_matrices[0]->Reset();
  this->BlockMatrix2D::rect_matrices[1]->Reset();
  
  // assemble
  Assemble2D_VectFE(2, &fe_spaces[0], 2, &sq_matrices[0], 2, &rect_matrices[0],
                    2, rhs_blocks, &fe_spaces[0], la, BoundaryConditions,
                    BoundaryValues);
  
  // copy Dirichlet values from rhs to solution vector (this is not really 
  // necessary in case of a direct solver)
  memcpy(sol+N_U_Active, rhs+N_U_Active, (N_U-N_U_Active)*SizeOfDouble);
} // void BlockMatrixDarcy2D::Assemble


void BlockMatrixDarcy2D::Solve(double *sol, double *rhs)
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
     DirectSolver(sq_matrices[0], sq_matrices[1], rect_matrices[0], 
                  rect_matrices[1], rhs, sol);
   break;
 
   default:
     OutPut("Unknown Solver\n");
     throw("Unknown Solver");
  }
}

void BlockMatrixDarcy2D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->BlockMatrix2D::sq_matrices[0]->GetN_Rows();
  // reset y
  memset(y, 0.0, n_total_rows * SizeOfDouble);
  this->apply_scaled_add(x, y, factor);
}

void BlockMatrixDarcy2D::apply_scaled_add(const double *x, double *y,
                                  double factor) const
{
  // number of velocity degrees of freedom
  unsigned int n_v = this->BlockMatrix2D::sq_matrices[0]->GetN_Rows();
  
  this->BlockMatrix2D::sq_matrices[0]->multiply(  x,     y,     factor);
  this->BlockMatrix2D::rect_matrices[0]->multiply(x+n_v, y,     factor);
  
  this->BlockMatrix2D::rect_matrices[1]->multiply(x,     y+n_v, factor);
  this->BlockMatrix2D::sq_matrices[1]->multiply(  x+n_v, y+n_v, factor);
}

unsigned int BlockMatrixDarcy2D::n_rows() const
{
  return 2;
}

unsigned int BlockMatrixDarcy2D::n_cols() const
{
  return 2;
}

unsigned int BlockMatrixDarcy2D::n_total_rows() const
{
  return this->BlockMatrix2D::sq_matrices[0]->GetN_Rows() 
       + this->BlockMatrix2D::rect_matrices[1]->GetN_Rows();
}

unsigned int BlockMatrixDarcy2D::n_total_cols() const
{
  return this->BlockMatrix2D::sq_matrices[0]->GetN_Columns() 
       + this->BlockMatrix2D::rect_matrices[0]->GetN_Columns();
}

unsigned int BlockMatrixDarcy2D::n_total_entries() const
{
  return this->BlockMatrix2D::sq_matrices[0]->GetN_Entries()
      + this->BlockMatrix2D::sq_matrices[1]->GetN_Entries()
      + this->BlockMatrix2D::rect_matrices[0]->GetN_Entries()
      + this->BlockMatrix2D::rect_matrices[1]->GetN_Entries();
}



