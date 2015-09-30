/** ************************************************************************ 
* @brief     source file for BlockMatrixCD2D
* @author    Ulrich Wilbrandt, 
* @date      11.09.15
* @History 
 ************************************************************************  */

#include <Database.h>
#include <BlockMatrixCD2D.h>
#include <SquareStructure2D.h>
#include <Assemble2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>

/** ************************************************************************ */
BlockMatrixCD2D::BlockMatrixCD2D(const TFESpace2D &fespace, 
                                 const BoundValueFunct2D *BoundValue)
 : BlockMatrix(Problem_type::ConvDiffReac, 2, 
               TDatabase::ParamDB->PROBLEM_TYPE != 2),
   boundary_values(BoundValue)
{
  // build matrices, first build matrix structure
  TSquareStructure2D* sqstructure = new TSquareStructure2D(&fespace);
  sqstructure->Sort();  // sort column numbers: numbers are in increasing order

  // the stiffness/system matrix for a convection diffusion problem
  this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(sqstructure));
  this->BlockMatrix::actives[0] = sqstructure->GetActiveBound();
  // create the mass matrix in case of the time dependent case
  if(TDatabase::ParamDB->PROBLEM_TYPE == 2)
  {
    this->BlockMatrix::blocks[1].reset(new TSquareMatrix2D(sqstructure));
    this->BlockMatrix::actives[1] = sqstructure->GetActiveBound();
  }
}

/** ************************************************************************ */
BlockMatrixCD2D::~BlockMatrixCD2D()
{
  // the combined matrix is this->get_matrix(), and will be deleted in 
  // ~BlockMatrix
  if(!this->combined_matrix)
    delete this->get_matrix()->GetStructure();
}

/** ************************************************************************ */
void BlockMatrixCD2D::Assemble(LocalAssembling2D& la, BlockVector& sol,
                               BlockVector& rhs)
{
  const TFESpace2D * fe_space = this->get_fe_space();
  
  // reset right hand side and matrix to zero
  rhs.reset();
  double * rhs_entries = rhs.get_entries();
  TSquareMatrix2D * matrix = this->get_matrix();
  matrix->Reset();
  
  BoundCondFunct2D * boundary_conditions = fe_space->GetBoundCondition();
  int N_Matrices = 1;
  // assemble
  Assemble2D(1, &fe_space, N_Matrices, &matrix, 0, NULL, 1, &rhs_entries, 
             &fe_space, &boundary_conditions, &boundary_values, la);
 
  // apply local projection stabilization method
  if(TDatabase::ParamDB->DISCTYPE==LOCAL_PROJECTION 
     && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
  {
    if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
    { 
      UltraLocalProjection((void *)&matrix, false);
    }
    else
    {
      OutPut("Check! LP_FULL_GRADIENT needs to be one to use LOCAL_PROJECTION" << endl);
      exit(4711);
    }
  }
  
  // copy Dirichlet values from rhs to solution vector (this is not really 
  // necessary in case of a direct solver)
  sol.copy_nonactive(rhs);
} // void BlockMatrixCD2D::Assemble

/** ************************************************************************ */
void BlockMatrixCD2D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->get_matrix()->GetN_Rows();
  // reset y
  memset(y, 0.0, n_total_rows*SizeOfDouble);
  this->apply_scaled_add(x, y, factor);
}

/** ************************************************************************ */
void BlockMatrixCD2D::apply_scaled_add(const double *x, double *y, 
                                          double factor) const
{
  this->get_matrix()->multiply(x, y, factor);
}





