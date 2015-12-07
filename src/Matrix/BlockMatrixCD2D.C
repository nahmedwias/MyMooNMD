/** ************************************************************************ 
* @brief     source file for BlockMatrixCD2D
* @author    Ulrich Wilbrandt, 
* @date      11.09.15
* @History 
 ************************************************************************  */

#include <Database.h>
#include <BlockMatrixCD2D.h>
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
                                 BoundValueFunct2D * const BoundValue,
                                 bool mass_matrix)
 : BlockMatrix(Problem_type::ConvDiffReac, 2, mass_matrix),
   boundary_values(BoundValue)
{
  // the stiffness/system matrix for a convection diffusion problem
  this->BlockMatrix::blocks[0].reset(new TSquareMatrix2D(&fespace));
  const TStructure& sqStructure = this->BlockMatrix::blocks[0]->GetStructure();
  // the number of active entries
  unsigned int n_active = sqStructure.GetN_Entries();
  // substract the number of non active entries (non active rows)
  n_active -= sqStructure.GetN_Rows() - sqStructure.GetActiveBound();
  this->BlockMatrix::actives[0] = n_active;
}

/** ************************************************************************ */
BlockMatrixCD2D::~BlockMatrixCD2D()
{
}

/** ************************************************************************ */
void BlockMatrixCD2D::Assemble(LocalAssembling2D& la, BlockVector& sol,
                               BlockVector& rhs)
{
  const TFESpace2D * fe_space = this->get_fe_space();
  BoundCondFunct2D * boundary_conditions = fe_space->GetBoundCondition();
  int N_Matrices = 1;
  double * rhs_entries = rhs.get_entries();
  TSquareMatrix2D * matrix = this->get_matrix();
  
  BoundValueFunct2D * non_const_bound_value = this->boundary_values;
  
  if(la.get_type() != TCD2D_Mass)
  {
    // reset right hand side and matrix to zero
    rhs.reset();
    matrix->reset();
    
    // assemble
    Assemble2D(1, &fe_space, N_Matrices, &matrix, 0, NULL, 1, &rhs_entries, 
               &fe_space, &boundary_conditions, &non_const_bound_value, la);
    
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
        ErrThrow("LP_FULL_GRADIENT needs to be one to use LOCAL_PROJECTION");
      }
    }
    
    // copy Dirichlet values from rhs to solution vector (this is not really 
    // necessary in case of a direct solver)
    sol.copy_nonactive(rhs);
  }
  else // assembling mass matrix
  //if(la.get_type() == TCD2D_Mass_Rhs_Galerkin)
  {
    // reset the matrix
    matrix->reset();
    Assemble2D(1, &fe_space, N_Matrices, &matrix, 0, NULL, 0, NULL, 
               NULL, &boundary_conditions, &non_const_bound_value, la);
  }
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

/** ************************************************************************ */
void BlockMatrixCD2D::apply_scaled_add_active(const double* x, double* y, 
					      double factor) const
{
  if(factor ==0)
    return;
  
  this->get_matrix()->multiplyActive(x,y,factor); 
  
}




