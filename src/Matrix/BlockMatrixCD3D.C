#include <Database.h>
#include <BlockMatrixCD3D.h>
#include <DiscreteForm3D.h>
#include <Assemble3D.h>
#include <AuxParam3D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
//#include <AssembleMat3D.h>

#include <stdlib.h>
#include <string.h>

#include <Solver.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>

#ifdef _MPI
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

BlockMatrixCD3D::BlockMatrixCD3D(const TFESpace3D &feSpace,
        						BoundValueFunct3D *BoundValue,
								bool massMatrix)
: BlockMatrix(Problem_type::ConvDiffReac, 3, massMatrix),
  boundaryValues_(BoundValue)
{
	// the stiffness/system matrix for a convection diffusion problem
	BlockMatrix::blocks[0].reset(new TSquareMatrix3D(&feSpace));

  const TStructure& sqStructure = this->BlockMatrix::blocks[0]->GetStructure();

	// START The following code is to determine the number of active entries
	unsigned int n_active = sqStructure.GetN_Entries();
	// substract the number of non active entries (non active rows)
	n_active -= sqStructure.GetN_Rows() - sqStructure.GetActiveBound();
	BlockMatrix::actives[0] = n_active;
	// END FIXME CB Is this correct - I'd rather like this not to be stored by the class!
}

/** ************************************************************************ */
//CB assemble method is copied from BlockMatrixCD2D and adapted.
void BlockMatrixCD3D::assemble(const LocalAssembling3D& la, BlockVector& sol,
                               BlockVector& rhs)
{

  const TFESpace3D * fe_space = this->get_fe_space();
  BoundCondFunct3D* boundary_conditions = fe_space->getBoundCondition();
  int N_Matrices = 1;
  double * rhs_entries = rhs.get_entries();
  TSquareMatrix3D * matrix = this->get_matrix();

  BoundValueFunct3D * non_const_bound_value = boundaryValues_;

//  if(la.get_type() != TCD2D_Mass_Rhs_Galerkin) //CB FIXME care for assembling type of the la-object
//  {
    // reset right hand side and matrix to zero
    rhs.reset();
    matrix->Reset();

    // assemble
    Assemble3D(1, &fe_space, N_Matrices, &matrix, 0, NULL, 1, &rhs_entries,
               &fe_space, &boundary_conditions, &non_const_bound_value, la);

    //FIXME CB which stabilization to begin with for 3D?

    // copy Dirichlet values from rhs to solution vector (this is not really
    // necessary in case of a direct solver)
    sol.copy_nonactive(rhs);

    //  }
//  else // assembling mass matrix
//  //if(la.get_type() == TCD2D_Mass_Rhs_Galerkin)
//  {
//    // reset the matrix
//    matrix->Reset();
//    Assemble2D(1, &fe_space, N_Matrices, &matrix, 0, NULL, 1, &rhs_entries,
//               &fe_space, &boundary_conditions, &non_const_bound_value, la);
//  }
}



/** ************************************************************************ */
void BlockMatrixCD3D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = get_matrix()->GetN_Rows();
  // reset y
  memset(y, 0.0, n_total_rows*SizeOfDouble);
  this->apply_scaled_add(x, y, factor);
}

/** ************************************************************************ */
void BlockMatrixCD3D::apply_scaled_add(const double *x, double *y,
                                          double factor) const
{
  get_matrix()->multiply(x, y, factor);
}
