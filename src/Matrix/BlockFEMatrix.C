#include <BlockFEMatrix.h>
#include <FEDatabase2D.h>
#include <BlockVector.h>
#include <Database.h>
#include <stdio.h>
#include <string.h>
#include <MooNMD_Io.h>
#include <LinAlg.h>
#include <limits>

/** ************************************************************************* */
BlockFEMatrix::BlockFEMatrix(const Problem_type type, unsigned int space_dimension,
                bool mass_matrix)
: BlockMatrix(type, space_dimension,  mass_matrix),
  actives(block_pattern->n_blocks(), std::numeric_limits<unsigned int>::max())
{
  ;
}


/** ************************************************************************* */

BlockFEMatrix::BlockFEMatrix(const TFESpace2D* testSpace, 
                             const TFESpace2D* ansatzSpace,
                             const TFESpace2D *pressureSpace)
 : BlockMatrix(3,1), TestSpace(testSpace), AnsatzSpace(ansatzSpace), 
   PressureSpace(pressureSpace)
{
  // checks
  if( (testSpace->GetCollection()) != (ansatzSpace->GetCollection() ))
  {
    Output::print<1>(testSpace->GetCollection()->GetN_Cells(), " ", 
                     ansatzSpace->GetCollection()->GetN_Cells());
    ErrThrow("grid is different for the test and ansatz sapces");
  }
  Output::print<2>("dof of testSpace  : ",   testSpace->GetN_DegreesOfFreedom());
  Output::print<2>("dof of ansatzSpace: ", ansatzSpace->GetN_DegreesOfFreedom());
  
  this->BlockMatrix::blocks[0].reset(new FEMatrix(testSpace, ansatzSpace));  
  this->BlockMatrix::blocks[1].reset(new FEMatrix(testSpace, ansatzSpace)); 
  std::shared_ptr<TStructure> structure(new 
                  TStructure(PressureSpace->GetN_DegreesOfFreedom(),
                             AnsatzSpace->GetN_DegreesOfFreedom()));
  this->BlockMatrix::blocks[2].reset(new TMatrix(structure));
  
  TCollection *coll = testSpace->GetCollection();
  unsigned int nCells = coll->GetN_Cells();
  // set normal orientation
  for(unsigned int i = 0; i<nCells; ++i)
    coll->GetCell(i)->SetNormalOrientation();
  // common declaration
  int nPoints;
  double *xi, *eta;
  
  for(unsigned int c=0; c<nCells; ++c)
  {
    TBaseCell *cell = coll->GetCell(c);
    
    TFE2D *elementAnsatz = TFEDatabase2D::GetFE2D(ansatzSpace->GetFE2D(c,cell));
    // basis funct obj
    TBaseFunct2D *baseFunctAnsatz = elementAnsatz->GetBaseFunct2D(); 
    // dimension of baseFunctAnsatz
    int baseVectDim = baseFunctAnsatz->GetBaseVectDim();
    // no of basis function
    int nDofAnsatz = baseFunctAnsatz->GetDimension(); 
    // compute the points on the reference element
    TNodalFunctional2D *nf = elementAnsatz->GetNodalFunctional2D();
    nf->GetPointsForAll(nPoints, xi, eta);
    
    // everything needed for the test space
    TFE2D *elementTest = TFEDatabase2D::GetFE2D(testSpace->GetFE2D(c,cell));
    TBaseFunct2D* baseFunctTest = elementTest->GetBaseFunct2D();
    int nDofTest = baseFunctTest->GetDimension();
        
    
    // the id of the reference transformation
    RefTrans2D refTransfID = elementAnsatz->GetRefTransID();
    TFEDatabase2D::SetCellForRefTrans(cell, refTransfID);
    
    // number of basis functions, this is the length of the array needed to 
    // evaluate the basis functions
    int nBaseFunct = nDofTest*baseVectDim;
    double uorig[nPoints][nBaseFunct];
    double AllPointValues[nDofTest];
    
    for(int i=0; i<nPoints; ++i)
    {
      baseFunctTest->GetDerivatives(D00, xi[i], eta[i], AllPointValues);
      TFEDatabase2D::GetOrigValues(refTransfID, xi[i], eta[i], baseFunctTest, 
                                   coll, (TGridCell *)cell, AllPointValues, 
                                   nullptr,nullptr, uorig[i], nullptr, nullptr );
    }
    
    double PointValuesx[nPoints* baseVectDim];
    double PointValuesy[nPoints* baseVectDim];
    memset(PointValuesx,0,nPoints* baseVectDim*sizeof(double));
    memset(PointValuesy,0,nPoints* baseVectDim*sizeof(double));
    double FunctionalValuesx[nDofAnsatz];
    double FunctionalValuesy[nDofAnsatz];
    
    for(int j=0; j<nDofTest; ++j)
    {
      if(testSpace->GetGlobalDOF(c)[j] >= testSpace->GetN_ActiveDegrees() )
        continue;
      for(int k=0; k<nPoints; ++k)
      {
        PointValuesx[k]         = uorig[k][j];
        PointValuesy[k+nPoints] = uorig[k][j];
      }
      nf->GetAllFunctionals(coll, cell, PointValuesx, FunctionalValuesx);
      nf->GetAllFunctionals(coll, cell, PointValuesy, FunctionalValuesy);
      
      for(int k=0; k<nDofAnsatz; k++)
      {
        // fill the correspoding blocks
        this->block(0)->set(testSpace->GetGlobalDOF(c)[j],
               ansatzSpace->GetGlobalDOF(c)[k], FunctionalValuesx[k]);
        this->block(1)->set(testSpace->GetGlobalDOF(c)[j],
                            ansatzSpace->GetGlobalDOF(c)[k], FunctionalValuesy[k]);
      }
    }
    // this->block(0)->PrintFull();
  }
  // this->get_combined_matrix()->PrintFull();
}
/** ************************************************************************* */
void BlockFEMatrix::add_scaled_active(const BlockMatrix& A, double factor)
{
  unsigned int n_blocks = A.n_blocks();
  if(this->n_blocks() != n_blocks)
  {
    ErrThrow("BlockFEMatrix::add_scaled_active : the two BlockMatrix objects do ",
             "not have the same number of blocks.");
  }

  for(unsigned int b = 0; b < n_blocks; b++)
  {
    if(this->actives[b] <= (unsigned int) this->blocks[b]->GetN_Entries())
    {
      // note: this could be a method of TMatrix as well, if TMatrix knew its
      // actives
      Daxpy(this->actives[b], factor, A.block(b)->GetEntries(),
            this->blocks[b]->GetEntries());
    }
    else
      ErrThrow("adding actives where the number of actives has not been set");
  }
}

/* ************************************************************************* */
void BlockFEMatrix::scale_active(double factor)
{
  for(unsigned int b = 0; b < this->n_blocks(); ++b)
  {
    if(this->actives[b] <= (unsigned int) this->blocks[b]->GetN_Entries())
    {
      // note: this could be a method of TMatrix as well, if TMatrix knew its
      // actives
      Dscal(this->actives[b], factor, this->blocks[b]->GetEntries());
    }
    else
      ErrThrow("scaling actives where the number of actives has not been set");
  }
}

/* ************************************************************************* */
size_t BlockFEMatrix::get_n_row_actives( size_t index ) const
{
  // check if there is as many block rows
  if(index > block_pattern->n_rows())
  {
    ErrThrow("Index out of bounds.");
  }
  return ((FEMatrix&)block( index , 0)).GetTestSpace()->GetActiveBound();
}

/* ************************************************************************* */
size_t BlockFEMatrix::get_n_column_actives( size_t index ) const
{
  // check if there is as many block rows
  if(index > block_pattern->n_cols())
  {
    ErrThrow("Index out of bounds.");
  }
  return ((FEMatrix&)block( 0, index)).GetTestSpace()->GetActiveBound();
}
