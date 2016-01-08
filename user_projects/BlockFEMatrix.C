/** ************************************************************************ 
* @brief     source file for BlockMatrixNSE2D
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History 
* ************************************************************************  */

#include <BlockFEMatrix.h>
#include <FEDatabase2D.h>
#include <stdio.h>
#include <string.h>

BlockFEMatrix::BlockFEMatrix(const TFESpace2D* testSpace, const TFESpace2D* ansatzSpace,
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
  std::shared_ptr<TStructure> structure(new TStructure(PressureSpace->GetN_DegreesOfFreedom(),
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
    TBaseFunct2D *baseFunctAnsatz = elementAnsatz->GetBaseFunct2D(); // basis funct obj
    int baseVectDim = baseFunctAnsatz->GetBaseVectDim(); // dimension of baseFunctAnsatz
    int nDofAnsatz = baseFunctAnsatz->GetDimension(); // no of basis function
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
        int e = elementAnsatz->GetFEDesc2D()->GetJointOfThisDOF(k);
        int n = 1; // normal orientation
        if(e != -1)
        {
          // n = cell->GetNormalOrientation(e);
        }
        // fill the correspoding blocks
        this->block(0)->set(testSpace->GetGlobalDOF(c)[j],
               ansatzSpace->GetGlobalDOF(c)[k], FunctionalValuesx[k]/**n*/);
        this->block(1)->set(testSpace->GetGlobalDOF(c)[j],
                            ansatzSpace->GetGlobalDOF(c)[k], FunctionalValuesy[k]/**n*/);
      }
    }
    // this->block(0)->PrintFull();
  }
  // this->get_combined_matrix()->PrintFull();
}

double BlockFEMatrix::addScaleActive(const double* x, double* y, double factor)
{  
}

