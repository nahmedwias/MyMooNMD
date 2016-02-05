#include <PrRobustNSE2D.h>
#include <Assemble2D.h>
#include <NSE2D.h>
#include <Domain.h>
#include <Database.h>
/** ************************************************************************ */
PrRobustNSE2D::SystemPerGrid::SystemPerGrid(const Example_NSE2D& example, 
TCollection& coll, unsigned int order, const TFESpace2D& space, 
                   const TFESpace2D& presSpace)
: projection_space(space.GetCollection(), (char*)"u", 
                   (char*)"projection velocity", example.get_bc(0),
                   TDatabase::ParamDB->PROJECTION_SPACE, nullptr), 
 ProjectionMatrix({&space}, {&projection_space, &projection_space, &presSpace}),
 rhsXh(projection_space.GetN_DegreesOfFreedom()) 
{
  // nothing to do here
  ProjectionMatrix.print_coloring_pattern("P");
  ProjectionMatrix = BlockFEMatrix::Projection_NSE2D(space, 
                                     projection_space, presSpace);
}
/** ************************************************************************ */
PrRobustNSE2D::PrRobustNSE2D(const TDomain& domain, 
       const Example_NSE2D& _example, unsigned int reference_id)
 : NSE2D(domain, _example, reference_id), Systems()
{
  // create the collection of cells from the domain 
  TCollection *coll = domain.GetCollection(It_Finest,0,reference_id);

  int proj_order = TDatabase::ParamDB->PROJECTION_METHOD;
  
  this->Systems.emplace_back(example, *coll, proj_order, 
                             this->NSE2D::get_velocity_space(), 
                             this->NSE2D::get_pressure_space());
  Output::print<1>("dof of RT/BDM : ",  setw(10), 
                   this->get_projection_space().GetN_DegreesOfFreedom());
  Output::print<1>("active dof    : ",  setw(10), 
                   this->get_projection_space().GetActiveBound());
}
/** ************************************************************************ */
void PrRobustNSE2D::assembleMatrixRhs()
{
  // TODO: give proper names to the functions
  // assembling the Navier-Stokes matrices and the right hand side
  this->NSE2D::assemble(); 
  
  // check if the reconstruction is not used then return
  if(TDatabase::ParamDB->DISCTYPE != RECONSTRUCTION)
    return;
  // in the case of DISCTYPE== RECONSTRUCTION
  // assemble the right hand side and projection matrix
  System_per_grid& s_base=this->NSE2D::systems.front();
  SystemPerGrid &s_derived = this->Systems.front();
  {
    TFEFunction2D *fe_functions[2] =  { s_base.u.GetComponent(0), 
                                        s_base.u.GetComponent(1) };
      // local assemble      
    LocalAssembling2D larhs(RECONSTR_GALERKIN_Rhs, fe_functions,
                              this->NSE2D::get_example().get_coeffs());
    // load the projection space 
    const TFESpace2D* space = &this->get_projection_space();
    // boundary conditions
    BoundCondFunct2D *bc[2] ={ space->GetBoundCondition(), 
                               space->GetBoundCondition() };
    // boundary values
    BoundValueFunct2D * const * const BoundValue 
       = this->NSE2D::get_example().get_bd();
    BoundValueFunct2D *bv[2] = {BoundValue[2], BoundValue[2] };
    // pointer to the space 
    const TFESpace2D * pointer_to_space[1] = { space };
    // reset the right hand side
    s_derived.rhsXh.reset(); 
    unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
    // rhs blocks 
    double *rhs_blocks[1] = { s_derived.rhsXh.get_entries() };
    
    // call the assemble function
    // assemble right hand side only
    Assemble2D_VectFE(1, pointer_to_space, 0, nullptr, 0, 
                      nullptr, 1, rhs_blocks, pointer_to_space, larhs, 
                      bc, bv);
    TDatabase::ParamDB->VELOCITY_SPACE = vel_space;    
  }
  // assemble the ProjectionMatrix  
  this->assembleProjectionMatrix();
  
  // update the right hand side with reconstrucion
  std::vector<std::shared_ptr<FEMatrix>> blocks 
         = s_derived.ProjectionMatrix.get_blocks_uniquely();
  // reset rhs 
  s_base.rhs.reset();
  // multiply the correspoding blocks with the right hand side 
  // and add it to the rhs of the base class
  blocks.at(0)->multiplyActive(s_derived.rhsXh.get_entries(), 
                               s_base.rhs.block(0));
  blocks.at(1)->multiplyActive(s_derived.rhsXh.get_entries(), 
                               s_base.rhs.block(1));
}

/** ************************************************************************ */
void PrRobustNSE2D::assembleProjectionMatrix()
{
  SystemPerGrid &s_derived = this->Systems.front();
  
  std::vector<std::shared_ptr<FEMatrix>> blocks 
         = s_derived.ProjectionMatrix.get_blocks_uniquely();
  std::vector<std::shared_ptr<FEMatrix>> matrices= {blocks.at(0), blocks.at(1)};
  const TFESpace2D *testSpace = &this->NSE2D::get_velocity_space();
  const TFESpace2D *ansatzSpace = &s_derived.projection_space;
  
  TCollection *coll = testSpace->GetCollection();
  unsigned int nCells = coll->GetN_Cells();

  for(unsigned int i=0; i<nCells; ++i)
    coll->GetCell(i)->SetNormalOrientation();
  
  int nPoints;
  double *xi, *eta;
  
  // loop over all cells
  for(unsigned int c=0; c<nCells; ++c)
  {
    TBaseCell *cell = coll->GetCell(c);
    TFE2D *eleAnsatz = 
      TFEDatabase2D::GetFE2D(ansatzSpace->GetFE2D(c,cell));
    TBaseFunct2D *baseFunctAnsatz = eleAnsatz->GetBaseFunct2D();
    int baseVectDim = baseFunctAnsatz->GetBaseVectDim();
    int nDofAnsatz = baseFunctAnsatz->GetDimension();
    // compute points on the reference element
    TNodalFunctional2D *nf = eleAnsatz->GetNodalFunctional2D();
    nf->GetPointsForAll(nPoints, xi, eta);
    
    // everything needed for the velocity space (test)
    TFE2D *eleTest = 
      TFEDatabase2D::GetFE2D(testSpace->GetFE2D(c,cell));
    // basis function for test space
    TBaseFunct2D *baseFunctTest = eleTest->GetBaseFunct2D();
    int nDofTest = baseFunctTest->GetDimension();
    
    // id for the reference transformation
    RefTrans2D refTransfID = eleAnsatz->GetRefTransID();
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
      if(testSpace->GetGlobalDOF(c)[j] >= 
         testSpace->GetN_ActiveDegrees() )
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
        matrices[0]->set(testSpace->GetGlobalDOF(c)[j],
               ansatzSpace->GetGlobalDOF(c)[k], FunctionalValuesx[k]);
        matrices[1]->set(testSpace->GetGlobalDOF(c)[j],
               ansatzSpace->GetGlobalDOF(c)[k], FunctionalValuesy[k]);
      }
    }// endfor j<nDofTest
  }// endfor c<nCells
  // delete coll;
}
/** ************************************************************************ */
void PrRobustNSE2D::solve()
{
  this->NSE2D::solve();
}

/** ************************************************************************ */
void PrRobustNSE2D::output(int i)
{
  this->NSE2D::output(i);
}
