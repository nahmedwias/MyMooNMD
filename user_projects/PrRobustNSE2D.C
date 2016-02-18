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
                   (char*)"projection velocity", example.get_bc(2),
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
void matrix_multiplication(double ***inputMat, int *nrowInput, int *ncolInput, 
                                double ***outputMat, int *nrowOut, int *ncolOut, 
                                double **inputrhs,  int *ndimInput, 
                                double **outputrhs, int *ndimOutput)
{
  /** matrix-vecotr manipulation*/
  auto matrix_vector_multiply = [] (std::pair<int,int> size, double ***matrix, 
                                    double *rhs, double **product)
  {
    unsigned int nrow = size.first;
    unsigned int ncol = size.second;
    
    for(unsigned int i=0; i<nrow; i++)
    {
      double temp = 0;
      double temp1 = 0;
      for(unsigned int j=0; j<ncol; j++)
      {
        temp += matrix[0][i][j] * rhs[j];
        temp1 += matrix[1][i][j] * rhs[j];
      }
      product[0][i] = temp;
      product[1][i] = temp1;
    }
    return;
  };
  // prepare the output right hand 
  std::pair<int,int> size (nrowInput[1], ncolInput[1]);
  
  double **matrix[2];
  matrix[0] = inputMat[1];
  matrix[1] = inputMat[2];  
  matrix_vector_multiply(size, matrix, inputrhs[0], outputrhs);
}
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
    TFEFunction2D *fe_functions[3] =  { s_base.u.GetComponent(0), 
                                        s_base.u.GetComponent(1), 
                                        &s_base.p };
      // local assemble      
    LocalAssembling2D larhs(RECONSTR_GALERKIN_Rhs, fe_functions,
                              this->NSE2D::get_example().get_coeffs());    
    // reset the right hand side
    s_derived.rhsXh.reset(); 
    unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
    // prepare everything which is needed for assembling matrices and 
    // rhs: 
    const int nFESpaces = 2;
    const TFESpace2D *prooject_space=&this->get_projection_space();
    const TFESpace2D *velocity_space=&this->NSE2D::get_velocity_space();
    const TFESpace2D *pressure_space=&this->NSE2D::get_pressure_space();
    // pointers to the fespace
    const TFESpace2D * pointer_to_space[2] = {velocity_space, prooject_space};
    
    const int nSqMatAssemble = 2;
    // const int nReMatAssemble = 0;
    
    const int nRhsAssemble = 1;
    //NOTE:for the assembling of matrices or right hand side, space numbers
    //are passed as array: this corresponds to the array "pointer_to_space"
    std::vector<int> rowSpace ={0,0,1}; // row space for assembling matrices
    std::vector<int> colSpace ={1,1,1}; // cols space for assembling matrices
    std::vector<int> rowSpaceRhs={1}; // row space for assembling rhs 
    
    // prepare everything for storing the matrices and rhs
    // this will be used latter inside the class
    const int nSqMatStored = 0;
    const int nReMatStored = 0;
    int nRhsStored = 2;
    const TFESpace2D * pointToRhsStored[3]={velocity_space,velocity_space,
                                                 pressure_space};
    s_base.rhs.reset();
    double *rhsStored[3] = {s_base.rhs.block(0), s_base.rhs.block(1), 
                            s_base.rhs.block(2)};
    // boundary conditions    
    BoundCondFunct2D *bc[2] ={ prooject_space->GetBoundCondition(), 
                             pressure_space->GetBoundCondition() };
    // boundary values
    BoundValueFunct2D * const * const BoundValue=this->NSE2D::get_example().get_bd();
    BoundValueFunct2D *bv[2] = {BoundValue[0], BoundValue[1] };
    // assemble the right hand side
    Assemble2D_VectFE(nFESpaces, pointer_to_space, 
                      nSqMatAssemble,rowSpace,colSpace, 
                      nRhsAssemble, rowSpaceRhs, 
                      nSqMatStored, nullptr, 
                      nReMatStored, nullptr, 
                      nRhsStored, rhsStored,pointToRhsStored, 
                      larhs, nullptr, MatVectMult,
                      projection_matrices, bc, bv);
    // reset the space
    TDatabase::ParamDB->VELOCITY_SPACE = vel_space;     
  }  
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
