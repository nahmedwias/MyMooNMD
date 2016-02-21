#include <PrRobustTime_NSE2D.h>
#include <Database.h>
#include <Assemble2D.h>
#include <DirectSolver.h>

PrRobustTime_NSE2D::SystemPerGrid::SystemPerGrid(const Example_NSE2D& example,
      TCollection& coll, const TFESpace2D& velocity_space, 
      const TFESpace2D& pressure_space, int order)
: projection_space(velocity_space.GetCollection(), (char*)"u", 
                    (char*)"projection space", example.get_bc(2), order, nullptr),
   ProjectionMatrix({&velocity_space}, {&projection_space, &projection_space, &pressure_space}),
   MassMatrix({&projection_space, &projection_space}),
   rhsXh(projection_space.GetN_DegreesOfFreedom()), 
   Modifed_Mass({&velocity_space, &velocity_space}),
   solutionVV(rhsXh),
   fefctVV(&projection_space, (char*)"uv", (char*)"uv", solutionVV.block(0), 
           solutionVV.length(0))
{
  ProjectionMatrix.print_coloring_pattern("P");
  Modifed_Mass.print_coloring_pattern("MM");
  MassMatrix.print_coloring_pattern("M");
  
  // projection matrix
  ProjectionMatrix = BlockFEMatrix::Projection_NSE2D(velocity_space, 
                                     projection_space, pressure_space);
  // MassMatrix of darcy type
  MassMatrix = BlockFEMatrix::Mass_NSE2D(projection_space);
  
  // modified mass matrix
  Modifed_Mass = BlockFEMatrix::Modified_MassMatrix(velocity_space);    
}

/**************************************************************************** */
PrRobustTime_NSE2D::PrRobustTime_NSE2D(const TDomain& domain, 
                                 Example_NSE2D& _example, int reference_id)
: Time_NSE2D(domain, _example, reference_id)
{
  TCollection *coll = domain.GetCollection(It_Finest,0,reference_id);
  
  int proj_order = TDatabase::ParamDB->PROJECTION_SPACE;
  
  this->Systems.emplace_back(_example, *coll, 
                             this->Time_NSE2D::get_velocity_space(), 
                             this->Time_NSE2D::get_pressure_space(),
                             proj_order);
  
  Output::print<1>("dof of RT/BDM : ",  setw(10), 
                   this->get_projection_space().GetN_DegreesOfFreedom());
  Output::print<1>("active dof    : ",  setw(10), 
                   this->get_projection_space().GetActiveBound());
  
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE != 5 
    || TDatabase::ParamDB->SOLVER_TYPE != 1)
  {
    Output::print<1>("Multigrid is not working correctly", 
                     "One can not easily mixed it up");
    return;
  }
}
void matrices_reconstruction(double ***inputMat, int *nrowInput, int *ncolInput, 
                                double ***outputMat, int *nrowOut, int *ncolOut, 
                                double **inputrhs,  int *ndimInput, 
                                double **outputrhs, int *ndimOutput)
{
  // using dgemm
  // inputMat[0] = M; inputMat[1]=P0; inputMat[2]=P1;
 /* double *Mass = new double[nrowInput[0]*ncolInput[0]];
  double *P0 =  new double[nrowInput[1]*ncolInput[1]];
  double *P1 =  new double[nrowInput[2]*ncolInput[2]];
  // copy the matrices
  for(int i=0; i<nrowInput[0]; ++i)
  {
    memcpy(Mass+i*ncolInput[0], inputMat[0][i], ncolInput[0]*SizeOfDouble);
  }
  for(int i=0; i<nrowInput[1]; ++i)
  {
    memcpy(P0+i*ncolInput[1], inputMat[1][i], ncolInput[1]*SizeOfDouble);
    memcpy(P1+i*ncolInput[2], inputMat[2][i], ncolInput[2]*SizeOfDouble);
  }
  int lda, ldb, ldc, m, n,k;
  char P0T[] = 'T';
  char P1T[] = 'T';
  
  lda = ncolInput[0];
  ldb = ncolInput[1];
  
  double *MP0T = new double[];
  */
 /** this lambda function compute the ABC^T
  * @param: matA projection matrix first velocity part (special case)
  * @param: matB mass matrix for BDM or RT which is multiplied in between
  * @param: matC projection matrix which multiplies as transpose from right 
  *         only for the off diagonal blocks (flag==1 & flag ==2)
  * @param: product the resulting matrix
  * @param: flag 0 and 3 are the cases when matA^T is multiplied from right
  *              1 and 2 when the matC^T is multiplied from right
  */
 auto compute_matrix_matrix = [](std::pair<int, int> rows, std::pair<int,int>cols, 
                                   double **matA, double **matB, double **matC, 
                                   double **&product, int flag)
 {
   unsigned int nRowA = rows.first;
   unsigned int nRowB = rows.second;
   unsigned int nColA = cols.first;
   unsigned int nColB = cols.second;
   
   if(nColB != nColA) // nColA is the nRowA of matA
   {
     ErrThrow("dimension mismatch ", nColB, " " , nColA);
   }
   // transpose of matrix A
   int nRowProduct = nColA; // column of A; 
   //int nColProduct = nRowA; // rows of A
   
   for(unsigned int i=0; i<nRowA; i++)
   {
     for(unsigned int j=0; j<nRowA; j++)
     {
       double temp = 0;
       for(unsigned int k=0; k<nRowB; k++)
       {
         for(unsigned int l=0; l<nRowProduct; l++)
         {
           if(flag == 0 || flag == 3)
             temp += matA[i][k] * matB[k][l] * matA[j][l];
           else
             temp += matA[i][k] * matB[k][l] * matC[j][l];
         }
       }
       product[i][j] = temp;
     }
   }
   return;
 };
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
  // prepare the outputMat[0] = P0 * M * P0^T
  std::pair<int, int> rowsArray(nrowInput[0], nrowInput[2]);
  std::pair<int, int> colsArray(ncolInput[0], ncolInput[2]);
  compute_matrix_matrix(rowsArray, colsArray, inputMat[0], inputMat[2], 
                          nullptr, outputMat[0], 0);
  // prepare the outputMat[0] = P0 * M * P1^T
  compute_matrix_matrix(rowsArray, colsArray, inputMat[0], inputMat[2], 
                          inputMat[1], outputMat[1], 1);
  // prepare the outputMat[0] = P1 * M * P0^T
  rowsArray.first=nrowInput[1]; rowsArray.second=nrowInput[2];
  colsArray.first=ncolInput[1]; rowsArray.second=ncolInput[2];
  compute_matrix_matrix(rowsArray, colsArray, inputMat[1], inputMat[2], 
                          inputMat[0], outputMat[2], 2);
  // prepare the outputMat[0] = P1 * M * P1^T
  compute_matrix_matrix(rowsArray, colsArray, inputMat[1], inputMat[2],
                          nullptr, outputMat[3], 3);
  
  // prepare the output right hand 
  std::pair<int,int> size (nrowInput[0], ncolInput[0]);
  double **matrix[2];
  matrix[0] = inputMat[0];
  matrix[1] = inputMat[1];
  matrix_vector_multiply(size, matrix, inputrhs[0], outputrhs);
  
  /** local matrices P0 and P1 needs to be globally assembled
   * so copy them to the outputMat here 
   */
  outputMat[4] = inputMat[0];
  outputMat[5] = inputMat[1];
  /*
  for(int i=0; i<size.first; i++)
    cout<<i<<"\t" << outputrhs[0][i]<<"\t" << outputrhs[1][i]<<"\n";  
    
  for(int i=0; i<nrowInput[1]; i++)
  {
    for(int j=0; j<nrowInput[1]; j++)
    {
      cout << outputMat[0][i][j] << "\t"<< outputMat[1][i][j]<<"\t"<<
              outputMat[2][i][j] << "\t"<< outputMat[3][i][j] << "\n";
    }
    cout <<endl;
  }
  exit(-4711);
  */  
}

/**************************************************************************** */
bool PrRobustTime_NSE2D::assembleProjMat()
{
  if(TDatabase::ParamDB->DISCTYPE != RECONSTRUCTION)
    return false;
  
  // setting the space for sign in GetSignOfThisDOF();
  unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
  TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
  // prepare everything to assemble the global projection matrix 
  SystemPerGrid& s_derived = this->Systems.front();
  
  const int nFESpaces = 2;
  const TFESpace2D *project_space  = &this->get_projection_space();
  const TFESpace2D *velocity_space = &this->Time_NSE2D::get_velocity_space();
  // pointers to the fespace
  const TFESpace2D * pointer_to_space[2] = { velocity_space, project_space };
  // assemble all square and rectangular matrices
  // nAllMatAssemble = nSquareMatrice + nRectMatrices;
  const int nAllMatAssemble = 2;
  //NOTE:for the assembling of matrices or right hand side, space numbers
  //are passed as array: this corresponds to the array "pointer_to_space"
  std::vector<int> rowSpace ={0, 0}; // row space for assembling matrices
  std::vector<int> colSpace ={1, 1}; // cols space for assembling matrices
  
  const int nRhsAssemble = 0;
  std::vector<int> rowSpaceRhs={};
  const int nSqMatStored = 0;
  TSquareMatrix2D** sqMatricesStored = nullptr;
  
  int nRhsStored = 0;
  double **rhsStored = nullptr;
  const TFESpace2D ** pointToRhsStored=nullptr;
  
  // rectangular matrices to be stored
  const int nReMatStored = 2;
  std::vector<std::shared_ptr<FEMatrix>> proj_matrices 
     = s_derived.ProjectionMatrix.get_blocks_uniquely();
  TMatrix2D *rectMatricesStored[2];
  rectMatricesStored[0]=reinterpret_cast<TMatrix2D*>(proj_matrices.at(0).get());
  rectMatricesStored[1]=reinterpret_cast<TMatrix2D*>(proj_matrices.at(1).get());
  rectMatricesStored[2]=reinterpret_cast<TMatrix2D*>(proj_matrices.at(1).get());
  for(int i=0; i<nReMatStored; i++)
      rectMatricesStored[i]->reset();
  BoundCondFunct2D * boundCond[2] = {
      project_space->GetBoundCondition(), project_space->GetBoundCondition()};
  
  std::array<BoundValueFunct2D*, 3> boundVal;
  boundVal[0] = example.get_bd()[0];
  boundVal[1] = example.get_bd()[1];
  boundVal[2] = example.get_bd()[2];
  
  TFEFunction2D ** fefct = {};
  LocalAssembling2D empty_constructor(NO_LOCAL_ASSEMBLE, fefct, 
                                      this->Time_NSE2D::get_example().get_coeffs());
  
  Assemble2D_VectFE(nFESpaces, pointer_to_space, 
                      nAllMatAssemble, rowSpace, colSpace, 
                      nRhsAssemble, rowSpaceRhs, 
                      nSqMatStored, sqMatricesStored, 
                      nReMatStored, rectMatricesStored, 
                      nRhsStored, rhsStored,pointToRhsStored, 
                      empty_constructor, nullptr, nullptr,
                      projection_matrices, boundCond, boundVal.data());
  rectMatricesStored[0]->Print("M");exit(0);
  
  if(TDatabase::ParamDB->PROBLEM_TYPE ==5)
    return true;
  else 
    return false;
}

/**************************************************************************** */
void PrRobustTime_NSE2D::assemble_initial_time()
{
  // assemble matrices and right hand side
  // till here no modification of the right hand side 
  // and the mass matrix is done  
  this->Time_NSE2D::assemble_initial_time(); 
  // if the DISCTYPE is RECONSTRUCTION proceed 
  // further
  if(TDatabase::ParamDB->DISCTYPE !=RECONSTRUCTION)
  {
    Output::print<1>("Standard Galerkin Discretization");
    return;
  }
  else
    Output::print<1>("PrRobustTime_NSE2D is considered");
  
  System_per_grid& s_base = this->Time_NSE2D::systems.front();
  // assemble mass matrix using vector space  
  for(SystemPerGrid &s_derived : this->Systems)
  {
    TFEFunction2D *fefct[3] = {this->Time_NSE2D::systems.front().u.GetComponent(0), 
                               this->Time_NSE2D::systems.front().u.GetComponent(1), 
                               &this->Time_NSE2D::systems.front().p};
    // local assemble routine
    LocalAssembling2D la_mass_rhs_vspace(RECONSTR_MASS, fefct, 
                                   this->Time_NSE2D::get_example().get_coeffs());    
    // assemble the right hand side separately  
    s_derived.rhsXh.reset();
    // setting the space for sign in GetSignOfThisDOF();
    unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
    // prepare everything which is needed for assembling matrices and 
    // the right hand side
    const int nFESpaces = 2;
    const TFESpace2D *project_space  = &this->get_projection_space();
    const TFESpace2D *velocity_space = &this->Time_NSE2D::get_velocity_space();
    const TFESpace2D *pressure_space = &this->Time_NSE2D::get_pressure_space();
    // pointers to the fespace
    const TFESpace2D * pointer_to_space[2] = { velocity_space, project_space };
    
    const int nAllMatAssemble=3;
    // const int nReMatAssemble = 0;
    const int nRhsAssemble = 1;
    //NOTE:for the assembling of matrices or right hand side, space numbers
    //are passed as array: this corresponds to the array "pointer_to_space"
    std::vector<int> rowSpace ={0, 0, 1}; // row space for assembling matrices
    std::vector<int> colSpace ={1, 1, 1}; // cols space for assembling matrices
    std::vector<int> rowSpaceRhs={1}; // row space for assembling rhs 
    
    // prepare everything for storing the matrices and rhs
    // this will be used latter inside the class
    const int nSqMatStored = 4;
    TSquareMatrix2D * sqMatricesStored[4];
    // 
    std::vector<std::shared_ptr<FEMatrix>> modified_mass
           = s_derived.Modifed_Mass.get_blocks_uniquely();
    //
    sqMatricesStored[0]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(0).get());
    sqMatricesStored[1]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(1).get());
    sqMatricesStored[2]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(2).get());
    sqMatricesStored[3]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(3).get());
    // reset the matrices
    for(int i=0; i<nSqMatStored; i++)
            sqMatricesStored[i]->reset();
    // rectangular matrices to be stored
    const int nReMatStored = 2;
    std::vector<std::shared_ptr<FEMatrix>> proj_matrices 
      = s_derived.ProjectionMatrix.get_blocks_uniquely();
    TMatrix2D *rectMatricesStored[3];
    rectMatricesStored[0]=reinterpret_cast<TMatrix2D*>(proj_matrices.at(0).get());
    rectMatricesStored[1]=reinterpret_cast<TMatrix2D*>(proj_matrices.at(1).get());
    rectMatricesStored[2]=reinterpret_cast<TMatrix2D*>(proj_matrices.at(2).get());
    for(int i=0; i<nReMatStored; i++)
      rectMatricesStored[i]->reset();
    // right hand vector to be stored
    int nRhsStored = 2;
    const TFESpace2D * pointToRhsStored[3]={velocity_space,velocity_space,
                                                 pressure_space};
    double *rhsStored[3] = {s_base.rhs.block(0), s_base.rhs.block(1), 
                            s_base.rhs.block(2)};
    s_base.rhs.reset();
    // boundary conditions    
    BoundCondFunct2D *bc[3] ={ velocity_space->GetBoundCondition(), 
                             velocity_space->GetBoundCondition(),
                             pressure_space->GetBoundCondition() };
    // boundary values
    BoundValueFunct2D * const * const BoundValue=this->Time_NSE2D::get_example().get_bd();
    BoundValueFunct2D *bv[3] = {BoundValue[0], BoundValue[1], BoundValue[2] };
    // assemble the correspoding matrices and right hand side
    Assemble2D_VectFE(nFESpaces, pointer_to_space, 
                      nAllMatAssemble, rowSpace, colSpace, 
                      nRhsAssemble, rowSpaceRhs, 
                      nSqMatStored, sqMatricesStored, 
                      nReMatStored, rectMatricesStored, 
                      nRhsStored, rhsStored,pointToRhsStored, 
                      la_mass_rhs_vspace, 
                      matrices_reconstruction, nullptr,
                      projection_matrices, bc, bv);

    // commit this line when the line below is tested
    TDatabase::ParamDB->VELOCITY_SPACE = vel_space;
    // sqMatricesStored[0]->Print("M");exit(0);
    // this piece of code is to check the old routines
    /*{
      double *rhs_blocks[1] = { s_derived.rhsXh.get_entries() };
      
      std::vector<std::shared_ptr<FEMatrix>> mass_block
           = s_derived.MassMatrix.get_blocks_uniquely();
          
      TSquareMatrix2D * sq_matrices[1];
      sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(mass_block.at(0).get());
      
      sq_matrices[0]->reset();
       Assemble2D_VectFE(1, pointer_to_space, 1, sq_matrices, 0,
                        nullptr, 1, rhs_blocks, pointer_to_space,
                        la_mass_rhs_vspace,
                        bc, bv);
      
      
      TDatabase::ParamDB->VELOCITY_SPACE = vel_space;
      // compute the projection matrix
      this->assembleProjectionMatrix();
      
      std::vector<std::shared_ptr<FEMatrix>> projBlocks
           = s_derived.ProjectionMatrix.get_blocks_uniquely();
      std::vector<std::shared_ptr<FEMatrix>> massBlocks
             = s_derived.MassMatrix.get_blocks_uniquely();
      
      std::shared_ptr<const TMatrix> M00
       (projBlocks.at(0)->multiply_with_transpose_from_right(*massBlocks.at(0)));
      
      // M00->Print("M");exit(0);
      std::shared_ptr<TMatrix> A2T(projBlocks.at(1)->GetTransposed());  
      std::shared_ptr<TMatrix> BA2T(massBlocks.at(0)->multiply(*A2T));
      std::shared_ptr<TMatrix> M01(projBlocks.at(0)->multiply(*BA2T));  
      
      std::shared_ptr<TMatrix> A1T(projBlocks.at(0)->GetTransposed());
      std::shared_ptr<TMatrix> BA1T(massBlocks.at(0)->multiply(*A1T));
      std::shared_ptr<TMatrix> M10(projBlocks.at(1)->multiply(*BA1T));
      
      std::shared_ptr<TMatrix> M11
        (projBlocks.at(1)->multiply_with_transpose_from_right(*massBlocks.at(0)));
      // projBlocks.at(1)->PrintFull("M");exit(0);
      // s_derived.rhsXh.print("r");exit(0);
    
    }*/
  }// endof assembling for all grids
  
  // copy the modified rhs to the old_rhs_modified for the right hand 
  // side at the initial time
  this->old_rhs_modified = this->Time_NSE2D::systems.front().rhs;
  this->old_solution = s_base.solution;
}

/**************************************************************************** */
void PrRobustTime_NSE2D::assembleProjectionMatrix()
{  
  SystemPerGrid &s_derived = this->Systems.front();
  
  std::vector<std::shared_ptr<FEMatrix>> blocks 
         = s_derived.ProjectionMatrix.get_blocks_uniquely();
  std::vector<std::shared_ptr<FEMatrix>> matrices= {blocks.at(0), blocks.at(1)};
  const TFESpace2D *testSpace = &this->Time_NSE2D::get_velocity_space();
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

/**************************************************************************** */
void PrRobustTime_NSE2D::assemble_rhs()
{
  double tau =TDatabase::TimeDB->TIMESTEPLENGTH;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;
  
  // System_per_grid for the member access fromt the Time_NSE2D class
  System_per_grid& s_base = this->Time_NSE2D::systems.front();
  // SystemPerGrid for the member access from the current class
  SystemPerGrid& s_derived= this->Systems.front();
  
  if(TDatabase::ParamDB->DISCTYPE == RECONSTRUCTION)
  { 
    TFEFunction2D *fefct[3] = {s_base.u.GetComponent(0), s_base.u.GetComponent(1), &s_base.p};
    // assembling the right hand side at current time step     
    LocalAssembling2D larhs(RECONSTR_GALERKIN_Rhs, fefct,
                          this->Time_NSE2D::get_example().get_coeffs());
    unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
    // assembling the right hand side
    // for this one have to assemble the projection matrices as well
    // FIXME: try to assemble the projection matrices separately in order 
    // to assemble the right hand side for the pressure robust method
    // prepare everything needed for the right hand side assembling
    s_derived.rhsXh.reset();
    s_base.rhs.reset();
    const int nFESpaces = 2;
    const TFESpace2D *prooject_space=&this->get_projection_space();
    const TFESpace2D *velocity_space=&this->Time_NSE2D::get_velocity_space();
    const TFESpace2D *pressure_space=&this->Time_NSE2D::get_pressure_space();
    // pointers to the fespace
    const TFESpace2D * pointer_to_space[2] = {velocity_space, prooject_space};
    // matrices to assemble in order to assemble the right hand 
    // side of the method
    const int nAllMatAssemble = 3;
    // const int nReMatAssemble = 0;
    
    const int nRhsAssemble = 1;    
    //NOTE:for the assembling of matrices or right hand side, space numbers
    //are passed as array: this corresponds to the array "pointer_to_space"
    std::vector<int> rowSpace ={0, 0, 1}; // row space for assembling matrices
    std::vector<int> colSpace ={1, 1, 1}; // cols space for assembling matrices
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
    BoundCondFunct2D *bc[3] ={ velocity_space->GetBoundCondition(), 
                             velocity_space->GetBoundCondition(), 
                             pressure_space->GetBoundCondition() };
    // boundary values
    BoundValueFunct2D * const * const BoundValue=this->Time_NSE2D::get_example().get_bd();
    BoundValueFunct2D *bv[3] = {BoundValue[0], BoundValue[1], BoundValue[2] };
    // assemble the right hand side
    Assemble2D_VectFE(nFESpaces, pointer_to_space, 
                      nAllMatAssemble,rowSpace,colSpace, 
                      nRhsAssemble, rowSpaceRhs, 
                      nSqMatStored, nullptr, 
                      nReMatStored, nullptr, 
                      nRhsStored, rhsStored,pointToRhsStored, 
                      larhs, nullptr, MatVectMult,
                      projection_matrices, bc, bv);
    TDatabase::ParamDB->VELOCITY_SPACE = vel_space;
    // copy nonactive
    s_base.solution.copy_nonactive(s_base.rhs);
    // assemble the right hand side vector 
    
    // scale the current right hand side with tau*theta4
    s_base.rhs.scaleActive(tau*theta4);
    // add right hand side from the previous time step, which 
    // is scaled by theta3*tau
    s_base.rhs.addScaledActive(this->old_rhs_modified, tau*theta3);
        
    // now it is this->systems[i].rhs = tau*theta3*f^{k-1} + tau*theta4*f^k
    // next we want to set old_rhs to f^k (to be used in the next time step)
    this->old_rhs_modified.addScaledActive(s_base.rhs, -1./(tau*theta3));
    this->old_rhs_modified.scaleActive(-theta3/theta4);
    // scale actives in each block
    double factor = -tau*theta2;
    s_base.matrix.scale_blocks_actives(factor, {{0,0},{0,1},{1,0},{1,1}});
    // add blocks from the mass matrix to the current matrix
    const FEMatrix& MFM00 = *s_derived.Modifed_Mass.get_blocks().at(0).get();    
    s_base.matrix.add_matrix_actives(MFM00, 1.0, {{0,0}}, {false});
    
    const FEMatrix& MFM01 = *s_derived.Modifed_Mass.get_blocks().at(1).get();    
    s_base.matrix.add_matrix_actives(MFM01, 1.0, {{0,1}}, {false});
    
    const FEMatrix& MFM10 = *s_derived.Modifed_Mass.get_blocks().at(2).get();    
    s_base.matrix.add_matrix_actives(MFM10, 1.0, {{1,0}}, {false});
    
    const FEMatrix& MFM11 = *s_derived.Modifed_Mass.get_blocks().at(3).get();    
    s_base.matrix.add_matrix_actives(MFM11, 1.0, {{1,1}}, {false});
    
    // multiplying matrix with the old_solution vector
    // FIXME FInd other solution than this submatrix method.
    s_base.matrix.apply_scaled_submatrix(old_solution, s_base.rhs, 2, 2, 1.0);
    
    // restore the matrix by subtracting the mass matrix and scale by 1./factor
    factor = 1./factor;
    s_base.matrix.add_matrix_actives(MFM00, -1.0, {{0,0}},{false});
    s_base.matrix.add_matrix_actives(MFM01, -1.0, {{0,1}},{false});
    s_base.matrix.add_matrix_actives(MFM10, -1.0, {{1,0}},{false});
    s_base.matrix.add_matrix_actives(MFM11, -1.0, {{1,1}},{false});
    // rescale now 
    s_base.matrix.scale_blocks_actives(factor, {{0,0}, {0,1}, {1, 0}, {1, 1}});
    
    // scaling the B blocks with current time step length
    // scale the BT blocks with time step length
    for(System_per_grid& s : this->Time_NSE2D::systems)
    {
      if(tau != oldtau)
      {
        // TODO: change the factor to be THETA1*tau;
        factor = /*TDatabase::TimeDB->THETA1**/tau;
        if(this->oldtau != 0.0)
        {
          factor /= this->oldtau;
          Output::print<1>("change in tau", this->oldtau, "->", tau);
        }
        // scale the BT transposed blocks with the current time step
        s.matrix.scale_blocks(factor, {{0,2}, {1,2}});      
        if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
        {
          s.matrix.scale_blocks(factor, {{2,0}, {2,1}});
        }
      }
    }
    this->Time_NSE2D::oldtau= tau;
  }
  else
  {
    this->Time_NSE2D::assemble_rhs();
  }
  s_base.rhs.copy_nonactive(s_base.solution);
  Output::print<3>("leaving the assemble_rhs PrRobustTime_NSE2D");
}

/**************************************************************************** */
void PrRobustTime_NSE2D::assemble_nonlinear()
{
  // check and return if the DISCTYPE==RECONSTRUCTION
  if(TDatabase::ParamDB->DISCTYPE != RECONSTRUCTION)
  {
    this->Time_NSE2D::assemble_nonlinear_term();
    return;
  }
  // assemble the nonlinear term  
  SystemPerGrid& s_derived = this->Systems.front();
  for(System_per_grid& s_base : this->Time_NSE2D::systems)
  {
    //NOTE: What to assemble:
    // pi_h u^old = projection_matrices*sol
    // assemble the Matrix A = (pi_h u^old grad u, pi_h v)
    // This will results in 
    // ModifiedA = (A00 A01)
    //             (A10 A11)
    // preparing pi_h u^old u in V_h, pi_h u in X_h
    // s_base.solution.print("r");
    s_derived.ProjectionMatrix.apply_scaled_transposed_submatrix(s_base.solution,
                               s_derived.solutionVV, 1, 2, 1.0);
    //s_derived.solutionVV.print("sv");exit(0);
    
  }
}

/**************************************************************************** */
void PrRobustTime_NSE2D::assemble_system_matrix()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;
  int i=0;
  if(TDatabase::ParamDB->DISCTYPE == RECONSTRUCTION)
  {
    for(System_per_grid& s_base : this->systems)
    {
      s_base.matrix.scale_blocks_actives(factor, {{0,0}, {0,1}, {1, 0}, {1, 1}});
      const FEMatrix& MFM00 = *this->Systems[i].Modifed_Mass.get_blocks().at(0).get();
      s_base.matrix.add_matrix_actives(MFM00, 1.0, {{0,0}}, {false});
      
      const FEMatrix& MFM01 = *this->Systems[i].Modifed_Mass.get_blocks().at(1).get();    
      s_base.matrix.add_matrix_actives(MFM01, 1.0, {{0,1}}, {false});
      
      const FEMatrix& MFM10 = *this->Systems[i].Modifed_Mass.get_blocks().at(2).get();    
      s_base.matrix.add_matrix_actives(MFM10, 1.0, {{1,0}}, {false});
      
      const FEMatrix& MFM11 = *this->Systems[i].Modifed_Mass.get_blocks().at(3).get();    
      s_base.matrix.add_matrix_actives(MFM11, 1.0, {{1,1}}, {false});
      i++;
    }
  }
  else
  {
    this->Time_NSE2D::assemble_system();
  }
  Output::print<3>("System matrix in PrRobustTime_NSE2D is assembled");
}

/**************************************************************************** */
void PrRobustTime_NSE2D::solve()
{
  System_per_grid& s_base = this->Time_NSE2D::systems.front();
  if(TDatabase::ParamDB->DISCTYPE == RECONSTRUCTION)
  {
    if((TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE !=5)
      || (TDatabase::ParamDB->SOLVER_TYPE !=1 ))
    {
      if(TDatabase::ParamDB->SOLVER_TYPE != 2)
        ErrThrow("only the direct solver is supported currently");
     
      /// @todo consider storing an object of DirectSolver in this class
      DirectSolver direct_solver(s_base.matrix, 
                                 DirectSolver::DirectSolverTypes::umfpack);
      direct_solver.solve(s_base.rhs, s_base.solution);
    }
    else
      ErrThrow("Multigrid is not tested yet");
  
    this->descaleMatrices();
  }
  else
    this->Time_NSE2D::solve();
  
  this->old_solution = s_base.solution;
  
  Output::print<3>("System matrix in PrRobustTime_NSE2D is solved");
}

/**************************************************************************** */
void PrRobustTime_NSE2D::descaleMatrices()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;  
  int i=0;
  for(System_per_grid& s : this->systems)
  {
    const FEMatrix& MFM00 = *this->Systems[i].Modifed_Mass.get_blocks().at(0).get();
    const FEMatrix& MFM01 = *this->Systems[i].Modifed_Mass.get_blocks().at(1).get();
    const FEMatrix& MFM10 = *this->Systems[i].Modifed_Mass.get_blocks().at(2).get();
    const FEMatrix& MFM11 = *this->Systems[i].Modifed_Mass.get_blocks().at(3).get();
    
    s.matrix.add_matrix_actives(MFM00, -1.0, {{0,0}}, {false});
    s.matrix.add_matrix_actives(MFM01, -1.0, {{0,1}}, {false});
    s.matrix.add_matrix_actives(MFM10, -1.0, {{1,0}}, {false});
    s.matrix.add_matrix_actives(MFM11, -1.0, {{1,1}}, {false});
    s.matrix.scale_blocks_actives(1./factor, {{0,0}, {0,1}, {1, 0}, {1, 1}});
    i++;
  }  
}

/**************************************************************************** */
void PrRobustTime_NSE2D::output(int m, int &image)
{
  this->Time_NSE2D::output(m, image);
}




















