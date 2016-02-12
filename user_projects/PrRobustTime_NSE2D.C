#include <PrRobustTime_NSE2D.h>
#include <Database.h>
#include <Assemble2D.h>

PrRobustTime_NSE2D::SystemPerGrid::SystemPerGrid(const Example_NSE2D& example,
      TCollection& coll, const TFESpace2D& velocity_space, 
      const TFESpace2D& pressure_space, int order)
: projection_space(velocity_space.GetCollection(), (char*)"u", 
                    (char*)"projection space", example.get_bc(0), order, nullptr),
   ProjectionMatrix({&velocity_space}, {&projection_space, &projection_space, &pressure_space}),
   MassMatrix({&projection_space, &projection_space}),
   rhsXh(projection_space.GetN_DegreesOfFreedom()), 
   Modifed_Mass({&velocity_space, &velocity_space})
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
 
  auto computeTransposed = [](int nrow, int ncol, 
                              double **mat)
  {
    double **matT= new double*[ncol];
    for(unsigned int i=0; i<ncol; ++i)
      matT[i] = new double[nrow];
    for(int i=0; i<nrow; i++)
    {
      for (int j=0; j<ncol; j++)
        matT[j][i] = mat[i][j];
    }
    return matT;
  };
  
  auto multiply = [](int *rows, int *cols, 
                     double **matA, double **matB)
  {
    int nrowA = rows[0];
    int ncolA = cols[0];
    int nrowB = rows[1];
    int ncolB = cols[1];
    
    double **product = new double*[nrowA];
    for(unsigned int i=0; i<nrowA; ++i)
      product[i] = new double[ncolB];
    
    for(unsigned int i=0; i<nrowA; ++i)
    {
      for(unsigned int j=0; j<ncolB; ++j)
      {
        product[i][j] = 0;
        for(unsigned int k=0; k<nrowB; ++k)
        {
          product[i][j] += matA[i][k] * matB[k][j];
        }
      }
    }
    return product;
  };
  // compute the transpose of the projection matrices
  double **matP0T = computeTransposed(nrowInput[1], ncolInput[1], inputMat[1]);
  double **matP1T = computeTransposed(nrowInput[2], ncolInput[2], inputMat[2]);
  // compute the product matrix-matrix-transpose, 
  int rows[2], cols[2];
  rows[0] = nrowInput[0]; rows[1] = ncolInput[1];
  cols[0] = ncolInput[0]; cols[1] = nrowInput[1];  
  double **matM_matP0T = multiply(rows, cols, inputMat[0], matP0T);// M*P0^T
  // compute the product matrix-matrix transpose
  rows[0] = nrowInput[0]; rows[1] = ncolInput[2];
  cols[0] = ncolInput[0]; cols[1] = nrowInput[2];  
  double **matM_matP1T = multiply(rows, cols, inputMat[0], matP1T); // M*P1^T
  
  // prepare the outputMat[0] = P0 * M * P0^T
  rows[0] = nrowInput[1]; rows[1] = nrowInput[0];
  cols[0] = ncolInput[1]; cols[1] = nrowInput[1];
  outputMat[0] = multiply(rows, cols, inputMat[1], matM_matP0T);
    
  // prepare the outputMat[1] = P0 * M * P1^T
  rows[0] = nrowInput[1]; rows[1] = nrowInput[0];
  cols[0] = ncolInput[1]; cols[1] = nrowInput[2];
  outputMat[1] = multiply(rows, cols, inputMat[1], matM_matP1T);
  
  // prepare the outputMat[2] = P1 * M * P0^T
  rows[0] = nrowInput[2]; rows[1] = nrowInput[0];
  cols[0] = ncolInput[2]; cols[1] = nrowInput[1];
  outputMat[2] = multiply(rows, cols, inputMat[2], matM_matP0T);
  
  // prepare the outputMat[2] = P1 * M * P1^T
  rows[0] = nrowInput[2]; rows[1] = nrowInput[0];
  cols[0] = ncolInput[2]; cols[1] = nrowInput[1];
  outputMat[3] = multiply(rows, cols, inputMat[2], matM_matP1T);
  
  /** matrix-vecotr manipulation*/
  auto matrix_vector_multiply = [] (double **matrix, double *rhs)
  {
    
  };
  
  /*
  cout<<rows[0] << "  " << cols[1] << endl;
  for(int i=0; i<rows[0]; i++)
  {
    for(int j=0; j<cols[1]; j++)
    {
      cout << outputMat[1][i][j] << "\t";
    }
    cout <<endl;
  }*/
}
void projection_matrices(int current_cell, const TFESpace2D* ansatzSpace, 
                         const TFESpace2D* testSpace, double ***locMatrix)
{
  // projection is the ansatzSpace
  // velocity is the test space
  int i,j, N_Rows, N_Columns;
  double **CurrentMatrix, *MatrixRow;
  
  TCollection *coll = testSpace->GetCollection();
  TBaseCell *cell = coll->GetCell(current_cell);
  
  TFE2D *eleAnsatz = 
    TFEDatabase2D::GetFE2D(ansatzSpace->GetFE2D(current_cell,cell));
  TBaseFunct2D *baseFunctAnsatz = eleAnsatz->GetBaseFunct2D();
  int baseVectDim = baseFunctAnsatz->GetBaseVectDim();
  int nDofAnsatz = baseFunctAnsatz->GetDimension();
  // compute points on the reference element
  TNodalFunctional2D *nf = eleAnsatz->GetNodalFunctional2D();
  int nPoints;
  double *xi, *eta;
  nf->GetPointsForAll(nPoints, xi, eta);
  
  // everything needed for the velocity space (test)
  TFE2D *eleTest = 
    TFEDatabase2D::GetFE2D(testSpace->GetFE2D(current_cell,cell));
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
  
  for(i=1; i<=2; ++i)
  {
    CurrentMatrix = locMatrix[i];
    N_Rows = nDofTest;
    N_Columns = nDofAnsatz;
    
    for(j=0;j<N_Rows;j++)
    {
      MatrixRow = CurrentMatrix[j];
      memset(MatrixRow, 0, SizeOfDouble*N_Columns);
    } 
  } 
  
  double **MatrixP0 = locMatrix[1];
  double **MatrixP1 = locMatrix[2];
  
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
    for(int k=0; k<nPoints; ++k)
    {
      PointValuesx[k]         = uorig[k][j];
      PointValuesy[k+nPoints] = uorig[k][j];
    }
    nf->GetAllFunctionals(coll, cell, PointValuesx, FunctionalValuesx);
    nf->GetAllFunctionals(coll, cell, PointValuesy, FunctionalValuesy);
    
    for(int k=0; k<nDofAnsatz; k++)
    {
      MatrixP0[j][k] = FunctionalValuesx[k];
      MatrixP1[j][k] = FunctionalValuesy[k];      
    }
  }
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
  for(auto &s_derived: this->Systems)
  {
    TFEFunction2D *fefct[3] = {this->Time_NSE2D::systems.front().u.GetComponent(0), 
                               this->Time_NSE2D::systems.front().u.GetComponent(1), 
                               &this->Time_NSE2D::systems.front().p};
    // local assemble routine
    LocalAssembling2D la_mass_rhs_vspace(RECONSTR_MASS, fefct, 
                                   this->Time_NSE2D::get_example().get_coeffs());
    // boundary conditions
    BoundCondFunct2D *bc[2] = 
       { this->Time_NSE2D::systems.front().velocity_space.GetBoundCondition(),
         this->get_projection_space().GetBoundCondition() };
    
    // change the velocity space to assign the correct sign 
    unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
    
    BoundValueFunct2D * const * const BoundValue 
       = this->Time_NSE2D::get_example().get_bd();
    BoundValueFunct2D *bv[2] = {BoundValue[2], BoundValue[2] };
      
    // assemble the right hand side separately  
    s_derived.rhsXh.reset();
    // setting the space for sign in GetSignOfThisDOF();
    const TFESpace2D *pr_space = &this->get_projection_space();
    const TFESpace2D *v_space = &this->Time_NSE2D::get_velocity_space();
    const TFESpace2D * pointer_to_space[2] = { pr_space, v_space };
    const TFESpace2D * pointer_to_space_rhs[2] = { v_space, v_space };
    
    double *rhs_blocks[1] = { s_derived.rhsXh.get_entries() };
    
    std::vector<std::shared_ptr<FEMatrix>> mass_block
         = s_derived.MassMatrix.get_blocks_uniquely();
    // BDM/RT mass matrix
    FEMatrix M(&s_derived.projection_space, 
               &s_derived.projection_space);
    
    FEMatrix P0(&this->Time_NSE2D::get_velocity_space(), 
                &s_derived.projection_space);
    
    FEMatrix P1(&this->Time_NSE2D::get_velocity_space(), 
                &s_derived.projection_space);
    
    TSquareMatrix2D * sq_matrices_assem[3];
    sq_matrices_assem[0] = reinterpret_cast<TSquareMatrix2D*>(&M);
    sq_matrices_assem[1] = reinterpret_cast<TSquareMatrix2D*>(&P0);
    sq_matrices_assem[2] = reinterpret_cast<TSquareMatrix2D*>(&P1);
    
    std::vector<std::shared_ptr<FEMatrix>> modified_mass
           = s_derived.Modifed_Mass.get_blocks_uniquely();
    
    TSquareMatrix2D * sq_matrices_stored[4];
    sq_matrices_stored[0]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(0).get());
    sq_matrices_stored[1]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(1).get());
    sq_matrices_stored[2]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(2).get());
    sq_matrices_stored[3]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(3).get());
    
    double *rhs_blocks_stored[2] = { s_base.rhs.block(0), 
                        s_base.rhs.block(1) };
    
    TSquareMatrix2D * sq_matrices[1];
    sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(mass_block.at(0).get());
    
    sq_matrices[0]->reset();
    std::vector<int> row_space={0, 1, 1};
    std::vector<int> col_space={0, 0, 0};
    std::vector<int> row_space_rhs={0};
    
    Assemble2D_VectFE(2, pointer_to_space, 3, row_space, col_space, 1, 
                      row_space_rhs, 4, sq_matrices_stored, 0, nullptr, 
                      2, rhs_blocks_stored, pointer_to_space_rhs, 
                      la_mass_rhs_vspace, matrices_reconstruction, 
                      projection_matrices, bc, bv);
    exit(0);
    // assemble right hand side only
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
    
    std::shared_ptr<TMatrix> A2T(projBlocks.at(1)->GetTransposed());  
    std::shared_ptr<TMatrix> BA2T(massBlocks.at(0)->multiply(*A2T));
    std::shared_ptr<TMatrix> M01(projBlocks.at(0)->multiply(*BA2T));  
    
    std::shared_ptr<TMatrix> A1T(projBlocks.at(0)->GetTransposed());
    std::shared_ptr<TMatrix> BA1T(massBlocks.at(0)->multiply(*A1T));
    std::shared_ptr<TMatrix> M10(projBlocks.at(1)->multiply(*BA1T));
     
    std::shared_ptr<TMatrix> M11
     (projBlocks.at(1)->multiply_with_transpose_from_right(*massBlocks.at(0)));
    
    s_base.matrix.add_matrix_actives(reinterpret_cast<const FEMatrix&>(*M00), 
                                     1.0, {{0,0}}, {false});
    // multiply the projection matrix with mass matrix
    // s_derived.Modifed_Mass.multiply(s_derived.MassMatrix, s_derived.ProjectionMatrix);
    
    exit(0);
    /*
    // std::shared_ptr<BlockMatrix> Modifed_Mass(s.projMat.multiply(s.MassMatrix));
    
    s_derived.projMat.apply(s_derived.rhsXh,this->Time_NSE2D::systems.front().rhs);
    */
  //  /**
  //   * checking the new block matrix
  //   */
  //  s_derived.Modifed_Mass.multiply(s_derived.MassMatrix, s_derived.projMat);
  }
  
  //Output::print<1>("leaving the assemble_initial_time PrRobustTime_NSE2D");
  // copy the modified rhs to the old_rhs_modified for the right hand 
  // side at the initial time
  //this->old_rhs_modified = this->Time_NSE2D::systems.front().rhs;
  
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
{/*
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
    LocalAssembling2D la_rhs(TNSE2D_Rhs, fefct,
                           this->example.get_coeffs());
    // assemble the standard right hand side
    this->Time_NSE2D::AssembleRhs(la_rhs, s_base.rhs);
    
    s_base.solution.copy_nonactive(s_base.rhs);
    // assemble the right hand side vector 
    this->AssembleRHS_Only();
    
    // modified Mass Matrix
//     std::shared_ptr<BlockMatrixNSE2D::BlockMatrix> 
//       Modifed_Mass(s_derived.projMat.multiply(s_derived.MassMatrix));
    
    s_derived.projMat.apply(s_derived.rhsXh, this->Time_NSE2D::systems.front().rhs);
  
    // scale the current right hand side with tau*theta4
    s_base.rhs.scaleActive(tau*theta4);
    // add right hand side from the previous time step, which 
    // is scaled by theta3*tau
    s_base.rhs.addScaledActive(this->old_rhs_modified, tau*theta3);
    
    // now it is this->systems[i].rhs = tau*theta3*f^{k-1} + tau*theta4*f^k
    // next we want to set old_rhs to f^k (to be used in the next time step)
    this->old_rhs_modified.addScaledActive(s_base.rhs, -1./(tau*theta3));
    this->old_rhs_modified.scaleActive(-theta3/theta4);
    
    // from now on: multiply the M and A blocks with previous solution
    // and add to the right hand side
    // One can use "applyScaledAddActive" if the fespaces are known
    s_derived.Modifed_Mass.BlockMatrix::apply_scaled_add(s_base.solution,
                                                s_base.rhs, 1.0);
    
    s_base.matrix.applyScaledAddActive(s_base.solution.get_entries(), 
                                  s_base.rhs.get_entries(), -tau*theta2);
    s_base.rhs.copy_nonactive(s_base.solution);
    // scaling the B blocks with current time step length
    for(System_per_grid &s : this->Time_NSE2D::systems)
    {
      if(tau != this->Time_NSE2D::oldtau)
      {
        double factor = tau;
        
        if(this->Time_NSE2D::oldtau != 0.0)
        {
          factor /= this->oldtau;
          Output::print<1>("change in tau", this->oldtau, "->", tau);
        }
        s.matrix.get_BT_block(0)->operator *=(factor);
        s.matrix.get_BT_block(1)->operator *=(factor);
        if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT >0 )
        {
          s.matrix.get_B_block(0)->operator*=(factor);
          s.matrix.get_B_block(1)->operator*=(factor);
        }
        Output::print<1>("NOTE: scaling in the B-blocks have to be checked",
                      "PrRobustTime_NSE2D " );
      }
    }
    this->Time_NSE2D::oldtau= 0.0;
  }
  else
  {
    this->Time_NSE2D::assemble_rhs();
  }
  Output::print<3>("leaving the assemble_rhs PrRobustTime_NSE2D");*/
}

/**************************************************************************** */
void PrRobustTime_NSE2D::AssembleRHS_Only()
{
  /*TFEFunction2D *fe_functions[2] = 
     { this->Time_NSE2D::systems.front().u.GetComponent(0), 
       this->Time_NSE2D::systems.front().u.GetComponent(1) };
  // local assemble      
  LocalAssembling2D larhs(RECONSTR_GALERKIN_Rhs, fe_functions,
                          this->Time_NSE2D::get_example().get_coeffs());
  
  BoundCondFunct2D *bc[2] = { this->get_projection_space().GetBoundCondition(),
                              this->get_projection_space().GetBoundCondition() };
  
  BoundValueFunct2D * const * const BoundValue 
       = this->Time_NSE2D::get_example().get_bd();
  BoundValueFunct2D *bv[2] = {BoundValue[2], BoundValue[2] };
  
  const TFESpace2D *v_space = &this->get_projection_space();
  const TFESpace2D * pointer_to_space[1] = { v_space };
  
  // assemble the right hand side separately  
  this->Systems.front().rhsXh.reset();
  
  // setting the space for sign in GetSignOfThisDOF();
  unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
  TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
  
  double *rhs_blocks[1] = { this->Systems.front().rhsXh.get_entries() };
  
  // assemble right hand side only
  Assemble2D_VectFE(1, pointer_to_space, 0, nullptr, 0, 
                    nullptr, 1, rhs_blocks, pointer_to_space, larhs, 
                    bc, bv);
  TDatabase::ParamDB->VELOCITY_SPACE = vel_space;*/
}

/**************************************************************************** */
void PrRobustTime_NSE2D::assemble_system_matrix()
{
 /* double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;
  System_per_grid& s_base = this->Time_NSE2D::systems.front();    
  if(TDatabase::ParamDB->DISCTYPE == RECONSTRUCTION)
  {
    SystemPerGrid& s_derived = this->Systems.front();
    s_base.matrix.scaleActive(factor);
    s_base.matrix.addScaledActive(s_derived.Modifed_Mass,1.0);
  }
  else
  {
    this->Time_NSE2D::assemble_system();
    // s_base.matrix.scaleActive(factor);
    // s_base.matrix.addScaledActive(s_base.Mass_Matrix,1.0);    
  }
  Output::print<3>("System matrix in PrRobustTime_NSE2D is assembled");*/
}

/**************************************************************************** */
void PrRobustTime_NSE2D::solve()
{/*
  System_per_grid& s_base = this->Time_NSE2D::systems.front();

  if((TDatabase::ParamDB->SC_PRE_SMOOTH_SADDLE != 5)
    || (TDatabase::ParamDB->SOLVER_TYPE !=1) )
  {
    // s_base.Solve(s_base.solution.get_entries(), 
    //                    s_base.rhs.get_entries());
    this->Time_NSE2D::solve();
  }*/
}

/**************************************************************************** */
void PrRobustTime_NSE2D::output(int m)
{
  /*int im=0;
  this->Time_NSE2D::output(m, im);*/
}




















