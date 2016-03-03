#include <PrRobustNSE2D.h>
#include <Assemble2D.h>
#include <NSE2D.h>
#include <Domain.h>
#include <Database.h>
#include <LinAlg.h>


extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *m, int *n, int *k, 
           double *alpha, double *A, int *lda, double *B, 
           int *ldb, double *beta, double *C, int *ldc);

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
void matrix_multiplication(double ***inputMat, int *nrowIn, int *ncolIn, 
                                double ***outputMat, int *nrowOut, int *ncolOut,
                                double **inputrhs,  int *ndimIn, 
                                double **outputrhs, int *ndimOut)
{
  /** inputMat[0] = A00,    inputMat[1] = A11;
   *  inputMat[2] = A00_nl  inputMat[3] = A11_nl
   *  inputMat[4] = P0      inputMat[5] = P1
   */
  // preparation for dgemm_
  double *A00 = new double[nrowIn[0]*ncolIn[0]];
  double *A11 = new double[nrowIn[1]*ncolIn[1]];  
  // A00 and A11 have the same size
  for(int i=0; i<nrowIn[0]; i++)
  {
    memcpy(A00+i*ncolIn[0], inputMat[0][i], ncolIn[0]*SizeOfDouble);
    memcpy(A11+i*ncolIn[1], inputMat[1][i], ncolIn[1]*SizeOfDouble);
  }
  // empty matrices filled with zeros
  double *A01 = new double[nrowIn[0]*ncolIn[0]];
  double *A10 = new double[nrowIn[1]*ncolIn[1]];
  memset(A01, 0, nrowIn[0]*ncolIn[0]*SizeOfDouble);
  memset(A10, 0, nrowIn[1]*ncolIn[1]*SizeOfDouble);
  
  // prepareing nonlinear matrices
  double *Anl00 = new double[nrowIn[2]*ncolIn[2]];
  double *Anl11 = new double[nrowIn[3]*ncolIn[3]];
  for(int i=0; i<nrowIn[2]; i++)
  {
    memcpy(Anl00+i*ncolIn[2], inputMat[2][i], ncolIn[2]*SizeOfDouble);
    memcpy(Anl11+i*ncolIn[3], inputMat[3][i], ncolIn[3]*SizeOfDouble);
  }
  // preparing P0 and p1
  double *P0 = new double[nrowIn[4]*ncolIn[4]];
  double *P1 = new double[nrowIn[5]*ncolIn[5]];
  for(int i=0; i<nrowIn[4]; i++)
  {
    memcpy(P0+i*ncolIn[4], inputMat[4][i], ncolIn[4]*SizeOfDouble);
    memcpy(P1+i*ncolIn[5], inputMat[5][i], ncolIn[5]*SizeOfDouble);
  }
  
  int m, n, k;//, lda, ldb, ldc;
  m = nrowIn[4];
  n = ncolIn[2];
  k = nrowIn[2]; // = ncolIn[4]
  
  double alpha = 1.;
  double beta  = 1.;
    
  dgemm_((char*)"N", (char*)"N", &n, &m, &k, &alpha, Anl00, &n, P0, &k, 
         &beta, A00, &n);
  
  dgemm_((char*)"N", (char*)"N", &n, &m, &k, &alpha, Anl00, &n, P1, &k, 
         &beta, A01, &n);
  
  dgemm_((char*)"N", (char*)"N", &n, &m, &k, &alpha, Anl11, &n, P0, &k, 
         &beta, A10, &n);
  
  dgemm_((char*)"N", (char*)"N", &n, &m, &k, &alpha, Anl11, &n, P1, &k, 
         &beta, A11, &n);
  
  for(int i=0; i<n; i++)
  {
    memcpy(outputMat[0][i], A00+i*n, n*SizeOfDouble);
    memcpy(outputMat[1][i], A01+i*n, n*SizeOfDouble);
    memcpy(outputMat[2][i], A10+i*n, n*SizeOfDouble);
    memcpy(outputMat[3][i], A11+i*n, n*SizeOfDouble);
  }
  delete [] A00;   delete [] A01; 
  delete [] A10;   delete [] A11;
  delete [] Anl00; delete Anl11;
  delete [] P0;    delete P1;  
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
    BoundCondFunct2D *bc[3] ={ velocity_space->GetBoundCondition(), 
                             velocity_space->GetBoundCondition(), 
                             pressure_space->GetBoundCondition() };
    // boundary values
    BoundValueFunct2D * const * const BoundValue=this->NSE2D::get_example().get_bd();
    BoundValueFunct2D *bv[3] = {BoundValue[0], BoundValue[1], BoundValue[2] };
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
  if(TDatabase::ParamDB->DISCTYPE != RECONSTRUCTION)
    return;
  
  unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
  TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
  // prepare everything to assemble the global projection matrix 
  
  for(SystemPerGrid& s_derived : this->Systems)
  {
    const int nFESpaces = 2;
    const TFESpace2D *project_space  = &this->get_projection_space();
    const TFESpace2D *velocity_space = &this->NSE2D::get_velocity_space();
    // pointers to the fespace
    const TFESpace2D * pointer_to_space[2] = { velocity_space, project_space };
    // assemble all square and rectangular matrices  
    const int nSqMatAssemble = 0;
    const int reqMatAssemble = 2;
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
    
    for(int i=0; i<nReMatStored; i++)
        rectMatricesStored[i]->reset();
    BoundCondFunct2D * boundCond[2] = {
        project_space->GetBoundCondition(), project_space->GetBoundCondition()};
    
    std::array<BoundValueFunct2D*, 3> boundVal;
    boundVal[0] = example.get_bd()[0];
    boundVal[1] = example.get_bd()[1];
    boundVal[2] = example.get_bd()[2];
    
    TFEFunction2D * fefct[3] = {this->NSE2D::systems.front().u.GetComponent(0), 
                                 this->NSE2D::systems.front().u.GetComponent(1), 
                                 &this->NSE2D::systems.front().p};
    // no local assembling routine is needed
    // the projection matrix uses alternate routine
    // namely "projection_matrices"
    LocalAssembling2D empty_constructor(NO_LOCAL_ASSEMBLE, fefct, 
                                        this->NSE2D::get_example().get_coeffs());
    
    Assemble2D_MixedFEM(nFESpaces, pointer_to_space, 
                        nSqMatAssemble, reqMatAssemble,
                        rowSpace, colSpace, 
                        nRhsAssemble, rowSpaceRhs, 
                        nSqMatStored, sqMatricesStored, 
                        nReMatStored, rectMatricesStored, 
                        nRhsStored, rhsStored,pointToRhsStored, 
                        empty_constructor, nullptr, nullptr,
                        projection_matrices, boundCond, boundVal.data());    
  }//
  TDatabase::ParamDB->VELOCITY_SPACE = vel_space;
}
/** ************************************************************************ */
void PrRobustNSE2D::assembleNL()
{
  // assembling of nonlinear term
  if(TDatabase::ParamDB->DISCTYPE != RECONSTRUCTION)
  {
    this->NSE2D::assemble_nonlinear_term();
    return;
  }
  
  System_per_grid& s_base = this->NSE2D::systems.front();
  for(SystemPerGrid& s_derived : this->Systems)
  {
    TFEFunction2D *fefct[3] = {s_base.u.GetComponent(0), s_base.u.GetComponent(1), 
                               &s_base.p};
    // local assembling object
    LocalAssembling2D la(RECONSTR_NLGALERKIN, fefct, 
                           this->NSE2D::get_example().get_coeffs());
    // setting the space for sign in GetSignOfThisDOF();
    unsigned int vel_space = TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->VELOCITY_SPACE = TDatabase::ParamDB->PROJECTION_SPACE;
    // prepare everything which is needed for assembling matrices and 
    // the right hand side
    const int nFESpaces = 2;
    const TFESpace2D *project_space  = &this->get_projection_space();
    const TFESpace2D *velocity_space = &this->NSE2D::get_velocity_space();
    const TFESpace2D *pressure_space = &this->NSE2D::get_pressure_space();
    // pointers to the fespace
    const TFESpace2D * pointer_to_space[2] = { velocity_space, project_space };
    
    const int nSqMatAssemble=2;
    const int nReMatAssemble=4;
    std::vector<int> rowSpace ={0, 0, 1, 1, 0, 0}; // row space for assembling matrices
    std::vector<int> colSpace ={0, 0, 0, 0, 1, 1}; // cols space for assembling matrices
    const int nRhsAssemble = 0;
    //NOTE:for the assembling of matrices or right hand side, space numbers
    //are passed as array: this corresponds to the array "pointer_to_space"
    std::vector<int> rowSpaceRhs={}; // row space for assembling rhs 
    
    const int nSqMatStored = 4;
    TSquareMatrix2D * sqMatricesStored[nSqMatStored]{nullptr}; // 4 matrices to be stored
    
    std::vector<std::shared_ptr<FEMatrix>> matrixA
           = s_base.matrix.get_blocks_uniquely();

    sqMatricesStored[0]=reinterpret_cast<TSquareMatrix2D*>(matrixA.at(0).get());
    sqMatricesStored[1]=reinterpret_cast<TSquareMatrix2D*>(matrixA.at(1).get());
    sqMatricesStored[2]=reinterpret_cast<TSquareMatrix2D*>(matrixA.at(3).get());
    sqMatricesStored[3]=reinterpret_cast<TSquareMatrix2D*>(matrixA.at(4).get());
    // reset the matrices
    for(int i=0; i<nSqMatStored; i++)
      sqMatricesStored[i]->reset();
    
    // rectangular matrices to be stored
    const int nReMatStored = 0;
    TMatrix2D ** rectMatricesStored = nullptr;
    
    // right hand vector to be stored
    int nRhsStored = 0;
    const TFESpace2D ** pointToRhsStored{nullptr};
    double **rhsStored{nullptr};    
    // boundary conditions    
    BoundCondFunct2D *boundCond[3] ={ velocity_space->GetBoundCondition(), 
                             velocity_space->GetBoundCondition(),
                             pressure_space->GetBoundCondition() };
    // boundary values
    std::array<BoundValueFunct2D*, 3> boundVal;
    boundVal[0] = example.get_bd()[0];
    boundVal[1] = example.get_bd()[1];
    boundVal[2] = example.get_bd()[2];
    
    Assemble2D_MixedFEM(nFESpaces, pointer_to_space, 
                        nSqMatAssemble, nReMatAssemble,
                        rowSpace, colSpace, 
                        nRhsAssemble, rowSpaceRhs, 
                        nSqMatStored, sqMatricesStored, 
                        nReMatStored, rectMatricesStored, 
                        nRhsStored, rhsStored,pointToRhsStored, 
                        la, matrix_multiplication, nullptr,
                        ProjectionMatricesNSE2D, boundCond, boundVal.data());    
    TDatabase::ParamDB->VELOCITY_SPACE=vel_space;
  }
}

/** ************************************************************************ */
bool PrRobustNSE2D::stopIteration(unsigned int it)
{
  System_per_grid& s = this->systems.front();
  unsigned int nuDof = s.solution.length(0);
  unsigned int npDof = s.solution.length(2);
  
  this->defect = s.rhs;   
  s.matrix.apply_scaled_add(s.solution, defect,-1.);  
  // 
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20FEFunction(&defect[2*nuDof], npDof, &this->get_pressure_space(),
                      TDatabase::ParamDB->VELOCITY_SPACE, 
                      TDatabase::ParamDB->PRESSURE_SPACE);
  double residual =  Ddot(2*nuDof+npDof, &this->defect[0], &this->defect[0]);
  double impulse_residual = Ddot(2*nuDof, &this->defect[0],
         &this->defect[0]);
  double mass_residual    = Ddot(npDof,&this->defect[2*nuDof],
         &this->defect[2*nuDof]);
  
  Output::print("nonlinear step  :  " , setw(3), it);
  Output::print("impulse_residual:  " , setw(3), impulse_residual);
  Output::print("mass_residual   :  " , setw(3), mass_residual);
  Output::print("residual        :  " , setw(3), sqrt(residual));
  
  if (it>0)
  {
    Output::print("rate:           :  " , setw(3), sqrt(residual)/oldResidual);
  }
   
  oldResidual = sqrt(residual);
  if(it == 0)
    initial_residual = sqrt(residual);

  int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  double limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  if (TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE)
  {
    limit *= sqrt(this->get_size());
    Output::print("stopping tolerance for nonlinear iteration ", limit);
  }
  if ((((sqrt(residual)<=limit)||(it==Max_It)))
    && (it >= TDatabase::ParamDB->SC_MINIT))
  {
    Output::print("ITE : ", setw(3), it, "  RES : ", sqrt(residual), 
                  " Reduction : ",  sqrt(residual)/initial_residual);
    return true;
  }
  else
    return false;
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
