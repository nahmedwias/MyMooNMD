#include <PrRobustTime_NSE2D.h>
#include <Database.h>
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <Output2D.h>

#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <LinAlg.h>

using namespace std;

extern "C" void dgemm_(char *TRANSA, char *TRANSB, int *m, int *n, int *k, 
           double *alpha, double *A, int *lda, double *B, 
           int *ldb, double *beta, double *C, int *ldc);

PrRobustTime_NSE2D::SystemPerGrid::SystemPerGrid(const Example_NSE2D& example,
      TCollection& coll, const TFESpace2D& velocity_space, 
      const TFESpace2D& pressure_space, int order)
: projection_space(velocity_space.GetCollection(), (char*)"u", 
                    (char*)"projection space", example.get_bc(2), order, nullptr),
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
/**************************************************************************** */


void matrices_reconstruction(double ***inputMat, int *nrowInput, int *ncolInput, 
                                double ***outputMat, int *nrowOut, int *ncolOut, 
                                double **inputrhs,  int *ndimInput, 
                                double **outputrhs, int *ndimOutput)
{
 /** this lambda function compute the ABC^T
  * @param: matA projection matrix first velocity part (special case)
  * @param: matB mass matrix for BDM or RT which is multiplied in between
  * @param: matC projection matrix which multiplies as transpose from right 
  *         only for the off diagonal blocks (flag==1 & flag ==2)
  * @param: product the resulting matrix
  * @param: flag 0 and 3 are the cases when matA^T is multiplied from right
  *              1 and 2 when the matC^T is multiplied from right
  */
 //FIXME slow try differently when everything works
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
  std::pair<int, int> rowsArray(nrowInput[5], nrowInput[2]);
  std::pair<int, int> colsArray(ncolInput[5], ncolInput[2]);
  compute_matrix_matrix(rowsArray, colsArray, inputMat[5], inputMat[2], 
                          nullptr, outputMat[0], 0);
  // prepare the outputMat[0] = P0 * M * P1^T
  compute_matrix_matrix(rowsArray, colsArray, inputMat[5], inputMat[2], 
                          inputMat[6], outputMat[1], 1);
  // prepare the outputMat[0] = P1 * M * P0^T
  rowsArray.first=nrowInput[6]; rowsArray.second=nrowInput[2];
  colsArray.first=ncolInput[6]; rowsArray.second=ncolInput[2];
  compute_matrix_matrix(rowsArray, colsArray, inputMat[6], inputMat[2], 
                          inputMat[5], outputMat[2], 2);
  // prepare the outputMat[0] = P1 * M * P1^T
  compute_matrix_matrix(rowsArray, colsArray, inputMat[6], inputMat[2],
                          nullptr, outputMat[3], 3);
  
  // prepare the output right hand 
  std::pair<int,int> size (nrowInput[0], ncolInput[0]);
  double **matrix[2];
  matrix[0] = inputMat[4];
  matrix[1] = inputMat[5];
  matrix_vector_multiply(size, matrix, inputrhs[0], outputrhs);
  //======================================================================
  // matrix-matrix multiplication using dgemm for the stiffness
  // matrix
  //======================================================================
  /**
   * A_NL * P^T + AL
   * 0 is loc matrix P0
   * 1 is loc matrix P1
   * 2 is loc matrix M
   * 3 is loc matrix A(grad grad)
   * 4 is loc matrix Anl0
   * 5 is loc matrix Anl1
   */
  double *P0 =  new double[nrowInput[5]*ncolInput[5]];
  double *P1 =  new double[nrowInput[6]*ncolInput[6]];
  // copy the matrices
  for(int i=0; i<nrowInput[5]; ++i)
  {
    memcpy(P0+i*ncolInput[5], inputMat[5][i], ncolInput[5]*SizeOfDouble);
    memcpy(P1+i*ncolInput[6], inputMat[6][i], ncolInput[6]*SizeOfDouble);
  }
  // nrowInput[0]=nrowInput[1]; and ncolInput[0] = ncolInput[1]
  double *A00= new double[nrowInput[0]*ncolInput[0]];
  double *A01= new double[nrowInput[0]*ncolInput[0]];
  double *A10= new double[nrowInput[0]*ncolInput[0]];
  double *A11= new double[nrowInput[1]*ncolInput[1]];
  
  memset(A10,0,nrowInput[0]*ncolInput[1]*SizeOfDouble);
  memset(A01,0,nrowInput[0]*ncolInput[1]*SizeOfDouble);
  
  for(int i=0; i<nrowInput[0]; ++i)
  {
    memcpy(A00+i*ncolInput[0], inputMat[0][i], ncolInput[0]*SizeOfDouble);
    memcpy(A11+i*ncolInput[1], inputMat[1][i], ncolInput[1]*SizeOfDouble);
  }
  double *Anl00 = new double[nrowInput[3]*ncolInput[3]];
  double *Anl11 = new double[nrowInput[4]*ncolInput[4]];
  for(int i=0; i<nrowInput[3]; ++i)
  {
    memcpy(Anl00+i*ncolInput[3], inputMat[3][i], ncolInput[3]*SizeOfDouble);
    memcpy(Anl11+i*ncolInput[4], inputMat[4][i], ncolInput[4]*SizeOfDouble);
  }
  int m, n, k;
  // local matrix A00 = P0*Anl0 + A(grad,grad)
  m = nrowInput[5];
  n = ncolInput[3];
  k = nrowInput[3];
  
  double alpha = 1.;
  double beta  = 1.;
  
  //======================================================================
  // normally dgemm uses fortran format where the entries are treated
  // as column storage format, but we can calkulate the matrices 
  // multiplication and addition in the row-storage format as follows:
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
    memcpy(outputMat[4][i], A00+i*n, n*SizeOfDouble);
    memcpy(outputMat[5][i], A01+i*n, n*SizeOfDouble);
    memcpy(outputMat[6][i], A10+i*n, n*SizeOfDouble);
    memcpy(outputMat[7][i], A11+i*n, n*SizeOfDouble);
  }

  delete [] P0;    delete [] P1;
  delete [] A00;   delete [] A01;
  delete [] A10;   delete [] A11;
  delete [] Anl00; delete [] Anl11;
}

/**************************************************************************** */
void nonlinear_term_reconstruct(double ***inputMat, int *nrowIn, int *ncolIn, 
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
  delete [] Anl00; delete [] Anl11;
  delete [] P0;    delete [] P1;  
}

/**************************************************************************** */
void PrRobustTime_NSE2D::assemble_initial_timeNS()
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
    //=========================
    // assemble all matrices and right hand side     
    // local assemble routine
    TFEFunction2D *fefct[3] = {s_base.u.GetComponent(0), s_base.u.GetComponent(1),&s_base.p};
    LocalAssembling2D la(RECONSTR_TNSE, fefct, 
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
     
    const int nSqMatAssemble=3;
    const int nReMatAssemble=4;
    std::vector<int> rowSpace = {0, 0, 1, 1, 1, 0, 0}; // row space for assembling matrices
    std::vector<int> colSpace = {0, 0, 1, 0, 0, 1, 1}; // cols space for assembling matrices
    
    const int nRhsAssemble = 1;
   //NOTE:for the assembling of matrices or right hand side, space numbers
   //are passed as array: this corresponds to the array "pointer_to_space"
   std::vector<int> rowSpaceRhs={1}; // row space for assembling rhs 

    // prepare everything for storing the matrices and rhs
    // this will be used latter inside the class
    const int nSqMatStored = 8;
    TSquareMatrix2D * sqMatricesStored[8]{nullptr}; // 8 matrices to be stored

    std::vector<std::shared_ptr<FEMatrix>> modified_mass
           = s_derived.Modifed_Mass.get_blocks_uniquely();
    std::vector<std::shared_ptr<FEMatrix>> matrixA
           = s_base.matrix.get_blocks_uniquely();
    sqMatricesStored[0]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(0).get());
    sqMatricesStored[1]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(1).get());
    sqMatricesStored[2]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(2).get());
    sqMatricesStored[3]=reinterpret_cast<TSquareMatrix2D*>(modified_mass.at(3).get());
    
    sqMatricesStored[4]=reinterpret_cast<TSquareMatrix2D*>(matrixA.at(0).get());
    sqMatricesStored[5]=reinterpret_cast<TSquareMatrix2D*>(matrixA.at(1).get());
    sqMatricesStored[6]=reinterpret_cast<TSquareMatrix2D*>(matrixA.at(3).get());
    sqMatricesStored[7]=reinterpret_cast<TSquareMatrix2D*>(matrixA.at(4).get());
    // reset the matrices
    for(int i=0; i<nSqMatStored; i++)
      sqMatricesStored[i]->reset();
    
    // rectangular matrices to be stored
    const int nReMatStored = 0;
    TMatrix2D ** rectMatricesStored = nullptr;
    
    // right hand vector to be stored
    int nRhsStored = 2;
    const TFESpace2D * pointToRhsStored[3]={velocity_space, velocity_space,
                                                 pressure_space};
    double *rhsStored[3] = {s_base.rhs.block(0), s_base.rhs.block(1), 
                            s_base.rhs.block(2)};
    s_base.rhs.reset();
    // boundary conditions    
    BoundCondFunct2D *boundCond[3] ={ velocity_space->GetBoundCondition(), 
                             velocity_space->GetBoundCondition(),
                             pressure_space->GetBoundCondition() };
    // boundary values
    std::array<BoundValueFunct2D*, 3> boundVal;
    boundVal[0] = example.get_bd()[0];
    boundVal[1] = example.get_bd()[1];
    boundVal[2] = example.get_bd()[2];
    // assemble the correspoding matrices and right hand side
    Assemble2D_MixedFEM(nFESpaces, pointer_to_space, 
                        nSqMatAssemble, nReMatAssemble,
                        rowSpace, colSpace, 
                        nRhsAssemble, rowSpaceRhs, 
                        nSqMatStored, sqMatricesStored, 
                        nReMatStored, rectMatricesStored, 
                        nRhsStored, rhsStored,pointToRhsStored, 
                        la, matrices_reconstruction, nullptr,
                        ProjectionMatricesTNSE2D, boundCond, boundVal.data());
    
    TDatabase::ParamDB->VELOCITY_SPACE = vel_space;
  }// endof assembling for all grids
  //s_base.matrix.get_blocks().at(4)->Print("A");
  //  cout<<"Its really my function: ???" << endl;exit(0);
  
  // copy the modified rhs to the old_rhs_modified for the right hand 
  // side at the initial time
  this->old_rhs_modified = this->Time_NSE2D::systems.front().rhs;
  this->old_solution = s_base.solution;
}

/**************************************************************************** */
void PrRobustTime_NSE2D::assemble_rhsNS()
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
    const int nSqMatAssemble = 0;
    const int nReMatAssemble = 2;
    
    const int nRhsAssemble = 1;    
    //NOTE:for the assembling of matrices or right hand side, space numbers
    //are passed as array: this corresponds to the array "pointer_to_space"
    std::vector<int> rowSpace ={0, 0}; // row space for assembling matrices
    std::vector<int> colSpace ={1, 1}; // cols space for assembling matrices
    std::vector<int> rowSpaceRhs={1}; // row space for assembling rhs 
    
    // prepare everything for storing the matrices and rhs
    // this will be used latter inside the class
    const int nSqMatStored = 0;
    TSquareMatrix2D **sqMatricesStored = nullptr;
    const int nReMatStored = 0;
   TMatrix2D **rectMatricesStored=nullptr;
    int nRhsStored = 2;
    const TFESpace2D * pointToRhsStored[3]={velocity_space,velocity_space,
                                                 pressure_space};
    s_base.rhs.reset();
    double *rhsStored[3] = {s_base.rhs.block(0), s_base.rhs.block(1), 
                            s_base.rhs.block(2)};
    // boundary conditions    
    BoundCondFunct2D * boundCond[2] = {
        velocity_space->GetBoundCondition(), velocity_space->GetBoundCondition()};
    
    std::array<BoundValueFunct2D*, 3> boundVal;
    boundVal[0] = example.get_bd()[0];
    boundVal[1] = example.get_bd()[1];
    boundVal[2] = example.get_bd()[2];
    
    // assemble the right hand side
    Assemble2D_MixedFEM(nFESpaces, pointer_to_space, 
                        nSqMatAssemble, nReMatAssemble,
                        rowSpace, colSpace, 
                        nRhsAssemble, rowSpaceRhs, 
                        nSqMatStored, sqMatricesStored, 
                        nReMatStored, rectMatricesStored, 
                        nRhsStored, rhsStored,pointToRhsStored, 
                        larhs, nullptr, MatVectMult,
                        projection_matrices, boundCond, boundVal.data());
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
    
    // FIXME FInd other solution than this submatrix method.
    // M u^{k-1}
    s_derived.Modifed_Mass.apply_scaled_submatrix(old_solution, s_base.rhs, 2, 2, 1.0);
    // -tau*theta2 * A u^{k-1}
    double factor = -tau*theta2;
    s_base.matrix.apply_scaled_submatrix(old_solution, s_base.rhs, 2, 2, factor);
    
    // scale actives in each block
    
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
void PrRobustTime_NSE2D::assemble_nonlinearNS()
{
  // check and return if the DISCTYPE==RECONSTRUCTION
  if(TDatabase::ParamDB->DISCTYPE != RECONSTRUCTION)
  {
    this->Time_NSE2D::assemble_nonlinear_term();
    return;
  }
  System_per_grid& s_base = this->Time_NSE2D::systems.front();
  // assemble mass matrix using vector space  
  for(SystemPerGrid &s_derived : this->Systems)
  {    
    //=========================
    // assemble all matrices and right hand side     
    // local assemble routine
    TFEFunction2D *fefct[3] = {s_base.u.GetComponent(0), s_base.u.GetComponent(1), 
                               &s_base.p};
    LocalAssembling2D la(RECONSTR_TNSENL, fefct, 
                           this->Time_NSE2D::get_example().get_coeffs());    
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
     
    const int nSqMatAssemble=2;
    const int nReMatAssemble=4;
    std::vector<int> rowSpace ={0, 0, 1, 1, 0, 0}; // row space for assembling matrices
    std::vector<int> colSpace ={0, 0, 0, 0, 1, 1}; // cols space for assembling matrices
    
    const int nRhsAssemble = 0;
   //NOTE:for the assembling of matrices or right hand side, space numbers
   //are passed as array: this corresponds to the array "pointer_to_space"
   std::vector<int> rowSpaceRhs={}; // row space for assembling rhs 

    // prepare everything for storing the matrices and rhs
    // this will be used latter inside the class
    const int nSqMatStored = 4;
    TSquareMatrix2D * sqMatricesStored[nSqMatStored]{nullptr}; // 8 matrices to be stored
   
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
    // assemble the correspoding matrices and right hand side
    Assemble2D_MixedFEM(nFESpaces, pointer_to_space, 
                        nSqMatAssemble, nReMatAssemble,
                        rowSpace, colSpace, 
                        nRhsAssemble, rowSpaceRhs, 
                        nSqMatStored, sqMatricesStored, 
                        nReMatStored, rectMatricesStored, 
                        nRhsStored, rhsStored,pointToRhsStored, 
                        la, nonlinear_term_reconstruct, nullptr,
                        ProjectionMatricesNSE2D, boundCond, boundVal.data());
    TDatabase::ParamDB->VELOCITY_SPACE=vel_space;
  }// endof assembling for all grids
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

/**************************************************************************** */
bool PrRobustTime_NSE2D::stopIte(unsigned int it_counter)
{
  if(TDatabase::ParamDB->DISCTYPE==RECONSTRUCTION)
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
    
    Output::print("nonlinear step  :  " , setw(3), it_counter , setw(14),
                   impulse_residual , setw(14) , mass_residual, 
                   setw(14) , sqrt(residual));
    
    if (it_counter>0)
    {
      Output::print("rate:           :  " , setw(3), sqrt(residual)/oldResidual);
    }
    
    oldResidual = sqrt(residual);
    if(it_counter == 0)
      initial_residual = sqrt(residual);
    
    int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
    double limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
    if (TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE)
    {
      limit *= sqrt(this->get_size());
      Output::print("stopping tolerance for nonlinear iteration ", limit);
    }
    
    if ((((sqrt(residual)<=limit)||(it_counter==Max_It)))
     && (it_counter>=TDatabase::ParamDB->SC_MINIT))
     {
       Output::print("ITE : ", setw(3), it_counter, "  RES : ", sqrt(residual), 
                     " Reduction : ",  sqrt(residual)/initial_residual);
       // descale the matrices, since only the diagonal A block will 
       // be reassembled in the next time step
       this->descaleMatrices();     
       return true;
     }
     else
       return false;
  }
  else
  {
    int ret;
    ret=this->Time_NSE2D::stopIte(it_counter);
    return ret;
  }
}

/**************************************************************************** */