#include "Extra_Routines_PR.h"

#include <Constants.h>
#include <Database.h>
#include <Enumerations.h>
#include <Convolution.h>
#include <LinAlg.h>
#include <ConvDiff.h>

#include <AllRefTrans.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
// #include <malloc.h>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>


#ifdef __2D__

int GetSignOfThisDOF(int N_DOF, int DOF)
{
  switch (N_DOF)
  {
  case 3:// Raviart-Thomas zeroth order, triangles
    return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[DOF];
    break;
    // note that RT2 on quadrilaterals and RT3 on triangles have 24 basis functions
  case 4:// Raviart-Thomas zeroth order, quadrilaterals
    return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[DOF];
    break;
  case 6: // BDM first order, triangles
    // degree of freedom on an edge, no inner degree of freedom
    return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%2)/2];
    break;
  case 8:
    // Raviart-Thomas first order, triangles
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1001)
    {
      if(DOF<6) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%2)/2];
      else // inner degree of freedom
        return 1;
    }
    else if (TDatabase::ParamDB->VELOCITY_SPACE == 1011) //BDM first order, quadrilaterals
    {
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%2)/2];
    }
    break;
  case 12:
    // Raviart-Thomas first order, quadrilaterals
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1001)
    {
    if(DOF<8) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%2)/2];
    else // inner degree of freedom
      return 1;
    }
    else if (TDatabase::ParamDB->VELOCITY_SPACE == 1012) //BDM second order, triangles
    {
      if(DOF<9) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%3)/3];
      else // inner degree of freedom
        return 1;
    }
    break;
  case 14://BDM second order, quadrilaterals
    if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%3)/3];
    else // inner degree of freedom
      return 1;
    break;
  case 15:// Raviart-Thomas second order, triangles
    if(DOF<9) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%3)/3];
    else // inner degree of freedom
      return 1;
    break;
  case 20:// BDM third order. triangles
    if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%4)/4];
    else // inner degree of freedom
      return 1;
    break;
  case 22://BDM third order, quadrilaterals
    if(DOF<16) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%4)/4];
    else // inner degree of freedom
      return 1;
    break;
  case 24:
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1002)
    {
      // Raviart-Thomas second order, quadrilaterals
      if(DOF<12) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%3)/3];
      else // inner degree of freedom
        return 1;
    }
    else if(TDatabase::ParamDB->VELOCITY_SPACE == 1003)
    {
      // Raviart-Thomas third order, triangles
      if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%4)/4];
      else // inner degree of freedom
        return 1;
    }
    else
    {
      ErrThrow("VELOCITY_SPACE has to be set to either 1002 or 1003\n");
    }
    break;
  case 40:// Raviart-Thomas third order, quadrilaterals
    if(DOF<16) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%4)/4];
    else
      return 1;
    break;
  default:
       ErrThrow("WARNING: Unknown Raviart-Thomas or BDM element !", N_DOF);
       break;
  }
  //dummy return, should not be reached
  return -4711;
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
  const double *xi, *eta;
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
  
  for(i=0; i<=1; ++i)
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
  
  double **MatrixP0 = locMatrix[0];
  double **MatrixP1 = locMatrix[1];
  
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

void ProjectionMatricesNSE2D(int current_cell, const TFESpace2D* ansatzSpace, 
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
  const double *xi, *eta;
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

  double **MatrixP0 = locMatrix[4];
  double **MatrixP1 = locMatrix[5];
  
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


void ProjectionMatricesTNSE2D(int current_cell, const TFESpace2D* ansatzSpace, 
                         const TFESpace2D* testSpace, double ***locMatrix)
{
  // projection is the ansatzSpace
  // velocity is the test space
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
  const double *xi, *eta;
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

  double **MatrixP0 = locMatrix[5];
  double **MatrixP1 = locMatrix[6];
  
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

void MatVectMult(double ***inputMat, std::pair<int,int>size, double *inputRhs, 
                 double **outputRhs)
{
  unsigned int nrow = size.first;
  unsigned int ncol = size.second;
  
  for(unsigned int i=0; i<nrow; i++)
  {
    double temp = 0;
    double temp1 = 0;
    for(unsigned int j=0; j<ncol; j++)
    { //FIXME the difference here is the number of 
      // matrices: currently even in the steady state case
      // three matrices are allocated but only two of them
      // are assembled namly the inputMat[1] and inputMat[2]
      // but in the time dependent case also the matrix inputMat[0]
      // is assembled. The thing which has to be fixed is the assembling
      // of matrices. Try to assemble only two matrices and then change 
      // the argument in the product case
      temp  += inputMat[0][i][j] * inputRhs[j];
      temp1 += inputMat[1][i][j] * inputRhs[j];//
    }
    outputRhs[0][i] = temp;
    outputRhs[1][i] = temp1;
  }
}
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
        // cout<<i + j*ncol<<"   " <<matrix[0][i][j] << "  " << matrix[1][i][j] << endl;
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
  // prepare the outputMat[1] = P0 * M * P1^T
  compute_matrix_matrix(rowsArray, colsArray, inputMat[5], inputMat[2],
                          inputMat[6], outputMat[1], 1);
  // prepare the outputMat[2] = P1 * M * P0^T
  rowsArray.first=nrowInput[6]; rowsArray.second=nrowInput[2];
  colsArray.first=ncolInput[6]; rowsArray.second=ncolInput[2];
  compute_matrix_matrix(rowsArray, colsArray, inputMat[6], inputMat[2],
                          inputMat[5], outputMat[2], 2);
  // prepare the outputMat[3] = P1 * M * P1^T
  compute_matrix_matrix(rowsArray, colsArray, inputMat[6], inputMat[2],
                          nullptr, outputMat[3], 3);

  // prepare the output right hand
  std::pair<int,int> size (nrowInput[5], ncolInput[6]);
  double **matrix[2];
  matrix[0] = inputMat[5];
  matrix[1] = inputMat[6];
  matrix_vector_multiply(size, matrix, inputrhs[0], outputrhs);
  //======================================================================
  // matrix-matrix multiplication using dgemm for the stiffness
  // matrix
  //======================================================================
  /**
   * A_NL * P^T + AL
   * 0 is loc matrix A11 linear part
   * 1 is loc matrix A22 linear part
   * 2 is loc matrix M   mass matrix vector valued
   * 3 is loc matrix A11 nonlinear (rectangular)
   * 4 is loc matrix A22 nonlinear (rectangular)
   * 5 is loc matrix P0 projection (rectangular)
   * 6 is loc matrix P1 projection (rectangular)
   */
  auto matrix_times_matrix_plus_matrix = [](unsigned int rows, unsigned int cols,
                                            double **matrix_left, double **matrix_right, 
                                            double **matrix_plus, double **&matrix_out)
  {
    for(unsigned int k=0; k<cols; ++k)
    {
      for(unsigned int i=0; i< rows; ++i)
      {
        double temp = matrix_left[i][k];
        for(unsigned int j=0; j< rows; ++j)
        {
          matrix_out[i][j] += temp*matrix_right[k][j];
        }
      }
    }
    for(unsigned int i=0; i<rows; ++i)
    {
      for(unsigned int j=0; j<rows; ++j)
        matrix_out[i][j] += matrix_plus[i][j];
    }
    return;
  };
  /*Output::print("row A11: ", nrowInput[0], " col A11: ", ncolInput[0]);
  Output::print("row A22: ", nrowInput[1], " col A22: ", ncolInput[1]);
  Output::print("row M: ", nrowInput[2], " col M: ", ncolInput[2]);
  Output::print("row A11nl: ", nrowInput[3], " col A11nl: ", ncolInput[3]);
  Output::print("row A22nl: ", nrowInput[4], " col A22nl: ", ncolInput[4]);
  Output::print("row P0: ", nrowInput[5], " col P0: ", ncolInput[5]);
  Output::print("row P1: ", nrowInput[6], " col P1: ", ncolInput[6]);*/
  // P0 * A00nl   P1 * A00nl
  // P0 * A11nl   P1 * A11nl
  // some checks
  if(nrowInput[5] != ncolInput[3])
  {
    ErrThrow("dimensions mismatch");
  }
  unsigned int r = nrowInput[5];
  unsigned int c = ncolInput[5];
  double **temp = new double*[nrowInput[0]];
  for(int i=0; i<nrowInput[0]; ++i)
    temp[i]=new double[ncolInput[0]];
  for(int i=0; i<nrowInput[0]; ++i)
    memset(temp[i], 0, ncolInput[0]*sizeof(double));
  // P0*A00nl + A00(linear);
  matrix_times_matrix_plus_matrix(r, c, inputMat[5], inputMat[3], inputMat[0], 
                                  outputMat[4]);
  // P1 * A00nl 
  matrix_times_matrix_plus_matrix(r, c, inputMat[5], inputMat[4], temp, 
                                  outputMat[5]);
  // P0*A11nl;
  matrix_times_matrix_plus_matrix(r, c, inputMat[6], inputMat[3], temp, 
                                  outputMat[6]);
  // P1 * A11nl + A11(linear) 
  matrix_times_matrix_plus_matrix(r, c, inputMat[6], inputMat[4], inputMat[1], 
                                  outputMat[7]);
  delete [] temp[0];
  delete [] temp;
}

/**************************************************************************** */
void nonlinear_term_reconstruct(double ***inputMat, int *nrowInput, int *ncolInput,
                                double ***outputMat, int *nrowOut, int *ncolOut,
                                double **inputrhs,  int *ndimIn,
                                double **outputrhs, int *ndimOut)
{
  /** inputMat[0] = A00,    inputMat[1] = A11;
   *  inputMat[2] = A00_nl  inputMat[3] = A11_nl
   *  inputMat[4] = P0      inputMat[5] = P1
   */
  auto matrix_times_matrix_plus_matrix = [](unsigned int rows, unsigned int cols,
                                            double **matrix_left, double **matrix_right, 
                                            double **matrix_plus, double **&matrix_out)
  {
    // rows of matrix_left are the columns of matrix_right
    for(unsigned int k=0; k<cols; ++k)
    {
      for(unsigned int i=0; i< rows; ++i)
      {
        double temp = matrix_left[i][k];
        for(unsigned int j=0; j< rows; ++j)
        {
          matrix_out[i][j] += temp*matrix_right[k][j];
        }
      }
    }
    for(unsigned int i=0; i<rows; ++i)
    {
      for(unsigned int j=0; j<rows; ++j)
        matrix_out[i][j] += matrix_plus[i][j];
    }
    return;
  };
  
  if(nrowInput[5] != ncolInput[3])
  {
    ErrThrow("dimensions mismatch");
  }
  unsigned int r = nrowInput[4];
  unsigned int c = ncolInput[4];
  double **temp = new double*[nrowInput[0]];
  for(int i=0; i<nrowInput[0]; ++i)
    temp[i]=new double[ncolInput[0]];
  for(int i=0; i<nrowInput[0]; ++i)
    memset(temp[i], 0, ncolInput[0]*sizeof(double));
  // P0*A00nl + A00(linear);
  matrix_times_matrix_plus_matrix(r, c, inputMat[4], inputMat[2], inputMat[0], 
                                  outputMat[0]);
  // P1 * A00nl 
  matrix_times_matrix_plus_matrix(r, c, inputMat[4], inputMat[3], temp, 
                                  outputMat[1]);
  // P0*A11nl;
  matrix_times_matrix_plus_matrix(r, c, inputMat[5], inputMat[2], temp, 
                                  outputMat[2]);
  // P1 * A11nl + A11(linear) 
  matrix_times_matrix_plus_matrix(r, c, inputMat[5], inputMat[3], inputMat[1], 
                                  outputMat[3]);
  delete [] temp[0];
  delete [] temp;
}


#endif
