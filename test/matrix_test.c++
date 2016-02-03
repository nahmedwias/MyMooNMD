#include <Matrix.h>
#include <MooNMD_Io.h>
#include <cmath>

bool equal(const double a, const double b)
{
  if(fabs(a) > 1e-12)
    return fabs((a-b)/a) < 1e-14;
  else
    return fabs(a-b) < 1e-14;
}

int main(int argc, char **argv) 
{
  /**
   * A = [ 0.5, 0,    0, 0  ; B = [ 4, 0, 0, 1; C = [ 7, 0, 8;
   *       0,   0.75, 0, 0  ;       0, 0, 5, 0;       9,10, 0;
   *       0,   0,    1, 0  ;       0, 6, 0, 1 ]      0,11,12;
   *       0,   0,    0, 1.5 ]                        13,0,14 ]
   *
   * D = [ 0.5 0   0  ;
   *       0   1   0.1;
   *       0   0.1 1.5 ]
   */
  
  //Matrix A
  int * RowA = new int[5];
  int * ColA = new int[4];
  ColA[0] = 0; ColA[1] = 1; ColA[2] =  2; ColA[3] = 3;
  RowA[0] = 0; RowA[1] = 1; RowA[2] = 2; RowA[3] = 3; RowA[4] = 4;
  std::shared_ptr<TStructure> structureA(new TStructure(4, 4, ColA, RowA));
  TMatrix matA(structureA);
  matA.GetEntries()[0] = 0.5;
  matA.GetEntries()[1] = 0.75;
  matA.GetEntries()[2] = 1;
  matA.GetEntries()[3] = 1.5;

  //Matrix B
  int * RowB = new int[4];
  int * ColB = new int[5];
  RowB[0] = 0; RowB[1] = 2; RowB[2] = 3; RowB[3] = 5;
  ColB[0] = 0; ColB[1] = 3; ColB[2] = 2; ColB[3] = 1; ColB[4] = 3;
  std::shared_ptr<TStructure> structureB(new TStructure(3, 4, 5, ColB, RowB));
  TMatrix matB(structureB);
  matB.GetEntries()[0] = 4;
  matB.GetEntries()[1] = 1;
  matB.GetEntries()[2] = 5;
  matB.GetEntries()[3] = 6;
  matB.GetEntries()[4] = 1;
  

  //Matrix C
  int * RowC = new int[5];
  int * ColC = new int[8];
  RowC[0] = 0; RowC[1] = 2; RowC[2] = 4; RowC[3] = 6; RowC[4] = 8;
  ColC[0] = 0; ColC[1] = 2; ColC[2] = 0; ColC[3] = 1; ColC[4] = 1; ColC[5] = 2;
  ColC[6] = 0; ColC[7] = 2;
  std::shared_ptr<TStructure> structureC(new TStructure(4, 3, 8, ColC, RowC));
  TMatrix matC(structureC);
  matC.GetEntries()[0] = 7;
  matC.GetEntries()[1] = 8;
  matC.GetEntries()[2] = 9;
  matC.GetEntries()[3] = 10;
  matC.GetEntries()[4] = 11;
  matC.GetEntries()[5] = 12;
  matC.GetEntries()[6] = 13;
  matC.GetEntries()[7] = 14;

  TMatrix* matCT = matC.GetTransposed();
  
  //Matrix D
  int * RowD = new int [4];
  int * ColD = new int [5];
  RowD[0] = 0; RowD[1] = 1; RowD[2] = 3; RowD[3] = 5;
  ColD[0] = 0; ColD[1] = 1; ColD[2] = 2; ColD[3] = 1; ColD[4] = 2;
  std::shared_ptr<TStructure> structureD(new TStructure(3, 3, 5, ColD, RowD));
  TMatrix matD(structureD);
  matD.GetEntries()[0] = 0.5;
  matD.GetEntries()[1] = 1;
  matD.GetEntries()[2] = 0.1;
  matD.GetEntries()[3] = 0.1;
  matD.GetEntries()[4] = 1.5;

  matA.PrintFull("matA", 4);
  matB.PrintFull("matB", 4);
  matC.PrintFull("matC", 4);
  matD.PrintFull("matD", 4);
  
  std::shared_ptr<TMatrix> matAC( matA.multiply(matC));
  std::shared_ptr<TMatrix> matBA( matB.multiply(matA));
  std::shared_ptr<TMatrix> matCB( matC.multiply(matB));
  
  matAC->PrintFull("matAC");
  matBA->PrintFull("matBA");
  matCB->PrintFull("matCB");
  
  
  // ##########################################################################
  // ##########################################################################
  // some tests
  {
    if(matAC->GetN_Rows() != 4 || matAC->GetN_Columns() != 3)
      ErrThrow("wrong dimension of product matrix AC");
    if(matBA->GetN_Rows() != 3 || matBA->GetN_Columns() != 4)
      ErrThrow("wrong dimension of product matrix BA");
    if(matCB->GetN_Rows() != 4 || matCB->GetN_Columns() != 4)
      ErrThrow("wrong dimension of product matrix CB");
    
    if(!equal(matA.GetNorm(-1), 1.5)
       || !equal(matA.GetNorm(-2), 2.015564437074637))
      ErrThrow("wrong matrix norm, matrix A");
    if(!equal(matB.GetNorm(-1), 7)
       || !equal(matB.GetNorm(-2), 8.888194417315589))
      ErrThrow("wrong matrix norm, matrix B");
    if(!equal(matC.GetNorm(-1), 27)
       || !equal(matC.GetNorm(-2), 30.39736830714133))
      ErrThrow("wrong matrix norm, matrix C");
    if(!equal(matCT->GetNorm(-1), 34)
       || !equal(matCT->GetNorm(-2), 30.39736830714133))
      ErrThrow("wrong matrix norm, matrix CT");
    if(!equal(matAC->GetNorm(-1), 40.5)
       || !equal(matAC->GetNorm(-2), 34.87567203653573))
      ErrThrow("wrong matrix norm, matrix AC");
    if(!equal(matBA->GetNorm(-1), 6) 
      || !equal(matBA->GetNorm(-2), 7.33143914930759))
      ErrThrow("wrong matrix norm, matrix BA");
    if(!equal(matCB->GetNorm(-1), 163)
       || !equal(matCB->GetNorm(-2), 161.3443522407896))
      ErrThrow("wrong matrix norm, matrix CB");
  }
  // ##########################################################################
  // ##########################################################################
  // checking the multiplication with a vector:
  {
    double * x = new double[4];
    double * y = new double[4];
    x[0] = 1.; x[1] = -2.5; x[2] = 11.34; x[3] = -0.004; // some data
    y[0] = 0.; y[1] = 0.; y[2] = 0.; y[3] = 0.;
    
    // y = -0.5*CB*x
    matCB->multiply(x, y, -0.5);
    if(!equal(y[0], 46.03) || !equal(y[1], -301.482) || !equal(y[2], -221.826)
       || !equal(y[3], 79.054))
      ErrThrow("wrong result in matrix vector multiplication, matrix CB");
    y[0] = 0.; y[1] = 0.; y[2] = 0.; y[3] = 0.; // reset y
    
    // y = 0.15*matC^Tx
    matC.transpose_multiply(x, y, 0.15); // last entry in y is not used
    if(!equal(y[0], -2.3328) || !equal(y[1], 14.961) || !equal(y[2], 21.6036)
       || !equal(y[3], 0.))
      ErrThrow("wrong result in (transposed) matrix vector multiplication, C");
    // y += -0.15*matCT*x
    matCT->multiply(x, y, -0.15); // y should be zero now
    if(!equal(y[0], 0.0) || !equal(y[1], 0.0) || !equal(y[2], 0.0)
       || !equal(y[3], 0.))
      ErrThrow("C and CT seem to be not transposed to each other");
    
    delete x; x = nullptr;
    delete y; y = nullptr;
  }
  // ##########################################################################
  // ##########################################################################
  // scaling
  {
    double norm = matD.GetNorm(-2); // Frobenius norm
    double scaling = 6.451;
    matD.scale(scaling); // scale all matrix entries
    if(!equal(matD.GetNorm(-2), norm * scaling))
      ErrMsg("wrong result after scaling a matrix, matrix D");
    
    TMatrix matB_copy(matB);
    
    double * x = new double[4];
    x[0] = 1.; x[1] = -2.; x[2] = 3.; x[3] = -4.;
    matB_copy.scale(x, true); // x[3] is not used
    if(!equal(matB_copy(0, 0), 4) || !equal(matB_copy(0, 3), 1)
       || !equal(matB_copy(1, 2), -10) || !equal(matB_copy(2, 1), 18)
       || !equal(matB_copy(2, 3), 3))
      ErrThrow("wrong results after scaling with vector, matrix B");
    matB_copy.scale(x, false);
    if(!equal(matB_copy(0, 0), 4) || !equal(matB_copy(0, 3), -4)
       || !equal(matB_copy(1, 2), -30)  || !equal(matB_copy(2, 1), -36)
       || !equal(matB_copy(2, 3), -12))
      ErrThrow("wrong results after scaling with vector, matrix B");
  }
  // ##########################################################################
  // ##########################################################################
  // matrix matrix^T multiplication
  {
    TMatrix *matBBT = matB.multiply_with_transpose_from_right();
    matBBT->PrintFull("BBT");
    if(!equal(matBBT->get(0, 0), 17) || !equal(matBBT->get(0, 2), 1)
       || !equal(matBBT->get(1, 1), 25) || !equal(matBBT->get(2, 0), 1)
       || !equal(matBBT->get(2, 2), 37))
      ErrThrow("wrong results after multiplication with transpose from right, ",
               "matrix B");
    // multiply with known structure
    std::vector<double> scaleB{1.0, -2.0, 3.0, -1.0};
    TMatrix * matBBT2 = matB.multiply_with_transpose_from_right(
      scaleB, matBBT->GetStructure());
    matBBT2->PrintFull("BBT2");
    if(!equal(matBBT2->get(0, 0), 15) || !equal(matBBT2->get(0, 2), -1)
       || !equal(matBBT2->get(1, 1), 75) || !equal(matBBT2->get(2, 0), -1)
       || !equal(matBBT2->get(2, 2), -73))
      ErrThrow("wrong results after multiplication with transpose from right ",
               "and known structure, matrix B");
    
    delete matBBT2;
    delete matBBT;
    
    TMatrix * matAAT = matA.multiply_with_transpose_from_right();
    matAAT->PrintFull("AAT");
    if(!equal(matAAT->get(0, 0), 0.25) || !equal(matAAT->get(1, 1), 0.5625)
       || !equal(matAAT->get(2, 2), 1) || !equal(matAAT->get(3, 3), 2.25) )
      ErrThrow("wrong results after multiplication with transpose from right, ",
               "matrix A");
    delete matAAT;
    
    std::vector<double> scale{1.0, 2.0, -3.0};
    TMatrix * matCCT = matC.multiply_with_transpose_from_right(scale);
    if(!equal(matAAT->get(0, 0), -143) || !equal(matAAT->get(0, 1), 63)
       || !equal(matAAT->get(0, 2), -288) || !equal(matAAT->get(0, 3), -245) )
      ErrThrow("wrong results after multiplication with transpose from right, ",
               "matrix C, first row");
    if(!equal(matAAT->get(1, 0), 63) || !equal(matAAT->get(1, 1), 281)
       || !equal(matAAT->get(1, 2), 220) || !equal(matAAT->get(1, 3), 117) )
      ErrThrow("wrong results after multiplication with transpose from right, ",
               "matrix C, second row");
    if(!equal(matAAT->get(2, 0), -288) || !equal(matAAT->get(2, 1), 220)
       || !equal(matAAT->get(2, 2), -190) || !equal(matAAT->get(2, 3), -504) )
      ErrThrow("wrong results after multiplication with transpose from right, ",
               "matrix C, third row");
    if(!equal(matAAT->get(3, 0), -245) || !equal(matAAT->get(3, 1), 117)
       || !equal(matAAT->get(3, 2), -504) || !equal(matAAT->get(3, 3), -419) )
      ErrThrow("wrong results after multiplication with transpose from right, ",
               "matrix C, fourth row");
    matCCT->PrintFull("CCT");
    delete matCCT;
  }
  // ##########################################################################
  // ##########################################################################    
  {
    //H =[1  0 0 0 ;
    //    0  3 0 4 ;
    //    0  0 0 6 ;
    int RowH[4] = {0, 1, 3, 4};
    int ColH[4] = {0, 1, 3, 3};
    std::shared_ptr<TStructure> structureE(new TStructure(3, 4, 4, ColH, RowH));
    TMatrix matH(structureE);
    matH.GetEntries()[0] = 1;
    matH.GetEntries()[1] = 3;
    matH.GetEntries()[2] = 4;
    matH.GetEntries()[3] = 6;
    
    matH.PrintFull("matH");
    
    // G =[    
    //      1 7 0 0 0; 
    //      0 2 0 0 0;
    //      5 0 0 9 0;
    //      0 0 0 0 4]
    //      
    
    int RowI[5] = {0, 2, 3, 5, 6};
    int ColI[6] = {0, 1, 1, 0, 3, 4} ;
    std::shared_ptr<TStructure> structureI(new TStructure(4, 4, 6, ColI, RowI));
    TMatrix matI(structureI);
    matI.GetEntries()[0] = 1;
    matI.GetEntries()[1] = 7;
    matI.GetEntries()[2] = 2;
    matI.GetEntries()[3] = 5;
    matI.GetEntries()[4] = 9;
    matI.GetEntries()[5] = 4;
    matI.PrintFull("I");
    std::shared_ptr<TMatrix> matHGHT(matH.multiply_with_transpose_from_right(matI));
    matHGHT->PrintFull("matHGHT");
    
    if(!equal(matHGHT->get(0, 0), 1) || !equal(matHGHT->get(0, 1), 21)
       || !equal(matHGHT->get(1, 1), 18))
      ErrThrow("wrong results after multiplication with transpose from right, ",
               "matrix H, and matrix I in between");
  }
  
  delete matCT;
  
  std::cout << "test successful\n";
  return 0;
}
