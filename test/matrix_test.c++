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
  double * EntriesA = new double[4];
  ColA[0] = 0; ColA[1] = 1; ColA[2] =  2; ColA[3] = 3;
  RowA[0] = 0; RowA[1] = 1; RowA[2] = 2; RowA[3] = 3; RowA[4] = 4;
  EntriesA[0] = 0.5; EntriesA[1] = 0.75; EntriesA[2] = 1; EntriesA[3] = 1.5;
  TStructure structureA(4, 4, ColA, RowA);
  TMatrix matA(&structureA, EntriesA);

  //Matrix B
  int * RowB = new int[4];
  int * ColB = new int[5];
  double * EntriesB = new double[5];
  RowB[0] = 0; RowB[1] = 2; RowB[2] = 3; RowB[3] = 5;
  ColB[0] = 0; ColB[1] = 3; ColB[2] = 2; ColB[3] = 1; ColB[4] = 3;
  EntriesB[0] = 4; EntriesB[1] = 1; EntriesB[2] = 5; EntriesB[3] = 6; 
  EntriesB[4] = 1;
  TStructure structureB(3, 4, 5, ColB, RowB);
  TMatrix matB(&structureB, EntriesB);

  //Matrix C
  int * RowC = new int[5];
  int * ColC = new int[8];
  double * EntriesC = new double[8];
  RowC[0] = 0; RowC[1] = 2; RowC[2] = 4; RowC[3] = 6; RowC[4] = 8;
  ColC[0] = 0; ColC[1] = 2; ColC[2] = 0; ColC[3] = 1; ColC[4] = 1; ColC[5] = 2;
  ColC[6] = 0; ColC[7] = 2;
  EntriesC[0] = 7; EntriesC[1] = 8; EntriesC[2] = 9; EntriesC[3] = 10;
  EntriesC[4] = 11; EntriesC[5] = 12; EntriesC[6] = 13; EntriesC[7] = 14;
  TStructure structureC(4, 3, 8, ColC, RowC);
  TMatrix matC(&structureC, EntriesC);

  TMatrix* matCT = matC.GetTransposed();
  
  //Matrix D
  int * RowD = new int [4];
  int * ColD = new int [5];
  double * EntriesD = new double[5];
  RowD[0] = 0; RowD[1] = 1; RowD[2] = 3; RowD[3] = 5;
  ColD[0] = 0; ColD[1] = 1; ColD[2] = 2; ColD[3] = 1; ColD[4] = 2;
  EntriesD[0] = 0.5; EntriesD[1] = 1; EntriesD[2] = 0.1; EntriesD[3] = 0.1; 
  EntriesD[4] = 1.5;
  TStructure structureD(3, 3, 5, ColD, RowD);
  TMatrix matD(&structureD, EntriesD);

  matA.PrintFull("matA", 4);
  matB.PrintFull("matB", 4);
  matC.PrintFull("matC", 4);
  matD.PrintFull("matD", 4);
  
  TMatrix* matAC = matA.multiply(&matC);
  TMatrix* matBA = matB.multiply(&matA);
  TMatrix* matCB = matC.multiply(&matB);
  
  matAC->PrintFull("matAC");
  matBA->PrintFull("matBA");
  matCB->PrintFull("matCB");
  
  // some tests
  if(matAC->GetN_Rows() != 4 || matAC->GetN_Columns() != 3)
    ErrThrow("wrong dimension of product matrix AC");
  if(matBA->GetN_Rows() != 3 || matBA->GetN_Columns() != 4)
    ErrThrow("wrong dimension of product matrix BA");
  if(matCB->GetN_Rows() != 4 || matCB->GetN_Columns() != 4)
    ErrThrow("wrong dimension of product matrix CB");
  
  if(!equal(matA.GetNorm(-1), 1.5) || !equal(matA.GetNorm(-2), 2.015564437074637))
    ErrThrow("wrong matrix norm, matrix A");
  if(!equal(matB.GetNorm(-1), 7) || !equal(matB.GetNorm(-2), 8.888194417315589))
    ErrThrow("wrong matrix norm, matrix B");
  if(!equal(matC.GetNorm(-1), 27) || !equal(matC.GetNorm(-2), 30.39736830714133))
    ErrThrow("wrong matrix norm, matrix C");
  if(!equal(matCT->GetNorm(-1), 34) || !equal(matCT->GetNorm(-2), 30.39736830714133))
    ErrThrow("wrong matrix norm, matrix CT");
  if(!equal(matAC->GetNorm(-1), 40.5) || !equal(matAC->GetNorm(-2), 34.87567203653573))
    ErrThrow("wrong matrix norm, matrix AC");
  if(!equal(matBA->GetNorm(-1), 6) || !equal(matBA->GetNorm(-2), 7.33143914930759))
    ErrThrow("wrong matrix norm, matrix BA");
  if(!equal(matCB->GetNorm(-1), 163) || !equal(matCB->GetNorm(-2), 161.3443522407896))
    ErrThrow("wrong matrix norm, matrix CB");
  
  // checking the multiplication with a vector:
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
  
  // scaling
  double norm = matD.GetNorm(-2); // Frobenius norm
  double scaling = 6.451;
  matD.scale(scaling); // scale all matrix entries
  if(!equal(matD.GetNorm(-2), norm * scaling))
    ErrMsg("wrong result after scaling a matrix, matrix D");
  
  x = new double[4];
  x[0] = 1.; x[1] = -2.; x[2] = 3.; x[3] = -4.;
  matB.scale(x, true); // x[3] is not used
  if(!equal(matB(0, 0), 4) || !equal(matB(0, 3), 1) || !equal(matB(1, 2), -10) 
     || !equal(matB(2, 1), 18) || !equal(matB(2, 3), 3))
    ErrThrow("wrong resulst after scaling with vector, matrix B");
  matB.scale(x, false);
  if(!equal(matB(0, 0), 4) || !equal(matB(0, 3), -4) || !equal(matB(1, 2), -30) 
     || !equal(matB(2, 1), -36) || !equal(matB(2, 3), -12))
    ErrThrow("wrong resulst after scaling with vector, matrix B");
  
  delete matCT->GetStructure();
  delete matCT;
  delete matAC->GetStructure();
  delete matAC;
  delete matBA->GetStructure();
  delete matBA;
  delete matCB->GetStructure();
  delete matCB;
  
  std::cout << "test successful\n";
  return 0;
}
