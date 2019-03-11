#include "LocalAssemblingPR.h"
#include "../include/General/Database.h"
#include "Extra_Routines_PR.h"
#include "../include/AssembleRoutines/LocalAssembling.h"
#include "string.h"

void press_robust_only_rhs(double Mult, double *coeff, double *param, double hK, 
double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  int N_U = N_BaseFuncts[1];
  double *Orig0 = OrigValues[1];

  double *Rhs = LocRhs[0];

  double c1 = coeff[1];
  double c2 = coeff[2];

  for(int i=0;i<N_U;i++)
  {
    int sign = GetSignOfThisDOF(N_U, i);
    double testx00 = sign*Orig0[i];
    double testy00 = sign*Orig0[N_U+i];

    Rhs[i] += Mult*( testx00*c1 + testy00*c2);
  }
}

void locAssemble_mass_laplacian(double Mult, double *coeff, double *param, double hK,
 double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0]; // sv
  double **MatrixA22 = LocMatrices[1];
  double **MatrixM   = LocMatrices[2]; // vv 
  // double **MatrixA11NL = LocMatrices[3]; // rectangular matrix
  // double **MatrixA22NL = LocMatrices[4]; // rectangular matrix
  
  double *Rhs = LocRhs[0];
  
  int scalar_nbf = N_BaseFuncts[0]; // nbasis for V_h
  int vector_nbf = N_BaseFuncts[1]; // nbasis for BDM
  
  double *Orig0  = OrigValues[0]; // u_x
  double *Orig1  = OrigValues[1]; // u_y
  // double *Orig2  = OrigValues[2]; // u
  double *OrigV0 = OrigValues[3]; // vector valued

  double nu = coeff[0];  // nu
  double f1 = coeff[1];  // f1
  double f2 = coeff[2];  // f2

  double ansatz10, ansatz01;
  double test10, test01, val;

  for(int i=0;i<scalar_nbf; i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];    

    for(int j=0;j<scalar_nbf;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];      

      val  = nu*(test10*ansatz10+test01*ansatz01);
      MatrixA11[i][j] += Mult * val;
      
      MatrixA22[i][j] += Mult * val;
    }
  }

  std::vector<int> sign(vector_nbf);
  for(int i=0;i<vector_nbf;i++)
    sign[i] = GetSignOfThisDOF(vector_nbf, i);
  
  double testx00, testy00;
  for(int i=0; i<vector_nbf; i++)
  {
    testx00 = sign[i] * OrigV0[i];
    testy00 = sign[i] * OrigV0[vector_nbf+i];
    
    // assmeble rhs 
    Rhs[i] += Mult*( testx00*f1 + testy00*f2);

    // mass matrix
    for(int j=0; j<vector_nbf; j++)
    {
      double ansatzx00 = sign[j]*OrigV0[j];
      double ansatzy00 = sign[j]*OrigV0[j+vector_nbf];
      
      val = testx00*ansatzx00 + testy00*ansatzy00;
      MatrixM[i][j] += Mult*val;
    }
  }  
}

LocalAssemblingPR::LocalAssemblingPR(ParameterDatabase db, CoeffFct2D coeffs, std::string name)
 : LocalAssembling<2>(db)
{
  this->N_Parameters = 0;
  this->N_ParamFct = 0;
  this->ParameterFct = {};
  this->N_FEValues = 0;
  this->FEValue_FctIndex = {};
  this->FEValue_MultiIndex = {};
  this->BeginParameter = {};
  
    
  this->Needs2ndDerivatives = new bool[3];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->Needs2ndDerivatives[2] = false;
  this->Manipulate = nullptr;
  
  if(name == std::string("pr_robust_mass_rhs"))
  {
    N_Rhs = 1;
    RhsSpace = {0};
    N_Terms = 2;
    Derivatives   = {D00, D00};
    FESpaceNumber = {0, 1};
    N_Matrices = 2;
    RowSpace    = { 0, 0, 1};
    ColumnSpace = { 1, 1, 1};
    using namespace std::placeholders;
    
    local_assemblings_routines.push_back(press_robust_only_rhs);
    
    AllOrigValues = new double** [N_Terms];
    OrigValues = new double* [N_Terms];
  }
  
  if(name == std::string("pr_robust_mass_stiff_rhs"))
  {
    
  }
  
}
