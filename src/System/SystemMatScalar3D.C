/** ************************************************************************ 
* @brief     source file for TSystemMatScalar3D
* @author    Sashikumaar Ganesan
* @date      23.01.15
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemMatScalar3D.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <Assemble3D.h>
#include <AuxParam3D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <AssembleMat3D.h>

#include <stdlib.h>
#include <string.h>

#include <Solver.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>

#ifdef _MPI
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

TSystemMatScalar3D::TSystemMatScalar3D(TFESpace3D *fespace)
 : SystemMat3D(1, 1, 0)
{
  //set number of multigrid levels
  int N_Levels = TDatabase::ParamDB->LEVELS;
  if(TDatabase::ParamDB->SC_MG_TYPE_SCALAR)
    ++N_Levels;
  
  //store the FEspace
  this->SystemMat3D::fe_spaces[0] = fespace;
  
  // build matrices
  // first build matrix structure
  TSquareStructure3D *sqstructure = new TSquareStructure3D(fespace);
  if(TDatabase::ParamDB->SOLVER_TYPE == DIRECT 
      || TDatabase::ParamDB->SOLVER_TYPE == GMG)
  //instance of the Assemble class
  AMatRhsAssemble = new TAssembleMat3D *[N_Levels];
  
  if(TDatabase::ParamDB->SOLVER_TYPE == DIRECT
      || TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    sqstructure->Sort();
  } // sort column numbers: numbers are in increasing order
  else if(TDatabase::ParamDB->SOLVER_TYPE == AMG_SOLVE)
  {
    sqstructure->SortDiagFirst();
  }
  this->SystemMat3D::sq_matrices[0] = new TSquareMatrix3D(sqstructure);
} //TSystemMatScalar3D::TSystemMatScalar3D

TSystemMatScalar3D::~TSystemMatScalar3D()
{
  delete this->SystemMat3D::sq_matrices[0]->GetStructure();
  delete this->SystemMat3D::sq_matrices[0];
}


void TSystemMatScalar3D::Init(CoeffFct3D *BilinearCoeffs,
                              BoundCondFunct3D *BoundCond,
                              BoundValueFunct3D *BoundValue)
{
  BoundaryConditions[0] = BoundCond;
  BoundaryValues[0] = BoundValue;
} // TSystemMatScalar3D::Init


void TSystemMatScalar3D::Assemble(CoeffFct3D *BilinearCoeffs, double *sol,
                                  double *rhs)
{
  TDiscreteForm3D *DiscreteFormGalerkin;
  InitializeDiscreteForms(DiscreteFormGalerkin, BilinearCoeffs);  
  TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
  
  Assemble3D(1, &this->fe_spaces[0], 1, &sq_matrices[0], 0, NULL, 1, &rhs,
             &fe_spaces[0], DiscreteFormGalerkin, BoundaryConditions,
             BoundaryValues, &aux);
  
  int N_DOF = this->SystemMat3D::fe_spaces[0]->GetN_DegreesOfFreedom();
  int N_Active = this->SystemMat3D::fe_spaces[0]->GetActiveBound();
  int N_DirichletDof = N_DOF - N_Active;
  // set rhs for Dirichlet nodes
  memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
}

// void TSystemMatScalar3D::Assemble()
// {
//   int i, N_DOF_low, N_Active;
// 
//    for(i=Start_Level;i<N_Levels;i++)
//     {    
//      N_DOF_low = FeSpaces[i]->GetN_DegreesOfFreedom();
//      N_Active =  FeSpaces[i]->GetActiveBound();
//     
//      // initialize matrices and rhs
//      AMatRhsAssemble[i]->Reset(); 
// 
//      // assemble
//      AMatRhsAssemble[i]->Assemble3D();
//   
//      // set rhs for Dirichlet nodes
//      memcpy(SolArray[i]+N_Active, RhsArray[i]+N_Active, (N_DOF_low - N_Active)*SizeOfDouble);   
//     } //  for(i=Start_Level;i<N_Levels;i++)    
// 
// //have to shift this in pardirectsolver    
// #ifdef _OMPONLY     
//     if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
//       DS->AssembleMatrix(sqmatrixA[N_Levels-1]);
// #endif
//     
// } // void TSystemMatScalar3D::Assemble(T


void TSystemMatScalar3D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->SystemMat3D::sq_matrices[0]->GetN_Rows();
  // reset y
  memset(y, 0.0, n_total_rows * SizeOfDouble);
  this->apply_scaled_add(x, y, factor);
}

void TSystemMatScalar3D::apply_scaled_add(const double *x, double *y,
                          double factor) const
{
  this->SystemMat3D::sq_matrices[0]->multiply(x, y, factor);
}
