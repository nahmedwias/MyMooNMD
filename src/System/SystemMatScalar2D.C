/** ************************************************************************ 
* @brief     source file for TSystemMatScalar2D
* @author    Sashikumaar Ganesan, 
* @date      08.08.14
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemMatScalar2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>

TSystemMatScalar2D::TSystemMatScalar2D(TFESpace2D *fespace)
 : SystemMat2D(1, 1, 0)
{
  //store the FEspace
  this->SystemMat2D::fe_spaces[0] = fespace;
  
  // build matrices
  // first build matrix structure
  TSquareStructure2D* sqstructure = new TSquareStructure2D(fespace);
  sqstructure->Sort();  // sort column numbers: numbers are in increasing order

  /** A is the stiffness/system matrix for a stationary convection diffusion 
   * problem */
  this->SystemMat2D::sq_matrices[0] = new TSquareMatrix2D(sqstructure);  
  this->SystemMat2D::defect = Defect_Scalar;
}

TSystemMatScalar2D::~TSystemMatScalar2D()
{
  delete this->SystemMat2D::sq_matrices[0]->GetStructure();
  delete this->SystemMat2D::sq_matrices[0];
}
  
  
void TSystemMatScalar2D::Init(BoundCondFunct2D *BoundCond,
                              BoundValueFunct2D *BoundValue)
{
  this->BoundaryConditions[0] = BoundCond;//should be the same as in fe_space[0]
  this->BoundaryValues[0] = BoundValue;
} // TSystemMatScalar2D::Init


void TSystemMatScalar2D::Assemble(LocalAssembling2D& la, double *sol,
                                  double *rhs)
{
  int N_DOF = this->SystemMat2D::fe_spaces[0]->GetN_DegreesOfFreedom();
  int N_Active = this->SystemMat2D::fe_spaces[0]->GetActiveBound();
  int N_DirichletDof = N_DOF - N_Active;
  
  // reset right hand side and matrix to zero
  memset(rhs, 0, N_DOF*SizeOfDouble);
  this->SystemMat2D::sq_matrices[0]->Reset();
  
  int N_Matrices = 1;
  // assemble
  Assemble2D(1, &fe_spaces[0], N_Matrices, &sq_matrices[0], 0, NULL, 1, &rhs, 
             &fe_spaces[0], BoundaryConditions, BoundaryValues, la);
 
  // apply local projection stabilization method
  if(TDatabase::ParamDB->DISCTYPE==LOCAL_PROJECTION 
     && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
  {
    if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
    { 
      UltraLocalProjection((void *)&sq_matrices[0], false);
    }
    else
    {
      OutPut("Check! LP_FULL_GRADIENT needs to be one to use LOCAL_PROJECTION" << endl);
      exit(4711);
    }
  }
  
  // copy Dirichlet values from rhs to solution vector (this is not really 
  // necessary in case of a direct solver)
  memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
} // void TSystemMatScalar2D::Assemble(T


void TSystemMatScalar2D::Solve(double *sol, double *rhs)
{
  switch(TDatabase::ParamDB->SOLVER_TYPE)
  {
   case AMG_SOLVE:
     cout << "AMG_SOLVE not yet implemented " <<endl;
   break;

   case GMG:
     cout << "GMG solver not yet implemented " <<endl;
   break;

   case DIRECT:
     DirectSolver(this->SystemMat2D::sq_matrices[0], rhs, sol);
   break;
 
   default:
     OutPut("Unknown Solver" << endl);
     exit(4711);;
  }
}

void TSystemMatScalar2D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->SystemMat2D::sq_matrices[0]->GetN_Rows();
  // reset y
  memset(y, 0.0, n_total_rows*SizeOfDouble);
  this->apply_scaled_add(x, y, factor);
}
    
void TSystemMatScalar2D::apply_scaled_add(const double *x, double *y, 
                                          double factor) const
{
  this->SystemMat2D::sq_matrices[0]->multiply(x, y, factor);
}







