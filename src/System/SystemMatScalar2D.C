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

TSystemMatScalar2D::TSystemMatScalar2D(TFESpace2D *fespace, int disctype, int solver)
{
  //store the FEspace
  FeSpace = fespace;
  
  //set the discretization type
  Disctype = disctype;
  
  //set the solver type
  SOLVER = solver;
  
  // build matrices
  // first build matrix structure
  sqstructure = new TSquareStructure2D(fespace);
  sqstructure->Sort();  // sort column numbers: numbers are in increasing order

  /** A is the stiffness/system matrix for a stationary convection diffusion 
   * problem */
  sqmatrixA = new TSquareMatrix2D(sqstructure);  
  N_Matrices = 1;
}

TSystemMatScalar2D::~TSystemMatScalar2D()
{
  delete sqstructure;
  delete sqmatrixA;
}
  
  
void TSystemMatScalar2D::Init(BoundCondFunct2D *BoundCond,
                              BoundValueFunct2D *BoundValue)
{
  BoundaryConditions[0] =  BoundCond;
  BoundaryValues[0] = BoundValue;
} // TSystemMatScalar2D::Init


void TSystemMatScalar2D::Assemble(LocalAssembling2D& la, double *sol, double *rhs)
{
  int N_DOF = FeSpace->GetN_DegreesOfFreedom();
  int N_Active =  FeSpace->GetActiveBound();
  int N_DirichletDof = N_DOF - N_Active;
  
  memset(rhs, 0, N_DOF*SizeOfDouble);
 
  // initialize matrices
  sqmatrixA->Reset();
  
  // assemble
  Assemble2D(1, &FeSpace, N_Matrices, &sqmatrixA, 0, NULL, 1, &rhs, &FeSpace,
             BoundaryConditions, BoundaryValues, la);
 
  // apply local projection stabilization method
  if(Disctype==LOCAL_PROJECTION && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
  {
    if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
    { 
      UltraLocalProjection((void *)sqmatrixA, false);
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
  switch(SOLVER)
  {
   case AMG_SOLVE:
     cout << "AMG_SOLVE not yet implemented " <<endl;
   break;

   case GMG:
     cout << "GMG solver not yet implemented " <<endl;
   break;

   case DIRECT:
     DirectSolver(sqmatrixA, rhs, sol);
   break;      
 
   default:
     OutPut("Unknown Solver" << endl);
     exit(4711);;
  }
}









