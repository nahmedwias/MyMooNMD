/** ************************************************************************ 
* @brief     source file for BlockMatrixCD2D
* @author    Sashikumaar Ganesan, 
* @date      08.08.14
* @History 
 ************************************************************************  */
#ifdef __2D__
#include <Database.h>
#include <BlockMatrixCD2D.h>
#include <SquareStructure2D.h>
#include <Assemble2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>

BlockMatrixCD2D::BlockMatrixCD2D(TFESpace2D *fespace)
 : BlockMatrix2D(1, 1, 0)
{
  //store the FEspace
  this->BlockMatrix2D::fe_spaces[0] = fespace;
  
  // build matrices
  // first build matrix structure
  TSquareStructure2D* sqstructure = new TSquareStructure2D(fespace);
  sqstructure->Sort();  // sort column numbers: numbers are in increasing order

  /** A is the stiffness/system matrix for a stationary convection diffusion 
   * problem */
  this->BlockMatrix2D::sq_matrices[0] = new TSquareMatrix2D(sqstructure);  
  this->BlockMatrix2D::defect = Defect_Scalar;
}

BlockMatrixCD2D::~BlockMatrixCD2D()
{
  delete this->BlockMatrix2D::sq_matrices[0]->GetStructure();
  delete this->BlockMatrix2D::sq_matrices[0];
}
  
  
void BlockMatrixCD2D::Init(BoundCondFunct2D *BoundCond,
                              BoundValueFunct2D *BoundValue)
{
  //TDiscreteForm2D *DiscreteFormHeatLine;
  this->BoundaryConditions[0] = BoundCond;//should be the same as in fe_space[0]
  this->BoundaryValues[0] = BoundValue;
} // BlockMatrixCD2D::Init


void BlockMatrixCD2D::Assemble(LocalAssembling2D& la, double *sol,
                                  double *rhs)
{
  int N_DOF = this->BlockMatrix2D::fe_spaces[0]->GetN_DegreesOfFreedom();
  int N_Active = this->BlockMatrix2D::fe_spaces[0]->GetActiveBound();
  int N_DirichletDof = N_DOF - N_Active;
  
  // reset right hand side and matrix to zero
  memset(rhs, 0, N_DOF*SizeOfDouble);
  this->BlockMatrix2D::sq_matrices[0]->Reset();
  
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
} // void BlockMatrixCD2D::Assemble


void BlockMatrixCD2D::Solve(double *sol, double *rhs)
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
     DirectSolver(this->BlockMatrix2D::sq_matrices[0], rhs, sol);
   break;
 
   default:
     OutPut("Unknown Solver" << endl);
     exit(4711);;
  }
}

void BlockMatrixCD2D::apply(const double *x, double *y, double factor) const
{
  unsigned int n_total_rows = this->BlockMatrix2D::sq_matrices[0]->GetN_Rows();
  // reset y
  memset(y, 0.0, n_total_rows*SizeOfDouble);
  this->apply_scaled_add(x, y, factor);
}
    
void BlockMatrixCD2D::apply_scaled_add(const double *x, double *y, 
                                          double factor) const
{
  this->BlockMatrix2D::sq_matrices[0]->multiply(x, y, factor);
}

#endif // #ifdef __2D__





