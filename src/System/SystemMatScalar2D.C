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

  /** A is the stiffness/system mat for stationary problem   */
  sqmatrixA = new TSquareMatrix2D(sqstructure);  
  N_Matrices = 1;

}

TSystemMatScalar2D::~TSystemMatScalar2D()
{
  delete sqstructure;
  delete sqmatrixA;
}
  
  
void TSystemMatScalar2D::Init(CoeffFct2D *BilinearCoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue)
{
  BoundaryConditions[0] =  BoundCond;
  BoundaryValues[0] = BoundValue;
    
  TDiscreteForm2D *DiscreteFormUpwind;  
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormGLS;  

  InitializeDiscreteForms_Stationary(DiscreteFormUpwind, DiscreteFormGalerkin, DiscreteFormSDFEM, DiscreteFormGLS,
                                     BilinearCoeffs);

    switch(Disctype)
     {
      case GALERKIN:
      case LOCAL_PROJECTION:
           DiscreteFormARhs = DiscreteFormGalerkin;
      break;

      case SUPG:
           DiscreteFormARhs = DiscreteFormSDFEM;
      break;

      case UPWIND:
           DiscreteFormARhs = DiscreteFormUpwind;
      break;      
      
      case GLS:
           DiscreteFormARhs = DiscreteFormGLS;
      break;
      
      default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
     }  
     
     
} // TSystemMatScalar2D::Init


void TSystemMatScalar2D::Assemble(TAuxParam2D *aux, double *sol, double *rhs)
{
  int N_DOF, N_Active, N_DirichletDof;
  double *RHSs[1];
 

  TFESpace2D *fesp[1], *ferhs[1];
//   TSquareMatrix2D *SQMATRICES[2];


   
    N_DOF = FeSpace->GetN_DegreesOfFreedom();
    N_Active =  FeSpace->GetActiveBound();
    N_DirichletDof = N_DOF - N_Active;
    
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);
  
    fesp[0] = FeSpace;
    ferhs[0] = FeSpace;
    
    // initialize matrices
    SQMATRICES[0] = sqmatrixA;
    SQMATRICES[0]->Reset(); 
    
    
    if(aux==NULL)
     { aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); }
    
    // assemble
    Assemble2D(1, fesp,
               N_Matrices, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteFormARhs,
               BoundaryConditions,
               BoundaryValues,
               aux);
 
     delete aux;
     
     // apply local projection stabilization method
     if(Disctype==LOCAL_PROJECTION && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
      {
       if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
        { 
         UltraLocalProjection((void *)SQMATRICES[0], FALSE);
        }
       else
        {
         OutPut("Check! LP_FULL_GRADIENT needs to be one to use LOCAL_PROJECTION" << endl);
         exit(4711);;
        }
      }
      
      
    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);     
      
//     cout << "Test Assemble " << endl;
    
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









