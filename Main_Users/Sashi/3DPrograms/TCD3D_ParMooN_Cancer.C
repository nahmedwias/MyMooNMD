// =======================================================================
//
// Purpose:     main program for time-dependent scalar equation with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 24.01.15

// =======================================================================
 
#include <Domain.h>
#include <Database.h>
#include <SystemMatTimeScalar3D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <LinAlg.h>
// #include <TCD3DErrorEstimator.h>
#include <MainUtilities.h>
#include <TimeUtilities.h>
#include <Solver.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

double bound = 0;
// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TCD_3D/Sin4.h"
// #include "../Examples/TCD_3D/ConstTSmooth.h"
//  #include "../Examples/TCD_3D/ConstT.h"
// #include "../Examples/TCD_3D/amc.h"
#include "../Examples_All/TCD_3D/cancer3D.h"

void AssembleCancerDensity(int Levels, TSquareMatrix3D **A, TFEFunction3D **u_FeFunction, TFEFunction3D **v_FeFunction)
{
  int ii, i, j, k, l, m, n, begin, end, N_Cells, N_LocalUsedElements;
  int N_Points, N_U_LocalDof, N_V_LocalDof, TestDOF;
  int *N_BaseFunct, *BeginIndex, *GlobalNumbers, *RowPtr, *KCol, *DOF, *VBeginIndex, *VGlobalNumbers, *VDOF;
  
  double v, *weights, *xi, *eta, *zeta, c0, lambda, *ValuesA;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double **uorig, **uxorig, **uyorig, **uzorig, *Orig, *x_Orig, *y_Orig, *z_Orig, U_value;
  double **vorig, **vxorig, **vyorig, **vzorig, *V_Orig, *Vx_Orig, *Vy_Orig, *Vz_Orig, V_value, Vx_value, Vy_value, Vz_value;
  double Mult, mass, *U_Sol, *V_Sol;
  double test000, test100, test010, test001, ansatz00;
  double uvfact, vgrad_test;
  double LocMatrixA[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  
  TFESpace3D *fespace, *V_fespace;
  TCollection *Coll;
  TBaseCell *Me;
  FE3D FEId, V_FEId;
  TFE3D *ele;
  FE3D LocalUsedElements[2];
  BaseFunct3D BaseFunct, V_BaseFunct, *BaseFuncts;
  
  bool *SecondDer;
  
  SecondDer = new bool[2];
  SecondDer[0] = FALSE;
  SecondDer[1] = FALSE;
  
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
  
//   U solution Space
  for(ii=0;ii<Levels; ii++)
   {
    fespace = A[ii]->GetFESpace();
    Coll = fespace->GetCollection();
    N_Cells = Coll->GetN_Cells();
    U_Sol = u_FeFunction[ii]->GetValues();
  
    BeginIndex = fespace->GetBeginIndex();
    GlobalNumbers = fespace->GetGlobalNumbers();
  
    RowPtr = A[ii]->GetRowPtr();
    KCol = A[ii]->GetKCol();
    ValuesA = A[ii]->GetEntries();
  
    //V solution Space
    V_Sol = v_FeFunction[ii]->GetValues();
    V_fespace = v_FeFunction[ii]->GetFESpace3D();
    VBeginIndex = V_fespace->GetBeginIndex();
    VGlobalNumbers = V_fespace->GetGlobalNumbers();
  
    N_LocalUsedElements = 2;
    c0 = TDatabase::ParamDB->REACTOR_P2; //chi
    lambda = TDatabase::ParamDB->REACTOR_P3;
  
    for(i=0;i<N_Cells;i++)
     {
      Me = Coll->GetCell(i);
      FEId = fespace->GetFE3D(i,Me);
      V_FEId = V_fespace->GetFE3D(i,Me);
    
      //Calculate values on original element
      LocalUsedElements[0]=FEId;
      LocalUsedElements[1]=V_FEId;
    
      TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                             Coll, Me, SecondDer, 
                             N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
    
      BaseFunct = BaseFuncts[FEId];
      N_U_LocalDof = N_BaseFunct[FEId];
    
      V_BaseFunct = BaseFuncts[V_FEId];
      N_V_LocalDof = N_BaseFunct[V_FEId];
    
      uorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
      uxorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
      uyorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
      uzorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
    
      vorig = TFEDatabase3D::GetOrigElementValues(V_BaseFunct, D000);
      vxorig = TFEDatabase3D::GetOrigElementValues(V_BaseFunct, D100);
      vyorig = TFEDatabase3D::GetOrigElementValues(V_BaseFunct, D010);
      vzorig = TFEDatabase3D::GetOrigElementValues(V_BaseFunct, D001);
    
      memset(LocMatrixA, 0, N_U_LocalDof*N_U_LocalDof*SizeOfDouble);
    
      DOF = GlobalNumbers + BeginIndex[i];
      VDOF = VGlobalNumbers + VBeginIndex[i];
    
      for(j=0;j<N_Points;j++)
       {
        Orig = uorig[j];
        x_Orig = uxorig[j];
        y_Orig = uyorig[j];
        z_Orig = uzorig[j];
      
        U_value = 0;
        for(l=0; l<N_U_LocalDof;l++)
         {
          U_value += U_Sol[DOF[l]]*Orig[l];
         }
      
        V_Orig = vorig[j];
        Vx_Orig = vxorig[j];
        Vy_Orig = vyorig[j];
        Vz_Orig = vzorig[j];
      
        V_value = 0.; Vx_value = 0.; Vy_value = 0.; Vz_value= 0.;
        for(l=0; l<N_V_LocalDof;l++)
         {
          v = V_Sol[VDOF[l]];
          V_value += v * V_Orig[l];
          Vx_value += v * Vx_Orig[l];
          Vy_value += v * Vy_Orig[l];
          Vz_value += v * Vz_Orig[l];
         }
      
        Mult = weights[j] * AbsDetjk[j];
      
       //assembling the local matrix
       for(k=0;k<N_U_LocalDof;k++)
        {
         test000 = Orig[k];
         test100 = x_Orig[k];
         test010 = y_Orig[k];
         test001 = z_Orig[k];

         vgrad_test = Mult * c0 * (Vx_value*test100 + Vy_value*test010 + Vz_value*test001);
         uvfact = Mult * lambda * (1.-U_value-V_value) * test000;

         for(l=0;l<N_U_LocalDof;l++)
          {
           ansatz00 = Orig[l];
           LocMatrixA[k*N_U_LocalDof + l] -= (vgrad_test + uvfact)*ansatz00;
          } //  for(l=0;l<N_U_LocalDof;l++
        } //  for(k=0;k<N_U_LocalDof;k++)
      } // for(j=0;j<N_Points;j++)
      
      //local to global matrix
     for(l=0;l<N_U_LocalDof;l++)
      {
       TestDOF = DOF[l];

       begin = RowPtr[TestDOF];
       end = RowPtr[TestDOF+1];

       for(n=begin;n<end;n++)
        {
         for(m=0;m<N_U_LocalDof;m++)
          {
           if(KCol[n] == DOF[m])
             {
              ValuesA[n] += LocMatrixA[l*N_U_LocalDof+m];      
              break;
             }
          }
         }
       } //  for(l=0;l<N_U_LocalDof;l++)
      } // for(i=0;i<N_Cells;i++)
    
//    cout << ii << " test assemble " << endl;

  } //  for(ii=0;ii<Levels; ii+
//      exit(0);
}

// void AssembleWArhs(TSquareMatrix3D *A, TFEFunction3D *u_FeFunction, double *RHS)
// {
//  int i, j, k, l, m, n, begin, end, N_Cells, N_LocalUsedElements;
//  int N_Points, N_U_LocalDof, N_V_LocalDof, TestDOF;
//  int *N_BaseFunct, *BeginIndex, *GlobalNumbers, *RowPtr, *KCol, *DOF, *UBeginIndex, *UGlobalNumbers, *UDOF;
//  
//  double *weights, *xi, *eta, *zeta, alpha, beta, *ValuesA;
//  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
//  double AbsDetjk[MaxN_QuadPoints_3D];
//  double **uorig, *Orig, *x_Orig, *y_Orig, *z_Orig, U_value;
//  double **vorig, *V_Orig, V_value;
//  double Mult, *U_Sol;
//  double test000, ansatz00;
//  double ufact, alphaU, alphaUbeta;
//  double LocMatrixA[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D], LocRhs[MaxN_BaseFunctions3D];
// 
//  TFESpace3D *fespace, *U_fespace;
//  TCollection *Coll;
//  TBaseCell *Me;
//  FE3D FEId, U_FEId;
//  TFE3D *ele;
//  FE3D LocalUsedElements[2];
//  BaseFunct3D BaseFunct, U_BaseFunct, *BaseFuncts;
// 
//  bool *SecondDer;
//  
//  SecondDer = new bool[2];
//  SecondDer[0] = FALSE;
//  SecondDer[1] = FALSE;
//  
//  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
//  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
//  
//  //W solution spaces
//  fespace = A->GetFESpace();
//  Coll = fespace->GetCollection();
//  N_Cells = Coll->GetN_Cells();
//  
//  BeginIndex = fespace->GetBeginIndex();
//  GlobalNumbers = fespace->GetGlobalNumbers();
//  
//  RowPtr = A->GetRowPtr();
//  KCol = A->GetKCol();
//  ValuesA = A->GetEntries();
//  
//  //U solution spaces
//  U_Sol = u_FeFunction->GetValues();
//  U_fespace = u_FeFunction->GetFESpace3D();
//  UBeginIndex = U_fespace->GetBeginIndex();
//  UGlobalNumbers = U_fespace->GetGlobalNumbers();
//  
//  N_LocalUsedElements = 2;
//  alpha = TDatabase::ParamDB->REACTOR_P5;
//  beta = TDatabase::ParamDB->REACTOR_P6;
//  
//  for(i=0;i<N_Cells;i++)
//  {
//    Me = Coll->GetCell(i);
//    
//    FEId = fespace->GetFE3D(i, Me);
//    U_FEId = U_fespace->GetFE3D(i, Me);
//    
//    LocalUsedElements[0]=FEId;
//    LocalUsedElements[1]=U_FEId;
//   
//    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
// 			   Coll, Me, SecondDer, 
// 			   N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
//    
//    BaseFunct = BaseFuncts[FEId];
//    N_V_LocalDof = N_BaseFunct[FEId];
//    
//    U_BaseFunct = BaseFuncts[U_FEId];
//    N_U_LocalDof = N_BaseFunct[FEId];
//    
//    vorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
//    uorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
//    
//    memset(LocMatrixA, 0, N_V_LocalDof*N_V_LocalDof*SizeOfDouble);
//    memset(LocRhs, 0, N_V_LocalDof*SizeOfDouble);
//    
//    DOF = GlobalNumbers + BeginIndex[i];
//    UDOF = UGlobalNumbers + UBeginIndex[i];
//    
//    for(j=0;j<N_Points;j++)
//    {
//      Orig = uorig[j];
//      
//      U_value = 0;
//      
//      for(l=0; l<N_U_LocalDof; l++)
//      {
//        U_value += U_Sol[UDOF[l]] * Orig[l];
//      }
//      
//      Mult = weights[j]*AbsDetjk[j];
//      alphaU = Mult * alpha * U_value;
//      alphaUbeta = Mult *(alpha* U_value + beta);
//      
//      //Assembling local matrix and rhs
//      
//      for(k=0;k<N_V_LocalDof;k++)
//      {
//        test000 = Orig[k];
//        LocRhs[k] += alphaU*test000;
//        
//        ufact = alphaUbeta*test000;
//        
//        for(l=0;l<N_V_LocalDof;l++)
//        {
// 	 ansatz00 = Orig[l];
// 	 LocMatrixA[k*N_V_LocalDof + l] += ufact*ansatz00;
//        }
//      }
//    }
//    
//    //local to global matrix
//    
//    for(l=0;l<N_V_LocalDof;l++)
//    {
//      TestDOF = DOF[l];
//      RHS[TestDOF] += LocRhs[l];
//      
//      begin = RowPtr[TestDOF];
//      end = RowPtr[TestDOF+1];
//      
//      for(n=begin;n<end;n++)
//      {
//        for(m=0;m<N_V_LocalDof;m++)
//        {
// 	 if(KCol[n] == DOF[m])
// 	 {
// 	   ValuesA[n] += LocMatrixA[l*N_V_LocalDof+m];
// 	   
// 	   break;
// 	 }
//        }
//      }
//    }
//  }
// }

int main(int argc, char* argv[])
{
// ======================================================================
// variable declaration
// ======================================================================
  int i, j, l, m, N_SubSteps, ORDER, LEVELS, mg_level, N_Cells, N_DOF, N_V_DOF, N_W_DOF, img=1;
  int N_Active, N_V_Active, N_W_Active, Max_It_scalar, Min_It_scalar;
  int mg_type;
  
  double *sol, *oldsol, *rhs, *oldrhs, t1, t2, errors[5], **Sol_array, **Rhs_array;
  double tau, end_time, *defect, olderror, olderror1, hmin, hmax;
  double *V_sol, *V_oldsol, *V_rhs, *V_oldrhs, *V_defect, **VSol_array, **VRhs_array;
  double *W_sol, *W_oldsol, *W_rhs, *W_oldrhs, *W_defect, **WSol_array, **WRhs_array;
  double residual_scalar, oldresidual_scalar, limit_scalar;
  double Parameters[10];
  
  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TCollection *coll;
  TFESpace3D **Scalar_FeSpaces, **Scalar_V_FeSpaces, **Scalar_W_FeSpaces, *fesp[1];
  TFEFunction3D *Scalar_FeFunction, *Scalar_V_FeFunction, *Scalar_W_FeFunction, 
                **Scalar_FeFunctions, **Scalar_V_FeFunctions, **Scalar_W_FeFunctions;
  TOutput3D *Output;
  TSystemMatTimeScalar3D *SystemMatrix, *V_SystemMatrix, *W_SystemMatrix;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = {D000, D100, D010, D001};
  TSquareMatrix3D *V_A_Matrix;
  
//   TSquareMatrix3D *sqmatrixA, *V_sqmatrixA, *W_sqmatrixA, **sqmatricesA, *SQMATRICES[3];;
//   TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
  char Name[] = "name";
  char Description[] = "description";
  char CString[] = "C";
  char VString[] = "V";
  char WString[] = "W";
  
  double Linfty=0;

  std::ostringstream os;
  os << " ";   
    
  mkdir(vtkdir, 0777);
  
// ======================================================================
// set the database values and generate mesh
// ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  
  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh  
  GEO = TDatabase::ParamDB->GEOFILE;
  Domain->Init(NULL, GEO);

  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
  
  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }
  
//=========================================================================
// set data for multigrid
//=========================================================================  
  LEVELS = TDatabase::ParamDB->LEVELS;

  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
     mg_type=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SCALAR = mg_type;
    }
  
  if(mg_type)
   {
    mg_level =  LEVELS + 1;
    ORDER = -1;
   }
  else
   {
    mg_level = LEVELS;
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
   }
   
  if(TDatabase::ParamDB->SOLVER_TYPE==GMG)
   {
    OutPut("=======================================================" << endl);
    OutPut("======           GEOMETRY  LEVEL ");
    OutPut(LEVELS << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level << "              ======" << endl);
    OutPut("=======================================================" << endl);   
   }
    
  Scalar_FeSpaces = new TFESpace3D*[mg_level];  
  Scalar_V_FeSpaces = new TFESpace3D*[mg_level]; 
  Scalar_W_FeSpaces = new TFESpace3D*[mg_level]; 
  
  Scalar_FeFunctions = new TFEFunction3D*[mg_level]; 
  Scalar_V_FeFunctions = new TFEFunction3D*[mg_level];
  Scalar_W_FeFunctions = new TFEFunction3D*[mg_level];
  
  Sol_array = new double*[mg_level];
  Rhs_array = new double*[mg_level];
  VSol_array = new double*[mg_level];
  VRhs_array = new double*[mg_level];
  WSol_array = new double*[mg_level];
  WRhs_array = new double*[mg_level];
  
//=========================================================================
// construct all finite element spaces
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
  for(i=0;i<LEVELS;i++)
   {   
    if(i)
     { Domain->RegRefineAll(); }

     coll=Domain->GetCollection(It_Finest, 0);
  
     // fespaces for scalar equation 
     Scalar_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);     
     Scalar_V_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, V_BoundCondition, ORDER); 
     Scalar_W_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, W_BoundCondition, ORDER); 
     
     //multilevel multigrid disc
     if(i==LEVELS-1 && mg_type==1) 
      {
       ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
       Scalar_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);
       Scalar_V_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, V_BoundCondition, ORDER);
       Scalar_W_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, W_BoundCondition, ORDER);
      } //  if(i==LEVELS-1 && i!=mg_level-1) 
     
//======================================================================
// construct all finite element functions
//======================================================================
    
    //===================================================================
    //u solution
    //===================================================================
    N_DOF = Scalar_FeSpaces[i]->GetN_DegreesOfFreedom();
    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    Sol_array[i] = sol;
    Rhs_array[i] = rhs;   
    
    Scalar_FeFunction  = new TFEFunction3D(Scalar_FeSpaces[i], CString, CString, sol, N_DOF);  
    Scalar_FeFunctions[i] = Scalar_FeFunction;
     
    if(i==LEVELS-1 && mg_type==1) 
     {  
      N_DOF = Scalar_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
      sol = new double[N_DOF];
      rhs = new double[N_DOF];
      Sol_array[mg_level-1] = sol;
      Rhs_array[mg_level-1] = rhs;

      Scalar_FeFunction = new TFEFunction3D(Scalar_FeSpaces[mg_level-1], CString, CString, sol, N_DOF);   
      Scalar_FeFunctions[mg_level-1] = Scalar_FeFunction;
     }//   if(i==LEVELS-1 && mg_type==1) 
     //========================================================================
     
     
     //========================================================================
   //v solution
   //========================================================================
    N_V_DOF = Scalar_V_FeSpaces[i]->GetN_DegreesOfFreedom();
    V_sol = new double[N_V_DOF];
    V_rhs = new double[N_V_DOF];
    VSol_array[i] = V_sol;
    VRhs_array[i] = V_rhs;   
    
    Scalar_V_FeFunction  = new TFEFunction3D(Scalar_V_FeSpaces[i], VString, VString, V_sol, N_V_DOF);  
    Scalar_V_FeFunctions[i] = Scalar_V_FeFunction;
     
    if(i==LEVELS-1 && mg_type==1) 
     {
      N_V_DOF = Scalar_V_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
      V_sol = new double[N_V_DOF];
      V_rhs = new double[N_V_DOF];
      VSol_array[mg_level-1] = V_sol;
      VRhs_array[mg_level-1] = V_rhs;

      Scalar_V_FeFunction = new TFEFunction3D(Scalar_V_FeSpaces[mg_level-1], VString, VString, V_sol, N_V_DOF);   
      Scalar_V_FeFunctions[mg_level-1] = Scalar_V_FeFunction;
     }//   if(i==LEVELS-1 && mg_type==1) 
  //  =============================================================================
    
    
   //  ========================================================================
  // w solution
   //========================================================================
    N_W_DOF = Scalar_W_FeSpaces[i]->GetN_DegreesOfFreedom();
    W_sol = new double[N_W_DOF];
    W_rhs = new double[N_W_DOF];
    WSol_array[i] = W_sol;
    WRhs_array[i] = W_rhs;   
    
    Scalar_W_FeFunction  = new TFEFunction3D(Scalar_W_FeSpaces[i], WString, WString, W_sol, N_W_DOF);  
    Scalar_W_FeFunctions[i] = Scalar_W_FeFunction;
     
    if(i==LEVELS-1 && mg_type==1) 
     {
      N_W_DOF = Scalar_W_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
      W_sol = new double[N_W_DOF];
      W_rhs = new double[N_W_DOF];
      WSol_array[mg_level-1] = W_sol;
      WRhs_array[mg_level-1] = W_rhs;

      Scalar_W_FeFunction = new TFEFunction3D(Scalar_W_FeSpaces[mg_level-1], WString, WString, W_sol, N_W_DOF);   
      Scalar_W_FeFunctions[mg_level-1] = Scalar_W_FeFunction;
     }//    // if(i==LEVELS-1 && mg_type==1) 
     //========================================================================
     
   }// for(i=0;i<LEVELS;i++)
   
   oldrhs = new double[N_DOF];
   oldsol = new double[N_DOF];
   
   V_oldrhs = new double[N_V_DOF];
   V_oldsol = new double[N_V_DOF];
   
   W_oldrhs = new double[N_W_DOF];
   W_oldsol = new double[N_W_DOF];
   
   N_Cells = coll->GetN_Cells();
   
   
   
   OutPut("N_Cells   : " << N_Cells <<endl);
   OutPut("Dof all   : " << N_DOF  << endl);  
   OutPut("Dof all   : " << N_V_DOF  << endl);
   OutPut("Dof all   : " << N_W_DOF  << endl);
   OutPut(endl); 
   
//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    /** interpolate the initial value */
    Scalar_FeFunction->Interpolate(InitialCondition);   
    Scalar_V_FeFunction->Interpolate(InitialCondition); 
    Scalar_W_FeFunction->Interpolate(InitialCondition); 
   
    
    // u system
    /** Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
     *  Solver: AMG_SOLVE (or) GMG  (or) DIRECT */
    SystemMatrix = new TSystemMatTimeScalar3D(mg_level, Scalar_FeSpaces, TDatabase::ParamDB->DISCTYPE, 
					      TDatabase::ParamDB->SOLVER_TYPE);
    
    /** initilize the system matrix with the functions defined in Example file */
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);
       
    /** assemble the system matrix with given aux, sol and rhs 
     *  aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
     *  otherwise, just pass it with NULL  */
    SystemMatrix->AssembleMRhs(NULL, Sol_array, Rhs_array); 
    
   /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
    memcpy(oldrhs, rhs, N_DOF*SizeOfDouble);   
    
    
    // v system
    V_SystemMatrix = new TSystemMatTimeScalar3D(mg_level, Scalar_V_FeSpaces, TDatabase::ParamDB->DISCTYPE, 
						TDatabase::ParamDB->SOLVER_TYPE);
    V_SystemMatrix->Init(V_BilinearCoeffs, V_BoundCondition, V_BoundValue);
    V_SystemMatrix->AssembleMRhs(NULL, VSol_array, VRhs_array);
    memcpy(V_oldrhs, V_rhs, N_V_DOF*SizeOfDouble);
    
    //w system
    W_SystemMatrix = new TSystemMatTimeScalar3D(mg_level, Scalar_W_FeSpaces, TDatabase::ParamDB->DISCTYPE, 
 						TDatabase::ParamDB->SOLVER_TYPE);
    W_SystemMatrix->Init(W_BilinearCoeffs, W_BoundCondition, W_BoundValue);
    W_SystemMatrix->AssembleMRhs(NULL, WSol_array, WRhs_array);
    memcpy(W_oldrhs, W_rhs, N_W_DOF*SizeOfDouble);
//======================================================================
// produce outout at t=0
//======================================================================
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput3D(2, 2, 1, 1, Domain);

    Output->AddFEFunction(Scalar_FeFunction);
    Output->AddFEFunction(Scalar_V_FeFunction);
    Output->AddFEFunction(Scalar_W_FeFunction);

//     Scalar_FeFunction->Interpolate(Exact);   
    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }   
// exit(0);
     
    /** measure errors to known solution */
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] = Scalar_FeSpaces[mg_level-1];       
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
     for(j=0;j<5;j++)
       errors[j] = 0;
     
      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);
     
      olderror = errors[0];
      olderror1 = errors[1]; 
     
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);     
//      Linfty=errors[0];
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  
       
       
       cout << "time " << TDatabase::TimeDB->CURRENTTIME << endl;
       
  coll->GetHminHmax(&hmin,&hmax);
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
  
  TDatabase::TimeDB->TIMESTEPLENGTH = hmax;
// exit(0);
//======================================================================
// time disc loop
//======================================================================    
   // parameters for time stepping scheme
   m = 0;
   N_SubSteps = GetN_SubSteps();
   end_time = TDatabase::TimeDB->ENDTIME; 
   
   Max_It_scalar = (int)TDatabase::ParamDB->REACTOR_P10;
   Min_It_scalar = TDatabase::ParamDB->REACTOR_P11;
   limit_scalar = TDatabase::ParamDB->REACTOR_P12;

   UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
   UpdateRhs = TRUE; //check BilinearCoeffs in example file
   ConvectionFirstTime=TRUE;
   
//    cout << "init " << Ddot(N_DOF, sol, sol)<< endl;
   /** time loop starts */
   while(TDatabase::TimeDB->CURRENTTIME< end_time)
    {
     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

     for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
      {
       SetTimeDiscParameters();

      if(m==1)
       {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
       }

      /** copy the sol to old sol */  
      memcpy(oldsol, sol, N_DOF*SizeOfDouble);  
      memcpy(V_oldsol, V_sol, N_V_DOF*SizeOfDouble);
      memcpy(W_oldsol, W_sol, N_W_DOF*SizeOfDouble);
       
      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
   
      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);  
      
      //Gauss-Seidal type linearization
      //exit(0);
//       sqmatrixA=sqmatricesA[mg_level-1];
//       SQMATRICES[0] = sqmatrixA;
      
      
      for(j=0;j<Max_It_scalar;j++)
      {
//===========================================================================================
// u system
//===========================================================================================	
	
      /** unless the stiffness matrix or rhs change in time, it is enough to assemble only once at the begning   */
      if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
       {  
        if(UpdateRhs)
         { SystemMatrix->AssembleARhs(NULL, Sol_array, Rhs_array); }
        else
         { SystemMatrix->AssembleARhs(NULL, Sol_array, Rhs_array); }
        
        AssembleCancerDensity(mg_level, SystemMatrix->GetAMatrix(), Scalar_FeFunctions, Scalar_V_FeFunctions);
        
        /**  M:= M + (tau*TDatabase::TimeDB->THETA1)*A
         *   rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs + [ M - (tau*THETA2)A]*oldsol **/
        SystemMatrix->AssembleSystMat(oldrhs, oldsol, rhs, sol);

        /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
        memcpy(oldrhs, rhs, N_DOF*SizeOfDouble); 
        
        ConvectionFirstTime = FALSE;
       }
       
       
      residual_scalar = SystemMatrix->GetResidual(sol);
        OutPut("Scalar nonlinear step " << setw(3) << j);
        OutPut(setw(14) << residual_scalar); // sqrt of residual_scalar is alread done in ScalarDefect

        if (j>0)
         { 
          if(fabs(oldresidual_scalar)>0)
          OutPut(setw(14) <<  residual_scalar/oldresidual_scalar ); 
          OutPut(endl);
         }
        else
         { OutPut(endl); 
	   
	}

        oldresidual_scalar = residual_scalar;
   
        if( ((residual_scalar<=limit_scalar)||(j==Max_It_scalar-1))  && (j>=Min_It_scalar) )
         { 
          if(UpdateStiffnessMat || UpdateRhs)
           {         
            SystemMatrix->RestoreMassMat();
           }   
          break;
         } 

      // solve the system matrix 
      SystemMatrix->Solve(sol);

      /** restore the mass matrix for the next time step unless the stiffness matrix 
       * or rhs change in time, it is not necessary to assemble the system matrix in every time step */
      if(UpdateStiffnessMat || UpdateRhs)
       {         
        SystemMatrix->RestoreMassMat();
       }
       
//=======================================================================       
//W SYSTEM
//======================================================================= 
//         if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
//          {  
//         if(UpdateRhs)
//          { W_SystemMatrix->AssembleARhs(NULL, WSol_array, WRhs_array); }
//         else
//          { W_SystemMatrix->AssembleARhs(NULL, WSol_array, WRhs_array); }
//         
//         AssembleWArhs(SystemMatrix->GetA_FineMatrix(), Scalar_FeFunction, W_rhs);
//         
//         /**  M:= M + (tau*TDatabase::TimeDB->THETA1)*A
//          *   rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs + [ M - (tau*THETA2)A]*oldsol **/
//         W_SystemMatrix->AssembleSystMat(W_oldrhs, W_oldsol, W_rhs, W_sol);
// 
//         /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
//         memcpy(W_oldrhs, W_rhs, N_W_DOF*SizeOfDouble); 
//         
// 	//?????????
// 	ConvectionFirstTime = FALSE;
//        }
//      
//         W_SystemMatrix->Solve(W_sol);
// 
//         if(UpdateStiffnessMat || UpdateRhs)
//          {         
//           W_SystemMatrix->RestoreMassMat();
//          }  
//=======================================================================       
// V SYSTEM
//======================================================================= 
//         no need a=to assemble A matrix
       /*  
         V_SystemMatrix->AssembleARhs(NULL, V_sol, V_rhs);
         
        V_A_Matrix = V_SystemMatrix->GetAMatrix();
        V_A_Matrix->Reset();
       // AssembleVMat(V_A_Matrix, Scalar_W_FeFunction, Scalar_FeFunction, Scalar_V_FeFunction);
        V_SystemMatrix->AssembleSystMat(V_oldrhs, V_oldsol, V_rhs, V_sol);
	
// 	exit(0);

        V_SystemMatrix->Solve(V_sol);

        if(UpdateStiffnessMat || UpdateRhs)
         {         
          V_SystemMatrix->RestoreMassMat();
         } */ 
//=======================================================================      
      
      } // for(j=o;j<Max_It_scalar;j++)
     } // for(l=0;l<N_SubSteps;l++) 
     
//======================================================================
// produce outout
//======================================================================
    if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      }
      
//======================================================================
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {      
      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);

      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);

      
      if(m>1)
       {      
        errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

        errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
       }
      
      olderror = errors[0];
      olderror1 = errors[1]; 
      
      if(Linfty<errors[0])
       Linfty=errors[0];

      OutPut( "Linfty " << Linfty << endl);      
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS)
//      exit(0); 
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

//======================================================================
// produce final outout
//======================================================================
  
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      }
      
  CloseFiles();
  
  return 0;
} // end main





