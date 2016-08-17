// =======================================================================
// @(#)Assemble2D.C        1.16 04/13/00
//
// Purpose:     bilinear form (discretized and stabilized assemble)
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================

#include <DefineParams.h>

#include <Assembler3.h>
//#include <Assemble2D.h>
#include <Enumerations.h>
#include <Matrix2D.h>
#include <AuxParam2D.h>
#include <DiscreteForm2D.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>
#include <FEDatabase2D.h>
#include <NodalFunctional2D.h>
#include <SquareMatrix2D.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <Convolution.h>

#include <string.h>
#include <stdlib.h>
#include <vector>



Assembler3::Assembler3(LocalAssembling2D_type type, TFEFunction2D **fefunctions2d,
                       CoeffFct2D *coeffs):la(type,fefunctions2d,coeffs)
{}

// HIER /////////////////////////////////////////////////////////////////////////////////////////////////////

void Assembler3::Assemble2D(int n_fespaces, const TFESpace2D** fespaces, int n_sqmatrices,
                            TSquareMatrix2D** sqmatrices, int n_matrices,
                            TMatrix2D** matrices, int n_rhs, double** rhs,
                            const TFESpace2D** ferhs,
                            //BoundCondFunct2D** BoundaryConditions,
                            //BoundValueFunct2D** BoundaryValues,
                            //LocalAssembling2D& la,
                            const Example2D& example,
                            int AssemblePhaseID)
{
#ifdef __3D__

    ErrThrow("Assembler not working in 3D");
#endif
    
    int **GlobalNumbers, **BeginIndex;
    int **TestGlobalNumbers, **TestBeginIndex;
    int **AnsatzGlobalNumbers, **AnsatzBeginIndex;
    int **RhsBeginIndex, **RhsGlobalNumbers;
    double *righthand;
    double **Matrices;
    double ***LocMatrices, **LocRhs;
    double **HangingEntries, **HangingRhs;
    double *Param[MaxN_QuadPoints_2D];
    for(size_t i=0;i<MaxN_QuadPoints_2D;++i) //initialize Param
    {
        // NOTE: The number 10 is magic here.
        // In the current setup it may well happen, that a Coefficient function
        // expects to evaluate parameters, but the local assembling object
        // will not make use of them (N_Parameters = 0) - nevertheless the array
        // must be initialized to something.
        Param[i]=new double[10]{0.0};
    }
    
    
  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  if(n_sqmatrices+n_matrices)
  {
    HangingEntries = new double *[n_sqmatrices+n_matrices];

    for(int i=0 ; i<n_sqmatrices ; i++)
    {
      int j = sqmatrices[i]->GetHangingN_Entries();
      HangingEntries[i] = new double [j];
      memset(HangingEntries[i], 0, SizeOfDouble*j);
    }

    for(int i=0;i<n_matrices;i++)
    {
      int j = matrices[i]->GetHangingN_Entries();
      HangingEntries[i+n_sqmatrices] = new double [j];
      memset(HangingEntries[i+n_sqmatrices], 0, SizeOfDouble*j);
    }
  }

  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    for(int i=0;i<n_sqmatrices;i++)
    {
      const TFESpace2D *fespace = sqmatrices[i]->GetFESpace2D();
      GlobalNumbers[i] = fespace->GetGlobalNumbers();
      BeginIndex[i] = fespace->GetBeginIndex();
    }                                             // endfor
  }                                               // endif n_sqmatrices

  if(n_matrices)
  {
    TestGlobalNumbers = new int* [n_matrices];
    AnsatzGlobalNumbers = new int* [n_matrices];
    TestBeginIndex = new int* [n_matrices];
    AnsatzBeginIndex = new int* [n_matrices];
    for(int i=0;i<n_matrices;i++)
    {
      const TFESpace2D *fespace = (TFESpace2D *) matrices[i]->GetTestSpace2D();
      TestGlobalNumbers[i] = fespace->GetGlobalNumbers();
      TestBeginIndex[i] = fespace->GetBeginIndex();

      fespace = (TFESpace2D *) matrices[i]->GetAnsatzSpace2D();
      AnsatzGlobalNumbers[i] = fespace->GetGlobalNumbers();
      AnsatzBeginIndex[i] = fespace->GetBeginIndex();
    }                                             // endfor
  }                                               // endif n_matrices

  if(n_rhs)
  {
    HangingRhs = new double* [n_rhs];
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    for(int i=0;i<n_rhs;i++)
    {
      const TFESpace2D *fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();

      int j = fespace->GetN_Hanging();
      HangingRhs[i] = new double [j];
      memset(HangingRhs[i], 0, SizeOfDouble*j);
    }                                             // endfor

    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(int i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

  //N_Parameters = Parameters->GetN_Parameters();
  int N_Parameters = la.GetN_Parameters();
  
  if(N_Parameters)
  {
    double *aux = new double [MaxN_QuadPoints_2D*N_Parameters];
    for(int j=0;j<MaxN_QuadPoints_2D;j++)
    {
      delete[] Param[j]; //clear away the default initialized array
      Param[j] = aux + j*N_Parameters;
    }
  };

  // 40 <= number of terms in bilinear form
  // DUE NOTE CHANGE BELOW 20 SINCE THE ENTRY 19 IS USED IN GetLocalForms
  double *aux = new double [MaxN_QuadPoints_2D*40];
  double *AuxArray[MaxN_QuadPoints_2D];
  for(int j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*40;

  int N_AllMatrices = n_sqmatrices+n_matrices;
  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
    for(int j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D;

    LocMatrices = new double** [N_AllMatrices];
    for(int i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
  }                                               // endif N_AllMatrices


    
    
  loop_over_cells(fespaces, AssemblePhaseID, n_fespaces,
                  //  la,
                  n_sqmatrices, LocMatrices,
                    sqmatrices,n_matrices,matrices,
                    TestGlobalNumbers, TestBeginIndex,
                    AnsatzGlobalNumbers, AnsatzBeginIndex,HangingEntries,HangingRhs,
                    righthand, RhsBeginIndex,RhsGlobalNumbers,
                    //BoundaryConditions,
                    //BoundaryValues,
                    N_AllMatrices, Matrices, Param,  n_rhs,  rhs,
                    AuxArray,LocRhs,ferhs,
                    N_BaseFunct,  BaseFuncts,GlobalNumbers, BeginIndex,example);

    
    
    
    
  // ####################################################################
  // modify matrix according to coupling
  // ####################################################################
  for(int j=0;j<n_sqmatrices;j++)
  {
    const TFESpace2D *fespace = sqmatrices[j]->GetFESpace2D();
    int N_ = fespace->GetN_Hanging();
    THangingNode **HangingNodes = fespace->GetHangingNodes();

    double *Entries = sqmatrices[j]->GetEntries();
    const int *RowPtr = sqmatrices[j]->GetRowPtr();
    const int *ColInd = sqmatrices[j]->GetKCol();

    double *CurrentHangingEntries = HangingEntries[j];
    const int *HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
    const int *HangingColInd = sqmatrices[j]->GetHangingKCol();

    int ActiveBound = fespace->GetActiveBound();

    for(int i=0;i<N_;i++)
    {
      THangingNode *hn = HangingNodes[i];
      HNDesc HNDescr = hn->GetType();
      THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      int k = HNDescr_Obj->GetN_Nodes();
      double *Coupling = HNDescr_Obj->GetCoeff();
      int *DOF = hn->GetDOF();

      int end = HangingRowPtr[i+1];
      for(int n=HangingRowPtr[i];n<end;n++)
      {
        double v = CurrentHangingEntries[n];
        int m = HangingColInd[n];
        for(int l=0;l<k;l++)
        {
          int l1 = DOF[l];
          if(l1<ActiveBound)
          {
            int last=RowPtr[l1+1];
            for(int l2=RowPtr[l1];l2<last;l2++)
            {
              if(ColInd[l2] == m)
              {
                Entries[l2] += Coupling[l] * v;
              }
            }                                     // endfor l2
          }                                       // endif
        }                                         // endfor l
      }                                           // endfor n
    }                                             // endfor i
  }                                               // endfor j

  for(int j=0;j<n_matrices;j++)
  {
    // hanging nodes in test space
    const TFESpace2D *fespace = (TFESpace2D *) (matrices[j]->GetTestSpace2D());
    int N_ = fespace->GetN_Hanging();
    THangingNode **HangingNodes = fespace->GetHangingNodes();

    double *Entries = matrices[j]->GetEntries();
    const int *RowPtr = matrices[j]->GetRowPtr();
    const int *ColInd = matrices[j]->GetKCol();

    double *CurrentHangingEntries = HangingEntries[j+n_sqmatrices];
    const int *HangingRowPtr = matrices[j]->GetHangingRowPtr();
    const int *HangingColInd = matrices[j]->GetHangingKCol();

    int ActiveBound = fespace->GetActiveBound();

    for(int i=0;i<N_;i++)
    {
      THangingNode *hn = HangingNodes[i];
      HNDesc HNDescr = hn->GetType();
      THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      int k = HNDescr_Obj->GetN_Nodes();
      double *Coupling = HNDescr_Obj->GetCoeff();
      int *DOF = hn->GetDOF();

      int end = HangingRowPtr[i+1];
      for(int n=HangingRowPtr[i];n<end;n++)
      {
        double v = CurrentHangingEntries[n];
        int m = HangingColInd[n];
        for(int l=0;l<k;l++)
        {
          int l1 = DOF[l];
          if(l1<ActiveBound)
          {
            int last=RowPtr[l1+1];
            for(int l2=RowPtr[l1];l2<last;l2++)
            {
              if(ColInd[l2] == m)
              {
                Entries[l2] += Coupling[l] * v;
              }
            }                                     // endfor l2
          }                                       // endif
        }                                         // endfor l
      }                                           // endfor n
    }                                             // endfor i

    // hanging nodes in ansatz space
    int N_Rows =  matrices[j]->GetN_Rows();
    fespace = (TFESpace2D *) (matrices[j]->GetAnsatzSpace2D());

    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();
     int AnsatzActiveBound = fespace->GetActiveBound();
     int AnsatzHangingBound = fespace->GetHangingBound();
    for(int i=0;i<N_Rows;i++)
    {
      int end = RowPtr[i+1];
      for(int k=RowPtr[i];k<end;k++)
      {
        int l = ColInd[k];
        if(l>=AnsatzActiveBound && l<AnsatzHangingBound)
        {
          // l is hanging node in ansatz space
          THangingNode *hn = HangingNodes[l-AnsatzActiveBound];
          HNDesc HNDescr = hn->GetType();
          THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
          int m = HNDescr_Obj->GetN_Nodes();
          double *Coupling = HNDescr_Obj->GetCoeff();
          int *DOF = hn->GetDOF();
          double v = Entries[k];
          for(int n=0;n<m;n++)
          {
            for(int l1=RowPtr[i];l1<end;l1++)
            {
              if(ColInd[l1] == DOF[n])
                Entries[l1] += v*Coupling[n];
            }
          }
          Entries[k] = 0;
        }                                         // endif l
      }                                           // endfor k
    }                                             // endfor i
  }                                               // endfor j

  for(int j=0;j<n_rhs;j++)
  {
    const TFESpace2D *fespace = ferhs[j];
    int N_Hanging = fespace->GetN_Hanging();
    THangingNode **HangingNodes = fespace->GetHangingNodes();

    double *RHS = rhs[j];
    double *CurrentHangingRhs = HangingRhs[j];

    int ActiveBound = fespace->GetActiveBound();

    for(int i=0;i<N_Hanging;i++)
    {
      THangingNode *hn = HangingNodes[i];
      HNDesc HNDescr = hn->GetType();
      THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      int N_ = HNDescr_Obj->GetN_Nodes();
      double *Coupling = HNDescr_Obj->GetCoeff();
      int *DOF = hn->GetDOF();

      for(int k=0;k<N_;k++)
      {
        int l = DOF[k];
        if(l<ActiveBound)
        {
          RHS[l] += Coupling[k] * CurrentHangingRhs[i];
        }
      }                                           // endfor k
    }                                             // endfor i
  }                                               // endfor j

  // ####################################################################
  // write coupling into matrix
  // ####################################################################
  for(int j=0;j<n_sqmatrices;j++)
  {
    const TFESpace2D *fespace = sqmatrices[j]->GetFESpace2D();
    int N_ = fespace->GetN_Hanging();
    THangingNode **HangingNodes = fespace->GetHangingNodes();

    double *Entries = sqmatrices[j]->GetEntries();
    const int *RowPtr = sqmatrices[j]->GetRowPtr();
    const int *ColInd = sqmatrices[j]->GetKCol();

    int ActiveBound = fespace->GetActiveBound();

    int n = RowPtr[ActiveBound];

    for(int i=0;i<N_;i++)
    {
      THangingNode *hn = HangingNodes[i];
      HNDesc HNDescr = hn->GetType();
      THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      int k = HNDescr_Obj->GetN_Nodes();
      double *Coupling = HNDescr_Obj->GetCoeff();
      int *DOF = hn->GetDOF();

      Entries[n] = 1.0;
      n++;

      for(int l=0;l<k;l++)
      {
        Entries[n] = - Coupling[l];
        n++;
      }                                           // endfor l
    }                                             // endfor i
  }

  if(n_sqmatrices)
  {
    delete [] GlobalNumbers;
    delete [] BeginIndex;
  }

  if(n_matrices)
  {
    delete [] AnsatzGlobalNumbers;
    delete [] AnsatzBeginIndex;
    delete [] TestGlobalNumbers;
    delete [] TestBeginIndex;
  }

  if(n_sqmatrices+n_matrices)
  {
    for(int i=0;i<n_sqmatrices+n_matrices;i++)
    delete [] HangingEntries[i];
    delete [] HangingEntries;
  }

  if(n_rhs)
  {
    for(int i=0;i<n_rhs;i++)
    delete [] HangingRhs[i];
    delete [] HangingRhs;
    delete [] righthand;
    delete [] LocRhs;
    delete [] RhsBeginIndex;
    delete [] RhsGlobalNumbers;
  }

  if(N_Parameters)
  {//Param was used and has to be deleted
    delete [] Param[0];
  }
  else
  { //we have to delete the default initialized thing
    for(int j=0;j<MaxN_QuadPoints_2D;j++)
    {
      delete[] Param[j];
    }
  }

  if(N_AllMatrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }

  delete [] AuxArray[0];
}                                                 // end of Assemble




    // ########################################################################
    // loop over all cells
    // ########################################################################


void Assembler3::loop_over_cells(const TFESpace2D** fespaces,
                                 int AssemblePhaseID,
                                 int n_fespaces,
                                 //LocalAssembling2D& la,
                                 int n_sqmatrices,
                                 double ***LocMatrices,
                                 TSquareMatrix2D **sqmatrices,
                                 int n_matrices,TMatrix2D **matrices,
                                 int **TestGlobalNumbers,
                                 int **TestBeginIndex,
                                 int **AnsatzGlobalNumbers,
                                 int **AnsatzBeginIndex,
                                 double **HangingEntries,
                                 double **HangingRhs,
                                 double *righthand,
                                 int **RhsBeginIndex,
                                 int **RhsGlobalNumbers,
                                 //BoundCondFunct2D** BoundaryConditions,
                                // BoundValueFunct2D** BoundaryValues,
                                 int N_AllMatrices,
                                 double **Matrices,
                                 double **Param,
                                 int n_rhs,
                                 double** rhs,
                                 double **AuxArray,
                                 double **LocRhs,
                                 const TFESpace2D** ferhs,
                                 int *N_BaseFunct,
                                 BaseFunct2D *BaseFuncts,
                                 int **GlobalNumbers,
                                 int **BeginIndex,
                                 const Example2D& example)
{
    TCollection *Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
    int N_Cells = Coll->GetN_Cells();
    //!!  for(int i=0;i<N_Cells;i++) // set the cell indices
    //  {
    //    TBaseCell *cell = Coll->GetCell(i);
    //    cell->SetCellIndex(i);
    //  }
    
    for(int i=0;i<N_Cells;i++)
    {
        TBaseCell *cell = Coll->GetCell(i);
        //!!    TDatabase::ParamDB->INTERNAL_CELL = i;
        
        // only for multiphase flows
        if(AssemblePhaseID >= 0)
        {
            if(AssemblePhaseID != cell->GetPhase_ID() )
                continue;
        }
        
        // ####################################################################
        // find local used elements on this cell
        // ####################################################################
        FE2D LocalUsedElements[N_FEs2D];
        int LocN_BF[N_BaseFuncts2D];
        BaseFunct2D LocBF[N_BaseFuncts2D];
        
        for(int j=0;j<n_fespaces;j++)
        {
            FE2D CurrentElement = fespaces[j]->GetFE2D(i, cell);
            LocalUsedElements[j] = CurrentElement;
            LocN_BF[j] = N_BaseFunct[CurrentElement];
            LocBF[j] = BaseFuncts[CurrentElement];
        }
        
        // ####################################################################
        // calculate values on original element
        // ####################################################################
        int N_Points;
        double *weights, *xi, *eta;
        double X[MaxN_QuadPoints_2D],Y[MaxN_QuadPoints_2D];
        double AbsDetjk[MaxN_QuadPoints_2D];
        
        //SecondDer = DiscreteForm->GetNeeds2ndDerivatives();
        bool *SecondDer = la.GetNeeds2ndDerivatives();
        
        TFEDatabase2D::GetOrig(n_fespaces, LocalUsedElements,
                               Coll, cell, SecondDer,
                               N_Points, xi, eta, weights, X, Y, AbsDetjk);
        
        //Parameters->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);
        la.GetParameters(N_Points, Coll, cell, i, X, Y, Param);
    
        
        if((TDatabase::ParamDB->DISCTYPE == SDFEM)
           || (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM)
           || (TDatabase::ParamDB->CELL_MEASURE == 4))
        {
            TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
            int N_Edges = cell->GetN_Edges();
            for (int ij=0; ij<N_Edges;ij++)
            {
                TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
                TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
            }
            if (N_Edges==3)
                TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
                TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
        }
        
        
        la.GetLocalForms(N_Points, weights, AbsDetjk, X, Y, LocN_BF, LocBF,
                         Param, AuxArray, cell, N_AllMatrices, n_rhs, LocMatrices,
                         LocRhs);
        
        int N_Joints = cell->GetN_Joints();
        // ####################################################################
        // add local/cellwise matrices to global matrices (ansatz == test)
        // ####################################################################
        for(int j=0;j<n_sqmatrices;j++)
        {
            // find space for this bilinear form
            const TFESpace2D *fespace = sqmatrices[j]->GetFESpace2D();
            FE2D CurrentElement = fespace->GetFE2D(i, cell);
            int N_ = N_BaseFunct[CurrentElement];
            
            double **Matrix = LocMatrices[j];
            double *Entries = sqmatrices[j]->GetEntries();
            const int *RowPtr = sqmatrices[j]->GetRowPtr();
            const int *ColInd = sqmatrices[j]->GetKCol();
            
            double *CurrentHangingEntries = HangingEntries[j];
            const int *HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
            const int *HangingColInd = sqmatrices[j]->GetHangingKCol();
            
            int ActiveBound = fespace->GetActiveBound();
            int DirichletBound = fespace->GetHangingBound();
            int *DOF = GlobalNumbers[j] + BeginIndex[j][i];
            
            // add local matrix to global
            for(int m=0;m<N_;m++)
            {
                int l=DOF[m];
                double *MatrixRow = Matrix[m];
                if(l < ActiveBound)
                {
                    for(int k=0;k<N_;k++)
                    {
                        // DOF[k] is the global index of the k-th local degree of freedom
                        // MatrixRow[k] is the assembled value corresponding to the m-th
                        // local test function and k-th local ansatz function. That means it
                        // corresponds to the l=DOF[m]-th global test function and the
                        // DOF[k]-th global ansatz function
                        sqmatrices[j]->add(l, DOF[k], MatrixRow[k]);
                    }
                }                                         // endif l
                else
                {
                    if(l<DirichletBound)
                    {
                        // hanging node
                        l -= ActiveBound;
                        int end = HangingRowPtr[l+1];
                        for(int n=HangingRowPtr[l];n<end;n++)
                        {
                            for(int k=0;k<N_;k++)
                            {
                                if(DOF[k] == HangingColInd[n])
                                {
                                    CurrentHangingEntries[n] += MatrixRow[k];
                                    break;
                                }                                 // endif
                            }                                   // endfor k
                        }                                     // endfor n
                    }
                    else
                    {// Dirichlet node
                        if(TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE) //this is needed e.g. when applying FEM-FCT
                        { // In that case, treat Dirichlet nodes like active dofs
                            // and add the calculated values to the square matrix
                            for(int k=0;k<N_;k++)
                            {
                                sqmatrices[j]->add(l, DOF[k], MatrixRow[k]);
                            }
                        } else {
                            //Dirichlet standard treatment
                            sqmatrices[j]->set(l,l,1.0);
                        }
                    }
                }
            }                                           // endfor m
        }                                             // endfor j
        
        // ################################################################################
        // add local/cellwise matrices to global matrices (A11,A12,B1,...) (ansatz != test)
        // ################################################################################
        for(int j=0;j<n_matrices;j++)
        {
            FE2D TestElement = ((TFESpace2D *) matrices[j]->GetTestSpace())
            ->GetFE2D(i, cell);
            FE2D AnsatzElement = ((TFESpace2D *) matrices[j]->GetAnsatzSpace())
            ->GetFE2D(i, cell);
            // cout << "non square matrix: " << j << endl;
            // cout << "TestElement: " << TestElement << endl;
            // cout << "AnsatzElement: " << AnsatzElement << endl;
            
            int N_Test = N_BaseFunct[TestElement];
            int N_Ansatz = N_BaseFunct[AnsatzElement];
            
            double **Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions2D;
            
            double *Entries = matrices[j]->GetEntries();
            const int *RowPtr = matrices[j]->GetRowPtr();
            const int *ColInd = matrices[j]->GetKCol();
            
            int *TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
            int *AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];
            
            const TFESpace2D *fespace = (TFESpace2D *)(matrices[j]->GetTestSpace2D());
            int ActiveBound = fespace->GetActiveBound();
            int DirichletBound = fespace->GetHangingBound();
            
            double *CurrentHangingEntries = HangingEntries[j+n_sqmatrices];
            const int *HangingRowPtr = matrices[j]->GetHangingRowPtr();
            const int *HangingColInd = matrices[j]->GetHangingKCol();
            
            // add local matrix to global
            for(int m=0;m<N_Test;m++)
            {
                int l=TestDOF[m];
                double *MatrixRow = Matrix[m];
                // cout << "DOF: " << l << endl;
                if(l < ActiveBound ||  l>=DirichletBound)
                {
                    int end=RowPtr[l+1];
                    for(int n=RowPtr[l];n<end;n++)
                    {
                        for(int k=0;k<N_Ansatz;k++)
                        {
                            if(AnsatzDOF[k] == ColInd[n])
                            {
                                // cout << m << "   " << k << endl << n << endl;
                                Entries[n] += MatrixRow[k];
                                break;
                            }                                   // endif
                        }                                     // endfor k
                    }                                       // endfor n
                }
                else
                {
                    // hanging node
                    l -= ActiveBound;
                    int end = HangingRowPtr[l+1];
                    for(int n=HangingRowPtr[l];n<end;n++)
                    {
                        // cout << l << " HangingColInd: " << HangingColInd[n] << endl;
                        for(int k=0;k<N_Ansatz;k++)
                        {
                            // cout << "AnsatzDOF: " << AnsatzDOF[k] << endl;
                            if(AnsatzDOF[k] == HangingColInd[n])
                            {
                                CurrentHangingEntries[n] += MatrixRow[k];
                                break;
                            }                                   // endif
                        }                                     // endfor k
                    }                                       // endfor n
                }
            }                                           // endfor m
        }                                             // endfor j  (n_matrices)
        
        // ####################################################################
        // add local right-hand sides to global right-hand side
        // ####################################################################
        for(int j=0;j<n_rhs;j++)
        {
            const TFESpace2D *fespace = ferhs[j];
            int ActiveBound = fespace->GetActiveBound();
            FE2D CurrentElement = fespace->GetFE2D(i, cell);
            BoundCond Cond0, Cond1;
            
            int N_ = N_BaseFunct[CurrentElement];
            
            double *local_rhs = righthand+j*MaxN_BaseFunctions2D;
            double *RHS = rhs[j];
            double *CurrentHangingRhs = HangingRhs[j];
            // find space for this linear form
            
            ActiveBound = fespace->GetActiveBound();
            int DirichletBound = fespace->GetHangingBound();
            
            // dof of the rhs nodes connected to this cell
            int *DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];
            
            // add local right-hand side to the global one
            for(int m=0;m<N_;m++)
            {
//                if (TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS)
//                {
//                    int l = 0;
//                    int N_Edges=cell->GetN_Edges();
//                    for(int jj=0;jj<N_Edges;jj++)                        // loop over all edges of cell
//                    {
//                        TVertex *ver0 = cell->GetVertex(jj);
//                        double t0 =  ver0->GetClipBoard();
//                        // vertex not on the boundary
//                        if (t0<-1e-8)
//                            continue;
//                        // component of boundary
//                        int comp = floor(t0+1e-8);
//                        // parameter
//                        t0 -= comp;
//                        // get boundary condition
//                        BoundaryConditions[0](comp, t0, Cond0);
//                        // Dirichlet
//                        if (Cond0== DIRICHLET)
//                        {
//                            l = -4711;
//                        }
//                    }
//                    if (l==-4711)
//                        break; // do nothing for this mesh cell
//                }
                
                int l=DOF[m];
                //cout << "DOF: " << l << endl;
                if(l<ActiveBound)
                {
                    // node l is inner or Neumann node
                    RHS[l] += local_rhs[m];
                    // cout << l << " " << RHS[l] << " " << local_rhs[m]<< " "<<endl;;
                }                                         // endif l
                else
                {
                    if(l<DirichletBound)
                    {
                        // hanging node
                        l -= ActiveBound;
                        CurrentHangingRhs[l] += local_rhs[m];
                    }
                }
            }                                           // endfor m
            
            BoundCondFunct2D *BoundaryCondition = fespace->GetBoundCondition();//BoundaryConditions[j];
            BoundValueFunct2D *BoundaryValue = example.get_bd()[j];//BoundaryValues[j];
            TFE2D *ele = TFEDatabase2D::GetFE2D(CurrentElement);
            double t0,t1;
            TBoundComp *BoundComp;
            int N_EdgePoints;
            double *EdgePoints;
            double eps=1e-4;
            
            //if ((ele >= D_P1_2D_Q_A)&&(ele<= D_P3_2D_Q_M))
            //  continue;
            
            TNodalFunctional2D *nf = ele->GetNodalFunctional2D();
            
            if(TDatabase::ParamDB->SUPERCONVERGENCE_ORDER)
            {
                /* Superconvergence boundary interpolation */
                if(nf->GetID() == NF_C_Q_Q2_2D)
                    nf = TFEDatabase2D::GetNodalFunctional2D(NF_S_Q_Q2_2D);
            }
            
            nf->GetPointsForEdge(N_EdgePoints, EdgePoints);
            
            TFEDesc2D *FEDesc_Obj = ele->GetFEDesc2D();
            int  N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
            // setting Dirichlet boundary condition
            N_Joints = cell->GetN_Edges();
            
            for(int m=0;m<N_Joints;m++)
            {
                TJoint *joint = cell->GetJoint(m);
                
                if(joint->GetType() == BoundaryEdge ||
                   joint->GetType() == IsoBoundEdge ||
                   joint->GetType() == InterfaceJoint)
                {
                    if(joint->GetType() == BoundaryEdge||
                       joint->GetType() == InterfaceJoint)
                    {
                        TBoundEdge *boundedge = (TBoundEdge *)joint;
                        BoundComp = boundedge->GetBoundComp();
                        boundedge->GetParameters(t0, t1);
                    }
                    else
                    {
                        TIsoBoundEdge *isoboundedge = (TIsoBoundEdge *)joint;
                        BoundComp = isoboundedge->GetBoundComp();
                        isoboundedge->GetParameters(t0, t1);
                    }
                    // get id of the boundary component
                    int comp=BoundComp->GetID();
                    // get type of the boundary condition at the beginning
                    // and at the end of the current edge
                    if (t0 < t1)
                    {
                        BoundaryCondition(comp, t0+eps, Cond0);
                        BoundaryCondition(comp, t1-eps, Cond1);
                    }
                    else
                    {
                        BoundaryCondition(comp, t0-eps, Cond0);
                        BoundaryCondition(comp, t1+eps, Cond1);
                    }
                    
                    // only one boundary condition per edge allowed
                    int l,k,l3;
                    QuadFormula1D LineQuadFormula;
                    TQuadFormula1D *qf1;
                    double s,t;
                    int *EdgeDOF;
                    double *LineWeights, *zeta;
                    double x0, x1, y0,y1;
                    double hE;
                    double **JointValues;
                    double *JointValue;
                    double FunctionalValues[MaxN_BaseFunctions2D];
                    double PointValues[MaxN_PointsForNodal2D];
                    int N_LinePoints;
#ifdef __3D__
                    double z0, z1;
#endif
                    
                    
                    if(Cond0 == Cond1)
                    {
                        switch(Cond0)
                        {
                            case DIRICHLET:
                                // if DG
                                if (N_EdgeDOF==0)
                                    break;
                                // read boundary values for each quadrature point
                                for(l=0;l<N_EdgePoints;l++)
                                {
                                    s = EdgePoints[l];
                                    t = 0.5*(t0*(1-s) + t1*(1+s));
                                    BoundaryValue(comp, t, PointValues[l]);
                                }                                 // endfor l
                                // compute boundary values for each dof on the
                                // boundary edge with the nodal functionals
                                
                                nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
                                                       FunctionalValues);
                                EdgeDOF = FEDesc_Obj->GetJointDOF(m);
                                // save boundary values of each dof on the boundary
                                // edge in the rhs
                                for( l=0;l<N_EdgeDOF;l++)
                                {
                                    RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
                                }
                                break;
                                
                            case NEUMANN:
                                // get polynomial degree of fe
                                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
                                // get a suitable line quadrature formula
                                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                                ->MakeRefElementData(LineQuadFormula);
                                JointValues=TFEDatabase2D::GetJointValues2D(
                                                                            BaseFuncts[CurrentElement], LineQuadFormula, m);
                                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                                ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                                // get vertices of boundary edge
                                
#ifdef __3D__
                                
                                cell->GetVertex(m)->GetCoords(x0, y0,z0);
                                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1,z1);
#else
                                cell->GetVertex(m)->GetCoords(x0, y0);
                                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                                // compute (half of the) length of the boundary edge
                                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                                // compute boundary integral
                                for( l=0;l<N_LinePoints;l++)
                                {
                                    // values of test functions in this quadrature point
                                    JointValue = JointValues[l];
                                    // get quadrature point on the boundary
                                    t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                                    // get value in this quadrature point (in s)
                                    BoundaryValue(comp, t, s);
                                    // multiply value with weights from quadrature formula
                                    // and determinant from integral transformation to the
                                    // unit edge (-1,1)
                                    s *= hE * LineWeights[l];
                                    // update rhs for all test functions
                                    for( k=0;k<N_;k++)
                                        if((l3 = DOF[k])<ActiveBound)
                                            RHS[l3] += s*JointValue[k];
                                }
                                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                                ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                                break;
                            case ROBIN:
                                // get polynomial degree of fe
                                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
                                // get a suitable line quadrature formula
                                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                                ->MakeRefElementData(LineQuadFormula);
                                JointValues=TFEDatabase2D::GetJointValues2D(
                                                                            BaseFuncts[CurrentElement], LineQuadFormula, m);
                                // get vertices of boundary edge
#ifdef __3D__
                                
                                cell->GetVertex(m)->GetCoords(x0, y0,z0);
                                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1,z1);
#else
                                cell->GetVertex(m)->GetCoords(x0, y0);
                                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif

                                // compute (half of the) length of the boundary edge
                                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                                // compute boundary integral
                                for(l=0;l<N_LinePoints;l++)
                                {
                                    // values of test functions in this quadrature point
                                    JointValue = JointValues[l];
                                    // get quadrature point on the boundary
                                    t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                                    // get value in this quadrature point (in s)
                                    BoundaryValue(comp, t, s);
                                    // multiply value with weights from quadrature formula
                                    // and determinant from integral transformation to the
                                    // unit edge (-1,1)
                                    s *= hE * LineWeights[l];
                                    // update rhs for all test functions
                                    for(k=0;k<N_;k++)
                                        if((l3 = DOF[k])<ActiveBound)
                                            RHS[l3] += s*JointValue[k];
                                }
                                break;
                            case SLIP:
                                OutPut("Use SLIP_FRICTION_PENETRATION_RESISTANCE boundary condition !"<< endl);
                                exit(4711);
                                break;
                            case SLIP_FRICTION_PENETRATION_RESISTANCE:
                                // do nothing here
                                // everything is done in Assemble2DSlipBC, see below
                                break;
                                
                            case FREESURF:
                            case INTERFACE:
                                // do nothing here
                                // everything is done in Freesurfint, see Freesurface2D.C
                                break;
                            case DIRICHLET_WEAK:
                                break;
                            default :
                                OutPut("Unknown boundary condition !"<< endl);
                                exit(4711);
                                
                        }                                     // endswitch Cond0
                    }                                       // endif (Cond0==Cond1)
                    else
                    {
                        OutPut("different boundary condition on one edge ");
                        OutPut("are not allowed!" << endl);
                        exit(4711);
                    }
                }                                         // endif (boundary joint)
            }                                           // endfor m (N_Joints)
        }                                             // endfor j (n_rhs)
    }                                               // endfor i (N_Cells)
    

}


Assembler3::~Assembler3(){};



