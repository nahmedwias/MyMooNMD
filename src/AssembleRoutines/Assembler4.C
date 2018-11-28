// =======================================================================
// @(#)Assembler4.C        1.16
//
// Purpose:     bilinear form (discretized and stabilized assemble)
//
// Author:      Alfonso Caiazzo & Laura Blank
//
// History:     start of implementation 10.08.2016
//
// =======================================================================

#include <Assembler4.h>
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



//================================================================================
Assembler4::Assembler4()
{
  Coll = nullptr;
  hangingEntries.resize(0);
  hangingRhs.resize(0);
  square_matrices.resize(0);
  rectangular_matrices.resize(0);
  rhs_blocks.resize(0);
  n_square_matrices = 0;
  n_rectangular_matrices = 0;
  n_rhs_blocks = 0;
  n_all_matrices = 0;
  maximum_number_base_functions = 0;
}

//================================================================================
void Assembler4::init(BlockFEMatrix &M,
    BlockVector &rhs,
    std::vector<const TFESpace2D*>& fespaces,
    std::vector<const TFESpace2D*>& ferhs)
{
  // set the collection
  this->Coll = fespaces[0]->GetCollection();

  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  this->maximum_number_base_functions = 0;

  /// this loop sets the CellIndex. Probably it should be done somewhere else?
  for (int i = 0; i < this->Coll->GetN_Cells(); i++)
  {
    TBaseCell *cell = this->Coll->GetCell(i);
    cell->SetCellIndex(i);

    // find the maximum number of basis functions used by the current collection
    for (unsigned int j = 0; j < fespaces.size(); j++)
    {
      FE2D CurrentElement = fespaces[j]->GetFE2D(i, this->Coll->GetCell(i));
      if (N_BaseFunct[CurrentElement] > this->maximum_number_base_functions)
      {
        this->maximum_number_base_functions = N_BaseFunct[CurrentElement];
      }
    }
  }
  // --------------------------------------------------------------------
  Output::print("Assembler: maximum_number_base_functions = ", maximum_number_base_functions);


  // set vector of blocks
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  // Note: This class is supposed to work only for NSType 14 at the moment
  if(blocks.size() == 9)
  {
    this->n_square_matrices = 5;
    square_matrices.resize(n_square_matrices);
    square_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    square_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
    square_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
    square_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
    square_matrices[4] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());

    this->n_rectangular_matrices = 4;
    rectangular_matrices.resize(n_rectangular_matrices);
    rectangular_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get());
    rectangular_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
    rectangular_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
    rectangular_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
  }
  else if(blocks.size() == 1)
  {
    this->n_square_matrices = 1;
    square_matrices.resize(n_square_matrices);
    square_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

    this->n_rectangular_matrices = 0;
    //rectangular_matrices.resize(n_rectangular_matrices);
  }
  else
  {
    Output::print("Assembler should be only used with NSType 14 Matrices");
    ErrThrow(" --> Wrong blocks.size() = ", blocks.size());
  }

  this->n_rhs_blocks = ferhs.size();
  rhs_blocks.resize(n_rhs_blocks);

  for (size_t k = 0; k < n_rhs_blocks; k++)
  {
    rhs_blocks[k] = rhs.block(k);
  }

  this->n_all_matrices = this->n_square_matrices + this->n_rectangular_matrices;

  // set the vectors for hanging nodes
  this->hangingEntries.resize(this->n_square_matrices + this->n_rectangular_matrices);

  for (int i = 0; i < this->n_square_matrices; i++)
  {
    int j = this->square_matrices[i]->GetHangingN_Entries();
    this->hangingEntries[i].resize(j);
  }

  for (int i = 0; i < this->n_rectangular_matrices; i++)
  {
    int j = rectangular_matrices[i]->GetHangingN_Entries();
    this->hangingEntries[i + this->n_square_matrices].resize(j);
  }

  this->hangingRhs.resize(rhs_blocks.size());

  for (size_t i = 0; i < this->n_rhs_blocks; i++)
  {
    int j = ferhs[i]->GetN_Hanging();
    this->hangingRhs[i].resize(j);
  }
}


//================================================================================
//================================================================================
void Assembler4::Assemble2D(BlockFEMatrix &M,
    BlockVector &b_rhs,
    std::vector<const TFESpace2D*>& fespaces,
    std::vector<const TFESpace2D*>& ferhs,
    const Example2D& example,
    std::vector< std::shared_ptr< LocalAssembling2D> > la_list,
    int AssemblePhaseID)
{
#ifdef __3D__
  ErrThrow("Assembler4::Assembler4() not yet working in 3D");
#endif

  // set matrices and rhs blocks
  this->init(M, b_rhs, fespaces, ferhs);

  // LocRhs: an array of pointers (size: number of rhs)
  double **LocRhs;
  // righthand: a big pointer
  double *righthand;
  if (n_rhs_blocks)
  {
    LocRhs = new double* [n_rhs_blocks];
    righthand = new double [n_rhs_blocks* maximum_number_base_functions];
    for (size_t i = 0; i < n_rhs_blocks; i++)
      LocRhs[i] = righthand + i * maximum_number_base_functions;
  }

  // set pointers to matrices
  double *aux;
  ///@todo all these pointers are used but might not be initialized
  double **Matrices;
  double ***LocMatrices;

  if (this->n_all_matrices)
  {
    aux = new double
        [this->n_all_matrices * maximum_number_base_functions * maximum_number_base_functions];
    Matrices = new double* [this->n_all_matrices*maximum_number_base_functions];

    for (int j = 0; j < (this->n_all_matrices * maximum_number_base_functions); j++)
      Matrices[j] = aux + j * maximum_number_base_functions;

    // LocMatrices[i] = a matrix of size NBaseFct x NBaseFct
    LocMatrices = new double** [this->n_all_matrices];

    for (int i = 0; i < this->n_all_matrices; i++)
      LocMatrices[i] = Matrices + i * maximum_number_base_functions;
  }

  //================================================================================
  // loop over all cells
  //================================================================================
  for (int i = 0; i < this->Coll->GetN_Cells(); i++)
  {
    TBaseCell *cell = this->Coll->GetCell(i);

    ///@attention the value of INTERNAL_CELL appears to be used somewhere later
    TDatabase::ParamDB->INTERNAL_CELL = i;

    // only for multiphase flows
    if ((AssemblePhaseID >= 0) &&
        (AssemblePhaseID != cell->GetPhase_ID()) )
      continue;
    for (size_t num_localAss = 0; num_localAss < la_list.size(); num_localAss++)
    {
      //assemble the selected local form on the i-th cell
      this->assemble_local_system(fespaces, i, LocMatrices, LocRhs, la_list[num_localAss]);

      // add local/cellwise matrices to global matrices
      //(ansatz == test: Aii, C; ansatz != test: A12,A12,B1,...)
      this->add_local_to_global_matrix(i, LocMatrices, Matrices);

      // add local/cellwise right-hand sides to global right-hand sides
      this->add_local_to_global_rhs(i, ferhs, LocRhs, example);

      // boundary condition part
      this->impose_boundary_conditions(i, ferhs, example);
    }
  }


  // modify matrix according to coupling of hanging nodes
  // this part is only relevant to the case with hanging nodes and
  // it could be put in a separate function
  this->handle_hanging_nodes(ferhs);

  // delete pointers
  if (n_rhs_blocks)
  {
    delete [] righthand;
    delete [] LocRhs;
  }

  if (this->n_all_matrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }
}


//================================================================================
// 
//================================================================================
void Assembler4::assemble_local_system(std::vector <const TFESpace2D*>& fespaces,
    int i, double ***LocMatrices, double **LocRhs,
    std::shared_ptr<LocalAssembling2D> la)
{
  BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  TBaseCell *cell = this->Coll->GetCell(i);

  // find local (on this cell) used elements
  FE2D LocalUsedElements[N_FEs2D];
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  int n_fespaces = fespaces.size();

  for (int j = 0; j < n_fespaces; j++)
  {
    FE2D CurrentElement = fespaces[j]->GetFE2D(i, cell);
    LocalUsedElements[j] = CurrentElement;
    LocN_BF[j] = N_BaseFunct[CurrentElement];
    LocBF[j] = BaseFuncts[CurrentElement];
  }

    // --------------------------------------------------------------------
    // calculate values on original element
    
    int N_Points;
    const double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D],Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];

  bool *SecondDer = la->GetNeeds2ndDerivatives();

  TFEDatabase2D::GetOrig(n_fespaces, LocalUsedElements,
      this->Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);

  bool is_sdfem = (la->get_disctype() == SDFEM);

  if (is_sdfem
      || (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM)
      || (TDatabase::ParamDB->CELL_MEASURE == 4))
  {
    TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
    int N_Edges = cell->GetN_Edges();
    for (int ij = 0; ij < N_Edges; ij++)
    {
      TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
      TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
    }
    if (N_Edges==3)
      TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
  }

  la->GetLocalForms(N_Points, weights, AbsDetjk, {{X, Y}}, LocN_BF, LocBF,
      cell, i, n_all_matrices, n_rhs_blocks, LocMatrices,
      LocRhs);
}

//================================================================================
// impose boundary conditions
//================================================================================
void Assembler4::impose_boundary_conditions(int i_cell,
    std::vector<const TFESpace2D*>& ferhs,
    const Example2D& example)
{
  TBaseCell *cell = this->Coll->GetCell(i_cell);
  for(unsigned int j=0;j<rhs_blocks.size();j++)
    {
      const TFESpace2D *fespace = ferhs[j];
      int *DOF = ferhs[j]->GetGlobalDOF(i_cell);
      double *RHS = this->rhs_blocks[j];
      
      FE2D CurrentElement = fespace->GetFE2D(i_cell, cell);
      BoundCondFunct2D *BoundaryCondition = fespace->get_boundary_condition();
      BoundValueFunct2D *BoundaryValue = example.get_bd()[j];
      TFE2D *ele = TFEDatabase2D::GetFE2D(CurrentElement);
      double t0,t1;
      const TBoundComp *BoundComp;
      int N_EdgePoints;
      const double *EdgePoints;
    double eps = 1e-4;

    TNodalFunctional2D *nf = ele->GetNodalFunctional2D();

    if (TDatabase::ParamDB->SUPERCONVERGENCE_ORDER)
    {
      /* Superconvergence boundary interpolation */
      if (nf->GetID() == NF_C_Q_Q2_2D)
        nf = TFEDatabase2D::GetNodalFunctional2D(NF_S_Q_Q2_2D);
    }

    nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

    TFEDesc2D *FEDesc_Obj = ele->GetFEDesc2D();
    int  N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
    // setting Dirichlet boundary condition
    BoundCond Cond0, Cond1;

    for (int m = 0; m < cell->GetN_Joints(); m++)
    {
      TJoint *joint = cell->GetJoint(m);

      if (joint->GetType() == BoundaryEdge ||
          joint->GetType() == IsoBoundEdge ||
          joint->GetType() == InterfaceJoint)
      {
        if (joint->GetType() == BoundaryEdge||
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
        int comp = BoundComp->GetID();
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

        double FunctionalValues[MaxN_BaseFunctions2D];
        double PointValues[MaxN_PointsForNodal2D];

        // only one boundary condition per edge allowed
        if(Cond0 == Cond1)
        {
          switch(Cond0)
          {
          case DIRICHLET:
          {
            // if DG
            if (N_EdgeDOF == 0)
              break;
            // read boundary values for each quadrature point
            for (int l = 0; l < N_EdgePoints; l++)
            {
              double s = EdgePoints[l];
              double t = 0.5 * (t0 * (1-s) + t1 * (1+s));
              BoundaryValue(comp, t, PointValues[l]);
            }

            // compute boundary values for each dof on the
            // boundary edge with the nodal functionals
            nf->GetEdgeFunctionals(this->Coll, cell, m, PointValues,
                FunctionalValues);
            int *EdgeDOF = FEDesc_Obj->GetJointDOF(m);
            // save boundary values of each dof on the boundary
            // edge in the rhs
            for (int l = 0; l < N_EdgeDOF; l++)
            {
              RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
            }
            break;
          }
          case NEUMANN:
          case ROBIN:
            Output::print<4>(" ** Assembler: WARNING: NEUMANN and ROBIN boundary conditions are not supported. ** ");
            Output::print<4>(" ** You can try to use the old Assemble2D(...) instead ** ");
            break;
          case DIRICHLET_WEAK:
            break;
          default :
            OutPut("Unknown boundary condition !"<< endl);
            exit(4711);
          }                                     // endswitch Cond0
        }                                       // end if (Cond0 == Cond1)
        else
        {
          OutPut("different boundary condition on one edge ");
          OutPut("are not allowed!" << endl);
          exit(4711);
        }
      }                                         // endif (boundary joint)
    }
  }
}


//===================================================================
// add local/cellwise matrices to global matrices
//(ansatz == test: Aii, C; ansatz != test: A12,A12,B1,...)
//===================================================================
void Assembler4::add_local_to_global_matrix(int i,
    double ***LocMatrices,
    double **Matrices)
{
  TBaseCell *cell = this->Coll->GetCell(i);
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // square matrices
  for (int j = 0; j < this->n_square_matrices; j++)
  {
    double **Matrix = LocMatrices[j];
    std::vector<double> CurrentHangingEntries = this->hangingEntries[j];
    const TFESpace2D *fespace = square_matrices[j]->GetFESpace2D();
    const int *HangingRowPtr = square_matrices[j]->GetHangingRowPtr();
    const int *HangingColInd = square_matrices[j]->GetHangingKCol();
    int ActiveBound = fespace->GetActiveBound();
    int DirichletBound = fespace->GetHangingBound();

    /**
         DOF[k] is the global index of the k-th local degree of freedom
         MatrixRow[k] is the assembled value corresponding to the m-th
         local test function and k-th local ansatz function. That means it
         corresponds to the l=DOF[m]-th global test function and the
         DOF[k]-th global ansatz function
     */
    //int *DOF = GlobalNumbers[j] + BeginIndex[j][i];
    ///@todo use a vector<int> dof = fespace.local_to_global_DOFs(i,cell)
    int *DOF = fespace->GetGlobalDOF(i);
    FE2D CurrentElement = fespace->GetFE2D(i, cell);
    /**
         The list of degrees of freedom consists of 3 parts:
         0,...,ActiveBounds-1: active nodes (no hanging, no Dirichlet)
         ActiveBounds, DirichletBounds-1: Hanging nodes
         DirichletBound,...end : the DOF with DIRICHLET condition
     */

    // add local matrix to global matrix
    for (int m = 0; m < N_BaseFunct[CurrentElement]; m++)
    {
      // active DOF
      if (DOF[m] < ActiveBound)
      {
        for (int k = 0; k < N_BaseFunct[CurrentElement]; k++)
          square_matrices[j]->add(DOF[m], DOF[k], Matrix[m][k]);
      }
      else
      {
        // non-active but not Dirichlet (hanging)
        if (DOF[m] < DirichletBound)
        {
          int l = DOF[m] - ActiveBound;
          int end = HangingRowPtr[l+1];

          for (int n = HangingRowPtr[l]; n < end; n++)
          {
            for (int k = 0; k < N_BaseFunct[CurrentElement]; k++)
            {
              if(DOF[k] == HangingColInd[n])
              {
                CurrentHangingEntries[n] += Matrix[m][k];
                break;
              }
            }
          }
        }
        else
        {
          // Dirichlet node
          if (TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE)
          {
            /*  When applying FEM-FCT, we treat Dirichlet nodes like active ones */
            for (int k = 0; k < N_BaseFunct[CurrentElement]; k++)
            {
              square_matrices[j]->add(DOF[m], DOF[k], Matrix[m][k]);
            }
          }
          else
          {
            //Dirichlet standard treatment: 1.0 on the diagonal
            square_matrices[j]->set(DOF[m], DOF[m], 1.0);
          }
        }
      }
    }                                           // endfor m
  }                                             // endfor j




  // rectangular matrices
  for (int j = 0; j < this->n_rectangular_matrices; j++)
  {
    FE2D TestElement = ((TFESpace2D *) rectangular_matrices[j]->GetTestSpace())
                ->GetFE2D(i, cell);
    FE2D AnsatzElement = ((TFESpace2D *) rectangular_matrices[j]->GetAnsatzSpace())
                ->GetFE2D(i, cell);

    int N_Test = N_BaseFunct[TestElement];
    int N_Ansatz = N_BaseFunct[AnsatzElement];

    double **Matrix = Matrices + (j+this->n_square_matrices) * maximum_number_base_functions;

    double *Entries = rectangular_matrices[j]->GetEntries();
    const int *RowPtr = rectangular_matrices[j]->GetRowPtr();
    const int *ColInd = rectangular_matrices[j]->GetKCol();

    //int *TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
    //int *AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];
    int *TestDOF =
        ((TFESpace2D *) rectangular_matrices[j]->GetTestSpace())->GetGlobalDOF(i);
    int *AnsatzDOF =
        ((TFESpace2D *) rectangular_matrices[j]->GetAnsatzSpace())->GetGlobalDOF(i);

    const TFESpace2D *fespace =
        (TFESpace2D *)(rectangular_matrices[j]->GetTestSpace2D());
    int ActiveBound = fespace->GetActiveBound();
    int DirichletBound = fespace->GetHangingBound();

    std::vector<double> CurrentHangingEntries = this->hangingEntries[j+this->n_square_matrices];
    const int *HangingRowPtr = rectangular_matrices[j]->GetHangingRowPtr();
    const int *HangingColInd = rectangular_matrices[j]->GetHangingKCol();

    // add local matrix to global
    for (int m = 0; m < N_Test; m++)
    {
      // the entries are added only for active nodes and for Dirichlet nodes
      if ( TestDOF[m] < ActiveBound || TestDOF[m]>=DirichletBound)
      {
        int end = RowPtr[ TestDOF[m] + 1];

        for (int n = RowPtr[ TestDOF[m] ]; n < end; n++)
        {
          for (int k = 0; k < N_Ansatz; k++)
          {
            if (AnsatzDOF[k] == ColInd[n])
            {
              Entries[n] += Matrix[m][k];
              break;
            }
          }
        }
      }
      else
      {
        // hanging node
        int l =  TestDOF[m] - ActiveBound;
        int end = HangingRowPtr[l+1];

        for (int n = HangingRowPtr[l]; n < end; n++)
        {
          for (int k = 0; k < N_Ansatz; k++)
          {
            if (AnsatzDOF[k] == HangingColInd[n])
            {
              CurrentHangingEntries[n] += Matrix[m][k];
              break;
            }
          }
        }
      }
    }
  }
}


//================================================================================
// add local right-hand sides to global right-hand side
//================================================================================
void Assembler4::add_local_to_global_rhs(int i,
    std::vector<const TFESpace2D*>& ferhs,
    double **LocRhs,
    const Example2D& example)
{
  TBaseCell *cell = this->Coll->GetCell(i);
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  for (unsigned int j = 0; j < rhs_blocks.size(); j++)
  {
    const TFESpace2D *fespace = ferhs[j];
    int ActiveBound = fespace->GetActiveBound();
    FE2D CurrentElement = fespace->GetFE2D(i, cell);

    int N_ = N_BaseFunct[CurrentElement];

    // get the pointer to j-th RHS from the global righthandside
    // local_rhs = LocRhs[j] ?
    //double *local_rhs = righthand+j*maximum_number_base_functions;
    double *RHS_block_j = rhs_blocks[j];
    std::vector<double> CurrentHangingRhs = this->hangingRhs[j];
    // find space for this linear form

    ActiveBound = fespace->GetActiveBound();
    int DirichletBound = fespace->GetHangingBound();

    // dof of the rhs nodes connected to this cell
    int *DOF = ferhs[j]->GetGlobalDOF(i);
    //RhsGlobalNumbers[j] + RhsBeginIndex[j][i];

    // add local right-hand side to the global one
    for (int m = 0; m < N_; m++)
    {
      /*
             /// Alfonso & Laura: this part was in the Assemble2D.
             if (TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS)
             {
             int l = 0;
             int N_Edges=cell->GetN_Edges();
             for(int jj=0;jj<N_Edges;jj++)                        // loop over all edges of cell
             {
             TVertex *ver0 = cell->GetVertex(jj);
             double t0 =  ver0->GetClipBoard();
             // vertex not on the boundary
             if (t0<-1e-8)
             continue;
             // component of boundary
             int comp = floor(t0+1e-8);
             // parameter
             t0 -= comp;
             // get boundary condition
             BoundaryConditions[0](comp, t0, Cond0);
             // Dirichlet
             if (Cond0== DIRICHLET)
             {
             l = -4711;
             }
             }
             if (l==-4711)
             break; // do nothing for this mesh cell
             }
       */
      int l = DOF[m];
      if (l < ActiveBound)
      {
        // node l is inner node
        RHS_block_j[l] += LocRhs[j][m];
      }
      else
      {
        if (l < DirichletBound)
        {
          // hanging node
          l -= ActiveBound;
          CurrentHangingRhs[l] += LocRhs[j][m];
        }
      }
    }
  }
}



//================================================================================
// hanging nodes
//================================================================================
void Assembler4::handle_hanging_nodes(std::vector<const TFESpace2D*>& ferhs)
{
  for (int j = 0; j < this->n_square_matrices; j++)
  {
    const TFESpace2D *fespace = square_matrices[j]->GetFESpace2D();
    THangingNode **HangingNodes = fespace->GetHangingNodes();

    double *Entries = square_matrices[j]->GetEntries();
    const int *RowPtr = square_matrices[j]->GetRowPtr();
    const int *ColInd = square_matrices[j]->GetKCol();

    std::vector<double> CurrentHangingEntries = this->hangingEntries[j];
    const int *HangingRowPtr = square_matrices[j]->GetHangingRowPtr();
    const int *HangingColInd = square_matrices[j]->GetHangingKCol();

    int ActiveBound = fespace->GetActiveBound();

    for (int i = 0; i < fespace->GetN_Hanging(); i++)
    {
      THangingNode *hn = HangingNodes[i];
      HNDesc HNDescr = hn->GetType();
      THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      int k = HNDescr_Obj->GetN_Nodes();
      double *Coupling = HNDescr_Obj->GetCoeff();
      int *DOF = hn->GetDOF();

      for (int n = HangingRowPtr[i]; n < HangingRowPtr[i+1]; n++)
      {
        //double v = CurrentHangingEntries[n];
        int m = HangingColInd[n];

        for (int l = 0; l < k; l++)
        {
          if (DOF[l] < ActiveBound)
          {
            for (int l2 = RowPtr[ DOF[l] ]; l2 < RowPtr[ DOF[l]+1 ]; l2++)
            {
              if (ColInd[l2] == m)
              {
                Entries[l2] += Coupling[l] * CurrentHangingEntries[n];
              }
            }
          }
        }
      }
    }
  }

  for (int j = 0; j < this->n_rectangular_matrices; j++)
  {
    // hanging nodes in test space
    const TFESpace2D *fespace = (TFESpace2D *) (rectangular_matrices[j]->GetTestSpace2D());
    THangingNode **HangingNodes = fespace->GetHangingNodes();

    double *Entries = rectangular_matrices[j]->GetEntries();
    const int *RowPtr = rectangular_matrices[j]->GetRowPtr();
    const int *ColInd = rectangular_matrices[j]->GetKCol();

    std::vector<double> CurrentHangingEntries = this->hangingEntries[j + this->n_square_matrices];

    const int *HangingRowPtr = rectangular_matrices[j]->GetHangingRowPtr();
    const int *HangingColInd = rectangular_matrices[j]->GetHangingKCol();

    int ActiveBound = fespace->GetActiveBound();

    for (int i = 0; i < fespace->GetN_Hanging(); i++)
    {
      THangingNode *hn = HangingNodes[i];
      HNDesc HNDescr = hn->GetType();
      THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      int k = HNDescr_Obj->GetN_Nodes();
      double *Coupling = HNDescr_Obj->GetCoeff();
      int *DOF = hn->GetDOF();

      for (int n = HangingRowPtr[i]; n < HangingRowPtr[i+1]; n++)
      {
        int m = HangingColInd[n];
        for (int l = 0; l < k; l++)
        {
          if (DOF[l] < ActiveBound)
          {
            for (int l2 = RowPtr[ DOF[l] ]; l2 < RowPtr[ DOF[l]+1 ]; l2++)
            {
              if (ColInd[l2] == m)
              {
                Entries[l2] += Coupling[l] * CurrentHangingEntries[n];
              }
            }
          }
        }
      }
    }

    // hanging nodes in ansatz space
    //int N_Rows =  matrices[j]->GetN_Rows();
    fespace = (TFESpace2D *) (rectangular_matrices[j]->GetAnsatzSpace2D());

    HangingNodes = fespace->GetHangingNodes();
    int AnsatzActiveBound = fespace->GetActiveBound();
    int AnsatzHangingBound = fespace->GetHangingBound();

    for (int i = 0; i < rectangular_matrices[j]->GetN_Rows(); i++)
    {
      //      int end = RowPtr[i+1];
      for (int k = RowPtr[i]; k < RowPtr[i+1]; k++)
      {
        int l = ColInd[k];

        if (l >= AnsatzActiveBound && l < AnsatzHangingBound)
        {
          // l is hanging node in ansatz space
          THangingNode *hn = HangingNodes[l - AnsatzActiveBound];
          HNDesc HNDescr = hn->GetType();
          THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
          int m = HNDescr_Obj->GetN_Nodes();
          double *Coupling = HNDescr_Obj->GetCoeff();
          int *DOF = hn->GetDOF();
          //double v = Entries[k];

          for (int n = 0; n < m; n++)
          {
            for (int l1 = RowPtr[i]; l1 < RowPtr[i+1]; l1++)
            {
              if (ColInd[l1] == DOF[n])
                Entries[l1] += Entries[k] * Coupling[n];
            }
          }
          Entries[k] = 0;
        }
      }
    }
  }

  for (size_t j = 0; j < n_rhs_blocks; j++)
  {
    const TFESpace2D *fespace = ferhs[j];
    THangingNode **HangingNodes = fespace->GetHangingNodes();

    double *RHS = rhs_blocks[j];
    std::vector<double> CurrentHangingRhs = this->hangingRhs[j];

    int ActiveBound = fespace->GetActiveBound();

    for (int i = 0; i < fespace->GetN_Hanging(); i++)
    {
      THangingNode *hn = HangingNodes[i];
      HNDesc HNDescr = hn->GetType();
      THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      int N_ = HNDescr_Obj->GetN_Nodes();
      double *Coupling = HNDescr_Obj->GetCoeff();
      int *DOF = hn->GetDOF();

      for (int k = 0; k < N_; k++)
      {
        if ( DOF[k] < ActiveBound )
        {
          RHS[ DOF[k] ] += Coupling[k] * CurrentHangingRhs[i];
        }
      }
    }
  }


  //================================================================================
  // write coupling into matrix
  //================================================================================
  for (int j = 0; j < this->n_square_matrices; j++)
  {
    const TFESpace2D *fespace = square_matrices[j]->GetFESpace2D();
    THangingNode **HangingNodes = fespace->GetHangingNodes();

    double *Entries = square_matrices[j]->GetEntries();
    const int *RowPtr = square_matrices[j]->GetRowPtr();
    int ActiveBound = fespace->GetActiveBound();
    int n = RowPtr[ActiveBound];

    for (int i = 0; i < fespace->GetN_Hanging(); i++)
    {
      THangingNode *hn = HangingNodes[i];
      HNDesc HNDescr = hn->GetType();
      THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      int k = HNDescr_Obj->GetN_Nodes();
      double *Coupling = HNDescr_Obj->GetCoeff();

      Entries[n] = 1.0;
      n++;

      for (int l = 0; l < k; l++)
      {
        Entries[n] = - Coupling[l];
        n++;
      }
    }
  }
}


//================================================================================
Assembler4::~Assembler4(){};




