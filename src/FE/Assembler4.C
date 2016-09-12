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

#include <DefineParams.h>

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
//================================================================================
Assembler4::Assembler4(LocalAssembling2D_type type,
                       TFEFunction2D **fefunctions2d,
                       CoeffFct2D *coeffs)
:la(type,fefunctions2d,coeffs)
{
    Coll = NULL;
    hangingEntries.resize(0);
    hangingRhs.resize(0);
    square_matrices.resize(0);
    rectangular_matrices.resize(0);
    rhs_blocks.resize(0);
    n_square_matrices=0;
    n_rectangular_matrices=0;
    n_rhs_blocks=0;
}

//================================================================================
//================================================================================
void Assembler4::init(BlockFEMatrix &M,
                      BlockVector &rhs,
                      std::vector<const TFESpace2D*>& fespaces,
                      std::vector<const TFESpace2D*>& ferhs)
{
    // set the collection
    this->Coll = fespaces[0]->GetCollection();
    
    // --------------------------------------------------------------------
    /// this loop set the CellIndex. Probably it should be done somewhere else?
    for(int i=0; i<this->Coll->GetN_Cells(); i++)
    {
        TBaseCell *cell = this->Coll->GetCell(i);
        cell->SetCellIndex(i);
    }
    // --------------------------------------------------------------------
    
    // set vector of blocks
    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
    // Note: This class is supposed to work only for NSType 14 at the moment
    if(blocks.size() != 9)
    {
        Output::print("Assembler should be only used with NSType 14 Matrices");
        ErrThrow(" --> Wrong blocks.size() = ", blocks.size());
    }
    
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
    
    this->n_rhs_blocks = ferhs.size();
    rhs_blocks.resize(n_rhs_blocks);
    for (size_t k=0; k<n_rhs_blocks;k++)
    {
    rhs_blocks[k] = rhs.block(k);
    }
    
    this->N_AllMatrices = this->n_square_matrices + this->n_rectangular_matrices;
    
    // --------------------------------------------------------------------
    // set the vectors for hanging nodes
    
    this->hangingEntries.resize(this->n_square_matrices + this->n_rectangular_matrices);
    
    for(int i=0 ; i< this->n_square_matrices; i++)
    {
        int j = this->square_matrices[i]->GetHangingN_Entries();
        this->hangingEntries[i].resize(j);
    }
    for(int i=0;i<this->n_rectangular_matrices;i++)
    {
        int j = rectangular_matrices[i]->GetHangingN_Entries();
        this->hangingEntries[i+this->n_square_matrices].resize(j);
    }
    this->hangingRhs.resize(rhs_blocks.size());
    for(int i=0;i<this->n_rhs_blocks;i++)
    {
        int j = ferhs[i]->GetN_Hanging();
        this->hangingRhs[i].resize(j);
    }
    // --------------------------------------------------------------------
}




//================================================================================
// impose boundary conditions
//================================================================================
void Assembler4::impose_boundary_conditions(const TFESpace2D *fespace,
                                          const Example2D& example,
                                          TBaseCell* cell,
                                          int i, int j,
                                          double *RHS,
                                          int *DOF)
{
    FE2D CurrentElement = fespace->GetFE2D(i, cell);
    BoundCondFunct2D *BoundaryCondition = fespace->GetBoundCondition();
    BoundValueFunct2D *BoundaryValue = example.get_bd()[j];
    TFE2D *ele = TFEDatabase2D::GetFE2D(CurrentElement);
    double t0,t1;
    TBoundComp *BoundComp;
    int N_EdgePoints;
    double *EdgePoints;
    double eps=1e-4;
    
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
    int N_Joints = cell->GetN_Joints();
    BoundCond Cond0, Cond1;
    
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
                        if (N_EdgeDOF==0)
                            break;
                        // read boundary values for each quadrature point
                        for(int l=0;l<N_EdgePoints;l++)
                        {
                            double s = EdgePoints[l];
                            double t = 0.5*(t0*(1-s) + t1*(1+s));
                            BoundaryValue(comp, t, PointValues[l]);
                        }                                 // endfor l
                        // compute boundary values for each dof on the
                        // boundary edge with the nodal functionals
                        
                        nf->GetEdgeFunctionals(this->Coll, cell, m, PointValues,
                                               FunctionalValues);
                        int *EdgeDOF = FEDesc_Obj->GetJointDOF(m);
                        // save boundary values of each dof on the boundary
                        // edge in the rhs
                        for(int l=0;l<N_EdgeDOF;l++)
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
            }                                       // endif (Cond0==Cond1)
            else
            {
                OutPut("different boundary condition on one edge ");
                OutPut("are not allowed!" << endl);
                exit(4711);
            }
        }                                         // endif (boundary joint)
    }
}


//================================================================================
// hanging nodes
//================================================================================
void Assembler4::handle_hanging_nodes(std::vector<const TFESpace2D*>& ferhs)
{
    for(int j=0;j< this->n_square_matrices;j++)
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
        
        for(int i=0;i< fespace->GetN_Hanging() ;i++)
        {
            THangingNode *hn = HangingNodes[i];
            HNDesc HNDescr = hn->GetType();
            THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
            int k = HNDescr_Obj->GetN_Nodes();
            double *Coupling = HNDescr_Obj->GetCoeff();
            int *DOF = hn->GetDOF();
            
            for(int n = HangingRowPtr[i]; n<HangingRowPtr[i+1]; n++)
            {
                //double v = CurrentHangingEntries[n];
                int m = HangingColInd[n];
                for(int l=0; l<k; l++) {
                    if(DOF[l] < ActiveBound) {
                        for(int l2=RowPtr[ DOF[l] ]; l2<RowPtr[ DOF[l]+1 ]; l2++) {
                            if(ColInd[l2] == m) {
                                Entries[l2] += Coupling[l] * CurrentHangingEntries[n];
                            }
                        }                                     // endfor l2
                    }                                       // endif
                }                                         // endfor l
            }                                           // endfor n
        }                                             // endfor i
        
    }                                               // endfor j
    
    for(int j=0;j<this->n_rectangular_matrices;j++)
    {
        // hanging nodes in test space
        const TFESpace2D *fespace = (TFESpace2D *) (rectangular_matrices[j]->GetTestSpace2D());
        THangingNode **HangingNodes = fespace->GetHangingNodes();
        
        double *Entries = rectangular_matrices[j]->GetEntries();
        const int *RowPtr = rectangular_matrices[j]->GetRowPtr();
        const int *ColInd = rectangular_matrices[j]->GetKCol();
        
        std::vector<double> CurrentHangingEntries = this->hangingEntries[j+this->n_square_matrices];
        
        const int *HangingRowPtr = rectangular_matrices[j]->GetHangingRowPtr();
        const int *HangingColInd = rectangular_matrices[j]->GetHangingKCol();
        
        int ActiveBound = fespace->GetActiveBound();
        
        for(int i=0;i<fespace->GetN_Hanging();i++)
        {
            THangingNode *hn = HangingNodes[i];
            HNDesc HNDescr = hn->GetType();
            THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
            int k = HNDescr_Obj->GetN_Nodes();
            double *Coupling = HNDescr_Obj->GetCoeff();
            int *DOF = hn->GetDOF();
            
            for(int n=HangingRowPtr[i]; n<HangingRowPtr[i+1]; n++) {
                int m = HangingColInd[n];
                for(int l=0;l<k;l++) {
                    if( DOF[l]<ActiveBound ) {
                        for(int l2=RowPtr[ DOF[l] ]; l2<RowPtr[ DOF[l]+1 ]; l2++) {
                            if(ColInd[l2] == m) {
                                Entries[l2] += Coupling[l] * CurrentHangingEntries[n];
                            }
                        }                                     // endfor l2
                    }                                       // endif
                }                                         // endfor l
            }                                           // endfor n
        }                                             // endfor i
        
        // hanging nodes in ansatz space
        //int N_Rows =  matrices[j]->GetN_Rows();
        fespace = (TFESpace2D *) (rectangular_matrices[j]->GetAnsatzSpace2D());
        
        HangingNodes = fespace->GetHangingNodes();
        int AnsatzActiveBound = fespace->GetActiveBound();
        int AnsatzHangingBound = fespace->GetHangingBound();
        for(int i=0; i<rectangular_matrices[j]->GetN_Rows(); i++)
        {
            //      int end = RowPtr[i+1];
            for(int k=RowPtr[i]; k<RowPtr[i+1];k++)
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
                    //double v = Entries[k];
                    
                    for(int n=0; n<m; n++)
                    {
                        for(int l1=RowPtr[i]; l1<RowPtr[i+1]; l1++)
                        {
                            if(ColInd[l1] == DOF[n])
                                Entries[l1] += Entries[k]*Coupling[n];
                        }
                    }
                    Entries[k] = 0;
                }                                         // endif l
            }                                           // endfor k
        }                                             // endfor i
    }                                               // endfor j
    
    for(int j=0;j<n_rhs_blocks;j++)
    {
        const TFESpace2D *fespace = ferhs[j];
        THangingNode **HangingNodes = fespace->GetHangingNodes();
        
        double *RHS = rhs_blocks[j];
        std::vector<double> CurrentHangingRhs = this->hangingRhs[j];
        
        int ActiveBound = fespace->GetActiveBound();
        
        for(int i=0; i<fespace->GetN_Hanging(); i++)
        {
            THangingNode *hn = HangingNodes[i];
            HNDesc HNDescr = hn->GetType();
            THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
            int N_ = HNDescr_Obj->GetN_Nodes();
            double *Coupling = HNDescr_Obj->GetCoeff();
            int *DOF = hn->GetDOF();
            
            for(int k=0;k<N_;k++)
            {
                if( DOF[k]<ActiveBound )
                {
                    RHS[ DOF[k] ] += Coupling[k] * CurrentHangingRhs[i];
                }
            }                                           // endfor k
        }                                             // endfor i
    }                                               // endfor j
    
    
    //================================================================================
    // write coupling into matrix
    //================================================================================
    for(int j=0;j<this->n_square_matrices;j++)
    {
        const TFESpace2D *fespace = square_matrices[j]->GetFESpace2D();
        THangingNode **HangingNodes = fespace->GetHangingNodes();
        
        double *Entries = square_matrices[j]->GetEntries();
        const int *RowPtr = square_matrices[j]->GetRowPtr();
        
        int ActiveBound = fespace->GetActiveBound();
        
        int n = RowPtr[ActiveBound];
        
        for(int i=0; i<fespace->GetN_Hanging(); i++)
        {
            THangingNode *hn = HangingNodes[i];
            HNDesc HNDescr = hn->GetType();
            THNDesc *HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
            int k = HNDescr_Obj->GetN_Nodes();
            double *Coupling = HNDescr_Obj->GetCoeff();
            
            Entries[n] = 1.0;
            n++;
            for(int l=0;l<k;l++)
            {
                Entries[n] = - Coupling[l];
                n++;
            }                                           // endfor l
        }                                             // endfor i
    }
}


//========================================================================================================
// add local/cellwise matrices to global matrices (ansatz == test: Aii, C; ansatz != test: A12,A12,B1,...)
//========================================================================================================
void Assembler4::add_local_to_global_matrix(int i,
                                        double ***LocMatrices,
                                        double **Matrices)
{
    TBaseCell *cell = this->Coll->GetCell(i);
    int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
    
    // --------------------------------------------------------------------
    // square matrices
    for(int j=0; j < this->n_square_matrices; j++)
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
        int N_ = N_BaseFunct[CurrentElement];
        /**
         The list of degrees of freedom consists of 3 parts:
         0,...,ActiveBounds-1: active nodes (no hanging, no Dirichlet)
         ActiveBounds,DirichletBounds-1: Hanging nodes
         DirichletBound,...end : the DOF with DIRICHLET condition
         */
        // add local matrix to global
        
        for(int m=0;m<N_;m++)
        {
            // active DOF
            if(DOF[m] < ActiveBound)
            {
                for(int k=0;k<N_;k++)
                    square_matrices[j]->add(DOF[m], DOF[k], Matrix[m][k]);
            }                                         // endif l
            else
            {
                // non-active but not Dirichlet (hanging)
                if( DOF[m] < DirichletBound)
                {
                    int l = DOF[m];
                    l -= ActiveBound;
                    int end = HangingRowPtr[l+1];
                    for(int n=HangingRowPtr[l];n<end;n++)
                    {
                        for(int k=0;k<N_;k++)
                        {
                            if(DOF[k] == HangingColInd[n])
                            {
                                CurrentHangingEntries[n] +=Matrix[m][k];
                                break;
                            }                                 // endif
                        }                                   // endfor k
                    }                                     // endfor n
                }
                else
                {
                    // Dirichlet node
                    if(TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE)
                    {
                        /*
                         When applying FEM-FCT, we treat Dirichlet nodes like active ones
                         */
                        for(int k=0;k<N_;k++)
                        {
                            square_matrices[j]->add(DOF[m], DOF[k], Matrix[m][k]);
                        }
                    } else { //Dirichlet standard treatment: 1.0 on the diagonal
                        square_matrices[j]->set(DOF[m], DOF[m], 1.0);
                    }
                }
            }
        }                                           // endfor m
    }                                             // endfor j
    // --------------------------------------------------------------------
    
    
    // --------------------------------------------------------------------
    // rectangular matrices
    for(int j=0 ; j< this->n_rectangular_matrices ; j++)
    {
        FE2D TestElement = ((TFESpace2D *) rectangular_matrices[j]->GetTestSpace())
        ->GetFE2D(i, cell);
        FE2D AnsatzElement = ((TFESpace2D *) rectangular_matrices[j]->GetAnsatzSpace())
        ->GetFE2D(i, cell);
        
        int N_Test = N_BaseFunct[TestElement];
        int N_Ansatz = N_BaseFunct[AnsatzElement];
        
        double **Matrix = Matrices+(j+this->n_square_matrices)*MaxN_BaseFunctions2D;
        
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
        for(int m=0;m<N_Test;m++)
        {
            // the entries are added only for active nodes and for Dirichlet nodes
            if( TestDOF[m] < ActiveBound  ||   TestDOF[m]>=DirichletBound)
            {
                int end=RowPtr[ TestDOF[m] + 1];
                for(int n=RowPtr[ TestDOF[m] ]; n<end; n++)
                {
                    for(int k=0; k<N_Ansatz; k++)
                    {
                        if(AnsatzDOF[k] == ColInd[n])
                        {
                            Entries[n] += Matrix[m][k];
                            break;
                        }                                   // endif
                    }                                     // endfor k
                }                                       // endfor n
            }
            else
            {
                // hanging node
                int l =  TestDOF[m]-ActiveBound;
                int end = HangingRowPtr[l+1];
                for(int n=HangingRowPtr[l];n<end;n++)
                {
                    for(int k=0;k<N_Ansatz;k++)
                    {
                        if(AnsatzDOF[k] == HangingColInd[n])
                        {
                            CurrentHangingEntries[n] += Matrix[m][k];
                            break;
                        }                                   // endif
                    }                                     // endfor k
                }                                       // endfor n
            }
        }                                           // endfor m
    }                                             // endfor j  (n_matrices)
}


//================================================================================
// add local right-hand sides to global right-hand side
//================================================================================
void Assembler4::add_local_to_global_rhs(int i,
                                     std::vector<const TFESpace2D*>& ferhs,
                                     double *righthand,
                                     const Example2D& example)
{
    TBaseCell *cell = this->Coll->GetCell(i);
    int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
    
    for(int j=0;j<rhs_blocks.size();j++)
    {
        const TFESpace2D *fespace = ferhs[j];
        int ActiveBound = fespace->GetActiveBound();
        FE2D CurrentElement = fespace->GetFE2D(i, cell);
        
        int N_ = N_BaseFunct[CurrentElement];
        
        double *local_rhs = righthand+j*MaxN_BaseFunctions2D;
        double *RHS = rhs_blocks[j];
        std::vector<double> CurrentHangingRhs = this->hangingRhs[j];
        // find space for this linear form
        
        ActiveBound = fespace->GetActiveBound();
        int DirichletBound = fespace->GetHangingBound();
        
        // dof of the rhs nodes connected to this cell
        int *DOF = ferhs[j]->GetGlobalDOF(i);
        //RhsGlobalNumbers[j] + RhsBeginIndex[j][i];
        
        // add local right-hand side to the global one
        for(int m=0;m<N_;m++)
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
            int l=DOF[m];
            if(l<ActiveBound)
            {
                // node l is inner node
                RHS[l] += local_rhs[m];
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
        // loop over all cell joints and impose corresponding BC on RHS
        // This version supports now only DIRICHLET BC type
        impose_boundary_conditions(fespace,example,cell,i,j,RHS,DOF);
    }                                             // endfor j (rhs_blocks.size())
}



//================================================================================
//================================================================================


void Assembler4::assemble_local_system(std::vector <const TFESpace2D*>& fespaces,
                                       int i,double ***LocMatrices,double **LocRhs)
{
    BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
    int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
    TBaseCell *cell = this->Coll->GetCell(i);

    double *aux;
    
    // --------------------------------------------------------------------
    // set the pointers of Parameters
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
    // note: if N_Parameters > 0, the default array is deleted
    int N_Parameters = la.GetN_Parameters();
    if(N_Parameters)
      {
	aux = new double [MaxN_QuadPoints_2D*N_Parameters];
	for(int j=0;j<MaxN_QuadPoints_2D;j++)
	  {
	    delete[] Param[j]; //clear away the default initialized array
	    Param[j] = aux + j*N_Parameters;
	  }
      }
    
    /*
    //HIER/////////////////////////////////////////////////////////////////////
    // --------------------------------------------------------------------
    // set the pointers of Parameters
    int N_Parameters = la.GetN_Parameters();
    std::vector<std::vector<double>> Parameters_new;
    Parameters_new.resize(MaxN_QuadPoints_2D);
    for (int t=0; t<MaxN_QuadPoints_2D; t++)
    {
        if(N_Parameters)
        {
	Parameters_new[t].resize(N_Parameters);
        }
        else
        {
        Parameters_new[t].resize(10 ,0.0);
        }
	}
    */
    
    // --------------------------------------------------------------------
    // 40 <= number of terms in bilinear form
    // Do not change below 20 since the entry 19 is used in GetLocalForms
    aux = new double [MaxN_QuadPoints_2D*40];
    double *AuxArray[MaxN_QuadPoints_2D];
    for(int j=0;j<MaxN_QuadPoints_2D;j++)
      AuxArray[j] = aux + j*40;
    
    /*
    //HIER/////////////////////////////////////////////////////////////////////
    // --------------------------------------------------------------------
    std::vector<std::vector<double>> AuxArray;
    AuxArray.resize(MaxN_QuadPoints_2D);
    for (int t=0; t<MaxN_QuadPoints_2D; t++)
    {
    AuxArray[t].resize(40);
    }
    
    */
    
    
    // --------------------------------------------------------------------
    // find local used elements on this cell
    
    FE2D LocalUsedElements[N_FEs2D];
    int LocN_BF[N_BaseFuncts2D];
    BaseFunct2D LocBF[N_BaseFuncts2D];
    int n_fespaces = fespaces.size();
    for(int j=0;j<n_fespaces;j++)
    {
        FE2D CurrentElement = fespaces[j]->GetFE2D(i,cell);
        LocalUsedElements[j] = CurrentElement;
        LocN_BF[j] = N_BaseFunct[CurrentElement];
        LocBF[j] = BaseFuncts[CurrentElement];
    }
    
    // --------------------------------------------------------------------
    // calculate values on original element
    
    int N_Points;
    double *weights, *xi, *eta;
    double X[MaxN_QuadPoints_2D],Y[MaxN_QuadPoints_2D];
    double AbsDetjk[MaxN_QuadPoints_2D];
    
    bool *SecondDer = la.GetNeeds2ndDerivatives();
    
    TFEDatabase2D::GetOrig(n_fespaces, LocalUsedElements,
                           this->Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // old version
    la.GetParameters(N_Points, this->Coll, cell, i, X, Y, Param);
    
    //HIER///////////////////////////////////////////////////////////////
    //la.compute_parameters(N_Points, this->Coll, i, X, Y, Parameters_new);
    
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

    // old version
    la.GetLocalForms(N_Points, weights, AbsDetjk, X, Y, LocN_BF, LocBF,
                     Param, AuxArray, cell, N_AllMatrices, n_rhs_blocks, LocMatrices,
                     LocRhs);

    
//Hier//////////////////////////////////////////////////////////////////////////////////
/*    la.get_local_forms(N_Points, weights, AbsDetjk, X, Y, LocN_BF, LocBF,
                       Parameters_new, AuxArray, cell, N_AllMatrices, n_rhs_blocks, LocMatrices,
                       LocRhs);
*/
}


//================================================================================
//================================================================================
Assembler4::~Assembler4(){};


//================================================================================
//================================================================================
void Assembler4::Assemble2D(BlockFEMatrix &M,
                            BlockVector &b_rhs,
                            std::vector<const TFESpace2D*>& fespaces,
                            std::vector<const TFESpace2D*>& ferhs,
                            const Example2D& example,
                            int AssemblePhaseID)
{
#ifdef __3D__
    ErrThrow("Assembler4::Assembler4() not yet working in 3D");
#endif
    
    // set matrices and rhs blocks
    this->init(M,b_rhs,fespaces,ferhs);
    
    ///@todo all these pointers are used but it might not be initialized
    double **Matrices;
    double ***LocMatrices;
    // LocRhs: an array of pointers (size: number of rhs)
    double **LocRhs;
    // righthand: a big pointer
    double *righthand;
    if(n_rhs_blocks)
    {
        LocRhs = new double* [n_rhs_blocks];
        righthand = new double [n_rhs_blocks*MaxN_BaseFunctions2D];
        for(int i=0;i<n_rhs_blocks;i++)
            LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
    }                                               // endif n_rhs_blocks
    

    
    // --------------------------------------------------------------------
    // set pointers to matrices
    double *aux;
    
    if(this->N_AllMatrices)
    {
        aux = new double
        [this->N_AllMatrices * MaxN_BaseFunctions2D * MaxN_BaseFunctions2D];
        Matrices = new double* [this->N_AllMatrices*MaxN_BaseFunctions2D];
        for(int j=0;j<this->N_AllMatrices*MaxN_BaseFunctions2D;j++)
            Matrices[j] = aux+j*MaxN_BaseFunctions2D;
        
        LocMatrices = new double** [this->N_AllMatrices];
        for(int i=0;i<this->N_AllMatrices;i++)
            LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
    }
    
    //================================================================================
    // loop over all cells
    //================================================================================
    
    for(int i=0 ; i < this->Coll->GetN_Cells() ; i++)
    {
        TBaseCell *cell = this->Coll->GetCell(i);
        
        ///@attention the value of INTERNAL_CELL appears to be used somewhere later
        TDatabase::ParamDB->INTERNAL_CELL = i;
        
        // only for multiphase flows
        if ((AssemblePhaseID >= 0) &&
            (AssemblePhaseID != cell->GetPhase_ID()) )
            continue;
        
        this->assemble_local_system(fespaces,
                                    i,LocMatrices, LocRhs);
        
        // add local/cellwise matrices to global matrices (ansatz == test: Aii, C; ansatz != test: A12,A12,B1,...)
        this->add_local_to_global_matrix(i, LocMatrices, Matrices);
        
        // add local/cellwise right-hand sides to global right-hand sides
        this->add_local_to_global_rhs(i,ferhs,righthand,example);
        
    }
    
    // --------------------------------------------------------------------
    // modify matrix according to coupling of hanging nodes
    // this part is only relevant to the case with hanging nodes and
    // it could be put in a separate function
    
    this->handle_hanging_nodes(ferhs);
    
    // --------------------------------------------------------------------
    // delete pointers
    if(n_rhs_blocks)
    {
        delete [] righthand;
        delete [] LocRhs;
    }
    
    /*
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
    */
    
    if(this->N_AllMatrices)
    {
        delete [] LocMatrices;
        delete [] Matrices[0];
        delete [] Matrices;
    }
    
    //delete [] AuxArray[0];
}                                                 // end of Assemble


