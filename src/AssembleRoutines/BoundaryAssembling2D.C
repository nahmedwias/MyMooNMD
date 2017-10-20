// ======================================================================
// @(#)BoundaryAssembling2D.C        05/18/16
//
// Functions for (external and internal) boundary integral
//
// ======================================================================
/*
 List of possible improvements
 * link each cell with its FE (in order not to use a global function to get FE properties)
 * create an FE class, in order not to recover all properties via FEDatabase::FunctioName(FeID)
 * rewrite the GetFormulaData using (e.g.) vector<> class
 * avoid to define double* with MaxN... (use vector instead?)
 * alternative: pass a list of joints (each with its own FE space on it)?
 */

#include <FEDatabase2D.h>
#include <TriaAffin.h>
#include <QuadAffin.h>
#include <BoundaryAssembling2D.h>
#include <BoundEdge.h>
#include <Collection.h>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ToDo: the edge list could be generated outside this class
void BoundaryAssembling2D::rhs_g_v_n(BlockVector &rhs,
                                     const TFESpace2D *U_Space,
                                     BoundValueFunct2D *given_boundary_data,
                                     int boundary_component_id,
                                     double mult)
{
    // Create a list of those boundary edges that are on the boundary component with given ID
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    
    rhs_g_v_n(rhs,U_Space, given_boundary_data, boundaryEdgeList, mult);
}

void BoundaryAssembling2D::rhs_g_v_n(BlockVector &rhs,
                                     const TFESpace2D *U_Space,
                                     BoundValueFunct2D *given_boundary_data,
                                     std::vector<TBoundEdge*> &boundaryEdgeList,
                                     double mult)
{
    // =========================================
    int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    
    // Go through all boundary edges on the current boundary component
    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
        
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        // get basis dimension and FE space data of the current cell
        FE2D FEId = U_Space->GetFE2D(0, cell);
        
        int BaseVectDim = 1; // we assume only scalar FE
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the argument of the integral
        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
        std::vector<double> quadWeights, quadPoints;
        get_quadrature_formula_data(quadPoints, quadWeights);
        
        // compute values of all basis functions and their first partial derivatives at all quadrature points
        std::vector< std::vector<double> > uorig, u_dx_orig ,u_dy_orig;
        get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, u_dx_orig, u_dy_orig);
        
        double x_0, x_1, y_0, y_1;
        boundedge->get_vertices(x_0, y_0, x_1, y_1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double n1, n2;
        boundedge->get_normal(n1, n2);
        
        // quadrature
        for(unsigned int k = 0; k < quadPoints.size(); k++)
        {
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            double x = x_0+(quadPoints[k]+1.)/2.*(x_1-x_0);
            double y = y_0+(quadPoints[k]+1.)/2.*(y_1-y_0);
            
            // Parametrization
            double T;
            boundedge->GetBoundComp()->GetTofXY(x, y, T);
            
            int BDComponent = boundedge->GetBoundComp()->GetID();
            
            // get the value of (\nabla u - pI) on the boundary component (here denoted by g) from the example, if possible
            double value;
            if(given_boundary_data != nullptr)
            {
                given_boundary_data(BDComponent, T, value);
            }
            else
            {
                value = 1;
            }
            
            // mapping from local (cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            for(unsigned int l = 0; l < uorig[k].size(); l++)
            {
                int global_dof_from_local = DOF[l];
                
                // if the DOF is Dirichlet, continue
                if(global_dof_from_local >= ActiveBound)
                    continue;
                
                // updating rhs: int_gamma rhsval v \cdot n
                double v1 = uorig[k][l]; // value of test function (vtest = vx = vy)
                double v2 = v1;
                // add for both components
                rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value * (v1 * n1) *
                (joint_length/reference_joint_length);
                rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value * (v2 * n2) *
                (joint_length/reference_joint_length);
            } //for(l=0;l<N_BaseFunct;l++)
        }
    } // endif
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::rhs_uD_v(BlockVector &rhs,
                                   const TFESpace2D *U_Space,
                                   BoundValueFunct2D *given_boundary_data1,
                                   BoundValueFunct2D *given_boundary_data2,
                                   int boundary_component_id,
                                   double mult,
                                   bool rescale_by_h)
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    rhs_uD_v(rhs,U_Space,given_boundary_data1,given_boundary_data2,boundaryEdgeList,mult,rescale_by_h);
}

void BoundaryAssembling2D::rhs_uD_v(BlockVector &rhs,
                                   const TFESpace2D *U_Space,
                                   BoundValueFunct2D *given_boundary_data1,
                                   BoundValueFunct2D *given_boundary_data2,
                                   std::vector<TBoundEdge*> &boundaryEdgeList,
                                   double mult,
                                   bool rescale_by_h)
{
    int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    
    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        
        // get basis dimension and FE space data of cell i
        FE2D FEId = U_Space->GetFE2D(0, cell);
        
        int BDComponent=boundedge->GetBoundComp()->GetID();
        
        int BaseVectDim = 1; // we assume only scalar FE // Only for BDM and RT elements \neq 1
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space (here exact to 2*fe_degree)
        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
        std::vector<double> quadWeights, quadPoints;
        get_quadrature_formula_data(quadPoints, quadWeights);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > uorig, u_dx_orig, u_dy_orig;
        get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, u_dx_orig, u_dy_orig);
        
        double x_0, x_1, y_0, y_1;
        boundedge->get_vertices(x_0, y_0, x_1, y_1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        
        // quadrature
        for(unsigned int k = 0; k < quadPoints.size(); k++)
        {
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            double x = x_0+(quadPoints[k]+1.)/2.*(x_1-x_0);
            double y = y_0+(quadPoints[k]+1.)/2.*(y_1-y_0);
            
            double T;
            boundedge->GetBoundComp()->GetTofXY(x, y, T);
            
            // get the boundary values of rhs
            double value1, value2;
            if(given_boundary_data1 != nullptr)
            {
                given_boundary_data1(BDComponent, T, value1);
            }
            else
            {
                value1 = 1.;
            }
            
            if(given_boundary_data2 != nullptr)
            {
                given_boundary_data2(BDComponent, T, value2);
            }
            else
            {
                value2 = 0.0;
            }
            


            // mapping from local (cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            for(unsigned int l = 0; l < uorig[k].size(); l++)
            {
                int global_dof_from_local = DOF[l];
                
                // if the DOF is Dirichlet, continue
                if(global_dof_from_local >= ActiveBound)
                    continue;
                
                // updating rhs: int_gamma rhsval[2] v
                double v1 = uorig[k][l]; // value of test function (vtest = vx = vy)
                double v2 = v1;
                
                // add for both components
                if (!rescale_by_h)
                {
                    rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value1 * v1 *
                    (joint_length/reference_joint_length);
                }
                else
                {
                    rhs.block(0)[global_dof_from_local] += (mult * quadWeights[k] * value1 * v1 *
                                                            (joint_length/reference_joint_length)) / joint_length;
                }
                if (!rescale_by_h)
                {
                    //rhs[0][global_dof_from_local+N_U] +=
                    rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value2 * v2 *
                    (joint_length/reference_joint_length);
                }
                else
                {
                    rhs.block(1)[global_dof_from_local] += (mult * quadWeights[k] * value2 * v2 *
                                                            (joint_length/reference_joint_length)) /joint_length;
                }
                
            }
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::rhs_gradv_n_uD(BlockVector &rhs,
                                   const TFESpace2D *U_Space,
                                   BoundValueFunct2D *given_boundary_data1,
                                   BoundValueFunct2D *given_boundary_data2,
                                   int boundary_component_id,
                                   double mult)
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    rhs_gradv_n_uD(rhs,U_Space,given_boundary_data1,given_boundary_data2,boundaryEdgeList,mult);
}

void BoundaryAssembling2D::rhs_gradv_n_uD(BlockVector &rhs,
                                   const TFESpace2D *U_Space,
                                   BoundValueFunct2D *given_boundary_data1,
                                   BoundValueFunct2D *given_boundary_data2,
                                   std::vector<TBoundEdge*> &boundaryEdgeList,
                                   double mult)
{
    int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    
    for(size_t m=0;m< boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        
        // get basis dimension and FE space data of cell i
        FE2D FEId = U_Space->GetFE2D(0,cell );
        
        int BDComponent=boundedge->GetBoundComp()->GetID();
        
        int BaseVectDim = 1; // we assume only scalar FE // Only for BDM and RT elements \neq 1
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space (here exact to 2*fe_degree)
        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
        std::vector<double> quadWeights,quadPoints;
        get_quadrature_formula_data(quadPoints,quadWeights);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > uorig,u_dx_orig,u_dy_orig;
        get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,u_dx_orig,u_dy_orig);
        
        double x_0, x_1, y_0, y_1;
        boundedge->get_vertices(x_0, y_0, x_1, y_1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double n1, n2;
        boundedge->get_normal(n1, n2);
        
        // quadrature
        for(unsigned int k=0;k<quadPoints.size();k++)
        {
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            double x = x_0+(quadPoints[k]+1.)/2.*(x_1-x_0);
            double y = y_0+(quadPoints[k]+1.)/2.*(y_1-y_0);
            
            double T;
            boundedge->GetBoundComp()->GetTofXY(x, y, T);
            
            // get the boundary values of rhs
            double value1, value2;
            if(given_boundary_data1 != nullptr)
            {
                given_boundary_data1(BDComponent, T, value1);
            }
            else
            {
                value1 = 1;
            }
            
            if(given_boundary_data2 != nullptr)
            {
                given_boundary_data2(BDComponent, T, value2);
            }
            else
            {
                value2 = 0.0;
            }
            
            
            // mapping from local (cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            for(unsigned int l = 0; l < uorig[k].size(); l++)
            {
                int global_dof_from_local = DOF[l];
                
                // if the DOF is Dirichlet, continue
                if(global_dof_from_local >= ActiveBound)
                    continue;
                
                // updating rhs: int_gamma rhsval[2] v
                double v1_dx = u_dx_orig[k][l]; // value of test function (vtest = (vx,vy) )
                double v1_dy = u_dy_orig[k][l];
                double v2_dy = v1_dy;
                double v2_dx = v1_dx;
                
//Falsch/alt?
//                    rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value1 * (vx_dx) * nx *
//                    (joint_length/reference_joint_length);
                rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value1 * (v1_dx * n1 + v1_dy * n2)*
                (joint_length/reference_joint_length);
                
                    //rhs[0][global_dof_from_local+N_U] +=
//Falsch/alt?
//                    rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value2 * (vy_dy) * ny *
//                    (joint_length/reference_joint_length);
                rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value2 * (v2_dx * n1 + v2_dy * n2) *
                (joint_length/reference_joint_length);
            }
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Löschen oder korrigieren
void BoundaryAssembling2D::matrix_u_n_v_n(BlockFEMatrix &M,
                                          const TFESpace2D *U_Space,
                                          int boundary_component_id,
                                          double mult
                                          )
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    
    matrix_u_n_v_n(M,U_Space,boundaryEdgeList,mult);
}

/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_u_n_v_n(BlockFEMatrix &M,
                                          const TFESpace2D *U_Space,
                                          std::vector<TBoundEdge*> &boundaryEdgeList,
                                          double mult)
{   int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
    /**
     * @todo: check if the matrix structure is correct:
     * we need 4 square matrices with the same FE spaces
     */
    
    for(size_t m=0;m< boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        // get basis dimension and FE space data of cell i
        FE2D FEId = U_Space->GetFE2D(0,cell );
        
        int BaseVectDim = 1; // we assume only scalar FE
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space
        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
        std::vector<double> quadWeights,quadPoints;
        get_quadrature_formula_data(quadPoints,quadWeights);
        
        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > uorig,uxorig,uyorig;
        get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,uxorig,uyorig);
        
        double x0, x1, y0, y1;
        boundedge->get_vertices(x0,  y0, x1, y1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double n1, n2;
        boundedge->get_normal(n1, n2);
        
        // quadrature
        for(unsigned int k=0;k<quadPoints.size();k++)
        {
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            
            // mapping from local(cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            // loop on test functions
            double scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
            for(unsigned int l1=0;l1< uorig[k].size();l1++)
            {
                int test_DOF = DOF[l1];
                
                // if the DOF is Dirichlet, continue
                if(test_DOF >= ActiveBound)
                    continue;
                
                double v1 = uorig[k][l1];
                double v2 = v1; // x and y component have the same FE space
                
                // loop on ansatz functions
                for(unsigned int l2=0; l2<uorig[k].size(); l2++)
                {
                    int ansatz_DOF = DOF[l2];
                    double u1 = uorig[k][l2];
                    double u2 = u1; // x and y component have the same FE space
                    
                    // (see the note about blocks at the beginning of the function)
                    blocks[0]->add(test_DOF, ansatz_DOF, scale_factor*(v1*n1)*(u1*n1) ); // A11
                    blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(v1*n1)*(u2*n2) ); // A12
                    blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(v2*n2)*(u1*n1) ); // A21
                    blocks[4]->add(test_DOF, ansatz_DOF, scale_factor*(v2*n2)*(u2*n2) ); // A22
                }
            } //for(l=0;l<N_BaseFunct;l++)
        }
    } // endif
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void BoundaryAssembling2D::matrix_gradu_n_v(BlockFEMatrix &M,
                                            const TFESpace2D *U_Space,
                                            int boundary_component_id,
                                            double mult
                                            )
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    matrix_gradu_n_v(M,U_Space,boundaryEdgeList,mult,boundary_component_id);
}

/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_gradu_n_v(BlockFEMatrix &M,
                                            const TFESpace2D *U_Space,
                                            std::vector<TBoundEdge*> &boundaryEdgeList,
                                            double mult,
                                            int boundary_component_id)
{   int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
    /**
     * @todo: check if the matrix structure is correct:
     * we need 4 square matrices with the same FE spaces
     */

////LB DEBUG start
//    Output::print("!!!!!!!!!!boundary_component_id: ", boundary_component_id);
//    if (boundary_component_id = 1)
//        blocks[0]->reset();
//    
//    blocks[0]->write("A11_before_gradunv.out");
////    std::vector<int> test_DOF_vec;
////    std::vector<int> ansatz_DOF_vec;
////    std::vector<std::vector<double>> A11_gradunv;
//    Output::print("boundaryEdgeList.size(): ",boundaryEdgeList.size());
////LB DEBUG end
    
    for(size_t m=0;m< boundaryEdgeList.size(); m++)
    {
        //LB Debug start
        //Output::print("boundary edge: ",m);
        //LB DEBUG end
        
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        // get basis dimension and FE space data of cell i
        FE2D FEId = U_Space->GetFE2D(0,cell );
        
        int BaseVectDim = 1; // we assume only scalar FE
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space
        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
        
        //LB DEBUG start
        //Output::print("fe_degree: ",fe_degree);
        //LB DEBUG end
        
        std::vector<double> quadWeights,quadPoints;
        get_quadrature_formula_data(quadPoints,quadWeights);
        
        //LB Debug start
        //Output::print("num_quadPoints_per_edge: ",quadPoints.size());
        //LB Debug end
        
        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > uorig,uxorig,uyorig;
        get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig);
        
//        //LB DEBUG start
//        Output::print("uorig[0].size(): ",uorig[0].size());
//        Output::print("uorig[1].size(): ",uorig[1].size());
//        Output::print("uorig[2].size(): ",uorig[2].size());
//        Output::print("uorig[3].size(): ",uorig[3].size());
//        
//        Output::print("uorig[0][0]: ",uorig[0][0]);
//        Output::print("uorig[0][1]: ",uorig[0][1]);
//        Output::print("uorig[0][2]: ",uorig[0][2]);
//        Output::print("uorig[0][3]: ",uorig[0][3]);
//        Output::print("uorig[0][4]: ",uorig[0][4]);
//        Output::print("uorig[0][5]: ",uorig[0][5]);
//        
//        Output::print("quadPoints[0]: ",quadPoints[0]);
//        Output::print("quadPoints[1]: ",quadPoints[1]);
//        Output::print("quadPoints[2]: ",quadPoints[2]);
//        //LB DEBUG end
        
        double x0, x1, y0, y1;
        boundedge->get_vertices(x0,  y0, x1, y1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double n1,n2;
        boundedge->get_normal(n1, n2);
        
        //LB Debug start
        //Output::print("Values on edge: ", m );
        //LB Debug end
        
        // quadrature
        for(unsigned int k=0;k<quadPoints.size();k++)
        {
            //LB Debug start
            //Output::print("quadpoint: ",k);
            //LB Debug end
            
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            
            // mapping from local(cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            // loop on test functions
            double scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
            for(unsigned int l1 = 0; l1 < uorig[k].size(); l1++)
            {
                int test_DOF = DOF[l1];
                
                // if the DOF is Dirichlet, continue
                if(test_DOF >= ActiveBound)
                    continue;
                
                double v1 = uorig[k][l1];
                double v2 = v1; // x and y component have the same FE space
                
                // loop on ansatz functions
                for(unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
                {
                    int ansatz_DOF = DOF[l2];
                    
                    double u1_dx = uxorig[k][l2];
                    double u1_dy = uyorig[k][l2];
                    double u2_dx = uxorig[k][l2]; // u1 and u2 component have the same FE space
                    double u2_dy = uyorig[k][l2]; // u1 and u2 component have the same FE space
                    
                    // (see the note about blocks at the beginning of the function)
                    blocks[0]->add(test_DOF, ansatz_DOF, scale_factor * v1 * u1_dx * n1 ); // A11
                    blocks[0]->add(test_DOF, ansatz_DOF, scale_factor * v1 * u1_dy * n2 ); // A11
                    //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(vx)*uyy*nx ); // A12
                    //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(vx)*uxx*ny ); // A21
                    blocks[4]->add(test_DOF, ansatz_DOF, scale_factor * v2 * u2_dy * n2 ); // A22
                    blocks[4]->add(test_DOF, ansatz_DOF, scale_factor * v2 * u2_dx * n1 ); // A22
                    
//LB DEBUG start
//
//                    Output::print("has value:",scale_factor * v1 * u1_dx * n1 + scale_factor * v1 * u1_dy * n2 );
//                    Output::print(test_DOF);
//                    test_DOF_vec[l1] = test_DOF;
//                    
//                     Output::print("HALLOOOO");
//                    ansatz_DOF_vec[l2] = ansatz_DOF;
//                     Output::print("HALLOOOO");
//                    A11_gradunv[l1][l2]=scale_factor * v1 * u1_dx * n1 + scale_factor * v1 * u1_dy * n2;
//LB DEBUG end
                    
                }
            } //for(l=0;l<N_BaseFunct;l++)
        }
    } // endif
    
////LB DEBUG start
//    for (int kbl=1; kbl<10; kbl++)
//    {
//    Output::print("test_DOF_vec",test_DOF_vec[kbl]);
//    }
//    std::ofstream matrixfile("A11values_gradunv.txt");
//    matrixfile << test_DOF_vec << " " << ansatz_DOF_vec <<  " " << A11_gradunv<<  endl;
//    matrixfile.close();
//    
//    blocks[0]->write("A11_gradunv.out");
//    
////LB DEBUG end
}




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::matrix_gradv_n_u(BlockFEMatrix &M,
                                            const TFESpace2D *U_Space,
                                            int boundary_component_id,
                                            double mult
                                            )
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    matrix_gradv_n_u(M,U_Space,boundaryEdgeList,mult);
}

/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_gradv_n_u(BlockFEMatrix &M,
                                            const TFESpace2D *U_Space,
                                            std::vector<TBoundEdge*> &boundaryEdgeList,
                                            double mult)
{   int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
    /**
     * @todo: check if the matrix structure is correct:
     * we need 4 square matrices with the same FE spaces
     */
    
    for(size_t m=0;m< boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        // get basis dimension and FE space data of cell i
        FE2D FEId = U_Space->GetFE2D(0,cell );
        
        int BaseVectDim = 1; // we assume only scalar FE
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space
        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
        std::vector<double> quadWeights, quadPoints;
        get_quadrature_formula_data(quadPoints, quadWeights);
        
        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > uorig,uxorig,uyorig;
        get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,uxorig,uyorig);
        
        double x0, x1, y0, y1;
        boundedge->get_vertices(x0, y0, x1, y1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double n1, n2;
        boundedge->get_normal(n1, n2);
        
        // quadrature
        for(unsigned int k = 0; k < quadPoints.size(); k++)
        {
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            
            // mapping from local(cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            // loop on test functions
            double scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
            for(unsigned int l1=0;l1< uorig[k].size();l1++)
            {
                int test_DOF = DOF[l1];
                
                // if the DOF is Dirichlet, continue
                if(test_DOF >= ActiveBound)
                    continue;
                
                
                double v1_dx = uxorig[k][l1];
                double v1_dy = uyorig[k][l1]; // x and y component have the same FE space
                double v2_dx = uxorig[k][l1];
                double v2_dy = uyorig[k][l1]; // x and y component have the same FE space
                
                // loop on ansatz functions
                for(unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
                {
                    int ansatz_DOF = DOF[l2];
                    
                   //falsch/alt?
                   // double u1 = uorig[k][l1];
                    double u1 = uorig[k][l2];
                    double u2 = u1; // x and y component have the same FE space
                    
                    // (see the note about blocks at the beginning of the function)
                    blocks[0]->add(test_DOF, ansatz_DOF, scale_factor * u1 * v1_dx * n1 ); // A11
                    blocks[0]->add(test_DOF, ansatz_DOF, scale_factor * u1 * v1_dy * n2 ); // A11
                    //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(ux)*vy_dy*nx ); // A12
                    //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(ux)*vx_dx*ny ); // A21
                    blocks[4]->add(test_DOF, ansatz_DOF, scale_factor * u2 * v2_dy * n2 ); // A22
                    blocks[4]->add(test_DOF, ansatz_DOF, scale_factor * u2 * v2_dx * n1 ); // A22
                }
            } //for(l=0;l<N_BaseFunct;l++)
        }
    } // endif
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::matrix_u_v(BlockFEMatrix &M,
                                      const TFESpace2D *U_Space,
                                      int boundary_component_id,
                                      double mult,
                                      bool rescale_by_h
                                      )
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    matrix_u_v(M,U_Space,boundaryEdgeList,mult,rescale_by_h);
}

/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_u_v(BlockFEMatrix &M,
                                      const TFESpace2D *U_Space,
                                      std::vector<TBoundEdge*> &boundaryEdgeList,
                                      double mult,
                                      bool rescale_by_h)
{
    int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
    /**
     * @todo: check if the matrix structure is correct:
     * we need 4 square matrices with the same FE spaces
     */
    
    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        // get basis dimension and FE space data of cell i
        FE2D FEId = U_Space->GetFE2D(0,cell);
        
        int BaseVectDim = 1; // we assume only scalar FE
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space
        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
        std::vector<double> quadWeights,quadPoints;
        get_quadrature_formula_data(quadPoints,quadWeights);
        
        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > uorig,uxorig,uyorig;
        get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig);
        
        double x0, x1, y0, y1;
        boundedge->get_vertices(x0, y0, x1, y1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        
        // quadrature
        for(unsigned int k = 0; k < quadPoints.size(); k++)
        {
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            
            // mapping from local(cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            // rescale local integral (Nitsche)
            double scale_factor;
            if (rescale_by_h)
            {
                scale_factor = (mult*quadWeights[k]*(joint_length/reference_joint_length)) /joint_length;
            }
            else
            {
                scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
            }
            
            // loop on test functions
            for(unsigned int l1=0;l1< uorig[k].size();l1++)
            {
                int test_DOF = DOF[l1];
                
                // if the DOF is Dirichlet, continue
                if(test_DOF >= ActiveBound)
                    continue;
                
                double v1 = uorig[k][l1];
                double v2 = v1; // x and y component have the same FE space
                
                // loop on ansatz functions
                for(unsigned int l2=0; l2<uorig[k].size(); l2++)
                {
                    int ansatz_DOF = DOF[l2];
                    double u1 = uorig[k][l2];
                    double u2 = u1; // x and y component have the same FE space
                    
                    // (see the note about blocks at the beginning of the function)
                    blocks[0]->add(test_DOF, ansatz_DOF, scale_factor * v1 * u1 ); // A11
                    //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(vx*nx)*uyy ); // A12
                    //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(vx*ny)*uxx ); // A21
                    blocks[4]->add(test_DOF, ansatz_DOF, scale_factor * v2 * u2 ); // A22
                }
            } //for(l=0;l<N_BaseFunct;l++)
        }
    } // endif
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// (the edge list can be generated outside this class)
void BoundaryAssembling2D::matrix_q_u_n(BlockFEMatrix &M,
                                        const TFESpace2D *U_Space,
                                        const TFESpace2D *P_Space,
                                        int boundary_component_id,
                                        double mult)
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    matrix_q_u_n(M,U_Space, P_Space,boundaryEdgeList,mult);
}

void BoundaryAssembling2D::matrix_q_u_n(BlockFEMatrix &M,
                                        const TFESpace2D *U_Space,
                                        const TFESpace2D *P_Space,
                                        std::vector<TBoundEdge*> &boundaryEdgeList,
                                        double mult)
{
    int ActiveBound = U_Space->GetActiveBound();
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
    /**
     * @todo: check if the matrix structure is correct:
     * we need 4 square matrices with the same FE spaces
     */
    
    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        // get basis dimension and FE space data of cell i
        FE2D FEId_U = U_Space->GetFE2D(0,cell);
        // get basis dimension and FE space data of cell i
        FE2D FEId_P = P_Space->GetFE2D(0,cell);
        
        int BaseVectDim = 1; // we assume only scalar FE
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space
        int fe_degree_U = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_U);
        int fe_degree_P = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_P);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(fe_degree_P+fe_degree_U);
        std::vector<double> quadWeights,quadPoints;
        get_quadrature_formula_data(quadPoints,quadWeights);
        
        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
        
        // compute values of all basis functions at all quadrature points
        std::vector<std::vector<double>> uorig, uxorig, uyorig;
        get_original_values(FEId_U, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig);
        
        int BaseVectDim_P = 1; // we assume only scalar FE; nur bei Raviart-Thomas & BDM \neq 1
       
        // compute values of all basis functions at all quadrature points
        std::vector<std::vector<double>> porig, pxorig, pyorig;
        get_original_values(FEId_P, joint_id, cell, quadPoints, BaseVectDim_P, porig, pxorig, pyorig);
        
        double x0, x1, y0, y1;
        boundedge->get_vertices(x0, y0, x1, y1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double n1, n2;
        boundedge->get_normal(n1, n2);
        
        // mapping from local(cell) DOF to global DOF
        // int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
        int *DOF_P=P_Space->GetGlobalDOF(cell->GetCellIndex());
        int *DOF_U=U_Space->GetGlobalDOF(cell->GetCellIndex());
        
        // quadrature
        for( unsigned int k = 0; k < quadPoints.size(); k++ )
        {
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            
            // loop on test functions
            double scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
            for( unsigned int l1 = 0; l1 < uorig[k].size(); l1++ )
            {
                int ansatz_DOF = DOF_U[l1];
                
                double u1 = uorig[k][l1];
                double u2 = u1; // x and y component have the same FE space

                
                // loop on ansatz functions
                for(unsigned int l2 = 0; l2 < porig[k].size(); l2++)
                {
                    int test_DOF = DOF_P[l2];
                    // if the DOF is Dirichlet, continue
                    if(test_DOF >= ActiveBound)
                        continue;
                    double q = porig[k][l2];
                                        //Output::print("HIIIEEERRRRRRR: ", quadPoints.size());
                    // (see the note about blocks at the beginning of the function)
                    blocks[6]->add(test_DOF, ansatz_DOF, scale_factor * q * u1 * n1 ); // B1
//                    Output::print("q * u1 * n1: ",q * u1 * n1);
//                    Output::print("scale_factor: ",scale_factor );
                    //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(vx*nx)*uyy ); // A12
                    //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(vx*ny)*uxx ); // A21
                    blocks[7]->add(test_DOF, ansatz_DOF, scale_factor * q * u2 * n2 ); // B2
//                    Output::print("q * u2 * n2: ",q * u2 * n2);
//                    Output::print("scale_factor: ",scale_factor );
//                    Output::print("q: ",q );
//                    Output::print("u1: ",u1 );
//                    Output::print("u2: ",u2 );
//                    Output::print("n1: ",n1 );
//                    Output::print("n2: ",n2 );
//                    Output::print("qun: scale_factor * q * u1 * n1 ",scale_factor * q * u1 * n1);
//                    Output::print("qun: scale_factor * q * u2 * n2", scale_factor * q * u2 * n2);

                } //for(l=0;l<N_BaseFunct;l++)
                
//                Output::print("HIIIEEERRRRRRR l1: ", l1);
//                int test_DOF = DOF_P[l1];
//                // if the DOF is Dirichlet, continue
//                if(test_DOF >= ActiveBound)
//                    continue;
//                double q = porig[k][l1];
//                
//                // loop on ansatz functions
//                for(unsigned int l2=0; l2<porig[k].size(); l2++)
//                {Output::print("HIIIEEERRRRRRR l2: ", l2);
//                    int ansatz_DOF = DOF_U[l2];
//                    
//                    double u1 = uorig[k][l2];
//                    double u2 = u1; // x and y component have the same FE space
//                    //Output::print("HIIIEEERRRRRRR: ", quadPoints.size());
//                    // (see the note about blocks at the beginning of the function)
//                    blocks[6]->add(test_DOF, ansatz_DOF, scale_factor*q*u1*nx ); // B1'
//                    //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(vx*nx)*uyy ); // A12
//                    //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(vx*ny)*uxx ); // A21
//                    blocks[7]->add(test_DOF, ansatz_DOF, scale_factor*q*u2*ny ); // B2'
//                } //for(l=0;l<N_BaseFunct;l++)

            }
        } // endif
    }
    
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// (the edge list can be generated outside this class)
void BoundaryAssembling2D::rhs_q_uD_n(BlockVector &rhs,
                                     const TFESpace2D *U_Space,
                                     const TFESpace2D *P_Space,
                                     BoundValueFunct2D *given_boundary_data1,
                                     BoundValueFunct2D *given_boundary_data2,
                                     int boundary_component_id,
                                     double mult)
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= P_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    
    rhs_q_uD_n(rhs, U_Space, P_Space, given_boundary_data1, given_boundary_data2, boundaryEdgeList, mult);
}

void BoundaryAssembling2D::rhs_q_uD_n(BlockVector &rhs,
                                     const TFESpace2D *U_Space,
                                     const TFESpace2D *P_Space,
                                     BoundValueFunct2D *given_boundary_data1,
                                     BoundValueFunct2D *given_boundary_data2,
                                     std::vector<TBoundEdge*> &boundaryEdgeList,
                                     double mult)
{
    // =========================================
    int *BeginIndex = P_Space->GetBeginIndex();
    int *GlobalNumbers = P_Space->GetGlobalNumbers();
    int ActiveBound = P_Space->GetActiveBound();
    
    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        
        int BaseVectDim_P = 1; // we assume only scalar FE
        int joint_id = boundedge->get_index_in_neighbour(cell);
        FE2D FEId_U = U_Space->GetFE2D(0,cell);
        // get basis dimension and FE space data of cell i
        FE2D FEId_P = P_Space->GetFE2D(0,cell);
        
        
        // get a quadrature formula good enough for the velocity FE space
        int fe_degree_U = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_U);
        int fe_degree_P = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_P);

        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(fe_degree_U+fe_degree_P);
        std::vector<double> quadWeights,quadPoints;
        get_quadrature_formula_data(quadPoints,quadWeights);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > porig,pxorig,pyorig;
        get_original_values(FEId_P, joint_id, cell, quadPoints, BaseVectDim_P, porig, pxorig, pyorig);
        
        double x_0, x_1, y_0, y_1;
        boundedge->get_vertices(x_0, y_0, x_1, y_1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double n1,n2; // n=(n_1, n_2)^T
        boundedge->get_normal(n1, n2);
        
        // quadrature
        for(unsigned int k = 0; k < quadPoints.size(); k++)
        {
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            double x = x_0 + (quadPoints[k]+1.)/2.*(x_1-x_0);
            double y = y_0 + (quadPoints[k]+1.)/2.*(y_1-y_0);
            
            double T;
            boundedge->GetBoundComp()->GetTofXY(x, y, T);
            
            int BDComponent=boundedge->GetBoundComp()->GetID();
            
            // get the boundary values of rhs
            double value1, value2;
            if(given_boundary_data1 != nullptr)
            {
                given_boundary_data1(BDComponent, T, value1);
            }
            else
            {
                value1 = 1.;
            }
            
            if(given_boundary_data2 != nullptr)
            {
                given_boundary_data2(BDComponent, T, value2);
            }
            else
            {
                value2 = 1.;//0.0;
            }
            
            // mapping from local (cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            for(unsigned int l = 0; l < porig[k].size(); l++)
            {
                int global_dof_from_local = DOF[l];
                
                // if the DOF is Dirichlet, continue
                if(global_dof_from_local >= ActiveBound)
                    continue;
                
                // updating rhs: int_gamma rhsval v \cdot n
                double q = porig[k][l]; // value of test function (vtest = vx = vy)
                // add for both components
                //Falsch/alt?
                //                rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * q * (value1 * n1) *
                //                (joint_length/reference_joint_length);
                rhs.block(2)[global_dof_from_local] += mult * quadWeights[k] * q * (value1 * n1) *
                (joint_length/reference_joint_length);
                //Falsch/alt?
                //                rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * q * (value2 * n2) *
                //                (joint_length/reference_joint_length);
                rhs.block(2)[global_dof_from_local] += mult * quadWeights[k] * q * (value2 * n2) *
                (joint_length/reference_joint_length);
            } //for(l=0;l<N_BaseFunct;l++)
        }
    } // endif
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// (the edge list can be generated outside this class)
void BoundaryAssembling2D::matrix_p_v_n(BlockFEMatrix &M,
                                        const TFESpace2D *U_Space,
                                        const TFESpace2D *P_Space,
                                        int boundary_component_id,
                                        double mult)
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    matrix_p_v_n(M,U_Space, P_Space,boundaryEdgeList,mult);
}

void BoundaryAssembling2D::matrix_p_v_n(BlockFEMatrix &M,
                                        const TFESpace2D *U_Space,
                                        const TFESpace2D *P_Space,
                                        std::vector<TBoundEdge*> &boundaryEdgeList,
                                        double mult)
{
    int ActiveBound = U_Space->GetActiveBound();
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
    /**
     * @todo: check if the matrix structure is correct:
     * we need 4 square matrices with the same FE spaces
     */
    
    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        // get basis dimension and FE space data of cell i
        FE2D FEId_U = U_Space->GetFE2D(0, cell);
        // get basis dimension and FE space data of cell i
        FE2D FEId_P = P_Space->GetFE2D(0, cell);
        
        int BaseVectDim = 1; // we assume only scalar FE
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space
        int fe_degree_U = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_U);
        int fe_degree_P = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_P);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(fe_degree_P*fe_degree_U);
        std::vector<double> quadWeights, quadPoints;
        get_quadrature_formula_data(quadPoints, quadWeights);
        
        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > uorig, uxorig, uyorig;
        get_original_values(FEId_U, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig);
        
        int BaseVectDim_P = 1; // we assume only scalar FE; nur bei Raviart-Thomas & BDM \neq 1
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double>> porig, pxorig, pyorig;
        get_original_values(FEId_P, joint_id, cell, quadPoints, BaseVectDim_P, porig, pxorig, pyorig);
        
        double x0, x1, y0, y1;
        boundedge->get_vertices(x0, y0, x1, y1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double n1, n2;
        boundedge->get_normal(n1, n2);
        
        // mapping from local(cell) DOF to global DOF
        // int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
        int *DOF_P=P_Space->GetGlobalDOF(cell->GetCellIndex());
        int *DOF_U=U_Space->GetGlobalDOF(cell->GetCellIndex());
        
        // quadrature
        for(unsigned int k=0;k<quadPoints.size();k++)
        {
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            
            // loop on test functions
            double scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
            for(unsigned int l1 = 0 ; l1 < uorig[k].size() ; l1++)
            {
                int test_DOF = DOF_U[l1];
                
                // if the DOF is Dirichlet, continue
                if(test_DOF >= ActiveBound)
                    continue;
                
                double v1 = uorig[k][l1];
                double v2 = v1; // x and y component have the same FE space
                
                // loop on ansatz functions
                for(unsigned int l2=0; l2 < porig[k].size(); l2++)
                {
                    int ansatz_DOF = DOF_P[l2];
                    double p = porig[k][l2];
                    
                    // (see the note about blocks at the beginning of the function)
                    blocks[2]->add(test_DOF, ansatz_DOF, scale_factor * p * v1 * n1 ); // B1
                    //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(v1*nx)*uyy ); // A12
                    //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(v1*ny)*uxx ); // A21
                    blocks[5]->add(test_DOF, ansatz_DOF, scale_factor * p * v2 * n2 ); // B2
//                    Output::print("pv n: scale_factor * p * v1 * n1 ",scale_factor * p * v1 * n1 );
//                    Output::print("pv n: scale_factor * p * v2 * n2", scale_factor * p * v2 * n2);
                } //for(l=0;l<N_BaseFunct;l++)
            }
        } // endif
    }
    
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::get_quadrature_formula_data(std::vector<double> &P,
                                                       std::vector<double> &W)
{
    // get the type of required quadrature (include/FE/Enumerations.h)
    // initialize points and weights of quadrature
    ///@attention LineQuadFormula must be set before calling the function GetQuadFormulaData
    TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(this->LineQuadFormula);
    
    int nQuadPoints;
    double *quadWeights, *quadPoints;
    qf1->GetFormulaData(nQuadPoints, quadWeights, quadPoints);
    P.resize(nQuadPoints);
    W.resize(nQuadPoints);
    
    for (int i=0; i<nQuadPoints; i++){
        P[i]=quadPoints[i];
        W[i]=quadWeights[i];
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void BoundaryAssembling2D::get_original_values(FE2D FEId, int joint_id,
                                               TBaseCell *cell,
                                               std::vector<double> quadPoints,
                                               int BaseVectDim,
                                               std::vector< std::vector<double> > &originalValues,
                                               std::vector< std::vector<double> > &originalValues_x,
                                               std::vector< std::vector<double> > &originalValues_y)
{
    BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
    BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    int *N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
    int N_BaseFunct = N_BaseFuncts[FEId];
    
    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
    
    originalValues.resize(quadPoints.size());
    originalValues_x.resize(quadPoints.size());
    originalValues_y.resize(quadPoints.size());
    
    for( unsigned int k=0;k<quadPoints.size();k++ )
    {
        // ----------------------------
        // for each quadrature point, get values and derivatives of basis functions
        // according to selected quadrature formula
        ///@todo write something like:
        // GetQuadratureValues(..., uorig, uxorig, uyorig)?
        // needs:
        // BaseFuncts[FEId],LineQuadFormula,joint_id, cell
        // RefElement, quadPoints[k], N_BaseFunct
        double **reference_values = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                                                                    LineQuadFormula,
                                                                    joint_id);
        double **derivative_xi = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                                                                      LineQuadFormula,
                                                                      joint_id, D10);
        double **derivative_eta = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                                                                       LineQuadFormula,
                                                                       joint_id, D01);
        
        double uorig[N_BaseFunct], uxorig[N_BaseFunct], uyorig[N_BaseFunct];
        RefTrans2D RefTrans;
        TRefTrans2D *F_K;
        switch(RefElement)
        {
            case BFUnitTriangle:
                RefTrans = TriaAffin;
                F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
                ((TTriaAffin*)F_K)->SetCell(cell);
                ((TTriaAffin*)F_K)->GetOrigValues(joint_id,quadPoints[k],N_BaseFunct,
                                                  reference_values[k],
                                                  derivative_xi[k], derivative_eta[k],
                                                  uorig, uxorig, uyorig,
                                                  BaseVectDim);
                break;
            case BFUnitSquare:
                RefTrans = QuadAffin;
                F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
                ((TQuadAffin*)F_K)->SetCell(cell);
                ((TQuadAffin*)F_K)->GetOrigValues(joint_id, quadPoints[k],N_BaseFunct,
                                                  reference_values[k], derivative_xi[k],
                                                  derivative_eta[k],
                                                  uorig, uxorig, uyorig,
                                                  BaseVectDim );
                break;
            default:
                Output::print("Unknown reference element in BoundaryAssembling2D.C");
                exit(0);
                break;
        }
        originalValues[k].resize(N_BaseFunct);
        originalValues_x[k].resize(N_BaseFunct);
        originalValues_y[k].resize(N_BaseFunct);
        for (int ib=0; ib<N_BaseFunct; ib++)
        {
            originalValues[k][ib] = uorig[ib];
            originalValues_x[k][ib] = uxorig[ib];
            originalValues_y[k][ib] = uyorig[ib];
        }
    }
}

