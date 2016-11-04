// ======================================================================
// @(#)BoundaryAssembling3D.C        28.10.16
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

#include <FEDatabase3D.h>
#include <TriaAffin.h>
#include <QuadAffin.h>
#include <BoundaryAssembling3D.h>
#include <Collection.h>
#include <BoundFace.h>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// (the edge list can be generated outside this class)
/*void BoundaryAssembling3D::rhs_g_v_n(BlockVector &rhs,
                                     const TFESpace3D *U_Space,
                                     BoundValueFunct3D *given_boundary_data,
                                     int boundary_component_id,
                                     double mult)
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    rhs_g_v_n(rhs,U_Space,given_boundary_data,boundaryEdgeList,mult);
}
*/

void BoundaryAssembling3D::rhs_g_v_n(BlockVector &rhs,
                                     const TFESpace3D *U_Space,
                                     BoundValueFunct3D *given_boundary_data,
                                     //---TBoundFace
                                     //std::vector<TBoundFace*> &boundaryFaceList,
                                     std::vector<TBaseCell*> &boundaryCells,
                                     double mult)
{
    // =========================================
    int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    TCollection* Coll =  U_Space->GetCollection();
    
    for(size_t i=0; i< boundaryCells.size(); i++) {
        
        TBaseCell* cell = boundaryCells[i];
        BaseFunct3D *BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
        // get basis dimension and FE space data of cell i
        FE3D FEId = U_Space->GetFE3D(i,cell );
        int fe_degree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);
        int* N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
        int nBasisFunctions = N_BaseFunct[FEId];
        for(size_t m=0;m< cell->GetN_Faces(); m++) {
            TJoint* joint = cell->GetJoint(m);
            
            if (joint->GetType() == BoundaryFace ||
                joint->GetType() == IsoBoundFace) {
                
                //----GIBT ES GETNEIGHBOUR????
                //TBaseCell *cell = boundface->GetNeighbour(0);
                
                //---joint_id is not available in 3d --> get_joint_values3d
                int joint_id_on_cell = m; //boundedge->get_index_in_neighbour(cell);
                const int *faceVertexMap, *faceVertexMapLength;
                int maxNVerticesPerFace;
                // get information of faces and local vertices
                cell->GetShapeDesc()->GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
                // number of vertices on face m
                size_t nFaceVertices = faceVertexMapLength[m];
                TBoundFace *boundface = (TBoundFace *)joint;
                // todo: use the class Point
                std::vector< std::vector<double> > faceVertices;
                faceVertices.resize(nFaceVertices);
                for (size_t l1=0; l1<nFaceVertices; l1++) {
                    faceVertices[l1].resize(3);
                    cell->GetVertex(faceVertexMap[m*maxNVerticesPerFace+l1])->
                    GetCoords(faceVertices[l1][0],faceVertices[l1][1],faceVertices[l1][2]);
                    
                }
                
                
                QuadFormula2D FaceQuadFormula; //=BaryCenterTria;
                switch(nFaceVertices) {
                    case 3:
                        // triangular face
                        FaceQuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*fe_degree);
                        FaceQuadFormula = Gauss3Tria;
                        break;
                    case 4:
                        // quadrilateral face
                        FaceQuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*fe_degree);
                        break;
                }
                
                // todo: use this code for setting quadrature formula
                //std::vector<double> quadWeights,quadPoints;
                //get_quadrature_formula_data(quadPoints,quadWeights);
                // compute values of all basis functions at all quadrature points
                //std::vector< std::vector<double> > uorig, u_dx_orig, u_dy_orig, u_dz_orig;
                //get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,u_dx_orig,u_dy_orig, ,u_dz_orig);

                int N_Points;
                double *t,*s;
                double* faceWeights;
                double* JointValue;
                double *RHS;
                // get a quadrature formula good enough for the velocity FE space
                TQuadFormula2D *qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
                qf2->GetFormulaData(N_Points, faceWeights, t, s);
                // generate data on reference mesh cell for the 2d face of 3d cell
                TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)
                ->MakeRefElementData(FaceQuadFormula);
                // values of base functions in all quadrature points on face
                double **JointValues = TFEDatabase3D::GetJointValues3D
                (BaseFuncts[FEId], FaceQuadFormula, m);
                TFEDatabase3D::GetBaseFunct3D(BaseFuncts[FEId])
                ->ChangeBF(Coll, cell, N_Points, JointValues);
                
                std::vector<double> normal;
                normal.resize(3);
                double transformationDeterminant;
                double xc1, yc1, zc1, xc2, yc2, zc2, xc3, yc3, zc3;
                switch(nFaceVertices) {
                    case 3:
                        // compute the 2 vectors that span the plane containing the current face
                        xc1 = faceVertices[1][0] - faceVertices[0][0];  // faceVertices[0][0]=X(0)
                        xc2 = faceVertices[2][0] - faceVertices[0][0];
                        
                        yc1 = faceVertices[1][1] - faceVertices[0][1]; // faceVertices[0][1]=Y(0)
                        yc2 = faceVertices[2][1] - faceVertices[0][1];
                        
                        zc1 = faceVertices[1][2] - faceVertices[0][2]; // faceVertices[0][2]=Z(0)
                        zc2 = faceVertices[2][2] - faceVertices[0][2];
                        
                        // normal vector on a plane is given by the cross product of the two spanning vectors (xc1, yc1, zc1) and (xc2, yc2, zc2) devided by the norm of this crossproduct
                        normal[0] = yc1*zc2 - zc1*yc2;
                        normal[1] = zc1*xc2 - xc1*zc2;
                        normal[2] = xc1*yc2 - yc1*xc2;
                        // determinant of reference trafo in order to get a normed normal vector
                        transformationDeterminant =
                        sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
                        normal[0] /= transformationDeterminant;
                        normal[1] /= transformationDeterminant;
                        normal[2] /= transformationDeterminant;
                        
                        
                for(size_t l=0;l<N_Points;l++)
                {
                    
                    /*double value;
                     if(given_boundary_data != nullptr)
                     {
                     given_boundary_data(BDComponent, T, value);
                     }
                     else
                     {
                     value = 1;
                     }*/
                    
                    // mapping from local (cell) DOF to global DOF
                    int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
                    
                    
                    // basis functions on Gauss point l
                    JointValue = JointValues[l];
                    
                    /*t0 = t[l];
                     t1 = s[l];
                     // cout << "t: " << t0 << " " << t1 << endl;
                     LinComb[0] = 1-t0-t1;
                     LinComb[1] = t0;
                     LinComb[2] = t1;
                     
                     xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                     +LinComb[2]*X[2];
                     yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                     +LinComb[2]*Y[2];
                     zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                     +LinComb[2]*Z[2];
                     BoundaryValue(xf, yf, zf, t0);*/
                    
                    /*
                     for(unsigned int l=0;l< uorig[k].size();l++)
                     {
                     int global_dof_from_local = DOF[l];
                     
                     // if the DOF is Dirichlet, continue
                     if(global_dof_from_local >= ActiveBound)
                     continue;
                     
                     // updating rhs: int_gamma rhsval v \cdot n
                     double v_x = uorig[k][l]; // value of test function (vtest = vx = vy =vz)
                     double v_y = v_x;
                     double v_z = v_x;
                     // add for both components
                     rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value * (v_x*n_x) *
                     (joint_length/reference_face_area);
                     rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value * (v_y*n_y) *
                     (joint_length/reference_face_area);
                     rhs.block(2)[global_dof_from_local] += mult * quadWeights[k] * value * (v_z*n_z) *
                     (joint_length/reference_face_area);
                     } //for(l=0;l<N_BaseFunct;l++)
                     */
                    
                    for(size_t k=0;k<nBasisFunctions;k++)
                        if(DOF[k]<ActiveBound)
                            RHS[ DOF[k] ] += faceWeights[l]*transformationDeterminant*mult*JointValue[k];
                } // endfor l
                        
                        break;
                        
                    case 4:
                        //Welche Einschränkung gibt es an die Vierecke?
                        xc1=(- faceVertices[0][0] + faceVertices[1][0] + faceVertices[2][0] - faceVertices[3][0] ) * 0.25;
                        xc2=(- faceVertices[0][0] - faceVertices[1][0] + faceVertices[2][0] + faceVertices[3][0] ) * 0.25;
                        xc3=(  faceVertices[0][0] - faceVertices[1][0] + faceVertices[2][0] - faceVertices[3][0] ) * 0.25;
                        
                        yc1=(- faceVertices[0][1] + faceVertices[1][1] + faceVertices[2][1] - faceVertices[3][1] ) * 0.25;
                        yc2=(- faceVertices[0][1] - faceVertices[1][1] + faceVertices[2][1] + faceVertices[3][1] ) * 0.25;
                        yc3=(  faceVertices[0][1] - faceVertices[1][1] + faceVertices[2][1] - faceVertices[3][1] ) * 0.25;
                        
                        zc1=(- faceVertices[0][2] + faceVertices[1][2] + faceVertices[2][2] - faceVertices[3][2] ) * 0.25;
                        zc2=(- faceVertices[0][2] - faceVertices[1][2] + faceVertices[2][2] + faceVertices[3][2] ) * 0.25;
                        zc3=(  faceVertices[0][2] - faceVertices[1][2] + faceVertices[2][2] - faceVertices[3][2] ) * 0.25;
                        
                        for(size_t l=0;l<N_Points;l++)
                        {
                            
                            // mapping from local (cell) DOF to global DOF
                            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
                            
                            
                            /*
                             JointValue = JointValues[l];
                            t0 = 0.5*(t[l]+1);
                            t1 = 0.5*(s[l]+1);
                            // cout << "t: " << t0 << " " << t1 << endl;
                            LinComb[0] = (1-t0)*(1-t1);
                            LinComb[1] = t0*(1-t1);
                            LinComb[2] = t0*t1;
                            LinComb[3] = (1-t0)*t1;
                            
                            xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                            +LinComb[2]*X[2] + LinComb[3]*X[3];
                            yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                            +LinComb[2]*Y[2] + LinComb[3]*Y[3];
                            zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                            +LinComb[2]*Z[2] + LinComb[3]*Z[3];
                            
                            // if(OuterBoundary)
                            // BoundComp->GetXYZandTS(TmpLen[m], LinComb,
                            // X, Y, Z, Param1, Param2,
                            // xf, yf, zf, t0, t1);
                            
                            // cout << xf << " " << yf << " " << zf << endl;
                            BoundaryValue(xf, yf, zf, t0);
                            // cout << "PV: " << t0 << endl;
                            nx = (yc1+s[l]*yc3)*(zc2+t[l]*zc3)
                            -(zc1+s[l]*zc3)*(yc2+t[l]*yc3);
                            ny = (zc1+s[l]*zc3)*(xc2+t[l]*xc3)
                            -(xc1+s[l]*xc3)*(zc2+t[l]*zc3);
                            nz = (xc1+s[l]*xc3)*(yc2+t[l]*yc3)
                            -(yc1+s[l]*yc3)*(xc2+t[l]*xc3);
                            t1 = nx*nx+ny*ny+nz*nz;
                            // cout << t1 << endl;
                            t0 *= faceWeights[l]*sqrt(t1);
                            */
                            for(size_t k=0;k<nBasisFunctions;k++)
                                if( DOF[k] < ActiveBound )
                                   RHS[ DOF[k] ] += JointValue[k]; //*t0
                        } // endfor l
                        break;
                }
                
            }
        } // endif
    }
}

////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//void BoundaryAssembling2D::rhs_g_v(BlockVector &rhs,
//                                   const TFESpace2D *U_Space,
//                                   BoundValueFunct2D *given_boundary_data1,
//                                   BoundValueFunct2D *given_boundary_data2,
//                                   int boundary_component_id,
//                                   double mult,
//                                   bool rescale_by_h)
//{
//    std::vector<TBoundEdge*> boundaryEdgeList;
//    TCollection *coll= U_Space->GetCollection();
//    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
//    rhs_g_v(rhs,U_Space,given_boundary_data1,given_boundary_data2,boundaryEdgeList,mult,rescale_by_h);
//}
//
//void BoundaryAssembling2D::rhs_g_v(BlockVector &rhs,
//                                   const TFESpace2D *U_Space,
//                                   BoundValueFunct2D *given_boundary_data1,
//                                   BoundValueFunct2D *given_boundary_data2,
//                                   std::vector<TBoundEdge*> &boundaryEdgeList,
//                                   double mult,
//                                   bool rescale_by_h)
//{
//    int *BeginIndex = U_Space->GetBeginIndex();
//    int *GlobalNumbers = U_Space->GetGlobalNumbers();
//    int ActiveBound = U_Space->GetActiveBound();
//    
//    for(size_t m=0;m< boundaryEdgeList.size(); m++)
//    {
//        TBoundEdge *boundedge = boundaryEdgeList[m];
//        TBaseCell *cell = boundedge->GetNeighbour(0);
//        
//        // get basis dimension and FE space data of cell i
//        FE2D FEId = U_Space->GetFE2D(0,cell );
//        
//        int BDComponent=boundedge->GetBoundComp()->GetID();
//        
//        int BaseVectDim = 1; // we assume only scalar FE // Only for BDM and RT elements \neq 1
//        int joint_id = boundedge->get_index_in_neighbour(cell);
//        
//        // get a quadrature formula good enough for the velocity FE space (here exact to 2*fe_degree)
//        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
//        std::vector<double> quadWeights,quadPoints;
//        get_quadrature_formula_data(quadPoints,quadWeights);
//        
//        // compute values of all basis functions at all quadrature points
//        std::vector< std::vector<double> > uorig,u_dx_orig,u_dy_orig;
//        get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,u_dx_orig,u_dy_orig);
//        
//        double x_0, x_1, y_0, y_1;
//        boundedge->get_vertices(x_0,  y_0, x_1, y_1);
//        // compute length of the edge
//        double joint_length = boundedge->get_length();
//        // normal vector to this boundary (normalized)
//        double n_x,n_y;
//        boundedge->get_normal(n_x, n_y);
//        
//        // quadrature
//        for(unsigned int k=0;k<quadPoints.size();k++)
//        {
//            ///@attention in 1D the reference joint is [-1,1] => length = 2
//            double reference_joint_length = 2;
//            double x = x_0+(quadPoints[k]+1.)/2.*(x_1-x_0);
//            double y = y_0+(quadPoints[k]+1.)/2.*(y_1-y_0);
//            
//            double T;
//            boundedge->GetBoundComp()->GetTofXY(x, y, T);
//            
//            // get the boundary values of rhs
//            double value1, value2;
//            if(given_boundary_data1 != nullptr)
//            {
//                given_boundary_data1(BDComponent, T, value1);
//            }
//            else
//            {
//                value1 = 1;
//            }
//            
//            if(given_boundary_data2 != nullptr)
//            {
//                given_boundary_data2(BDComponent, T, value2);
//            }
//            else
//            {
//                value2 = 0.0;
//            }
//            
//            // mapping from local (cell) DOF to global DOF
//            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
//            
//            for(unsigned int l=0;l< uorig[k].size();l++)
//            {
//                int global_dof_from_local = DOF[l];
//                
//                // if the DOF is Dirichlet, continue
//                if(global_dof_from_local >= ActiveBound)
//                    continue;
//                
//                // updating rhs: int_gamma rhsval[2] v
//                double v_x = uorig[k][l]; // value of test function (vtest = vx = vy)
//                double v_y=v_x;
//                
//                // add for both components
//                if (!rescale_by_h)
//                {
//                    rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value1 * (v_x) *
//                    (joint_length/reference_joint_length);
//                }
//                else
//                {
//                    rhs.block(0)[global_dof_from_local] += (mult * quadWeights[k] * value1 * (v_x) *
//                                                            (joint_length/reference_joint_length)) /joint_length;
//                }
//                if (!rescale_by_h)
//                {
//                    //rhs[0][global_dof_from_local+N_U] +=
//                    rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value2 * (v_y) *
//                    (joint_length/reference_joint_length);
//                }
//                else
//                {
//                    rhs.block(1)[global_dof_from_local] += (mult * quadWeights[k] * value2 * (v_y) *
//                                                            (joint_length/reference_joint_length)) /joint_length;
//                }
//                
//            }
//        }
//    }
//}
//
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//void BoundaryAssembling2D::matrix_v_n_v_n(BlockFEMatrix &M,
//                                          const TFESpace2D *U_Space,
//                                          int boundary_component_id,
//                                          double mult
//                                          )
//{
//    std::vector<TBoundEdge*> boundaryEdgeList;
//    TCollection *coll= U_Space->GetCollection();
//    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
//    
//    matrix_v_n_v_n(M,U_Space,boundaryEdgeList,mult);
//}
//
///**
// * @attention this functions assumes implicitely Matrix Type 14.
// * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
// * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
// * for A11, A12, A21, A22
// **/
//void BoundaryAssembling2D::matrix_v_n_v_n(BlockFEMatrix &M,
//                                          const TFESpace2D *U_Space,
//                                          std::vector<TBoundEdge*> &boundaryEdgeList,
//                                          double mult)
//{   int *BeginIndex = U_Space->GetBeginIndex();
//    int *GlobalNumbers = U_Space->GetGlobalNumbers();
//    int ActiveBound = U_Space->GetActiveBound();
//    
//    
//    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
//    /**
//     * @todo: check if the matrix structure is correct:
//     * we need 4 square matrices with the same FE spaces
//     */
//    
//    for(size_t m=0;m< boundaryEdgeList.size(); m++)
//    {
//        TBoundEdge *boundedge = boundaryEdgeList[m];
//        TBaseCell *cell = boundedge->GetNeighbour(0);
//        // get basis dimension and FE space data of cell i
//        FE2D FEId = U_Space->GetFE2D(0,cell );
//        
//        int BaseVectDim = 1; // we assume only scalar FE
//        int joint_id = boundedge->get_index_in_neighbour(cell);
//        
//        // get a quadrature formula good enough for the velocity FE space
//        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
//        std::vector<double> quadWeights,quadPoints;
//        get_quadrature_formula_data(quadPoints,quadWeights);
//        
//        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
//        
//        // compute values of all basis functions at all quadrature points
//        std::vector< std::vector<double> > uorig,uxorig,uyorig;
//        get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,uxorig,uyorig);
//        
//        double x0, x1, y0, y1;
//        boundedge->get_vertices(x0,  y0, x1, y1);
//        // compute length of the edge
//        double joint_length = boundedge->get_length();
//        // normal vector to this boundary (normalized)
//        double nx,ny;
//        boundedge->get_normal(nx, ny);
//        
//        // quadrature
//        for(unsigned int k=0;k<quadPoints.size();k++)
//        {
//            ///@attention in 1D the reference joint is [-1,1] => length = 2
//            double reference_joint_length = 2;
//            
//            // mapping from local(cell) DOF to global DOF
//            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
//            
//            // loop on test functions
//            double scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
//            for(unsigned int l1=0;l1< uorig[k].size();l1++)
//            {
//                int test_DOF = DOF[l1];
//                
//                // if the DOF is Dirichlet, continue
//                if(test_DOF >= ActiveBound)
//                    continue;
//                
//                double vx = uorig[k][l1];
//                double vy = vx; // x and y component have the same FE space
//                
//                // loop on ansatz functions
//                for(unsigned int l2=0; l2<uorig[k].size(); l2++)
//                {
//                    int ansatz_DOF = DOF[l2];
//                    double ux = uorig[k][l2];
//                    double uy = ux; // x and y component have the same FE space
//                    
//                    // (see the note about blocks at the beginning of the function)
//                    blocks[0]->add(test_DOF, ansatz_DOF, scale_factor*(vx*nx)*(ux*nx) ); // A11
//                    blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(vx*nx)*(uy*ny) ); // A12
//                    blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(vy*ny)*(ux*nx) ); // A21
//                    blocks[4]->add(test_DOF, ansatz_DOF, scale_factor*(vy*ny)*(uy*ny) ); // A22
//                }
//            } //for(l=0;l<N_BaseFunct;l++)
//        }
//    } // endif
//}
//
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//void BoundaryAssembling2D::matrix_gradv_n_v(BlockFEMatrix &M,
//                                            const TFESpace2D *U_Space,
//                                            int boundary_component_id,
//                                            double mult
//                                            )
//{
//    std::vector<TBoundEdge*> boundaryEdgeList;
//    TCollection *coll= U_Space->GetCollection();
//    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
//    matrix_gradv_n_v(M,U_Space,boundaryEdgeList,mult);
//}
//
///**
// * @attention this functions assumes implicitely Matrix Type 14.
// * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
// * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
// * for A11, A12, A21, A22
// **/
//void BoundaryAssembling2D::matrix_gradv_n_v(BlockFEMatrix &M,
//                                            const TFESpace2D *U_Space,
//                                            std::vector<TBoundEdge*> &boundaryEdgeList,
//                                            double mult)
//{   int *BeginIndex = U_Space->GetBeginIndex();
//    int *GlobalNumbers = U_Space->GetGlobalNumbers();
//    int ActiveBound = U_Space->GetActiveBound();
//    
//    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
//    /**
//     * @todo: check if the matrix structure is correct:
//     * we need 4 square matrices with the same FE spaces
//     */
//    
//    for(size_t m=0;m< boundaryEdgeList.size(); m++)
//    {
//        TBoundEdge *boundedge = boundaryEdgeList[m];
//        TBaseCell *cell = boundedge->GetNeighbour(0);
//        // get basis dimension and FE space data of cell i
//        FE2D FEId = U_Space->GetFE2D(0,cell );
//        
//        int BaseVectDim = 1; // we assume only scalar FE
//        int joint_id = boundedge->get_index_in_neighbour(cell);
//        
//        // get a quadrature formula good enough for the velocity FE space
//        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
//        std::vector<double> quadWeights,quadPoints;
//        get_quadrature_formula_data(quadPoints,quadWeights);
//        
//        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
//        
//        // compute values of all basis functions at all quadrature points
//        std::vector< std::vector<double> > uorig,uxorig,uyorig;
//        get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,uxorig,uyorig);
//        
//        double x0, x1, y0, y1;
//        boundedge->get_vertices(x0,  y0, x1, y1);
//        // compute length of the edge
//        double joint_length = boundedge->get_length();
//        // normal vector to this boundary (normalized)
//        double nx,ny;
//        boundedge->get_normal(nx, ny);
//        
//        // quadrature
//        for(unsigned int k=0;k<quadPoints.size();k++)
//        {
//            ///@attention in 1D the reference joint is [-1,1] => length = 2
//            double reference_joint_length = 2;
//            
//            // mapping from local(cell) DOF to global DOF
//            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
//            
//            // loop on test functions
//            double scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
//            for(unsigned int l1=0;l1< uorig[k].size();l1++)
//            {
//                int test_DOF = DOF[l1];
//                
//                // if the DOF is Dirichlet, continue
//                if(test_DOF >= ActiveBound)
//                    continue;
//                
//                double vx = uorig[k][l1];
//                double vy = vx; // x and y component have the same FE space
//                
//                // loop on ansatz functions
//                for(unsigned int l2=0; l2<uorig[k].size(); l2++)
//                {
//                    int ansatz_DOF = DOF[l2];
//                    
//                    double uxx = uxorig[k][l2];
//                    double uyy = uyorig[k][l2]; // x and y component have the same FE space
//                    
//                    // (see the note about blocks at the beginning of the function)
//                    blocks[0]->add(test_DOF, ansatz_DOF, -scale_factor*(vx)*uxx*nx ); // A11
//                    blocks[0]->add(test_DOF, ansatz_DOF, -scale_factor*(vx)*uyy*ny ); // A11
//                    //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(vx)*uyy*nx ); // A12
//                    //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(vx)*uxx*ny ); // A21
//                    blocks[4]->add(test_DOF, ansatz_DOF, -scale_factor*(vy)*uyy*ny ); // A22
//                    blocks[4]->add(test_DOF, ansatz_DOF, -scale_factor*(vy)*uxx*nx ); // A22
//                }
//            } //for(l=0;l<N_BaseFunct;l++)
//        }
//    } // endif
//}
//
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//void BoundaryAssembling2D::matrix_u_v(BlockFEMatrix &M,
//                                      const TFESpace2D *U_Space,
//                                      int boundary_component_id,
//                                      double mult,
//                                      bool rescale_by_h
//                                      )
//{
//    std::vector<TBoundEdge*> boundaryEdgeList;
//    TCollection *coll= U_Space->GetCollection();
//    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
//    matrix_u_v(M,U_Space,boundaryEdgeList,mult,rescale_by_h);
//}
//
///**
// * @attention this functions assumes implicitely Matrix Type 14.
// * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
// * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
// * for A11, A12, A21, A22
// **/
//void BoundaryAssembling2D::matrix_u_v(BlockFEMatrix &M,
//                                      const TFESpace2D *U_Space,
//                                      std::vector<TBoundEdge*> &boundaryEdgeList,
//                                      double mult,
//                                      bool rescale_by_h)
//{
//    int *BeginIndex = U_Space->GetBeginIndex();
//    int *GlobalNumbers = U_Space->GetGlobalNumbers();
//    int ActiveBound = U_Space->GetActiveBound();
//    
//    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
//    /**
//     * @todo: check if the matrix structure is correct:
//     * we need 4 square matrices with the same FE spaces
//     */
//    
//    for(size_t m=0;m< boundaryEdgeList.size(); m++)
//    {
//        TBoundEdge *boundedge = boundaryEdgeList[m];
//        TBaseCell *cell = boundedge->GetNeighbour(0);
//        // get basis dimension and FE space data of cell i
//        FE2D FEId = U_Space->GetFE2D(0,cell );
//        
//        int BaseVectDim = 1; // we assume only scalar FE
//        int joint_id = boundedge->get_index_in_neighbour(cell);
//        
//        // get a quadrature formula good enough for the velocity FE space
//        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
//        std::vector<double> quadWeights,quadPoints;
//        get_quadrature_formula_data(quadPoints,quadWeights);
//        
//        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
//        
//        // compute values of all basis functions at all quadrature points
//        std::vector< std::vector<double> > uorig,uxorig,uyorig;
//        get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,uxorig,uyorig);
//        
//        double x0, x1, y0, y1;
//        boundedge->get_vertices(x0,  y0, x1, y1);
//        // compute length of the edge
//        double joint_length = boundedge->get_length();
//        // normal vector to this boundary (normalized)
//        double nx,ny;
//        boundedge->get_normal(nx, ny);
//        
//        // quadrature
//        for(unsigned int k=0;k<quadPoints.size();k++)
//        {
//            ///@attention in 1D the reference joint is [-1,1] => length = 2
//            double reference_joint_length = 2;
//            
//            // mapping from local(cell) DOF to global DOF
//            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
//            
//            // rescale local integral (Nitsche)
//            double scale_factor;
//            if (rescale_by_h)
//            {
//                scale_factor = (mult*quadWeights[k]*(joint_length/reference_joint_length)) /joint_length;
//            }
//            else
//            {
//                scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
//            }
//            
//            // loop on test functions
//            for(unsigned int l1=0;l1< uorig[k].size();l1++)
//            {
//                int test_DOF = DOF[l1];
//                
//                // if the DOF is Dirichlet, continue
//                if(test_DOF >= ActiveBound)
//                    continue;
//                
//                double vx = uorig[k][l1];
//                double vy = vx; // x and y component have the same FE space
//                
//                // loop on ansatz functions
//                for(unsigned int l2=0; l2<uorig[k].size(); l2++)
//                {
//                    int ansatz_DOF = DOF[l2];
//                    double ux = uorig[k][l2];
//                    double uy = ux; // x and y component have the same FE space
//                    
//                    // (see the note about blocks at the beginning of the function)
//                    blocks[0]->add(test_DOF, ansatz_DOF, scale_factor*vx*ux ); // A11
//                    //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(vx*nx)*uyy ); // A12
//                    //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(vx*ny)*uxx ); // A21
//                    blocks[4]->add(test_DOF, ansatz_DOF, scale_factor*vy*uy ); // A22
//                }
//            } //for(l=0;l<N_BaseFunct;l++)
//        }
//    } // endif
//}
//
//
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//// (the edge list can be generated outside this class)
//void BoundaryAssembling2D::matrix_p_v_n(BlockFEMatrix &M,
//                                        const TFESpace2D *U_Space,
//                                        const TFESpace2D *P_Space,
//                                        int boundary_component_id,
//                                        double mult)
//{
//    std::vector<TBoundEdge*> boundaryEdgeList;
//    TCollection *coll= U_Space->GetCollection();
//    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
//    matrix_p_v_n(M,U_Space, P_Space,boundaryEdgeList,mult);
//}
//
//void BoundaryAssembling2D::matrix_p_v_n(BlockFEMatrix &M,
//                                        const TFESpace2D *U_Space,
//                                        const TFESpace2D *P_Space,
//                                        std::vector<TBoundEdge*> &boundaryEdgeList,
//                                        double mult)
//{
//    int ActiveBound = U_Space->GetActiveBound();
//    
//    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
//    /**
//     * @todo: check if the matrix structure is correct:
//     * we need 4 square matrices with the same FE spaces
//     */
//    
//    for(size_t m=0;m< boundaryEdgeList.size(); m++)
//    {
//        TBoundEdge *boundedge = boundaryEdgeList[m];
//        TBaseCell *cell = boundedge->GetNeighbour(0);
//        // get basis dimension and FE space data of cell i
//        FE2D FEId_U = U_Space->GetFE2D(0,cell );
//        // get basis dimension and FE space data of cell i
//        FE2D FEId_P = P_Space->GetFE2D(0,cell );
//        
//        int BaseVectDim = 1; // we assume only scalar FE
//        int joint_id = boundedge->get_index_in_neighbour(cell);
//        
//        // get a quadrature formula good enough for the velocity FE space
//        int fe_degree_U = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_U);
//        int fe_degree_P = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_P);
//        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(fe_degree_P*fe_degree_U);
//        std::vector<double> quadWeights,quadPoints;
//        get_quadrature_formula_data(quadPoints,quadWeights);
//        
//        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
//        
//        // compute values of all basis functions at all quadrature points
//        std::vector< std::vector<double> > uorig,uxorig,uyorig;
//        get_original_values(FEId_U, joint_id, cell, quadPoints, BaseVectDim, uorig,uxorig,uyorig);
//        
//        int BaseVectDim_P = 1; // we assume only scalar FE; nur bei Raviart-Thomas & BDM \neq 1
//        
//        // compute values of all basis functions at all quadrature points
//        std::vector< std::vector<double>> porig, pxorig, pyorig;
//        get_original_values(FEId_P, joint_id, cell, quadPoints, BaseVectDim_P, porig, pxorig, pyorig);
//        
//        double x0, x1, y0, y1;
//        boundedge->get_vertices(x0,  y0, x1, y1);
//        // compute length of the edge
//        double joint_length = boundedge->get_length();
//        // normal vector to this boundary (normalized)
//        double nx,ny;
//        boundedge->get_normal(nx, ny);
//        
//        // mapping from local(cell) DOF to global DOF
//        // int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
//        int *DOF_P=P_Space->GetGlobalDOF(cell->GetCellIndex());
//        int *DOF_U=U_Space->GetGlobalDOF(cell->GetCellIndex());
//        
//        // quadrature
//        for(unsigned int k=0;k<quadPoints.size();k++)
//        {
//            ///@attention in 1D the reference joint is [-1,1] => length = 2
//            double reference_joint_length = 2;
//            
//            // loop on test functions
//            double scale_factor = mult*quadWeights[k]*(joint_length/reference_joint_length);
//            for( unsigned int l1=0 ; l1 < uorig[k].size() ; l1++ )
//            {
//                int test_DOF = DOF_U[l1];
//                
//                // if the DOF is Dirichlet, continue
//                if(test_DOF >= ActiveBound)
//                    continue;
//                
//                double vx = uorig[k][l1];
//                double vy = vx; // x and y component have the same FE space
//                
//                // loop on ansatz functions
//                for(unsigned int l2=0; l2<porig[k].size(); l2++)
//                {
//                    int ansatz_DOF = DOF_P[l2];
//                    double p = porig[k][l2];
//                    
//                    // (see the note about blocks at the beginning of the function)
//                    blocks[2]->add(test_DOF, ansatz_DOF, scale_factor*p*vx*nx ); // B1
//                    //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(vx*nx)*uyy ); // A12
//                    //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(vx*ny)*uxx ); // A21
//                    blocks[5]->add(test_DOF, ansatz_DOF, scale_factor*p*vy*ny ); // B2
//                } //for(l=0;l<N_BaseFunct;l++)
//            }
//        } // endif
//    }
//    
//}
//
//
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//void BoundaryAssembling2D::get_quadrature_formula_data(std::vector<double> &P,
//                                                       std::vector<double> &W)
//{
//    // get the type of required quadrature (include/FE/Enumerations.h)
//    // initialize points and weights of quadrature
//    ///@attention LineQuadFormula must be set before calling the function GetQuadFormulaData
//    TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(this->LineQuadFormula);
//    
//    int nQuadPoints;
//    double *quadWeights, *quadPoints;
//    qf1->GetFormulaData(nQuadPoints, quadWeights, quadPoints);
//    P.resize(nQuadPoints);
//    W.resize(nQuadPoints);
//    
//    for (int i=0; i<nQuadPoints; i++){
//        P[i]=quadPoints[i];
//        W[i]=quadWeights[i];
//    }
//}
//
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//void BoundaryAssembling2D::get_original_values(FE2D FEId, int joint_id,
//                                               TBaseCell *cell,
//                                               std::vector<double> quadPoints,
//                                               int BaseVectDim,
//                                               std::vector< std::vector<double> > &originalValues,
//                                               std::vector< std::vector<double> > &originalValues_x,
//                                               std::vector< std::vector<double> > &originalValues_y)
//{
//    BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
//    BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
//    int *N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
//    int N_BaseFunct = N_BaseFuncts[FEId];
//    
//    TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
//    
//    originalValues.resize(quadPoints.size());
//    originalValues_x.resize(quadPoints.size());
//    originalValues_y.resize(quadPoints.size());
//    
//    for( unsigned int k=0;k<quadPoints.size();k++ )
//    {
//        // ----------------------------
//        // for each quadrature point, get values and derivatives of basis functions
//        // according to selected quadrature formula
//        ///@todo write something like:
//        // GetQuadratureValues(..., uorig, uxorig, uyorig)?
//        // needs:
//        // BaseFuncts[FEId],LineQuadFormula,joint_id, cell
//        // RefElement, quadPoints[k], N_BaseFunct
//        double **reference_values = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
//                                                                    LineQuadFormula,
//                                                                    joint_id);
//        double **derivative_xi = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
//                                                                      LineQuadFormula,
//                                                                      joint_id, D10);
//        double **derivative_eta = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
//                                                                       LineQuadFormula,
//                                                                       joint_id, D01);
//        
//        double uorig[N_BaseFunct], uxorig[N_BaseFunct], uyorig[N_BaseFunct];
//        RefTrans2D RefTrans;
//        TRefTrans2D *F_K;
//        switch(RefElement)
//        {
//            case BFUnitTriangle:
//                RefTrans = TriaAffin;
//                F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
//                ((TTriaAffin*)F_K)->SetCell(cell);
//                ((TTriaAffin*)F_K)->GetOrigValues(joint_id,quadPoints[k],N_BaseFunct,
//                                                  reference_values[k],
//                                                  derivative_xi[k], derivative_eta[k],
//                                                  uorig, uxorig, uyorig,
//                                                  BaseVectDim);
//                break;
//            case BFUnitSquare:
//                RefTrans = QuadAffin;
//                F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
//                ((TQuadAffin*)F_K)->SetCell(cell);
//                ((TQuadAffin*)F_K)->GetOrigValues(joint_id, quadPoints[k],N_BaseFunct,
//                                                  reference_values[k], derivative_xi[k],
//                                                  derivative_eta[k],
//                                                  uorig, uxorig, uyorig,
//                                                  BaseVectDim );
//                break;
//            default:
//                Output::print("Unknown reference element in BoundaryAssembling2D.C");
//                exit(0);
//                break;
//        }
//        originalValues[k].resize(N_BaseFunct);
//        originalValues_x[k].resize(N_BaseFunct);
//        originalValues_y[k].resize(N_BaseFunct);
//        for (int ib=0; ib<N_BaseFunct; ib++)
//        {
//            originalValues[k][ib] = uorig[ib];
//            originalValues_x[k][ib] = uxorig[ib];
//            originalValues_y[k][ib] = uyorig[ib];
//        }
//    }
//}
//
