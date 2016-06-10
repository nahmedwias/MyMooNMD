// ======================================================================
// @(#)BoundaryAssembling2D.C        05/18/16
//
// Functions for (external and internal) boundary integral
//
// ======================================================================
/**
 List of possible improvements
 * link each cell with its FE (in order not to use a global function to get FE properties)
 * create a FE class, in order not to recover all properties via FEDatabase::FunctioName(FeID)
 * rewrite the GetFormulaData using (e.g.) vector<> class
 * avoid to define double* with MaxN... (use vector instead?)
 * alternative: pass a list of joints (each with its own FE space on it)?
 */

#include <FEDatabase2D.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>
#include <TriaAffin.h>
#include <QuadAffin.h>
#include <BoundaryAssembling2D.h>
#include <BoundEdge.h>
#include <Collection.h>
#include <Database.h>




void BoundaryAssembling2D::BoundaryAssemble_on_rhs_g_v_n(double **rhs,
                                                         const TFESpace2D *U_Space,
                                                         TFEFunction2D *given_data,
                                                         int boundary_component_id,
                                                         double mult
                                                         )
{

    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    
    BoundaryAssemble_on_rhs_g_v_n(rhs,U_Space,given_data,boundaryEdgeList,mult);
    
}


void BoundaryAssembling2D::BoundaryAssemble_on_rhs_g_v_n(double **rhs,
                                                         const TFESpace2D *U_Space,
                                                         TFEFunction2D *given_data,
                                                         std::vector<TBoundEdge*> &boundaryEdgeList,
                                                         double mult)
{   // =========================================
    int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    
    for(int m=0;m< boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        // get basis dimension and FE space data of cell i
        FE2D FEId = U_Space->GetFE2D(0,cell );
        
        int BaseVectDim = 1; // we assume only scalar FE
        // ---------------------------------------------------------------
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space
        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
        std::vector<double> quadWeights,quadPoints;
        get_quadrature_formula_data(fe_degree,quadPoints,quadWeights);
        
        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > uorig,uxorig,uyorig;
        get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,uxorig,uyorig);
        
        double x0, x1, y0, y1;
        boundedge->get_vertices(x0,  y0, x1, y1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double nx,ny;
        boundedge->get_normal(nx, ny);
        
        // -------------------------------------------
        // quadrature
        for(unsigned int k=0;k<quadPoints.size();k++)
        {
            
            //double uorig[N_BaseFunct], uxorig[N_BaseFunct];
            //double uyorig[N_BaseFunct];
            //get_original_values(FEId, joint_id, cell, quadPoints[k],BaseVectDim, uorig,uxorig,uyorig);
            // ----------------------------
            
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            double x = x0+(quadPoints[k]+1.)/2.*(x1-x0);
            double y = y0+(quadPoints[k]+1.)/2.*(y1-y0);
            // given_data_value[0] = G(x,y),
            // given_data_value[1] = dG/dx, given_data_value[2] = dG/dy
            double given_data_value[3];
            // get the value of rhs
            if(given_data)
                given_data->FindGradientLocal(cell, cell->GetCellIndex(), x, y, given_data_value);
            else
            {
                given_data_value[0] = 1.0;
                given_data_value[1] = 0.0;
                given_data_value[2] = 0.0;
            }
            // mapping from local(cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            for(unsigned int l=0;l< uorig[k].size();l++)
            {
                int global_dof_from_local = DOF[l];
                
                // if the DOF is Dirichlet, continue
                if(global_dof_from_local >= ActiveBound)
                    continue;
                
                // updating rhs: int_gamma rhsval v \cdot n
                double vtest = uorig[k][l]; // value of test function (vtest = vx = vy)
                // add for both components
                rhs[0][global_dof_from_local] +=
                mult * quadWeights[k] * given_data_value[0] * (vtest*nx) *
                (joint_length/reference_joint_length);
                //rhs[0][global_dof_from_local+N_U] +=
                rhs[1][global_dof_from_local] +=
                mult * quadWeights[k] * given_data_value[0] * (vtest*ny) *
                (joint_length/reference_joint_length);
            } //for(l=0;l<N_BaseFunct;l++) 
        }
    } // endif
}
// ---------------------------------------------------------------


//%// Mod for Joints instead of cells END



//    // loop over all cells
//    for(int i=0;i< U_Space->GetCollection()->GetN_Cells(); i++)
//    {
//
//        // get cell from Joint (not from i)
//        TBaseCell *cell = U_Space->GetColction()->GetCell(i);
//
//        // get basis dimension and FE space data of cell i
//        FE2D FEId = U_Space->GetFE2D(0, cell);
//
//        int BaseVectDim = 1; // we assume only scalar FE
//
//        // get number of basis functions
//        int *N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
//        int N_BaseFunct = N_BaseFuncts[FEId];
//
//        int nEdges = cell->GetN_Edges();
//        // ---------------------------------------------------------------
//        for(int joint_id=0; joint_id<nEdges; joint_id++)
//        {
//            if(!(cell->GetJoint(joint_id)->InnerJoint()))
//            {
//                TJoint *joint = cell->GetJoint(joint_id);
//                TBoundEdge *boundedge = (TBoundEdge *)joint;
//                TBoundComp *BoundComp = boundedge->GetBoundComp();
//                if (BoundComp->GetID() == boundary_component_id) {
//                    
//                    
//                    // get a quadrature formula good enough for the velocity FE space
//                    int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//                    this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
//                    std::vector<double> quadWeights,quadPoints;
//                    get_quadrature_formula_data(fe_degree,quadPoints,quadWeights);
//                    
//                    
//                    //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
//                    
//                    // compute values of all basis functions at all quadrature points
//                    std::vector< std::vector<double> > uorig,uxorig,uyorig;
//                    get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,uxorig,uyorig);
//                    
//                    double x0, x1, y0, y1;
//                    boundedge->get_vertices(x0,  y0, x1, y1);
//                    // compute length of the edge
//                    double joint_length = boundedge->get_length();
//                    // normal vector to this boundary (normalized)
//                    double nx,ny;
//                    boundedge->get_normal(nx, ny);
//                    
//
//                    // -------------------------------------------
//                    // quadrature
//                    for(unsigned int k=0;k<quadPoints.size();k++)
//                    {
//                        
//                        //double uorig[N_BaseFunct], uxorig[N_BaseFunct];
//                        //double uyorig[N_BaseFunct];
//                        //get_original_values(FEId, joint_id, cell, quadPoints[k],BaseVectDim, uorig,uxorig,uyorig);
//                        
//                        
//                        // ----------------------------
//                        
//                        ///@attention in 1D the reference joint is [-1,1] => length = 2
//                        double reference_joint_length = 2;
//                        double x = x0+(quadPoints[k]+1.)/2.*(x1-x0);
//                        double y = y0+(quadPoints[k]+1.)/2.*(y1-y0);
//                        // given_data_value[0] = G(x,y),
//                        // given_data_value[1] = dG/dx, given_data_value[2] = dG/dy
//                        double given_data_value[3];
//                        // get the value of rhs
//                        if(given_data)
//                            given_data->FindGradientLocal(cell,i,x, y, given_data_value);
//                        else
//                        {
//                            given_data_value[0] = 1.0;
//                            given_data_value[1] = 0.0;
//                            given_data_value[2] = 0.0;
//                        }
//                        // mapping from local(cell) DOF to global DOF
//                        int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
//                        
//                        for(unsigned int l=0;l< uorig[k].size();l++)
//                        {
//                            int global_dof_from_local = DOF[l];
//                            
//                            // if the DOF is Dirichlet, continue
//                            if(global_dof_from_local >= ActiveBound)
//                                continue;
//                            
//                            // updating rhs: int_gamma rhsval v \cdot n
//                            double vtest = uorig[k][l]; // value of test function (vtest = vx = vy)
//                            // add for both components
//                            rhs[0][global_dof_from_local] +=
//                            mult * quadWeights[k] * given_data_value[0] * (vtest*nx) *
//                            (joint_length/reference_joint_length);
//                            //rhs[0][global_dof_from_local+N_U] +=
//                            rhs[1][global_dof_from_local] +=
//                            mult * quadWeights[k] * given_data_value[0] * (vtest*ny) *
//                            (joint_length/reference_joint_length);
//                            
//                            
//                        } //for(l=0;l<N_BaseFunct;l++)
//                    }
//                    
//                }
//            } // endif
//        }//  for(j=0;j<N_Edges;j++)
//        // ---------------------------------------------------------------
//    } //  for(k=0;k<N_LinePoints;k++
//} //  for(i=0;i<N_Cells;i++)



void BoundaryAssembling2D::get_quadrature_formula_data(int fe_degree,
					      std::vector<double> &P,
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
  
  for (unsigned int i=0; i<nQuadPoints; i++){
    P[i]=quadPoints[i];
    W[i]=quadWeights[i];
  }
}


void BoundaryAssembling2D::get_original_values(FE2D FEId, int joint_id, TBaseCell *cell, std::vector<double> quadPoints, int BaseVectDim,
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
    
    for(unsigned int k=0;k<quadPoints.size();k++)
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
                
        } // endswitch
        
        originalValues[k].resize(N_BaseFunct);
        originalValues_x[k].resize(N_BaseFunct);
        originalValues_y[k].resize(N_BaseFunct);
        for (int ib=0; ib<N_BaseFunct; ib++) {
            originalValues[k][ib] = uorig[ib];
            originalValues_x[k][ib] = uxorig[ib];
            originalValues_y[k][ib] = uyorig[ib];
        }
        
    }
}

// quadrature over the joints on the selected boundary
//    for(unsigned int j=0;j<jointsOnComponent.size();j++) {

//      int joint_id = jointsOnComponent[j];






//      // get geometrical properties of the boundary joint
//      ///@todo GetJointData(Joint)
//      double x0,x1,y0,y1;
//      double nx,ny;
//      double joint_length;
//
//      x0 = cell->GetVertex(joint_id)->GetX();
//      x1 = cell->GetVertex((joint_id+1)%nEdges)->GetX();
//      y0 = cell->GetVertex(joint_id)->GetY();
//      y1 = cell->GetVertex((joint_id+1)%nEdges)->GetY();
//
//      // compute length of the edge
//      joint_length = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
//
//      // normal vector to this boundary (normalized)
//      nx = (y1-y0)/joint_length;
//      ny = (x0-x1)/joint_length;

//      // ---------------------------------------------------------------
//      ///@todo GetFormulaData(FEId,quadWeights, quadPoints)
//      // get a quadrature formula good enough for the velocity FE space
//      int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//      // get the type of required quadrature (include/FE/Enumerations.h)
//      QuadFormula1D LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
//      // initialize points and weights of quadrature
//      TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
//
//      ///@todo rewrite the GetFormulaData using (e.g.) vector<> class
//      int nQuadPoints;
//      double *quadWeights, *quadPoints;
//      qf1->GetFormulaData(nQuadPoints, quadWeights, quadPoints);
//      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
//      // ---------------------------------------------------------------
//      

//}

/*
// assemble
// int_E v.n * (u.n) for a given function p and Dp/Dn = u.n
 
void  Assemble_RHS_v_n_u_n(double *rhs, 
                           TFESpace2D *U_Space,
                           TFEFunction2D *p,
                           int N_InnerInterfaceJoints,
                           TInnerInterfaceJoint **innerInterfaceJoints,
                           double mult, bool NitscheBC)
{
  cout << " Assemble_RHS_v_n_u_n on inner interface " << endl;
  //double tx, ty;
  double nx, ny;
  double x0,y0,vec_x,vec_y;
  int N_BaseFunct, *N_BaseFuncts;
  int ActiveBound, N_Unkonwns;
  TInnerInterfaceJoint *IJoint;
  TBaseCell *cell;
  BaseFunct2D *BaseFuncts;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  int fe_degree;
  int N_LinePoints;
  int BaseVectDim=1;
  int *BeginIndex, *GlobalNumbers;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans2D *F_K;
  FE2D FEId;
  double *LineWeights, *zeta;
  double hE;
  double **uref, **uxiref, **uetaref;
  double xonEdge,yonEdge;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double length_1D_reference;
  double pvalue[3];
  double u1val, u2val;
  double u_n;
  double v,value, localfactor;
  int *DOF;
  int TestDOF;
  
  BeginIndex = U_Space->GetBeginIndex();
  GlobalNumbers = U_Space->GetGlobalNumbers();
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  ActiveBound = U_Space->GetActiveBound();
  N_Unkonwns = U_Space->GetN_DegreesOfFreedom();
  // all cells in this subdomain have the same id
  int ID = U_Space->GetCollection()->GetCell(0)->GetReference_ID();
  
  for(int j=0; j<N_InnerInterfaceJoints; j++) 
  {
    IJoint = innerInterfaceJoints[j];
    // get physical coordinates on edge, one vertex is at (x0,y0), the other 
    // vertex is at ( x0+vec_x, y0+vec_y )
    IJoint->GetParams(x0, y0, vec_x,vec_y);
    // compute length of the edge
    hE = IJoint->GetLength();
    // tangential normal vector to this boundary (normalized)
    IJoint->GetNormal(nx,ny);
    //IJoint->GetTangent(tx,ty);
    
    localfactor = (NitscheBC) ? (1.0/hE) : (1.0);
    
    // loop over the two neighboring cells
    for (int i=0; i<2; i++)
    {
      cell = IJoint->GetNeighbour(i);
      if(cell->GetReference_ID() != ID) 
        // assemble only on one side of the interface.
        continue;
      FEId = U_Space->GetFE2D(0, cell);
      fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(
                                                              LineQuadFormula);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
      N_BaseFunct = N_BaseFuncts[FEId];
      
      RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEId);
      F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
      
      DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()];
      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], 
                            LineQuadFormula, IJoint->GetIndexInNeighbor(cell));
      uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], 
                       LineQuadFormula, IJoint->GetIndexInNeighbor(cell), D10);
      uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], 
                       LineQuadFormula, IJoint->GetIndexInNeighbor(cell), D01);
      F_K->SetCell(cell);
      // loop of quadrature points
      for(int k=0;k<N_LinePoints;k++)
      {
        switch(RefTrans) 
        {
          case TriaAffin:
            ((TTriaAffin*)F_K)->GetOrigValues(IJoint->GetIndexInNeighbor(cell),
                                              zeta[k],N_BaseFunct,
                                              uref[k], uxiref[k], uetaref[k],
                                              uorig, uxorig, uyorig, 
                                              BaseVectDim);
          break;
          
          case QuadAffin:
            ((TQuadAffin*)F_K)->GetOrigValues(IJoint->GetIndexInNeighbor(cell),
                                              zeta[k], N_BaseFunct,
                                              uref[k], uxiref[k], uetaref[k],
                                              uorig, uxorig, uyorig,
                                              BaseVectDim);
          break;
          
          default: cout << "Quad or Tria only allowed, check VectFEBoundInt() " << endl;
            exit(0);
          break;
        } // endswitch
        
        //zeta \in [-1,1]
        length_1D_reference = 2.;
        xonEdge = x0+((zeta[k]+1)/2)*(vec_x);
        yonEdge = y0+((zeta[k]+1)/2)*(vec_y);
        
        // get the values of rhs
        p->FindGradient(xonEdge, yonEdge, pvalue);
        //p->FindGradientLocal(cell, cell->GetCellIndex(), xonEdge, yonEdge, pvalue);
        // fctvalue: u,dx_u,dy_u
        u1val = pvalue[1]; // dx_p
        u2val = pvalue[2]; // dy_p
        
        u_n = u1val*nx + u2val*ny;
        // exact for solution from ana_sol_Stokes_Darcy.h
        //u_n = -sin(yonEdge*Pi/2)*cos(xonEdge*Pi/2) + 1 - xonEdge;
        // exact for solution from rectangle.h
        //OutPut("x " << xonEdge << " \t u_n " << u_n << ", p " << pvalue[0] 
        //            << ", u1 " << u1val << ", u2 " << u2val 
        //            << ", n(" << nx << "," << ny << ")" << endl);
        //u_n = -4*xonEdge*(1-xonEdge);
        for(int l=0;l<N_BaseFunct;l++)
        {
          TestDOF = DOF[l];
          if(TestDOF>=ActiveBound)
            continue;
          // updating rhs: int_gamma u\cdot n v \cdot n 
          v = uorig[l]; // value of test function
          value  = u_n*LineWeights[k]*(hE/2);
          value *= v*mult*localfactor;
          // add to right hand side
          rhs[TestDOF] += value*nx;
          rhs[TestDOF+N_Unkonwns]+= value*ny;
        } //for(l=0;l<N_BaseFunct;l++)
      }
    }
  }
}

*/
//__________________________________________________________________
//v_n_u_n
//__________________________________________________________________

void BoundaryAssembling2D::BoundaryAssemble_on_Matrix_u_n_v_n(TMatrix2D **MatrixA11,
                                                              TMatrix2D **MatrixA12,
                                                              TMatrix2D **MatrixA21,
                                                              TMatrix2D **MatrixA22,
                                                              const TFESpace2D *U_Space,
                                                              int boundary_component_id,
                                                              double mult
                                                              )
{
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll= U_Space->GetCollection();
    coll->get_edge_list_on_component(boundary_component_id,boundaryEdgeList);
    
    
    BoundaryAssemble_on_Matrix_u_n_v_n(MatrixA11,
                                       MatrixA12,
                                       MatrixA21,
                                       MatrixA22,
                                       U_Space,
                                       boundaryEdgeList,
                                       mult);
}

void BoundaryAssembling2D::BoundaryAssemble_on_Matrix_u_n_v_n(TMatrix2D **MatrixA11,
                                                              TMatrix2D **MatrixA12,
                                                              TMatrix2D **MatrixA21,
                                                              TMatrix2D **MatrixA22,
                                                              const TFESpace2D *U_Space,
                                                              std::vector<TBoundEdge*> &boundaryEdgeList,
                                                              double mult)
{   int *BeginIndex = U_Space->GetBeginIndex();
    int *GlobalNumbers = U_Space->GetGlobalNumbers();
    int ActiveBound = U_Space->GetActiveBound();
    
    for(int m=0;m< boundaryEdgeList.size(); m++)
    {
        TBoundEdge *boundedge = boundaryEdgeList[m];
        TBaseCell *cell = boundedge->GetNeighbour(0);
        // get basis dimension and FE space data of cell i
        FE2D FEId = U_Space->GetFE2D(0,cell );
        
        int BaseVectDim = 1; // we assume only scalar FE
        // ---------------------------------------------------------------
        int joint_id = boundedge->get_index_in_neighbour(cell);
        
        // get a quadrature formula good enough for the velocity FE space
        int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        this->LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
        std::vector<double> quadWeights,quadPoints;
        get_quadrature_formula_data(fe_degree,quadPoints,quadWeights);
        
        //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);
        
        // compute values of all basis functions at all quadrature points
        std::vector< std::vector<double> > uorig,uxorig,uyorig;
        get_original_values(FEId, joint_id, cell, quadPoints,BaseVectDim, uorig,uxorig,uyorig);
        
        double x0, x1, y0, y1;
        boundedge->get_vertices(x0,  y0, x1, y1);
        // compute length of the edge
        double joint_length = boundedge->get_length();
        // normal vector to this boundary (normalized)
        double nx,ny;
        boundedge->get_normal(nx, ny);
        
        // -------------------------------------------
        // quadrature
        for(unsigned int k=0;k<quadPoints.size();k++)
        {
            
            //double uorig[N_BaseFunct], uxorig[N_BaseFunct];
            //double uyorig[N_BaseFunct];
            //get_original_values(FEId, joint_id, cell, quadPoints[k],BaseVectDim, uorig,uxorig,uyorig);
            // ----------------------------
            
            ///@attention in 1D the reference joint is [-1,1] => length = 2
            double reference_joint_length = 2;
            double x = x0+(quadPoints[k]+1.)/2.*(x1-x0);
            double y = y0+(quadPoints[k]+1.)/2.*(y1-y0);
            
            // mapping from local(cell) DOF to global DOF
            int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
            
            for(unsigned int l=0;l< uorig[k].size();l++)
            {
                int global_dof_from_local = DOF[l];
                
                // if the DOF is Dirichlet, continue
                if(global_dof_from_local >= ActiveBound)
                    continue;
                
                // updating Matrix: int_gamma u \cdot n v \cdot n
                double vtest = uorig[k][l]; // value of test function (vtest = vx = vy)
                double u = uorig[k][l]; // value of test function (u = ux = uy)
                
                double value11 = mult*quadWeights[k]*(vtest*nx)*(u*nx)*(joint_length/reference_joint_length);
                MatrixA11[k][l] += value11;			//->add(global_dof_from_local ,AnsatzDOF,value);
                double value12 = mult*quadWeights[k]*(vtest*nx)*(u*ny)*(joint_length/reference_joint_length);
                MatrixA12[k][l] += value12; 			//->add(global_dof_from_local ,AnsatzDOF,value);
                double value21 = mult*quadWeights[k]*(vtest*ny)*(u*nx)*(joint_length/reference_joint_length);
                MatrixA21[k][l] += value21; 			//->add(global_dof_from_local ,AnsatzDOF,value);
                double value22 = mult*quadWeights[k]*(vtest*ny)*(u*ny)*(joint_length/reference_joint_length);
                MatrixA22[k][l] += value22; 			//->add(global_dof_from_local ,AnsatzDOF,value);
                
            } //for(l=0;l<N_BaseFunct;l++)
        }
    } // endif
}


