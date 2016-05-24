// ======================================================================
// @(#)BoundaryAssembling2D.C        05/18/16
//
// Functions for (external and internal) boundary integral
//
// ======================================================================
#include <BoundaryAssembling2D.h>

#include <FEDatabase2D.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>
#include <TriaAffin.h>
#include <QuadAffin.h>

void  BoundaryAssemble_on_rhs_g_v_n(double **rhs, 
				    const TFESpace2D *U_Space, 
				    TFEFunction2D *given_data,
				    int boundary_component_id,
				    double mult
				    )
{
    
  // =========================================
  int *BeginIndex = U_Space->GetBeginIndex();
  int *GlobalNumbers = U_Space->GetGlobalNumbers();
  int N_U = U_Space->GetN_DegreesOfFreedom();
  int ActiveBound = U_Space->GetActiveBound();
  
  
  BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  int *N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
    
  // loop over all cells
  int N_Cells = U_Space->GetCollection()->GetN_Cells();
  for(int i=0;i<N_Cells;i++)
  {
    TBaseCell *cell = U_Space->GetCollection()->GetCell(i);

    // get basis dimension and FE space data of cell i
    FE2D FEId = U_Space->GetFE2D(0, cell);
    ///@todo check how to extend to RT elements
    int BaseVectDim = 1; // we assume only scalar FE
    
    // get reference element type and number of basis functions
    BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    int N_BaseFunct = N_BaseFuncts[FEId];
  
    
    // for each cell, we store the list of joints belonging to
    // component boundary_component_id
    std::vector<int> jointsOnComponent;
    jointsOnComponent.clear();

    int nEdges = cell->GetN_Edges();
    // ---------------------------------------------------------------
      for(int j=0;j<nEdges;j++)
    {
      if(!(cell->GetJoint(j)->InnerJoint())) 
	{
	  TJoint *joint = cell->GetJoint(j);
	  TBoundEdge *boundedge = (TBoundEdge *)joint;
	  TBoundComp *BoundComp = boundedge->GetBoundComp();
	  if (BoundComp->GetID() == boundary_component_id) {
	    jointsOnComponent.push_back(j);
	  }
	} // endif
    }//  for(j=0;j<N_Edges;j++)
    // ---------------------------------------------------------------
    
    // quadrature over the joints on the selected boundary
    for(unsigned int j=0;j<jointsOnComponent.size();j++) {

      int joint_id = jointsOnComponent[j];

      // get geometrical properties of the boundary joint
      double x0,x1,y0,y1;
      double nx,ny;
      double joint_length;
      
      x0 = cell->GetVertex(joint_id)->GetX();
      x1 = cell->GetVertex((joint_id+1)%nEdges)->GetX();
      y0 = cell->GetVertex(joint_id)->GetY();
      y1 = cell->GetVertex((joint_id+1)%nEdges)->GetY();
      
      // compute length of the edge
      joint_length = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));

      // normal vector to this boundary (normalized)
      nx = (y1-y0)/joint_length;
      ny = (x0-x1)/joint_length;

      // ---------------------------------------------------------------
      // get a quadrature formula good enough for the velocity FE space
      int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      // get the type of required quadrature (include/FE/Enumerations.h)
      QuadFormula1D LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
      // initialize points and weights of quadrature
      TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);

      ///@todo rewrite the GetFormulaData using (e.g.) vector<> class
      int nQuadPoints;
      double *quadWeights, *quadPoints;
      qf1->GetFormulaData(nQuadPoints, quadWeights, quadPoints);
      // ---------------------------------------------------------------
      
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);

      // -------------------------------------------
      // quadrature
      for(int k=0;k<nQuadPoints;k++)
      {

	// get values and derivatives of basis functions
	// according to selected quadrature formula
	double **reference_values = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], 
                                               LineQuadFormula, joint_id);
        double **derivative_xi = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], 
                                                      LineQuadFormula, 
                                                      joint_id, D10);
        double **derivative_eta = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], 
                                                        LineQuadFormula, 
                                                        joint_id, D01);

	
	double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
	double uyorig[MaxN_BaseFunctions2D];

	RefTrans2D RefTrans;
	TRefTrans2D *F_K;
        switch(RefElement)
	  {
	  case BFUnitTriangle:
            RefTrans = TriaAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaAffin*)F_K)->SetCell(cell);
            ((TTriaAffin*)F_K)->GetOrigValues(joint_id,  quadPoints[k],N_BaseFunct, 
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

	///@attention in 1D the reference joint is [-1,1] => length = 2
        double reference_joint_length = 2;
        double x = x0+(quadPoints[k]+1.)/2.*(x1-x0);
        double y = y0+(quadPoints[k]+1.)/2.*(y1-y0);
	// given_data_value[0] = G(x,y),
	// given_data_value[1] = dG/dx, given_data_value[2] = dG/dy
	double given_data_value[3];
        // get the value of rhs
        if(given_data)
          given_data->FindGradient(x, y, given_data_value);
        else
          given_data_value[0] = 1.0;

	// mapping from local(cell) DOF to global DOF
	int *DOF = GlobalNumbers + BeginIndex[i];

        for(int l=0;l<N_BaseFunct;l++)
        {
          int global_dof_from_local = DOF[l];

	  // if the DOF is Dirichlet, continue
          if(global_dof_from_local >= ActiveBound)
            continue;
	  
          // updating rhs: int_gamma rhsval v \cdot n
          double vx,vy;
	  vx = uorig[l]; // value of test function
          vy = uorig[l]; // value of test function
          // add for both components
          rhs[0][global_dof_from_local] +=
	    mult * quadWeights[k] * given_data_value[0] * (vx*nx) *
	    (joint_length/reference_joint_length);
          //rhs[0][global_dof_from_local+N_U] +=
	  rhs[1][global_dof_from_local] +=
	    mult * quadWeights[k] * given_data_value[0] * (vy*ny) *
	    (joint_length/reference_joint_length);

	  
        } //for(l=0;l<N_BaseFunct;l++)
      } //  for(k=0;k<N_LinePoints;k++
    } // for j=0,..joint.size()
  } //  for(i=0;i<N_Cells;i++)
}
