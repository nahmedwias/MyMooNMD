// ======================================================================
// @(#)BoundaryAssembling2D.C        05/18/16
//
// Functions for (external and internal) boundary integral
//
// authors: Alfonso Caiazzo and Laura Blank
// ======================================================================
/*
   List of possible improvements
 * link each cell with its FE (in order not to use a global function to get FE properties)
 * create an FE class, in order not to recover all properties via FEDatabase::FunctionName(FeID)
 * rewrite the GetFormulaData using (e.g.) vector<> class
 * avoid to define double* with MaxN... (use vector instead?)
 * alternative: pass a list of joints (each with its own FE space on it)?
 * the edge list could be generated outside this class
 */

#include <FEDatabase2D.h>
#include <TriaAffin.h>
#include <QuadAffin.h>
#include <BoundaryAssembling2D.h>
#include <BoundEdge.h>
#include <Collection.h>


//========================================================================
void BoundaryAssembling2D::rhs_g_v_n(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data,
    int boundary_component_id,
    double mult)
{
  // Create a list of those boundary edges that are on the boundary component with given ID
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

  rhs_g_v_n(rhs, U_Space, given_boundary_data, boundaryEdgeList, mult);
}

void BoundaryAssembling2D::rhs_g_v_n(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
{
 // Go through all boundary edges on the current boundary component
  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    // current edge
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // mapping from local (cell) DOF to global DOF
    int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    // --------------------------------------------------------------
    ///@todo put the following part into FESpace2D::getEdgeQuadratureData
    // compute values of all basis functions
    // and their first partial derivatives at all quadrature points
    // get basis dimension and FE space data of the current cell
    FE2D FEId = U_Space->GetFE2D(0, cell);
    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);
    // get a quadrature formula good enough for the argument of the integral
    int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2 * fe_degree);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    std::vector< std::vector<double> > uorig, u_dx_orig, u_dy_orig;
    get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim,
        uorig, u_dx_orig, u_dy_orig,	LineQuadFormula);

    // get normal and edge length
    double x_0, x_1, y_0, y_1;
    boundedge->get_vertices(x_0, y_0, x_1, y_1);
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);
    double reference_joint_length = 2;
    double transformationDeterminant = joint_length/reference_joint_length;

    // ***** COMPUTE INTEGRAL *****
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      // get the value of (\nabla u - pI) on the boundary component (here denoted by g) 
      double value;

      if (given_boundary_data != nullptr)
      {
        double x = x_0 + (quadPoints[k]+1.)/2. * (x_1-x_0);
        double y = y_0 + (quadPoints[k]+1.)/2. * (y_1-y_0);
        double T;
        boundedge->GetBoundComp()->GetTofXY(x, y, T);
        given_boundary_data(boundedge->GetBoundComp()->GetID(), T, value);
      }
      else
      {
        value = 1;
      }

      double commonFactor = mult * quadWeights[k] * transformationDeterminant;

      for (unsigned int l = 0; l < uorig[k].size(); l++)
      {
        int global_dof_from_local = DOF[l];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local < U_Space->GetActiveBound())
        {
          // updating rhs: int_gamma rhsval v \cdot n
          double v1 = uorig[k][l]; // value of test function (vtest = v1 = v2)
          double v2 = v1;

          // add for both components
          rhs.block(0)[global_dof_from_local] += commonFactor * value * v1 * n1;

          if (rhs.n_blocks() == 2)
          {
            rhs.block(1)[global_dof_from_local] += commonFactor * value * v2 * n2;
          }
        }
      }
    }
  }
}

//=================================================================================
void BoundaryAssembling2D::rhs_uD_v(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data1,
    BoundValueFunct2D *given_boundary_data2,
    int boundary_component_id,
    double mult,
    bool rescale_by_h)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  rhs_uD_v(rhs, U_Space, given_boundary_data1, given_boundary_data2,
      boundaryEdgeList, mult, rescale_by_h);
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

  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);

    // get basis dimension and FE space data of cell i
    FE2D FEId = U_Space->GetFE2D(0, cell);

    int BDComponent = boundedge->GetBoundComp()->GetID();

    int BaseVectDim = 1; // we assume only scalar FE // Only for BDM and RT elements \neq 1
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space (here exact to 2*fe_degree)
    int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, u_dx_orig, u_dy_orig;
    get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, u_dx_orig, u_dy_orig, LineQuadFormula);

    double x_0, x_1, y_0, y_1;
    boundedge->get_vertices(x_0, y_0, x_1, y_1);
    // compute length of the edge
    double joint_length = boundedge->get_length();

    // quadrature
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;
      double x = x_0+(quadPoints[k]+1.)/2.*(x_1-x_0);
      double y = y_0+(quadPoints[k]+1.)/2.*(y_1-y_0);

      double T;
      boundedge->GetBoundComp()->GetTofXY(x, y, T);

      // get the boundary values of rhs
      double value1, value2;
      if (given_boundary_data1 != nullptr)
      {
        given_boundary_data1(BDComponent, T, value1);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_uD_v(), a default value for "
            "given_boundary_data1 is used, which might not fit your example.");
        value1 = 1.;
      }

      if(given_boundary_data2 != nullptr)
      {
        given_boundary_data2(BDComponent, T, value2);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_uD_v(), a default value for "
            "given_boundary_data2 is used, which might not fit your example.");
        value2 = 0.0;
      }

      // mapping from local (cell) DOF to global DOF
      int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];

      for (unsigned int l = 0; l < uorig[k].size(); l++)
      {
        int global_dof_from_local = DOF[l];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local >= ActiveBound)
          continue;

        // updating rhs: int_gamma rhsval[2] v
        double v1 = uorig[k][l]; // value of test function (vtest = vx = vy)
        double v2 = v1;

        // add for both components
        if (!rescale_by_h)
        {
          rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value1 * v1 *
            (joint_length/reference_joint_length);

          rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value2 * v2 *
            (joint_length/reference_joint_length);
        }
        else
        {
          rhs.block(0)[global_dof_from_local] += (mult * quadWeights[k] * value1 * v1 *
              (joint_length/reference_joint_length)) / joint_length;

          rhs.block(1)[global_dof_from_local] += (mult * quadWeights[k] * value2 * v2 *
              (joint_length/reference_joint_length)) /joint_length;
        }
      }
    }
  }
}

//========================================================================
void BoundaryAssembling2D::rhs_gradv_n_uD(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data1,
    BoundValueFunct2D *given_boundary_data2,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  rhs_gradv_n_uD(rhs, U_Space, given_boundary_data1, given_boundary_data2, boundaryEdgeList, mult);
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

  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);

    // get basis dimension and FE space data of cell i
    FE2D FEId = U_Space->GetFE2D(0, cell);

    int BDComponent = boundedge->GetBoundComp()->GetID();

    int BaseVectDim = 1; // we assume only scalar FE // Only for BDM and RT elements \neq 1
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space (here exact to 2*fe_degree)
    int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, u_dx_orig, u_dy_orig;
    get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, u_dx_orig, u_dy_orig, LineQuadFormula);

    double x_0, x_1, y_0, y_1;
    boundedge->get_vertices(x_0, y_0, x_1, y_1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // quadrature
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;
      double x = x_0 + (quadPoints[k]+1.)/2. * (x_1-x_0);
      double y = y_0 + (quadPoints[k]+1.)/2. * (y_1-y_0);

      double T;
      boundedge->GetBoundComp()->GetTofXY(x, y, T);

      // get the boundary values of rhs
      double value1, value2;
      if (given_boundary_data1 != nullptr)
      {
        given_boundary_data1(BDComponent, T, value1);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_gradv_n_uD(), a default value for "
            "given_boundary_data1 is used, which might not fit your example.");
        value1 = 1;
      }

      if(given_boundary_data2 != nullptr)
      {
        given_boundary_data2(BDComponent, T, value2);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_gradv_n_uD(), a default value for "
            "given_boundary_data2 is used, which might not fit your example.");
        value2 = 0.0;
      }

      // mapping from local (cell) DOF to global DOF
      int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];

      for (unsigned int l = 0; l < uorig[k].size(); l++)
      {
        int global_dof_from_local = DOF[l];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local >= ActiveBound)
          continue;

        // updating rhs: int_gamma rhsval[2] v
        double v1_dx = u_dx_orig[k][l]; // value of test function (vtest = (vx,vy) )
        double v1_dy = u_dy_orig[k][l];
        double v2_dy = v1_dy;
        double v2_dx = v1_dx;

        rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value1 *
            (v1_dx * n1 + v1_dy * n2) *  (joint_length/reference_joint_length);

        rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value2 *
            (v2_dx * n1 + v2_dy * n2) * (joint_length/reference_joint_length);
      }
    }
  }
}

//======================================================================
void BoundaryAssembling2D::matrix_u_n_v_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    int boundary_component_id,
    double mult,
    bool rescale_by_h)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

  matrix_u_n_v_n(M, U_Space, boundaryEdgeList, mult, rescale_by_h);
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

  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // get basis dimension and FE space data of cell i
    FE2D FEId = U_Space->GetFE2D(0, cell);

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig, LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // quadrature
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // mapping from local(cell) DOF to global DOF
      int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];

      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      if (rescale_by_h)
      {
        scale_factor = scale_factor/joint_length;
      }

      // loop on test functions
      for (unsigned int l1 = 0; l1 < uorig[k].size(); l1++)
      {
        int test_DOF = DOF[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1 = uorig[k][l1];
        double v2 = v1; // x and y component have the same FE space

        // loop on ansatz functions
        for (unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
        {
          int ansatz_DOF = DOF[l2];
          double u1 = uorig[k][l2];
          double u2 = u1; // x and y component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * (v1 * n1) * (u1 * n1) ); // A11
          blocks[1]->add( test_DOF, ansatz_DOF, scale_factor * (v1 * n1) * (u2 * n2) ); // A12
          blocks[3]->add( test_DOF, ansatz_DOF, scale_factor * (v2 * n2) * (u1 * n1) ); // A21
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * (v2 * n2) * (u2 * n2) ); // A22
        }
      }
    }
  }
}
// ToDo: implement the rhs of this term

//===============================================================================
void BoundaryAssembling2D::matrix_cornerjump_u_n_cornerjump_v_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    //int boundary_component_id,
    int nBoundaryParts,
    double mult)
{ 
  //Check if there is only one connected boundary (i.e., e.g., no holes)
  if (nBoundaryParts != 1)
  {
    Output::print("WARNING BoundaryAssembling2D::matrix_cornerjump_u_n_cornerjump_v_n() is defined for one single Boundary Part only. In your example you have ", nBoundaryParts , " many boundary parts. Note that nBoundaryParts refers to the number of connected boundary paths. ");
  }
  else
  { 
    std::vector<TBoundEdge*> boundaryEdgeList;
    TCollection *coll = U_Space->GetCollection();
    coll->get_boundary_edge_list(boundaryEdgeList);
    matrix_cornerjump_u_n_cornerjump_v_n(M, U_Space, boundaryEdgeList, mult);
  }
}

/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/

void BoundaryAssembling2D::matrix_cornerjump_u_n_cornerjump_v_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
{
  int *BeginIndex = U_Space->GetBeginIndex();
  int *GlobalNumbers = U_Space->GetGlobalNumbers();
  //int ActiveBound = U_Space->GetActiveBound();

  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  // find a 'pair' of corner edges
  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge_1 = boundaryEdgeList[m];

    // Check if boundary edge is on boundary components of Nitsche type
    for (int km = 0; km < TDatabase::ParamDB->n_nitsche_boundary; km++)           
    {
      if (boundedge_1->GetBoundComp()->GetID() == TDatabase::ParamDB->nitsche_boundary_id[km])
      {
        TBoundEdge *boundedge_2;
        double n1_E1, n1_E2, n2_E1, n2_E2;
        boundedge_1->get_normal(n1_E1, n1_E2);

        for (size_t k = m+1; k < boundaryEdgeList.size(); k++)
        { 
          boundedge_2 = boundaryEdgeList[k]; 

          for (int kl = 0; kl < TDatabase::ParamDB->n_nitsche_boundary; kl++)           
          {
            // Check if boundary edge is on boundary components of Nitsche type
            if (boundedge_2->GetBoundComp()->GetID() == TDatabase::ParamDB->nitsche_boundary_id[kl])
            {
              boundedge_2->get_normal(n2_E1, n2_E2);

              if ( fabs(n1_E1 * n2_E1 + n1_E2 * n2_E2) != 1 )  // a.b=|a||b|cos(alpha)
              {
                // get coordinates of relevant (Nitsche) corners xc,yc
                double x0_E1, y0_E1, x1_E1, y1_E1, x0_E2, y0_E2, x1_E2, y1_E2;
                boundedge_1->get_vertices(x0_E1, y0_E1, x1_E1, y1_E1);
                boundedge_2->get_vertices(x0_E2, y0_E2, x1_E2, y1_E2);

                //Check if the edges E1 and E2 share a vertex
                std::vector<double> xc, yc;
                if ( ( (x0_E1 == x0_E2) & (y0_E1 == y0_E2)) || (( x0_E1 == x1_E2) & (y0_E1 == y1_E2) ) ) 
                {
                  xc.push_back(x0_E1);
                  yc.push_back(y0_E1);
                  //Output::print("HIER!!!!!!");
                  //Output::print("lenght xc: ", xc.size());
                  //Output::print("(xc,yc)= ", xc[0], ", ", yc[0]);
                  //Output::print("Detected a pair of corner edges.");
                  //Output::print("num boundary edges: ", boundaryEdgeList.size());
                  //Output::print("OK, (m,k) = ", m ,",",k);

                  int locdof_corner_1, locdof_corner_2;
                  find_cornerDofs_in_boundarycells(xc, yc, U_Space, boundedge_1, boundedge_2, locdof_corner_1, locdof_corner_2);

                  //Output::print("locdof_corner_1, locdof_corner_2: ", locdof_corner_1, locdof_corner_2);

                  // mapping from local(cell) DOF to global DOF
                  TBaseCell *cell_1 = boundedge_1->GetNeighbour(0);
                  int *DOF = GlobalNumbers + BeginIndex[cell_1->GetCellIndex()];
                  int test_DOF = DOF[locdof_corner_1];
                  int ansatz_DOF = test_DOF;

                  //Output::print("test_Dof: ", test_DOF);
                  //see the note about blocks at the beginning of the function)
                  // In each corner point, there is exactly one basis function which is nonzero at (xc,yc).
                  // In fact it has the value (1,1) in (xc,yc).
                  blocks[0]->add( test_DOF, ansatz_DOF, mult * (n1_E1 - n1_E2) * (n1_E1 - n1_E2) ); // (1,0)^T * (n_E1-n_E2) --> A11
                  blocks[1]->add( test_DOF, ansatz_DOF, mult * (n2_E1 - n2_E2) * (n1_E1 - n1_E2) ); // A12
                  blocks[3]->add( test_DOF, ansatz_DOF, mult * (n1_E1 - n1_E2) * (n2_E1 - n2_E2) ); // A21
                  blocks[4]->add( test_DOF, ansatz_DOF, mult * (n2_E1 - n2_E2) * (n2_E1 - n2_E2)  ); // (0,1)^T * (n_E1-n_E2) --> A22
                  Output::print("Here, something was added due to corner stabilization.");
                }
                else if ( ( (x1_E1 == x1_E2) & (y1_E1 == y1_E2)) || (( x1_E1 == x0_E2) & (y1_E1 == y0_E2) ) )
                {
                  xc.push_back(x1_E1);
                  yc.push_back(y1_E1);
                  //Output::print("HIER!!!!!!");
                  //Output::print("lenght xc: ", xc.size());
                  //Output::print("(xc,yc)= ", xc[0], ", ", yc[0]);
                  //Output::print("Detected a pair of corner edges.");
                  //Output::print("num boundary edges: ", boundaryEdgeList.size());
                  //Output::print("OK, (m,k) = ", m ,",",k);

                  int locdof_corner_1, locdof_corner_2;
                  find_cornerDofs_in_boundarycells(xc, yc, U_Space, boundedge_1, boundedge_2, locdof_corner_1, locdof_corner_2);

                  //Output::print("locdof_corner_1, locdof_corner_2: ", locdof_corner_1, locdof_corner_2);

                  // mapping from local(cell) DOF to global DOF
                  TBaseCell *cell_1 = boundedge_1->GetNeighbour(0);
                  int *DOF = GlobalNumbers + BeginIndex[cell_1->GetCellIndex()];
                  int test_DOF = DOF[locdof_corner_1];
                  int ansatz_DOF = test_DOF;

                  //Output::print("test_Dof: ", test_DOF);
                  //see the note about blocks at the beginning of the function)
                  // In each corner point, there is exactly one basis function which is nonzero at (xc,yc).
                  // In fact it has the value (1,1) in (xc,yc).
                  blocks[0]->add( test_DOF, ansatz_DOF, mult * (n1_E1 - n1_E2) * (n1_E1 - n1_E2) ); // (1,0)^T * (n_E1-n_E2) --> A11
                  blocks[1]->add( test_DOF, ansatz_DOF, mult * (n2_E1 - n2_E2) * (n1_E1 - n1_E2) ); // A12
                  blocks[3]->add( test_DOF, ansatz_DOF, mult * (n1_E1 - n1_E2) * (n2_E1 - n2_E2) ); // A21
                  blocks[4]->add( test_DOF, ansatz_DOF, mult * (n2_E1 - n2_E2) * (n2_E1 - n2_E2) ); // (0,1)^T * (n_E1-n_E2) --> A22
                  Output::print("Here, something was added due to corner stabilization.");
                }
                else
                { 
                  //Output::print("x0_E1, y0_E1, x1_E1, y1_E1: ", x0_E1 ,",", y0_E1,",", x1_E1,",", y1_E1);
                  //Output::print("x0_E2, y0_E2, x1_E2, y1_E2: ", x0_E2 ,",", y0_E2,",", x1_E2,",", y1_E2);
                }
              }
            }
          }
        }
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::find_cornerDofs_in_boundarycells(std::vector<double> xc, std::vector<double> yc, 
    const TFESpace2D *U_Space,
    TBoundEdge *boundedge_1, TBoundEdge *boundedge_2, 
    int &locdof_corner_1, int &locdof_corner_2)
{
  int *BeginIndex = U_Space->GetBeginIndex();
  int *GlobalNumbers = U_Space->GetGlobalNumbers();
  TBaseCell *cell_1 = boundedge_1->GetNeighbour(0);
  TBaseCell *cell_2 = boundedge_2->GetNeighbour(0);

  FE2D FEId = U_Space->GetFE2D(0, cell_1);
  //FE2D FEId_2 = U_Space->GetFE2D(0, cell_2);

  // BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  // BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
  int *N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  int N_BaseFunct = N_BaseFuncts[FEId];

  for (int locdof = 0; locdof < N_BaseFunct; locdof++)
  {
    for (unsigned int i = 0; i < xc.size(); i++)
    {
      int *DOF = GlobalNumbers + BeginIndex[cell_1->GetCellIndex()];
      int globaldof = DOF[locdof];
      double x, y; 
      U_Space->GetDOFPosition(globaldof, x, y);
      if ( (x == xc[i]) & (y == yc[i]) )
      { 
        locdof_corner_1 = locdof;
        break;
      }
    }
  }
  //Output::print("number of cell_1: ", cell_1->GetCellIndex());
  //Output::print("number of cell_2: ", cell_2->GetCellIndex());

  if (cell_1->GetCellIndex()  == cell_2->GetCellIndex() )
  {
    //Output::print("Detected a single corner cell.");
    locdof_corner_2 = locdof_corner_1;
  }
  else
  {
    FE2D FEId = U_Space->GetFE2D(0, cell_2);

    //    BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
    //    BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    int *N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
    int N_BaseFunct = N_BaseFuncts[FEId];

    for (int locdof = 0; locdof < N_BaseFunct; locdof++)
    {
      for (unsigned int i = 0; i < xc.size(); i++)
      {
        int *DOF = GlobalNumbers + BeginIndex[cell_2->GetCellIndex()];
        int globaldof = DOF[locdof];
        double x, y;
        U_Space->GetDOFPosition(globaldof, x, y);
        if ( (x == xc[i]) & (y == yc[i]))
        { 
          locdof_corner_2 = locdof;
          break;
        }
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// <d(u.n)/dtau, d(v.n)/dtau> at corners of the domain
void BoundaryAssembling2D::matrix_gradu_n_t_gradv_n_t(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

  matrix_gradu_n_t_gradv_n_t(M, U_Space, boundaryEdgeList, mult);
}
/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_gradu_n_t_gradv_n_t(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
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
    FE2D FEId = U_Space->GetFE2D(0, cell);

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig, LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double t1, t2, n1, n2;
    boundedge->get_tangent(t1, t2);
    boundedge->get_normal(n1, n2);

    // quadrature
    for(unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // mapping from local (cell) DOF to global DOF
      int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];

      // loop on test functions
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      for(unsigned int l1 = 0; l1 < uorig[k].size(); l1++)
      {
        int test_DOF = DOF[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1_dx = uxorig[k][l1];
        double v1_dy = uyorig[k][l1];
        double v2_dx = v1_dx; // v1 and v2 component have the same FE space
        double v2_dy = v1_dy; // v1 and v2 component have the same FE space

        // loop over ansatz functions
        for(unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
        {
          int ansatz_DOF = DOF[l2];
          double u1_dx = uxorig[k][l2];
          double u1_dy = uyorig[k][l2];
          double u2_dx = u1_dx; // u1 and u2 component have the same FE space
          double u2_dy = u1_dy; // u1 and u2 component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * ( (u1_dx * n1 * t1) * (v1_dx * n1 * t1) 
                +(u1_dx * n1 * t1) * (v1_dy * n1 * t2)
                +(u1_dy * n1 * t2) * (v1_dx * n1 * t1)
                +(u1_dy * n1 * t2) * (v1_dy * n1 * t2)) ); // A11
          blocks[1]->add( test_DOF, ansatz_DOF, scale_factor * ( (u2_dx * n2 * t1) * (v1_dx * n1 * t1) 
                +(u2_dx * n2 * t1) * (v1_dy * n1 * t2)
                +(u2_dy * n2 * t2) * (v1_dx * n1 * t1)
                +(u2_dy * n2 * t2) * (v1_dy * n1 * t2)) ); // A12
          blocks[3]->add( test_DOF, ansatz_DOF, scale_factor * ( (u1_dx * n1 * t1) * (v2_dx * n2 * t1) 
                +(u1_dx * n1 * t1) * (v2_dy * n2 * t2)
                +(u1_dy * n1 * t2) * (v2_dx * n2 * t1)
                +(u1_dy * n1 * t2) * (v2_dy * n2 * t2)) );// A21
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * (  (u2_dx * n2 * t1) * (v2_dx * n2 * t1) 
                +(u2_dx * n2 * t1) * (v2_dy * n2 * t2)
                +(u2_dy * n2 * t2) * (v2_dx * n2 * t1)
                +(u2_dy * n2 * t2) * (v2_dy * n2 * t2)) ); // A22
        }
      } 
    }
  } 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void BoundaryAssembling2D::matrix_gradu_n_v(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_gradu_n_v(M, U_Space, boundaryEdgeList, mult, boundary_component_id);
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
{
  int *BeginIndex = U_Space->GetBeginIndex();
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

  for(size_t m = 0; m < boundaryEdgeList.size(); m++)
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
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);

    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    //LB Debug start
    //Output::print("num_quadPoints_per_edge: ",quadPoints.size());
    //LB Debug end

    //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig, LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0,  y0, x1, y1);
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
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
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
          double u2_dx = u1_dx; // u1 and u2 component have the same FE space
          double u2_dy = u1_dy; // u1 and u2 component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * v1 * u1_dx * n1 ); // A11
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * v1 * u1_dy * n2 ); // A11
          //blocks[1]->add(test_DOF, ansatz_DOF, scale_factor*(vx)*uyy*nx ); // A12
          //blocks[3]->add(test_DOF, ansatz_DOF, scale_factor*(vx)*uxx*ny ); // A21
          blocks[4]->add(test_DOF, ansatz_DOF, scale_factor * v2 * u2_dy * n2 ); // A22
          blocks[4]->add(test_DOF, ansatz_DOF, scale_factor * v2 * u2_dx * n1 ); // A22

          //LB DEBUG start
          //                    Output::print("has value:",scale_factor * v1 * u1_dx * n1 + scale_factor * v1 * u1_dy * n2 );
          //                    Output::print(test_DOF);
          //                    test_DOF_vec[l1] = test_DOF;
          //                    ansatz_DOF_vec[l2] = ansatz_DOF;
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
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_gradv_n_u(M, U_Space, boundaryEdgeList, mult);
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
    FE2D FEId = U_Space->GetFE2D(0, cell);

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig, LineQuadFormula);

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
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      for (unsigned int l1 = 0; l1 < uorig[k].size(); l1++)
      {
        int test_DOF = DOF[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1_dx = uxorig[k][l1];
        double v1_dy = uyorig[k][l1];
        double v2_dx = v1_dx;
        double v2_dy = v1_dy; // x and y component have the same FE space

        // loop on ansatz functions
        for(unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
        {
          int ansatz_DOF = DOF[l2];

          double u1 = uorig[k][l2];
          double u2 = u1; // x and y component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * u1 * v1_dx * n1 ); // A11
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * u1 * v1_dy * n2 ); // A11
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * u2 * v2_dx * n1 ); // A22
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * u2 * v2_dy * n2 ); // A22
        }
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::matrix_u_v(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    int boundary_component_id,
    double mult,
    bool rescale_by_h)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_u_v(M, U_Space, boundaryEdgeList, mult, rescale_by_h);
}
/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11, A12, B1T, A21, A22, B2T, B1, B2, C)
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
    FE2D FEId = U_Space->GetFE2D(0, cell);

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*fe_degree);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(FEId, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig, LineQuadFormula);

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
        scale_factor = ( mult * quadWeights[k] * (joint_length/reference_joint_length)) /joint_length;
      }
      else
      {
        scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      }

      // loop on test functions
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
          double u1 = uorig[k][l2];
          double u2 = u1; // x and y component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * u1 * v1 ); // A11
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * u2 * v2 ); // A22
        }
      }
    }
  }
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
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_q_u_n(M, U_Space, P_Space, boundaryEdgeList, mult);
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
    FE2D FEId_U = U_Space->GetFE2D(0, cell);
    // get basis dimension and FE space data of cell i
    FE2D FEId_P = P_Space->GetFE2D(0, cell);

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree_U = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_U);
    int fe_degree_P = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_P);
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(fe_degree_P+fe_degree_U);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector<std::vector<double>> uorig, uxorig, uyorig;
    get_original_values(FEId_U, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig, LineQuadFormula);

    int BaseVectDim_P = 1; // we assume only scalar FE; nur bei Raviart-Thomas & BDM \neq 1

    // compute values of all basis functions at all quadrature points
    std::vector<std::vector<double>> porig, pxorig, pyorig;
    get_original_values(FEId_P, joint_id, cell, quadPoints, BaseVectDim_P, porig, pxorig, pyorig, LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // mapping from local(cell) DOF to global DOF
    // int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
    int *DOF_P = P_Space->GetGlobalDOF(cell->GetCellIndex());
    int *DOF_U = U_Space->GetGlobalDOF(cell->GetCellIndex());

    // quadrature
    for( unsigned int k = 0; k < quadPoints.size(); k++ )
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // loop on test functions
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);

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
          blocks[6]->add( test_DOF, ansatz_DOF, scale_factor * q * u1 * n1 ); // B1
          blocks[7]->add( test_DOF, ansatz_DOF, scale_factor * q * u2 * n2 ); // B2
        }
      }
    } 
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::rhs_q_uD_n(BlockVector &rhs,
    const TFESpace2D *U_Space,
    const TFESpace2D *P_Space,
    BoundValueFunct2D *given_boundary_data1,
    BoundValueFunct2D *given_boundary_data2,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = P_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

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
    FE2D FEId_U = U_Space->GetFE2D(0, cell);
    // get basis dimension and FE space data of cell i
    FE2D FEId_P = P_Space->GetFE2D(0, cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree_U = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_U);
    int fe_degree_P = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_P);

    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(fe_degree_U + fe_degree_P);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > porig, pxorig, pyorig;
    get_original_values(FEId_P, joint_id, cell, quadPoints, BaseVectDim_P, porig, pxorig, pyorig, LineQuadFormula);

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

      int BDComponent = boundedge->GetBoundComp()->GetID();

      // get the boundary values of rhs
      double value1, value2;
      if(given_boundary_data1 != nullptr)
      {
        given_boundary_data1(BDComponent, T, value1);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_q_uD_n(), a default value for "
            "given_boundary_data1 is used, which might not fit your example.");
        value1 = 1.;
      }

      if(given_boundary_data2 != nullptr)
      {
        given_boundary_data2(BDComponent, T, value2);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_q_uD_n(), a default value for "
            "given_boundary_data2 is used, which might not fit your example.");
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
        rhs.block(2)[global_dof_from_local] += mult * quadWeights[k] * q * (value1 * n1) *
          (joint_length/reference_joint_length);
        rhs.block(2)[global_dof_from_local] += mult * quadWeights[k] * q * (value2 * n2) *
          (joint_length/reference_joint_length);
      }   
    }
  } 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::matrix_p_v_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    const TFESpace2D *P_Space,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_p_v_n(M,U_Space, P_Space, boundaryEdgeList, mult);
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
    QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(fe_degree_P*fe_degree_U);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);

    //TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(FEId_U, joint_id, cell, quadPoints, BaseVectDim, uorig, uxorig, uyorig, LineQuadFormula);

    int BaseVectDim_P = 1; // we assume only scalar FE; nur bei Raviart-Thomas & BDM \neq 1

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double>> porig, pxorig, pyorig;
    get_original_values(FEId_P, joint_id, cell, quadPoints, BaseVectDim_P, porig, pxorig, pyorig, LineQuadFormula);

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
    for(unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // loop on test functions
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      for(unsigned int l1 = 0 ; l1 < uorig[k].size() ; l1++)
      {
        int test_DOF = DOF_U[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1 = uorig[k][l1];
        double v2 = v1; // x and y component have the same FE space

        // loop on ansatz functions
        for(unsigned int l2 = 0; l2 < porig[k].size(); l2++)
        {
          int ansatz_DOF = DOF_P[l2];
          double p = porig[k][l2];
          // ToDo: Check if block reference is correct:
          // (see the note about blocks at the beginning of the function)
          blocks[2]->add( test_DOF, ansatz_DOF, scale_factor * p * v1 * n1 ); // B1
          blocks[5]->add( test_DOF, ansatz_DOF, scale_factor * p * v2 * n2 ); // B2
        }  
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::get_quadrature_formula_data(std::vector<double> &P,
    std::vector<double> &W, const QuadFormula1D& LineQuadFormula) 
{
  // get the type of required quadrature (include/FE/Enumerations.h)
  // initialize points and weights of quadrature
  ///@attention LineQuadFormula has to be set before calling the function GetQuadFormulaData
  TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);

  int nQuadPoints;
  double *quadWeights, *quadPoints;
  qf1->GetFormulaData(nQuadPoints, quadWeights, quadPoints);
  P.resize(nQuadPoints);
  W.resize(nQuadPoints);

  for (int i = 0; i < nQuadPoints; i++){
    P[i] = quadPoints[i];
    W[i] = quadWeights[i];
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
    std::vector< std::vector<double> > &originalValues_y,
    const QuadFormula1D& LineQuadFormula)
{
  BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
  int *N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  int N_BaseFunct = N_BaseFuncts[FEId];

  TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);

  originalValues.resize(quadPoints.size());
  originalValues_x.resize(quadPoints.size());
  originalValues_y.resize(quadPoints.size());

  for( unsigned int k = 0; k < quadPoints.size(); k++ )
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


// =================================================================================
void BoundaryAssembling2D::nitsche_bc(BlockFEMatrix &s_matrix,BlockVector &s_rhs,
				      const TFESpace2D * v_space, const TFESpace2D *p_space,
				      BoundValueFunct2D * U1, BoundValueFunct2D *U2,
				      int bd_comp, double gamma, double mu,
				      int sym_u, int sym_p)
{
  //============================== PENALTY TERMS ===================================
  // mueff * gamma/h * (u,v)
  matrix_u_v(s_matrix, v_space, bd_comp, gamma*mu, true);  // rescale local integral by edge values

  // mueff * gamma/h * (uD,v) [rhs]
  rhs_uD_v(s_rhs, v_space, U1, U2, bd_comp, gamma*mu, true);   // rescale local integral by edge values

 /* // sigma *L_0^2 * gamma/h (u.n,v.n)
  matrix_u_n_v_n(s_matrix, v_space,  bd_comp,
      gamma* sigma * TDatabase::ParamDB->L_0 * TDatabase::ParamDB->L_0, true); // rescale local integral by edge values

   //todo: rhs for matrix_u_n_v_n()
 */

  //=========================== PENALTY-FREE TERMS =================================
  // - (mu grad(u)n,v)
  matrix_gradu_n_v(s_matrix, v_space, bd_comp, -1. * mu); 

  // - sign_u * (u,mu grad(v)n) [sign_u=1: symmetrix, -1: skew-symmetric]
  matrix_gradv_n_u(s_matrix, v_space, bd_comp, (-1) * sym_u * mu);
  
  // - sign_u * (uD,mu grad(v)n) [rhs]
  rhs_gradv_n_uD(s_rhs, v_space, U1, U2, bd_comp, (-1) * sym_u * mu );

  // (pn,v)
  matrix_p_v_n(s_matrix, v_space, p_space, bd_comp, 1.);

  // sign_p * (u,qn)
  matrix_q_u_n(s_matrix, v_space, p_space, bd_comp, sym_p);

  // sign_p * (uD,qn) [rhs]
  rhs_q_uD_n(s_rhs, v_space, p_space, U1, U2, bd_comp, sym_p);

  //todo:  define corner stab in 3D

}
