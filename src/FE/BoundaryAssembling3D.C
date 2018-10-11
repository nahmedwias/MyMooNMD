// ============================================================================ //
// @(#)BoundaryAssembling3D.C        28.10.16                                   //
//                                                                              //
// Functions for (external and internal) boundary integral                      //
//                                                                              //
// ============================================================================ //

#include <FEDatabase3D.h>
#include <TriaAffin.h>
#include <QuadAffin.h>
#include <BoundaryAssembling3D.h>
#include <Collection.h>
#include <BoundFace.h>
#include <BaseCell.h>

// ===========================================================================

// int_{Gamma} mult*given_boundary_data(x,y,z)*<v,normal>
void BoundaryAssembling3D::rhs_g_v_n(BlockVector &rhs,
    const TFESpace3D *U_Space,
    BoundValueFunct3D *given_boundary_data,
    std::vector<TBaseCell*> &boundaryCells,
    int componentID,
    double mult)
{
  /// @todo: create a boundary cell list for each boundary component

  // now we loop always over the whole collection: inefficient!
  TCollection* coll = U_Space->GetCollection();

  for(size_t i = 0; i < coll->GetN_Cells(); i++)
  {
    TBaseCell* cell = coll->GetCell(i); //boundaryCells[i];

    // mapping from local (cell) DOF to global DOF
    int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    for(size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace) {

        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        if (boundface->GetBoundComp()->get_physical_id() == componentID)
        {
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights, qPointsT, qPointsS;
          std::vector< std::vector<double> > basisFunctionsValues;
          //this->getQuadratureData(U_Space,cell,joint_id,qWeights,
          //			  qPointsT,qPointsS,basisFunctionsValues);
          U_Space->getFaceQuadratureData(cell, joint_id,
              qWeights, qPointsT, qPointsS,
              basisFunctionsValues);

          std::vector<double> normal;
          double transformationDeterminant;
          cell->computeNormalAndTransformationData(joint_id, normal,
              transformationDeterminant);

          double x, y, z;
          double value;
          normal.resize(3);
          boundface->GetXYZofTS(qPointsT[0], qPointsS[0], x, y, z);
          boundface->get_normal_vector(x, y, z, normal[0], normal[1], normal[2]);

          // loop over Gauss points
          for (size_t l = 0; l < qWeights.size(); l++)
          {
            // NEW LB 11.10.18
            // get the boundary values of rhs
            boundface->GetBoundComp()->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);

            if(given_boundary_data != nullptr)
              /// NEW LB 11.10.18
              given_boundary_data(x, y, z, value);
            // OLD value = 1.;
            else
              value = 1.;

            double scale_factor = mult * qWeights[l] * transformationDeterminant;

            for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
            {
              int global_dof_from_local = DOF[k];

              if (global_dof_from_local < U_Space->GetActiveBound())
              {
                double v1 = basisFunctionsValues[l][k]; // value of test function (vtest = vx = vy =vz)
                double v2 = v1;
                double v3 = v1;

                rhs.block(0)[global_dof_from_local] += scale_factor * value * (v1*normal[0]);
                rhs.block(1)[global_dof_from_local] += scale_factor * value * (v2*normal[1]);
                rhs.block(2)[global_dof_from_local] += scale_factor * value * (v3*normal[2]);
              }
            }
          }
        }
      }
    }
  }
}


// ===========================================================================
void BoundaryAssembling3D::getQuadratureData(const TFESpace3D *fespace,TBaseCell *cell, int m,
    std::vector<double>& qWeights,std::vector<double>& qPointsT,
    std::vector<double>& qPointsS,
    std::vector< std::vector<double> >& basisFunctionsValues)
{
  int nFaceVertices = cell->getNumberOfFaceVertices(m);

  // set quadrature formula and compute quadrature info
  FE3D FEId = fespace->GetFE3D(cell->GetCellIndex(),cell);
  int fe_degree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);

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

  int N_Points;
  double* faceWeights;
  double *t,*s;
  // get a quadrature formula good enough for the velocity FE space
  TQuadFormula2D *qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
  qf2->GetFormulaData(N_Points, faceWeights, t, s);

  // ====================================
  // generate data on reference mesh cell for the 2d face of 3d cell
  TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)
  ->MakeRefElementData(FaceQuadFormula);

  BaseFunct3D *BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  int* N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // values of base functions in all quadrature points on face
  double **JointValues = TFEDatabase3D::GetJointValues3D
      (BaseFuncts[FEId], FaceQuadFormula, m);

  TFEDatabase3D::GetBaseFunct3D(BaseFuncts[FEId])->ChangeBF(fespace->GetCollection(),
      cell, N_Points, JointValues);

  // ====================================
  // convert the double* to vectors

  qWeights.resize(N_Points);
  qPointsT.resize(N_Points);
  qPointsS.resize(N_Points);
  for (size_t k=0; k<(size_t)N_Points; k++) {
    qWeights[k] = faceWeights[k];
    qPointsT[k] = t[k];
    qPointsS[k] = s[k];
  }


  basisFunctionsValues.resize(qWeights.size());                       
  for (unsigned int l=0; l<qWeights.size(); l++) {
    basisFunctionsValues[l].resize(N_BaseFunct[FEId]);                          
    for (unsigned int k=0; k<basisFunctionsValues[l].size(); k++) {
      basisFunctionsValues[l][k]=JointValues[l][k];
    }
  }
  // ====================================

}




// ------------------------------
// NOT TESTED //

// ===========================================================================
// int_{Gamma} p v n
void BoundaryAssembling3D::matrix_p_v_n(BlockFEMatrix &M,
    const TFESpace3D *U_Space,
    const TFESpace3D *P_Space,
    std::vector<TBaseCell*> &boundaryCells,
    int componentID,
    double mult)
{
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */
  ///@todo create a boundary cell list
  // loop over all boundary cells (at the moment it consists of all cells!!)
  for (size_t i = 0; i< boundaryCells.size(); i++)
  {
    TBaseCell* cell = boundaryCells[i];
    // mapping from local (cell) DOF to global DOF
    int *DOF_u = U_Space->GetGlobalDOF(cell->GetCellIndex());
    // mapping from local (cell) DOF to global DOF
    int *DOF_p = P_Space->GetGlobalDOF(cell->GetCellIndex());

    // loop over the faces (joints)
    for (size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      // select boundary faces
      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        //Output::print(boundface->GetBoundComp()->get_physical_id()," comp, ",componentID);

        if (boundface->GetBoundComp()->get_physical_id() == componentID)
        {
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights_u, qPointsT_u, qPointsS_u, qWeights_p, qPointsT_p, qPointsS_p; // S,T for parametrization in 3D
          std::vector< std::vector<double> > basisFunctionsValues_u, basisFunctionsValues_p;
          this->getQuadratureData(U_Space, cell, joint_id,
              qWeights_u, qPointsT_u, qPointsS_u, basisFunctionsValues_u);
          this->getQuadratureData(P_Space, cell, joint_id,
              qWeights_p, qPointsT_p, qPointsS_p, basisFunctionsValues_p);

          std::vector<double> normal;
          double transformationDeterminant;
          this->computeNormalAndTransformationData(cell, joint_id, normal,
              transformationDeterminant);

          // rescale local integral by mesh-size (important for Nitsche boundary)
          for (size_t l = 0; l < qWeights_u.size(); l++)
          {
            double scale_factor  = mult * qWeights_u[l] * transformationDeterminant;

            for (size_t k1 = 0; k1 < basisFunctionsValues_u[l].size(); k1++)
            {
              int global_dof_from_local_u = DOF_u[k1]; // Test-DOF

              if (global_dof_from_local_u < U_Space->GetActiveBound())
              {
                double v1 = basisFunctionsValues_u[l][k1]; // value of test function (vtest = vx = vy =vz)
                double v2 = v1;
                double v3 = v1;

                for (size_t k2 = 0; k2 < basisFunctionsValues_p[l].size(); k2++)
                {
                  int global_dof_from_local_p = DOF_p[k2]; // Ansatz-DOF

                  if (global_dof_from_local_p < P_Space->GetActiveBound())
                  {
                    double p = basisFunctionsValues_p[l][k2]; // value of ansatz function

                    double  n1 = normal[0];
                    double  n2 = normal[1];
                    double  n3 = normal[2];

                    // add for all three components
                    blocks[3]->add(DOF_u[k1],  DOF_p[k2], scale_factor * p * v1 * n1 ); // B1
                    blocks[7]->add(DOF_u[k1],  DOF_p[k2], scale_factor * p * v2 * n2 ); // B2
                    blocks[11]->add(DOF_u[k1], DOF_p[k2], scale_factor * p * v3 * n3 ); // B3
                  }
                } //for(l=0;l<N_BaseFunct;l++)
              }
            } // endif
          }
        }
      }
    }
  }
}

// ===========================================================================
// int_{Gamma} q*u*n
void BoundaryAssembling3D::matrix_q_u_n(BlockFEMatrix &M,
    const TFESpace3D *U_Space,
    const TFESpace3D *P_Space,
    std::vector<TBaseCell*> &boundaryCells,
    int componentID,
    double mult)
{
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */
  ///@todo create a boundary cell list
  // loop over all boundary cells (at the moment it consists of all cells!!)
  for (size_t i = 0; i < boundaryCells.size(); i++)
  {
    TBaseCell* cell = boundaryCells[i];
    // mapping from local (cell) DOF to global DOF
    int *DOF_u = U_Space->GetGlobalDOF(cell->GetCellIndex());
    // mapping from local (cell) DOF to global DOF
    int *DOF_p = P_Space->GetGlobalDOF(cell->GetCellIndex());

    // loop over the faces (joints)
    for (size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      // select boundary faces
      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        //Output::print(boundface->GetBoundComp()->get_physical_id()," comp, ",componentID);

        if (boundface->GetBoundComp()->get_physical_id()==componentID)
        {
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights_u, qPointsT_u, qPointsS_u, qWeights_p, qPointsT_p, qPointsS_p; // S,T for parametrization in 3D
          std::vector< std::vector<double> > basisFunctionsValues_u, basisFunctionsValues_p;
          this->getQuadratureData(U_Space, cell,joint_id,
              qWeights_u, qPointsT_u, qPointsS_u, basisFunctionsValues_u);
          this->getQuadratureData(P_Space, cell,joint_id,
              qWeights_p, qPointsT_p, qPointsS_p, basisFunctionsValues_p);

          std::vector<double> normal;
          double transformationDeterminant;
          this->computeNormalAndTransformationData(cell, joint_id, normal,
              transformationDeterminant);

          // rescale local integral by mesh-size (important for Nitsche boundary)
          for (size_t l = 0; l < qWeights_u.size(); l++)
          {
            double scale_factor  = mult * qWeights_u[l] * transformationDeterminant;

            for (size_t k1 = 0; k1 < basisFunctionsValues_p[l].size(); k1++)
            {
              int global_dof_from_local_p = DOF_p[k1]; // Test-DOF

              if (global_dof_from_local_p < U_Space->GetActiveBound())
              {
                double q = basisFunctionsValues_p[l][k1]; // value of test function (vtest = vx = vy =vz)

                for (size_t k2 = 0; k2 < basisFunctionsValues_u[l].size(); k2++)
                {
                  int global_dof_from_local_u = DOF_u[k2]; // Ansatz-DOF

                  if (global_dof_from_local_u < U_Space->GetActiveBound())
                  {
                    double u1 = basisFunctionsValues_u[l][k2]; // value of ansatz function
                    double u2 = u1; // value of ansatz function
                    double u3 = u1; // value of ansatz function

                    double  n1 = normal[0];
                    double  n2 = normal[1];
                    double  n3 = normal[2];

                    // add for all three components
                    blocks[12]->add(DOF_u[k1],  DOF_p[k2], scale_factor * q * u1 * n1 ); // B1
                    blocks[13]->add(DOF_u[k1],  DOF_p[k2], scale_factor * q * u2 * n2 ); // B2
                    blocks[14]->add(DOF_u[k1], DOF_p[k2], scale_factor * q * u3 * n3 ); // B3
                  }
                } //for(l=0;l<N_BaseFunct;l++)
              }
            } // endif
          }
        }
      }
    }
  }
}


// ===========================================================================
// int_{Gamma} q*u*n
void BoundaryAssembling3D::rhs_q_uD_n(BlockVector &rhs,
    const TFESpace3D *U_Space,
    const TFESpace3D *P_Space,
    BoundValueFunct3D *given_boundary_data0,
    BoundValueFunct3D *given_boundary_data1,
    BoundValueFunct3D *given_boundary_data2,
    std::vector<TBaseCell*> &boundaryCells,
    int componentID,
    double mult)
{
  //    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */
  ///@todo create a boundary cell list
  // loop over all boundary cells (at the moment it consists of all cells!!)
  for (size_t i = 0; i< boundaryCells.size(); i++)
  {
    TBaseCell* cell = boundaryCells[i];
    //        // mapping from local (cell) DOF to global DOF
    //        int *DOF_u = U_Space->GetGlobalDOF(cell->GetCellIndex());
    // mapping from local (cell) DOF to global DOF
    int *DOF_p = P_Space->GetGlobalDOF(cell->GetCellIndex());

    // loop over the faces (joints)
    for (size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      // select boundary faces
      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        //Output::print(boundface->GetBoundComp()->get_physical_id()," comp, ",componentID);

        if (boundface->GetBoundComp()->get_physical_id() == componentID)
        {
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights_u, qPointsT_u, qPointsS_u, qWeights_p,qPointsT_p,qPointsS_p; // S,T for parametrization in 3D
          std::vector< std::vector<double> > basisFunctionsValues_u, basisFunctionsValues_p;
          this->getQuadratureData(U_Space, cell, joint_id,
              qWeights_u, qPointsT_u, qPointsS_u, basisFunctionsValues_u);
          this->getQuadratureData(P_Space, cell, joint_id,
              qWeights_p, qPointsT_p, qPointsS_p, basisFunctionsValues_p);

          std::vector<double> normal;
          double transformationDeterminant;
          this->computeNormalAndTransformationData(cell, joint_id, normal,
              transformationDeterminant);

          double x, y, z;
          std::vector<double> uDirichlet(3);

          // rescale local integral by mesh-size (important for Nitsche boundary)
                        for (size_t l = 0; l < qWeights_u.size(); l++)
                        {
                          // NEW LB 11.10.18
                          // get the boundary values of rhs
                          boundface->GetBoundComp()->GetXYZofTS(qPointsT_u[l], qPointsS_u[l], x, y, z);

                          if(given_boundary_data0 != nullptr)
                            given_boundary_data0(x, y, z, uDirichlet[0]);
                          else
                            uDirichlet[0] = 0.;

                          if(given_boundary_data1 != nullptr)
                            given_boundary_data1(x, y, z, uDirichlet[1]);
                          else
                            uDirichlet[1] = 0.;

                          if(given_boundary_data2 != nullptr)
                            given_boundary_data2(x, y, z, uDirichlet[2]);
                          else
                            uDirichlet[2] = 0.;

                          double scale_factor  = mult * qWeights_u[l] * transformationDeterminant;

                          for (size_t k1 = 0; k1 < basisFunctionsValues_p[l].size(); k1++)
                          {
                            int global_dof_from_local_p = DOF_p[k1]; // Test-DOF

                            if (global_dof_from_local_p < P_Space->GetActiveBound())
                            {
                              double q = basisFunctionsValues_p[l][k1]; // value of test function (vtest = vx = vy =vz)

                              double  n1 = normal[0];
                              double  n2 = normal[1];
                              double  n3 = normal[2];

                              // add for all three components
                              // rhs.block(0)[global_dof_from_local_p] += scale_factor * q * uDirichlet[0] * n1;
                              // rhs.block(1)[global_dof_from_local_p] += scale_factor * q * uDirichlet[1] * n2;
                              // rhs.block(2)[global_dof_from_local_p] += scale_factor * q * uDirichlet[2] * n3;
                              rhs.block(2)[global_dof_from_local_p] += scale_factor * q * uDirichlet[0] * n1;
                              rhs.block(2)[global_dof_from_local_p] += scale_factor * q * uDirichlet[1] * n2;
                              rhs.block(2)[global_dof_from_local_p] += scale_factor * q * uDirichlet[2] * n3;

                              //todo: check if scaling with transformation size is necessary (see Boundary_Assembling_2D.C)
                            } //for(l=0;l<N_BaseFunct;l++)
                          }
                        } // endif
        }
      }
    }
  }
}



// ===========================================================================
// int_{Gamma} mult * ( (grad u).n, v )
void BoundaryAssembling3D::matrix_gradu_n_v(BlockFEMatrix &M,
    const TFESpace3D *U_Space,
    std::vector<TBaseCell*> &boundaryCells,
    int componentID,
    double mult)

{
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

  ///@todo create a boundary cell list
  // loop over all boundary cells (at the moment it consists of all cells!!)
  for (size_t i = 0; i < boundaryCells.size(); i++)
  {
    TBaseCell* cell = boundaryCells[i];

    // mapping from local (cell) DOF to global DOF
    int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    // loop over the faces (joints)
    for (size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      // select boundary faces
      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        //Output::print(boundface->GetBoundComp()->get_physical_id()," comp, ",componentID);

        if (boundface->GetBoundComp()->get_physical_id() == componentID)
        {
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights, qPointsT, qPointsS; // S,T for parametrization in 3D
          std::vector< std::vector<double> > basisFunctionsValues, basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z;
          this->getQuadratureDataIncludingFirstDerivatives(U_Space, cell,joint_id,
              qWeights, qPointsT, qPointsS,basisFunctionsValues,
              basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z );
          std::vector<double> normal;
          double transformationDeterminant;
          this->computeNormalAndTransformationData(cell, joint_id, normal,
              transformationDeterminant);

          // rescale local integral by mesh-size (important for Nitsche boundary)
                        for (size_t l = 0; l < qWeights.size(); l++)
                        {
                          double scale_factor = mult * qWeights[l] * transformationDeterminant;

                          for (size_t k1 = 0; k1 < basisFunctionsValues[l].size(); k1++)
                          {
                            int global_dof_from_local = DOF[k1]; // Test-DOF

                            if (global_dof_from_local < U_Space->GetActiveBound())
                            {
                              double v1 = basisFunctionsValues[l][k1]; // value of test function (vtest = vx = vy =vz)
                              double v2 = v1;
                              double v3 = v1;

                              for (size_t k2 = 0; k2 < basisFunctionsValues[l].size(); k2++)
                              {
                                int global_dof_from_local = DOF[k2]; // Ansatz-DOF

                                if (global_dof_from_local < U_Space->GetActiveBound())
                                {
                                  //double u_x = basisFunctionsValues[l][k2]; // value of test function (vtest = vx = vy =vz)
                                  //double u_y = u_x;
                                  //double u_z = u_x;

                                  double  n1=normal[0];
                                  double  n2=normal[1];
                                  double  n3=normal[2];

                                  double u_dx = basisFunctionsValues_derivative_x[l][k2];
                                  double u_dy = basisFunctionsValues_derivative_x[l][k2];
                                  double u_dz = basisFunctionsValues_derivative_x[l][k2];

                                  // (see the note about blocks at the beginning of the function)
                                  blocks[0]->add(DOF[k1], DOF[k2],  scale_factor * v1 * u_dx * n1 ); // A11
                                  blocks[0]->add(DOF[k1], DOF[k2],  scale_factor * v1 * u_dy * n2 ); // A11
                                  blocks[0]->add(DOF[k1], DOF[k2],  scale_factor * v1 * u_dz * n3 ); // A11

                                  blocks[5]->add(DOF[k1], DOF[k2],  scale_factor * v2 * u_dx * n1 ); // A22
                                  blocks[5]->add(DOF[k1], DOF[k2],  scale_factor * v2 * u_dy * n2 ); // A22
                                  blocks[5]->add(DOF[k1], DOF[k2],  scale_factor * v2 * u_dz * n3 ); // A22

                                  blocks[10]->add(DOF[k1], DOF[k2], scale_factor * v3 * u_dx * n1 ); // A33
                                  blocks[10]->add(DOF[k1], DOF[k2], scale_factor * v3 * u_dy * n2 ); // A33
                                  blocks[10]->add(DOF[k1], DOF[k2], scale_factor * v3 * u_dz * n3 ); // A33
                                }
                              } //for(l=0;l<N_BaseFunct;l++)
                            }
                          } // endif
                        }
        }
      }
    }
  }
}

// ===========================================================================
// int_{Gamma} mult*grad v*n*u
void BoundaryAssembling3D::matrix_gradv_n_u(BlockFEMatrix &M,
    const TFESpace3D *U_Space,
    std::vector<TBaseCell*> &boundaryCells,
    int componentID,
    double mult)

{
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

  ///@todo create a boundary cell list
  // loop over all boundary cells (at the moment it consists of all cells!!)
  for (size_t i = 0; i< boundaryCells.size(); i++)
  {
    TBaseCell* cell = boundaryCells[i];

    // mapping from local (cell) DOF to global DOF
    int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());
    // loop over the faces (joints)
    for (size_t joint_id = 0; joint_id< (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      // select boundary faces
      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        //Output::print(boundface->GetBoundComp()->get_physical_id()," comp, ",componentID);

        if (boundface->GetBoundComp()->get_physical_id() == componentID)
        {
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights, qPointsT, qPointsS; // S,T for parametrization in 3D
          std::vector< std::vector<double> > basisFunctionsValues, basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z ;
          this->getQuadratureDataIncludingFirstDerivatives(U_Space, cell,joint_id,
              qWeights, qPointsT, qPointsS,
              basisFunctionsValues,
              basisFunctionsValues_derivative_x,
              basisFunctionsValues_derivative_y,
              basisFunctionsValues_derivative_z );
          std::vector<double> normal;
          double transformationDeterminant;
          this->computeNormalAndTransformationData(cell, joint_id, normal,
              transformationDeterminant);

          // rescale local integral by mesh-size (important for Nitsche boundary)
                        for (size_t l = 0 ; l < qWeights.size(); l++)
                        {
                          double scale_factor = mult * qWeights[l] * transformationDeterminant;

                          for (size_t k1 = 0; k1 < basisFunctionsValues[l].size(); k1++)
                          {
                            int global_dof_from_local = DOF[k1]; // Test-DOF

                            if (global_dof_from_local < U_Space->GetActiveBound())
                            {
                              //double v_x = basisFunctionsValues[l][k1]; // value of test function (vtest = vx = vy =vz)
                              //double v_y = v_x;
                              //double v_z = v_x;

                              double v_dx = basisFunctionsValues_derivative_x[l][k1];
                              double v_dy = basisFunctionsValues_derivative_x[l][k1];
                              double v_dz = basisFunctionsValues_derivative_x[l][k1];

                              for(size_t k2=0;k2<basisFunctionsValues[l].size();k2++) {
                                int global_dof_from_local = DOF[k2]; // Ansatz-DOF

                                if(global_dof_from_local < U_Space->GetActiveBound()) {
                                  double u1 = basisFunctionsValues[l][k2]; // value of test function (vtest = vx = vy =vz)
                                  double u2 = u1;
                                  double u3 = u1;

                                  double  n1=normal[0];
                                  double  n2=normal[1];
                                  double  n3=normal[2];

                                  // (see the note about blocks at the beginning of the function)
                                  blocks[0]->add(DOF[k1], DOF[k2], scale_factor * u1 * v_dx * n1 ); // A11
                                  blocks[0]->add(DOF[k1], DOF[k2], scale_factor * u1 * v_dy * n2); // A11
                                  blocks[0]->add(DOF[k1], DOF[k2], scale_factor * u1 * v_dz * n3 ); // A11

                                  blocks[5]->add(DOF[k1], DOF[k2], scale_factor * u2 * v_dx * n1 ); // A22
                                  blocks[5]->add(DOF[k1], DOF[k2], scale_factor * u2 * v_dy * n2 ); // A22
                                  blocks[5]->add(DOF[k1], DOF[k2], scale_factor * u2 * v_dz * n3 ); // A22

                                  blocks[10]->add(DOF[k1], DOF[k2], scale_factor * u3 * v_dx * n1 ); // A33
                                  blocks[10]->add(DOF[k1], DOF[k2], scale_factor * u3 * v_dy * n2 ); // A33
                                  blocks[10]->add(DOF[k1], DOF[k2], scale_factor * u3 * v_dz * n3 ); // A33
                                }
                              } //for(l=0;l<N_BaseFunct;l++)
                            }
                          } // endif
                        }
        }
      }
    }
  }
}


// ===========================================================================
// int_{Gamma} mult*grad v*n*uD
void BoundaryAssembling3D::rhs_gradv_n_uD(BlockVector &rhs,
    const TFESpace3D *U_Space,
    BoundValueFunct3D *given_boundary_data0,
    BoundValueFunct3D *given_boundary_data1,
    BoundValueFunct3D *given_boundary_data2,
    std::vector<TBaseCell*> &boundaryCells,
    int componentID,
    double mult)

{
  //    std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

  ///@todo create a boundary cell list
  // loop over all boundary cells (at the moment it consists of all cells!!)
  for (size_t i = 0; i < boundaryCells.size(); i++)
  {
    TBaseCell* cell = boundaryCells[i];

    // mapping from local (cell) DOF to global DOF
    int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());
    // loop over the faces (joints)
    for (size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      // select boundary faces
      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        if (boundface->GetBoundComp()->get_physical_id() == componentID)
        {
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights, qPointsT, qPointsS; // S,T for parametrization in 3D
          std::vector< std::vector<double> > basisFunctionsValues, basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z ;
          this->getQuadratureDataIncludingFirstDerivatives(U_Space, cell, joint_id,
              qWeights, qPointsT, qPointsS,
              basisFunctionsValues, basisFunctionsValues_derivative_x, basisFunctionsValues_derivative_y, basisFunctionsValues_derivative_z );
          std::vector<double> normal;
          double transformationDeterminant;
          this->computeNormalAndTransformationData(cell, joint_id, normal,
              transformationDeterminant);

          std::vector<double> uDirichlet(3);
          double x, y, z;

          // loop over Gauss points
          for (size_t l = 0; l < qWeights.size(); l++)
          {
            // get the boundary values of rhs
            boundface->GetBoundComp()->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);

            // NEW LB 11.10.18
            if(given_boundary_data0 != nullptr)
              given_boundary_data0(x, y, z, uDirichlet[0]);
            else
              uDirichlet[0] = 0.;

            if(given_boundary_data1 != nullptr)
              given_boundary_data1(x, y, z, uDirichlet[1]);
            else
              uDirichlet[1] = 0.;

            if(given_boundary_data2 != nullptr)
              given_boundary_data2(x, y, z, uDirichlet[2]);
            else
              uDirichlet[2] = 0.;

            // rescale local integral by mesh-size (important for Nitsche boundary)
            double scale_factor = mult * qWeights[l] * transformationDeterminant;

            for (size_t k1 = 0; k1 < basisFunctionsValues[l].size(); k1++)
            {
              int global_dof_from_local = DOF[k1]; // Test-DOF

              if (global_dof_from_local < U_Space->GetActiveBound())
              {
                //double v_x = basisFunctionsValues[l][k1]; // value of test function (vtest = vx = vy =vz)
                //double v_y = v_x;
                //double v_z = v_x;

                double v_dx = basisFunctionsValues_derivative_x[l][k1];
                double v_dy = basisFunctionsValues_derivative_x[l][k1];
                double v_dz = basisFunctionsValues_derivative_x[l][k1];

                double  n1 = normal[0];
                double  n2 = normal[1];
                double  n3 = normal[2];

                rhs.block(0)[global_dof_from_local] += scale_factor * uDirichlet[0] * v_dx * n1;
                rhs.block(1)[global_dof_from_local] += scale_factor * uDirichlet[1] * v_dy * n2;
                rhs.block(2)[global_dof_from_local] += scale_factor * uDirichlet[2] * v_dz * n3;
              }
            } //for(l=0;l<N_BaseFunct;l++)
          }
        } // endif
      }
    }
  }
}



// ===========================================================================
// int_{Gamma} mult * (u. v)
void BoundaryAssembling3D::matrix_u_v(BlockFEMatrix &M,
    const TFESpace3D *U_Space,
    std::vector<TBaseCell*> &boundaryCells,
    int componentID,
    double mult,
    bool rescale_by_h)
{
  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

  ///@todo create a boundary cell list
  // loop over all boundary cells (at the moment it consists of all cells!!)
  for (size_t i = 0; i < boundaryCells.size(); i++)
  {
    TBaseCell* cell = boundaryCells[i];

    // mapping from local (cell) DOF to global DOF
    int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    // loop over the faces (joints)
    for (size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      // select boundary faces
      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        if (boundface->GetBoundComp()->get_physical_id() == componentID)
        {
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights, qPointsT, qPointsS; // S,T for parametrization in 3D
          std::vector< std::vector<double> > basisFunctionsValues;
          this->getQuadratureData(U_Space, cell, joint_id,
              qWeights, qPointsT, qPointsS, basisFunctionsValues);
          std::vector<double> normal;
          double transformationDeterminant;
          this->computeNormalAndTransformationData(cell, joint_id, normal,
              transformationDeterminant);

          // rescale local integral by mesh-size (important for Nitsche boundary)
                        for (size_t l = 0; l < qWeights.size(); l++)
                        {
                          double scale_factor = mult * qWeights[l] * transformationDeterminant;

                          if (rescale_by_h)
                          {
                            double h = cell->Get_hK(0);
                            scale_factor = scale_factor / h;
                          }

                          for (size_t k1 = 0; k1 < basisFunctionsValues[l].size(); k1++)
                          {
                            int global_dof_from_local = DOF[k1]; // Test-DOF

                            if (global_dof_from_local < U_Space->GetActiveBound())
                            {
                              double v1 = basisFunctionsValues[l][k1]; // value of test function (vtest = vx = vy =vz)
                              double v2 = v1;
                              double v3 = v1;

                              for (size_t k2 = 0; k2 < basisFunctionsValues[l].size(); k2++)
                              {
                                int global_dof_from_local = DOF[k2]; // Ansatz-DOF

                                if (global_dof_from_local < U_Space->GetActiveBound())
                                {
                                  double u1 = basisFunctionsValues[l][k2]; // value of test function (vtest = vx = vy =vz)
                                  double u2 = u1;
                                  double u3 = u2;

                                  // add for all three components
                                  blocks[0]->add(DOF[k1], DOF[k2], scale_factor * u1 * v1 ); // A11
                                  blocks[5]->add(DOF[k1], DOF[k2], scale_factor * u2 * v2 ); // A22
                                  blocks[10]->add(DOF[k1], DOF[k2], scale_factor * u3 * v3 ); // A33
                                }
                              }
                            }
                          }
                        }
        }
      }
    }
  }
}


// ===========================================================================
// int_{Gamma} mult*given_boundary_data(x,y,z)*v
void BoundaryAssembling3D::rhs_uD_v(BlockVector &rhs,
    const TFESpace3D *U_Space,
    BoundValueFunct3D *given_boundary_data0,
    BoundValueFunct3D *given_boundary_data1,
    BoundValueFunct3D *given_boundary_data2,
    std::vector<TBaseCell*> &boundaryCells,
    int componentID,
    double mult,
    bool rescale_by_h)
{
  // Todo: create a boundary cell list
  // loop over all boundary cells (at the moment it consists of all cells!!)
  for (size_t i = 0; i < boundaryCells.size(); i++)
  {
    TBaseCell* cell = boundaryCells[i];

    // mapping from local (cell) DOF to global DOF
    int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    for (size_t joint_id = 0; joint_id < (size_t) cell->GetN_Faces(); joint_id++)
    {
      TJoint* joint = cell->GetJoint(joint_id);

      if (joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // convert the joint to an object of BoundFace type
        TBoundFace *boundface = (TBoundFace *)joint;

        if (boundface->GetBoundComp()->get_physical_id() == componentID)
        {
          // get all data necessary for computing the integral:
          // quadrature weights, points, functions values, normal, determinant
          std:: vector<double> qWeights, qPointsT, qPointsS;
          std::vector< std::vector<double> > basisFunctionsValues;
          this->getQuadratureData(U_Space, cell, joint_id, qWeights,
              qPointsT, qPointsS, basisFunctionsValues);

          std::vector<double> normal;
          double transformationDeterminant;
          this->computeNormalAndTransformationData(cell, joint_id, normal,
              transformationDeterminant);

          double x, y, z;
          std::vector<double> uDirichlet(3);

          // loop over Gauss points
          for (size_t l = 0; l < qWeights.size(); l++)
          {
            // NEW LB 11.10.18
            // get the boundary values of rhs
            boundface->GetBoundComp()->GetXYZofTS(qPointsT[l], qPointsS[l], x, y, z);

            if(given_boundary_data0 != nullptr)
              given_boundary_data0(x, y, z, uDirichlet[0]);
            else
              uDirichlet[0] = 0.;

            if(given_boundary_data1 != nullptr)
              given_boundary_data1(x, y, z, uDirichlet[1]);
            else
              uDirichlet[1] = 0.;

            if(given_boundary_data2 != nullptr)
              given_boundary_data2(x, y, z, uDirichlet[2]);
            else
              uDirichlet[2] = 0.;


            double scale_factor = mult * qWeights[l] * transformationDeterminant;

            for (size_t k = 0; k < basisFunctionsValues[l].size(); k++)
            {
              int global_dof_from_local = DOF[k];

              if (global_dof_from_local < U_Space->GetActiveBound())
              {
                double v1 = basisFunctionsValues[l][k]; // value of test function (vtest = vx = vy =vz)
                double v2 = v1;
                double v3 = v1;

                // add for all three components
                if (!rescale_by_h)
                {
                  rhs.block(0)[global_dof_from_local] += scale_factor * uDirichlet[0] * v1;
                  //rhs[0][global_dof_from_local+N_U] +=
                  rhs.block(1)[global_dof_from_local] += scale_factor * uDirichlet[1] * v2;
                  //rhs[0][global_dof_from_local+N_U] +=
                  rhs.block(2)[global_dof_from_local] += scale_factor * uDirichlet[2] * v3;
                }
                else
                {
                  double h = cell->Get_hK(0);

                  rhs.block(0)[global_dof_from_local] += ( scale_factor * uDirichlet[0] * v1) /h;
                  rhs.block(1)[global_dof_from_local] += ( scale_factor * uDirichlet[1] * v2) /h;
                  rhs.block(2)[global_dof_from_local] += ( scale_factor * uDirichlet[2] * v3) /h;
                }
              }
            }
          }
        }
      }
    }
  }
}





// ===========================================================================
void BoundaryAssembling3D::getQuadratureDataIncludingFirstDerivatives(
    const TFESpace3D *fespace,TBaseCell *cell, int m,
    std::vector<double>& qWeights,std::vector<double>& qPointsT,
    std::vector<double>& qPointsS,
    std::vector< std::vector<double> >& basisFunctionsValues,
    std::vector< std::vector<double> >& basisFunctionsValues_derivative_x,
    std::vector< std::vector<double> >& basisFunctionsValues_derivative_y,
    std::vector< std::vector<double> >& basisFunctionsValues_derivative_z)
{
  int nFaceVertices = getNumberOfFaceVertices(cell,m);
  // set quadrature formula and compute quadrature info
  FE3D FEId = fespace->GetFE3D(cell->GetCellIndex(),cell);
  int fe_degree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);

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
  int N_Points;
  double* faceWeights;
  double *t,*s;
  // get a quadrature formula good enough for the velocity FE space
  TQuadFormula2D *qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
  qf2->GetFormulaData(N_Points, faceWeights, t, s);
  // ====================================
  qWeights.resize(N_Points);
  qPointsT.resize(N_Points);
  qPointsS.resize(N_Points);
  for (size_t k=0; k<(size_t)N_Points; k++) {
    qWeights[k] = faceWeights[k];
    qPointsT[k] = t[k];
    qPointsS[k] = s[k];
  }

  // ====================================
  // generate data on reference mesh cell for the 2d face of 3d cell
  TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)
  ->MakeRefElementData(FaceQuadFormula);

  BaseFunct3D *BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  int* N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();


  // values of base functions in all quadrature points on face
  double **JointValues = TFEDatabase3D::GetJointValues3D
      (BaseFuncts[FEId], FaceQuadFormula, m);



  double **JointValues_derivative_xi = TFEDatabase3D::GetJointDerivatives3D
      (BaseFuncts[FEId], FaceQuadFormula, m, D100);

  double **JointValues_derivative_eta = TFEDatabase3D::GetJointDerivatives3D
      (BaseFuncts[FEId], FaceQuadFormula, m, D010);

  double **JointValues_derivative_rho = TFEDatabase3D::GetJointDerivatives3D
      (BaseFuncts[FEId], FaceQuadFormula, m, D001);



  TFEDatabase3D::GetBaseFunct3D(BaseFuncts[FEId])->ChangeBF(fespace->GetCollection(), cell, qWeights.size(), JointValues);
  TFEDatabase3D::GetBaseFunct3D(BaseFuncts[FEId])->ChangeBF(fespace->GetCollection(), cell, qWeights.size(), JointValues_derivative_xi);
  TFEDatabase3D::GetBaseFunct3D(BaseFuncts[FEId])->ChangeBF(fespace->GetCollection(), cell, qWeights.size(), JointValues_derivative_eta);
  TFEDatabase3D::GetBaseFunct3D(BaseFuncts[FEId])->ChangeBF(fespace->GetCollection(), cell, qWeights.size(), JointValues_derivative_rho);

  // convert the double** to a vector
  basisFunctionsValues.resize(qWeights.size());
  basisFunctionsValues_derivative_x.resize(qWeights.size());
  basisFunctionsValues_derivative_y.resize(qWeights.size());
  basisFunctionsValues_derivative_z.resize(qWeights.size());

  for (unsigned int l=0; l<qWeights.size(); l++) {
    basisFunctionsValues[l].resize(N_BaseFunct[FEId]);
    basisFunctionsValues_derivative_x[l].resize(N_BaseFunct[FEId]);
    basisFunctionsValues_derivative_y[l].resize(N_BaseFunct[FEId]);
    basisFunctionsValues_derivative_z[l].resize(N_BaseFunct[FEId]);

    for (unsigned int k=0; k<basisFunctionsValues[l].size(); k++) {
      basisFunctionsValues[l][k]=JointValues[l][k];
      basisFunctionsValues_derivative_x[l][k]=JointValues_derivative_xi[l][k];
      basisFunctionsValues_derivative_y[l][k]=JointValues_derivative_eta[l][k];
      basisFunctionsValues_derivative_z[l][k]=JointValues_derivative_rho[l][k];
    }
  }
}

// ===========================================================================
void BoundaryAssembling3D::computeNormalAndTransformationData(TBaseCell *cell, int m,
    std::vector<double>& normal,
    double &transformationDeterminant)
{
  const int *faceVertexMap, *faceVertexMapLength;
  int maxNVerticesPerFace;
  // For the current cell, get information of faces and local vertices
  // faceVertexMap should be seen as an array of arrays, e.g.
  // faceVertexMap = { {a,b,c},{b,c,d},{a,c,d},{a,b,d}}
  // where faceVertexMap[i] contains the id of vertices defining face i
  // faceVertexMapLength is an array specifying the length of each list
  // note: in the case that faces of an element have differennt number of
  // vertices (e.g. a pyramid), the faceVertexMap lists have all lenght equal to
  // maxNVerticesPerFace, and these are filled with 0 for the faces with less vertices
  cell->GetShapeDesc()->GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
  // simplify: number of vertices on face m (m=joint_id)
  size_t nFaceVertices = faceVertexMapLength[ m ];
  std::vector< Point > faceVertices(nFaceVertices,Point((unsigned int) 3));
  for (size_t l1=0; l1<nFaceVertices; l1++) {
    double _x,_y,_z;
    cell->GetVertex(faceVertexMap[ m*maxNVerticesPerFace+l1 ])->GetCoords(_x,_y,_z);
    faceVertices[l1].x() = _x;
    faceVertices[l1].y() = _y;
    faceVertices[l1].z() = _z;
  }

  normal.resize(3);
  double xc1, yc1, zc1, xc2, yc2, zc2;
  switch(faceVertices.size()) {
  case 3:
    // compute the 2 vectors that span the plane containing the current face
    xc1 = faceVertices[1].x() - faceVertices[0].x();
    xc2 = faceVertices[2].x() - faceVertices[0].x();

    yc1 = faceVertices[1].y() - faceVertices[0].y();
    yc2 = faceVertices[2].y() - faceVertices[0].y();

    zc1 = faceVertices[1].z() - faceVertices[0].z();
    zc2 = faceVertices[2].z() - faceVertices[0].z();

    // plane spanned by vectors v1=(xc1, yc1, zc1) and v2=(xc2, yc2, zc2)
    // Area of the triangle: 0.5*||v1 x v2||
    // normed Normal vector = v1 x v2/||v1 x v2||
    // Area of reference triangle (0,0)-(0,1)-(1,0): 1/2*g*h=0.5
    // Determinant of tranform.: A(triangle)/A(ref triangle) = ||v1 x v2||
    normal[0] = yc1*zc2 - zc1*yc2;
    normal[1] = zc1*xc2 - xc1*zc2;
    normal[2] = xc1*yc2 - yc1*xc2;
    // determinant of reference trafo in order to get a normed normal vector
    transformationDeterminant =
        sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    normal[0] /= transformationDeterminant;
    normal[1] /= transformationDeterminant;
    normal[2] /= transformationDeterminant;

    break;

  case 4:
    // We consider a quadrilateral (P0,P1,P2,P3) as composed by 2 triangles
    // T1: P0,P1,P2
    // T2: P2,P3,P0
    // and we do the same as above (twice)
    // normed normal: ( (P1-P0) x (P2-P0) ) / || (P1-P0) x (P2-P0) ||
    // area: || (P1-P0) x (P2-P0) || / 2 + || (P3-P2) x (P0-P2) || / 2
    // area reference element [-1,1]x[-1,1]: 4
    // first triangle
    xc1 = faceVertices[1].x() - faceVertices[0].x();
    xc2 = faceVertices[2].x() - faceVertices[0].x();

    yc1 = faceVertices[1].y() - faceVertices[0].y();
    yc2 = faceVertices[2].y() - faceVertices[0].y();

    zc1 = faceVertices[1].z() - faceVertices[0].z();
    zc2 = faceVertices[2].z() - faceVertices[0].z();

    // normal vector (the same (except for length) for T1 and T2)
    normal[0] = yc1*zc2 - zc1*yc2;
    normal[1] = zc1*xc2 - xc1*zc2;
    normal[2] = xc1*yc2 - yc1*xc2;

    // determinant of reference transformation in order to get a normed normal vector
    double areaT1 =
        sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    areaT1 /= 2.0;
    // second triangle
    xc1 = faceVertices[3].x() - faceVertices[2].x();
    xc2 = faceVertices[0].x() - faceVertices[2].x();

    yc1 = faceVertices[3].y() - faceVertices[2].y();
    yc2 = faceVertices[0].y() - faceVertices[2].y();

    zc1 = faceVertices[3].z() - faceVertices[2].z();
    zc2 = faceVertices[0].z() - faceVertices[2].z();


    // normal vector (the same (except for length) for T1 and T2)
    normal[0] = yc1*zc2 - zc1*yc2;
    normal[1] = zc1*xc2 - xc1*zc2;
    normal[2] = xc1*yc2 - yc1*xc2;

    // determinant of reference trasformation in order to get a normed normal vector
    double areaT2 =
        sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    normal[0] /= areaT2;
    normal[1] /= areaT2;
    normal[2] /= areaT2;

    areaT2 /= 2.0;

    // note: the reference element is [-1,1] x [-1,1]
    transformationDeterminant = (areaT1+areaT2)/4.;

    break;

  } // tria or quads

}

//// ===========================================================================
/////@todo This function should be a member of Boundface.C
//void BoundaryAssembling3D::compute_h(TBaseCell *cell,
//                                     int m,
//                                     double &h)
//{
//    const int *faceVertexMap, *faceVertexMapLength;
//    int maxNVerticesPerFace;
//    // For the current cell, get information of faces and local vertices
//    // faceVertexMap should be seen as an array of arrays, e.g.
//    // faceVertexMap = { {a,b,c},{b,c,d},{a,c,d},{a,b,d}}
//    // where faceVertexMap[i] contains the id of vertices defining face i
//    // faceVertexMapLength is an array specifying the length of each list
//    // note: in the case that faces of an element have differennt number of
//    // vertices (e.g. a pyramid), the faceVertexMap lists have all lenght equal to
//    // maxNVerticesPerFace, and these are filled with 0 for the faces with less vertices
//    cell->GetShapeDesc()->GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
//    // simplify: number of vertices on face m (m=joint_id)
//    size_t nFaceVertices = faceVertexMapLength[ m ];
//    std::vector< Point > faceVertices(nFaceVertices,Point((unsigned int) 3));
//    for (size_t l1=0; l1<nFaceVertices; l1++)
//    {
//        double _x,_y,_z;
//        cell->GetVertex(faceVertexMap[ m*maxNVerticesPerFace+l1 ])->GetCoords(_x,_y,_z);
//        faceVertices[l1].x() = _x;
//        faceVertices[l1].y() = _y;
//        faceVertices[l1].z() = _z;
//    }
//    
//    std::vector<double> edge_lengths;
//    double min_edge_length, max_edge_length;
//    double edge_length1, edge_length2, edge_length3, edge_length4;
//    switch(faceVertices.size()) {
//       case 3:
//            // we consider a triangle (P0,P1,P2)
//            // and compute the edge lengths ||P1-P0||_2, ||P2-P0||_2, ||P1-P2||_2
//            edge_lengths[0] = sqrt  ((faceVertices[1].x() - faceVertices[0].x())*(faceVertices[1].x() - faceVertices[0].x())+
//                                     (faceVertices[1].y() - faceVertices[0].y())*(faceVertices[1].y() - faceVertices[0].y())+
//                                     (faceVertices[1].z() - faceVertices[0].z())*(faceVertices[1].z() - faceVertices[0].z())  );
//            
//            edge_lengths[1]= sqrt  ((faceVertices[2].x() - faceVertices[0].x())*(faceVertices[2].x() - faceVertices[0].x())+
//                                    (faceVertices[2].y() - faceVertices[0].y())*(faceVertices[2].y() - faceVertices[0].y())+
//                                    (faceVertices[2].z() - faceVertices[0].z())*(faceVertices[2].z() - faceVertices[0].z())  );
//            
//            edge_lengths[2] = sqrt  ((faceVertices[1].x() - faceVertices[2].x())*(faceVertices[1].x() - faceVertices[2].x())+
//                                     (faceVertices[1].y() - faceVertices[2].y())*(faceVertices[1].y() - faceVertices[2].y())+
//                                     (faceVertices[1].z() - faceVertices[2].z())*(faceVertices[1].z() - faceVertices[2].z())  );
//            
//            // compute the maximum and minimum edge length
//            min_edge_length = *min_element(edge_lengths.begin(), edge_lengths.end());
//            max_edge_length = *max_element(edge_lengths.begin(), edge_lengths.end());
//            
//            // set h to be the maximum edge length
//            h=max_edge_length;
//            
//            break;
//            
//        case 4:
//            // we consider a quadrilateral (P0,P1,P2,P3)
//            // and compute the edge lengths ||P1-P0||_2, ||P2-P1||_2, ||P3-P2||_2, ||P3-P0||_2
//            edge_lengths[0] = sqrt  ((faceVertices[1].x() - faceVertices[0].x())*(faceVertices[1].x() - faceVertices[0].x())+
//                                     (faceVertices[1].y() - faceVertices[0].y())*(faceVertices[1].y() - faceVertices[0].y())+
//                                     (faceVertices[1].z() - faceVertices[0].z())*(faceVertices[1].z() - faceVertices[0].z())  );
//            
//            edge_lengths[1] = sqrt  ((faceVertices[2].x() - faceVertices[1].x())*(faceVertices[2].x() - faceVertices[1].x())+
//                                     (faceVertices[2].y() - faceVertices[1].y())*(faceVertices[2].y() - faceVertices[1].y())+
//                                     (faceVertices[2].z() - faceVertices[1].z())*(faceVertices[2].z() - faceVertices[1].z())  );
//            
//            edge_lengths[2] = sqrt  ((faceVertices[3].x() - faceVertices[2].x())*(faceVertices[3].x() - faceVertices[2].x())+
//                                     (faceVertices[3].y() - faceVertices[2].y())*(faceVertices[3].y() - faceVertices[2].y())+
//                                     (faceVertices[3].z() - faceVertices[2].z())*(faceVertices[3].z() - faceVertices[2].z())  );
//            
//            edge_lengths[3] = sqrt  ((faceVertices[3].x() - faceVertices[0].x())*(faceVertices[3].x() - faceVertices[0].x())+
//                                     (faceVertices[3].y() - faceVertices[0].y())*(faceVertices[3].y() - faceVertices[0].y())+
//                                     (faceVertices[3].z() - faceVertices[0].z())*(faceVertices[3].z() - faceVertices[0].z())  );
//            
//            // compute the maximum and minimum edge length
//            min_edge_length = *min_element(edge_lengths.begin(), edge_lengths.end());
//            max_edge_length = *max_element(edge_lengths.begin(), edge_lengths.end());
//            
//            // set h to be the maximum edge length
//            h=max_edge_length;
//            
//            break;
//            
//    } // triangles or quadrilaterals
//    
//}


// ===========================================================================
///@todo this functions should belong to TBaseCell
int BoundaryAssembling3D::getNumberOfFaceVertices(TBaseCell *cell, int m)
{
  const int *faceVertexMap, *faceVertexMapLength;
  int maxNVerticesPerFace;
  cell->GetShapeDesc()->GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
  // simplify: number of vertices on face m
  return faceVertexMapLength[ m ];
}

// ===========================================================================


void BoundaryAssembling3D::nitsche_bc(BlockFEMatrix &s_matrix, BlockVector &s_rhs,
    const TFESpace3D * v_space, const TFESpace3D *p_space,
    BoundValueFunct3D * U1, BoundValueFunct3D *U2, BoundValueFunct3D *U3,
    std::vector<TBaseCell*> &boundaryCells,
    int bd_comp, double gamma, double mu,
    int sym_u, int sym_p)
{
  // gamma/h (u,v)
  matrix_u_v(s_matrix, v_space, boundaryCells, bd_comp, gamma*mu, true);  // rescale local integral by edge values

  // gamma/h (uD,v) [rhs]
  rhs_uD_v(s_rhs, v_space, U1, U2, U3, boundaryCells, bd_comp, gamma*mu, true);   // rescale local integral by edge values

  // - (mu grad(u)n,v)
  matrix_gradu_n_v(s_matrix, v_space, boundaryCells, bd_comp, -1. * mu);  // OK

  // - sign_u * (u, mu grad(v)n) [sign_u=1: symmetrix, -1: skew-symmetric]
  matrix_gradv_n_u(s_matrix, v_space, boundaryCells, bd_comp, (-1) * sym_u * mu);

  // - sign_u * (uD,mu grad(v)n) [rhs]
  rhs_gradv_n_uD(s_rhs, v_space, U1, U2, U3, boundaryCells, bd_comp, (-1) * sym_u * mu );

  // (pn,v)
  matrix_p_v_n(s_matrix, v_space, p_space, boundaryCells, bd_comp, 1.);

  // sign_p * (u,qn)
  matrix_q_u_n(s_matrix, v_space, p_space, boundaryCells, bd_comp, sym_p);

  // sign_p * (uD,qn) [rhs]
  rhs_q_uD_n(s_rhs, v_space, p_space, U1, U2, U3, boundaryCells, bd_comp, sym_p);
}



