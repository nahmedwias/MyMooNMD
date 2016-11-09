// ======================================================================
// @(#)BoundaryAssembling3D.C        28.10.16
//
// Functions for (external and internal) boundary integral
//
// ======================================================================

#include <FEDatabase3D.h>
#include <TriaAffin.h>
#include <QuadAffin.h>
#include <BoundaryAssembling3D.h>
#include <Collection.h>
#include <BoundFace.h>


// int_{Gamma} mult*given_boundary_data(x,y,z)*<v,normal>
void BoundaryAssembling3D::rhs_g_v_n(BlockVector &rhs,
                                     const TFESpace3D *U_Space,
                                     BoundValueFunct3D *given_boundary_data,
                                     std::vector<TBaseCell*> &boundaryCells,
                                     double mult)
{
  for(size_t i=0; i< boundaryCells.size(); i++) {
    
    TBaseCell* cell = boundaryCells[i];
    // mapping from local (cell) DOF to global DOF
    int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex()); 
    
    for(size_t joint_id=0; joint_id< cell->GetN_Faces(); joint_id++) {
      TJoint* joint = cell->GetJoint(joint_id);
      
      if (joint->GetType() == BoundaryFace ||
	  joint->GetType() == IsoBoundFace) {
	
	// convert the joint to an object of BoundFace type
	//TBoundFace *boundface = (TBoundFace *)joint;
	//TBaseCell *cell =  boundface->GetNeighbour(0);
        
	// ====================================
	// get all data necessary for computing the integral:
	// quadrature weights, points, functions values, normal, determinant
	std:: vector<double> qWeights,qPointsT,qPointsS;
	std::vector< std::vector<double> > basisFunctionsValues;
  	this->getQuadratureData(U_Space, cell,joint_id,
				qWeights,qPointsT,qPointsS,basisFunctionsValues);
	
	
	std::vector<double> normal;
	double transformationDeterminant;
	this->computeNormalAndTransformationData(cell,joint_id,
						 normal,transformationDeterminant);
	// ====================================

	// loop over Gauss points
	for(size_t l=0;l<qWeights.size();l++) {
	  double value;
	  if(given_boundary_data != nullptr)
	    value = 1.;//given_boundary_data(...);
	  else
	    value = 1;
	  
	  double commonFactor = mult * qWeights[l] * transformationDeterminant;
	  for(size_t k=0;k<basisFunctionsValues[l].size();k++) {
	    int global_dof_from_local = DOF[k];
	    
	    if(global_dof_from_local < U_Space->GetActiveBound()) {
	      double v_x = basisFunctionsValues[l][k]; // value of test function (vtest = vx = vy =vz)
	      double v_y = v_x;
	      double v_z = v_x;
	      // add for both components
	      rhs.block(0)[global_dof_from_local] += commonFactor * value * (v_x*normal[0]);
	      rhs.block(1)[global_dof_from_local] += commonFactor * value * (v_y*normal[1]);
	      rhs.block(2)[global_dof_from_local] += commonFactor * value * (v_z*normal[2]);

	    }
	    
	  } // endfor k (n. basis fcts)
	} // endfor l (Gauss points)
	
      } // if boundary face
      
    } // for n. of cells
  } // loop over cells
}

void BoundaryAssembling3D::getQuadratureData(const TFESpace3D *fespace,TBaseCell *cell, int m,
					     std::vector<double>& qWeights,std::vector<double>& qPointsT,
					     std::vector<double>& qPointsS,
					     std::vector< std::vector<double> >& basisFunctionsValues)
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
  for (size_t k=0; k<N_Points; k++) {
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
  TFEDatabase3D::GetBaseFunct3D(BaseFuncts[FEId])
    ->ChangeBF(fespace->GetCollection(), cell, qWeights.size(), JointValues);
  // convert the double** to a vector
  basisFunctionsValues.resize(qWeights.size());
  for (unsigned int l=0; l<qWeights.size(); l++) {
    basisFunctionsValues[l].resize(N_BaseFunct[FEId]);
    
    for (unsigned int k=0; k<basisFunctionsValues[l].size(); k++) {
      basisFunctionsValues[l][k]=JointValues[l][k];
    }
  }
  
  
}

void BoundaryAssembling3D::computeNormalAndTransformationData(TBaseCell *cell, int m,
							      std::vector<double>& normal,
							      double &transformationDeterminant)
{
  const int *faceVertexMap, *faceVertexMapLength;
  int maxNVerticesPerFace;
  // get information of faces and local vertices
  // faceVertexMap should be seen as an array of arrays, e.g.
  // faceVertexMap = { {a,b,c},{b,c,a},{a,c,b}}
  // where faceVertexMap[i] contains the id of vertices defining face i
  // faceVertexMapLength is an array specifying the length of each list
  // note: in the case that faces of an element have differennt number of
  // vertices (e.g. a pyramid), the faceVertexMap lists have all lenght equal to
  // maxNVerticesPerFace, and these are filled with 0 for the faces with less vertices
  cell->GetShapeDesc()->
    GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
  // simplify: number of vertices on face m
  size_t nFaceVertices = faceVertexMapLength[ m ];
  std::vector< Point > faceVertices(nFaceVertices,Point((unsigned int) 3));
  for (size_t l1=0; l1<nFaceVertices; l1++) {
    double _x,_y,_z;
    cell->GetVertex(faceVertexMap[ m*maxNVerticesPerFace+l1 ])->
      GetCoords(_x,_y,_z);
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
    // Normal vector = v1 x v2/||v1 x v2||
    // Area of reference triangle: 0.5
    // Determinant of tranform: A(triangle)/A(ref triangle) = ||v1 x v2||
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
    // normal: (P1-P0) x (P2-P0)/ || ... ||
    // area: || (P1-P0) x (P2-P0) || / 2 + || (P3-P2) x (P0-P2) || / 2
    // area reference element: 4
    // first triangle
    xc1 = faceVertices[1].x() - faceVertices[0].x(); 
    xc2 = faceVertices[2].x() - faceVertices[0].x();
    
    yc1 = faceVertices[1].y() - faceVertices[0].y(); 
    yc2 = faceVertices[2].y() - faceVertices[0].y();
    
    zc1 = faceVertices[1].z() - faceVertices[0].z(); 
    zc2 = faceVertices[2].z() - faceVertices[0].z();
    
    // normal vector (the same for T1 and T2)
    normal[0] = yc1*zc2 - zc1*yc2;
    normal[1] = zc1*xc2 - xc1*zc2;
    normal[2] = xc1*yc2 - yc1*xc2;
    
    // determinant of reference trasformation in order to get a normed normal vector
    double areaT1 =
      sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2])/2.;
    
    // second triangle
    xc1 = faceVertices[3].x() - faceVertices[2].x();  
    xc2 = faceVertices[0].x() - faceVertices[2].x();
    
    yc1 = faceVertices[3].y() - faceVertices[2].y();
    yc2 = faceVertices[0].y() - faceVertices[2].y();
    
    zc1 = faceVertices[3].z() - faceVertices[2].z(); 
    zc2 = faceVertices[0].z() - faceVertices[2].z();
    
    
    // normal vector (the same for T1 and T2)
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
///@todo this functions should belong to TBaseCell
int BoundaryAssembling3D::getNumberOfFaceVertices(TBaseCell *cell, int m)
{
  const int *faceVertexMap, *faceVertexMapLength;
  int maxNVerticesPerFace;
  cell->GetShapeDesc()->
    GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
  // simplify: number of vertices on face m
  return faceVertexMapLength[ m ];
}
