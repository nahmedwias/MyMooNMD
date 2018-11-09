
#ifdef _MPI
#  include "mpi.h"
#endif

#include <BaseCell.h>
#include <Joint.h>
#include <stdlib.h>

#include <BoundComp3D.h>
#include <BoundEdge.h>
#include <BoundFace.h>
#include <Edge.h>
#include <Point.h>


// Constructor
TBaseCell::TBaseCell(TRefDesc *refdesc)
{
 int i, N_;

  RefDesc = refdesc;

  N_ = RefDesc->GetShapeDesc()->GetN_Joints();  
  Joints = new TJoint*[N_];

  for (i=0;i<N_;i++)
   Joints[i] = nullptr;

  ClipBoard = 0;
  Phase_ID = 0;
  Reference_ID = 0;
  CellIndex = -1;
  region = 0;
  LayerCell = 0;  
  
  normalOrientation = nullptr;

#ifdef __3D__   
  N_ = RefDesc->GetN_OrigEdges();
  Edges = new TEdge*[N_];
  for (i=0;i<N_;i++)
   Edges[i] = nullptr;
#endif

#ifdef  _MPI
  ClipBoard_Par  = -1;
  SubDomainNumber = 0;
  GlobalCellNo = -1;
  SubDomainLocalCellNo = -1;

  OwnCell=false;
  HaloCell=false;
  SubDomainInterfaceCell=false;
  DependentCell=false;
  N_NeibProcesses = 0;
  NeibProcessesIds = nullptr;
#endif
}

// Destructor
TBaseCell::~TBaseCell()
{
  int i, N_;
  TJoint *CurrJoint;

#ifdef __2D__
  N_ = RefDesc->GetN_OrigEdges();
#else
  N_ = RefDesc->GetN_OrigFaces();
#endif
  for (i=0;i<N_;i++)
  {
    CurrJoint = Joints[i];
    switch (CurrJoint->GetType())
    {
      case JointEqN:
      case InterfaceJoint:
      case IsoInterfaceJoint:
      case InnerInterfaceJoint:
#ifdef __3D__
      case InterfaceJoint3D:
      case IsoInterfaceJoint3D:
#endif
        if (CurrJoint->GetNeighbour(this))
	{ CurrJoint->Delete(this);}
        else
	{ delete CurrJoint;}
        break;

      default:
        delete CurrJoint;
        break;
    }
  }

  delete[] Joints;
}

double TBaseCell::Get_hK(int cell_measure) const
{
  switch (cell_measure)
  {
    case 0:                                     // diameter
      return this->GetDiameter();
      //case 1: // with reference map
      //OutPut("cell measure " << endl);
      //return this->GetLengthWithReferenceMap();
    case 2:                                     // shortest edge
      return this->GetShortestEdge();
      break;
    case 1:                                     // with reference map
    case 3:                                     // measure
      return sqrt(this->GetMeasure());
      break;
    case 4:                                     // mesh size in convection direction, this is just a dummy
      return this->GetDiameter();
      break;
    case 5:                                     // take value from an array
      // this is in general not the diameter but a pw constant value
      // which is needed for some reasons
      return this->GetDiameter();
      break;
    default:                                    // diameter
      OutPut("CELL_MEASURE " << cell_measure << " not available!!!" << endl);
      return this->GetDiameter();
      break;
  }
}

// added 25.04.2010 for fixing refinement problem
#ifdef __3D__
void TBaseCell::CorrectBoundaryVertices(TVertex **NewVertices, TJoint **NewJoints)
{
  int N_NewFaces = RefDesc->GetN_Faces();
  const int *FaceVertex;
  int MaxN_VpF;
  
  TBoundComp3D *BoundComp;
  TVertex *Vertex;
  double x, y, z, t, s;

  for (int i=0;i<N_NewFaces;++i)
  {
    if ( NewJoints[i]->GetType() == BoundaryFace )
    {
      BoundComp = ((TBoundFace*) NewJoints[i])->GetBoundComp();
      
      if ( BoundComp->GetType() == Plane ) continue;
      
      RefDesc->GetFaceVertex(FaceVertex, MaxN_VpF);
      
      for (int j=0;j<MaxN_VpF;++j)
      {
	Vertex = NewVertices[FaceVertex[MaxN_VpF*i+j]];
	
	Vertex->GetCoords(x ,y, z);
	
	BoundComp->GetTSofXYZ(x, y, z, t, s);
	BoundComp->GetXYZofTS(t, s, x, y, z);
	
	Vertex->SetCoords(x, y, z);
      }     
    }
    else continue;
  }
  
}
#endif

// on each joint, decide whether the 
// global normal is outgoing (+1) or ingoing (-1)
void TBaseCell::SetNormalOrientation()
{
  TBaseCell *neighbCell;
  int nEdges = RefDesc->GetN_OrigEdges();
  #ifdef __3D__
  nEdges = RefDesc->GetN_OrigFaces();
  #endif
  if(normalOrientation != nullptr) // nothing more to do
    return;
  normalOrientation = new int[nEdges];
  for (int i=0; i<nEdges;i++)
  {
    normalOrientation[i] = 1;
    TJoint *joint = Joints[i];
    if(joint->InnerJoint())
    {
      neighbCell = joint->GetNeighbour(this);
      if(neighbCell->GetCellIndex() < GetCellIndex()&& 
         Reference_ID == neighbCell->GetReference_ID()) 
        normalOrientation[i] = -1;
    }
  }
}

//LB ====================================================
bool TBaseCell::IsBoundaryCell( int BoundComp_id ) const
{   int num_inner_joints = 0;
    for(int j = 0;  j < this->GetN_Joints(); j++)
    {   const TJoint *joint = this->GetJoint(j);
        if( joint->InnerJoint())
        {
            num_inner_joints=num_inner_joints+1;
            continue;
        }
        else if(!(joint->InnerJoint()) && BoundComp_id == -4711)
            {
                return true;  
            }
        else if(!(joint->InnerJoint()) && BoundComp_id != -4711)
            {
#ifdef __2D__
               const TBoundEdge *boundedge = (const TBoundEdge *)joint;
                const TBoundComp *BoundComp = boundedge->GetBoundComp();
#elif __3D__
               const TBoundFace *boundface = (const TBoundFace *)joint;
                const TBoundComp *BoundComp = boundface->GetBoundComp();
#endif //__3D__
                if (BoundComp->GetID() == BoundComp_id)
                {
                    return true;
                }
            }
        if(this->GetN_Joints()-1  == num_inner_joints )
            return false;
    }
return false;   
}
//LB ====================================================

#ifdef __3D__
void TBaseCell::computeNormalAndTransformationData(int m,
					std::vector<double>& normal,
					double &transformationDeterminant) const
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
    this->RefDesc->GetShapeDesc()->GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
    // simplify: number of vertices on face m (m=joint_id)
    size_t nFaceVertices = faceVertexMapLength[ m ];
    std::vector< Point > faceVertices(nFaceVertices,Point((unsigned int) 3));
    for (size_t l1=0; l1<nFaceVertices; l1++) {
        double _x,_y,_z;
        this->GetVertex(faceVertexMap[ m*maxNVerticesPerFace+l1 ])->GetCoords(_x,_y,_z);
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
#endif

// Methods
#ifdef  _MPI
void TBaseCell::SetNeibProcessesIds(int *Neiblist)
{
 int i;

 if(N_NeibProcesses==0)
  {
   printf(" Set the N_NeibProcesses for this cell first !!!! \n");
   MPI_Finalize();
   exit(0);
  }

 if(NeibProcessesIds) delete [] NeibProcessesIds;

 NeibProcessesIds = new int[N_NeibProcesses];

 for(i=0;i<N_NeibProcesses;i++)
  NeibProcessesIds[i] = Neiblist[i];

}
#endif

double norm(double x,double y, double z)
{
	return std::sqrt(x*x+y*y+z*z);
}
#ifdef __3D__
//Returns normed triple product
double triple_product(const TVertex * vert0, const TVertex * vert1, const TVertex * vert2, const TVertex * vert3)
{
	double x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3;
	vert0->GetCoords(x0,y0,z0);
	vert1->GetCoords(x1,y1,z1);
	vert2->GetCoords(x2,y2,z2);
	vert3->GetCoords(x3,y3,z3);

	//calculate the triple product (volume of the hexahedron)
	//TODO: Normiere den Spass, damit sich ablesen laesst, ob Zelle entartet ist.
	double product =((x2-x0)*(y3-y0)-(x3-x0)*(y2-y0))*(z1-z0)
			+((x3-x0)*(y1-y0)-(x1-x0)*(y3-y0))*(z2-z0)
			+((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0))*(z3-z0);
	//double norm = norm(x1-x0,y1-0,z1-z0)*norm(x2-x0,y2-y0,z2-z0)*norm(x3-x0,y3-y0,z3-z0);
	//if (norm==0)
	//	ErrThrow("Two vertices");
	return product;
}

//Check if the edge between vert0 and vert1 and the edge between vert2 and vert3 cross
//Solve LGS
/*{ (x0-x1)  (x2-x3)	{ s		(x3-x1)
 *  (y0-y1)  (y2-y3)	  t}  =	(y3-y1)
 *  (z0-z1)  (z2-z3) }	  		(z3-z1)
 * */
//This method may be numerical instable.
bool find_cross_point(const TVertex * vert0, const TVertex * vert1,
		const TVertex * vert2, const TVertex * vert3)
		//double s,double t
{
	//Input check: Are the vertices pairwise different
	const TVertex* Vertices[4]={vert0 , vert1 , vert2 , vert3};
	for (int iter1 =0; iter1 <4; iter1++)
	{
		for (int iter2=iter1+1; iter2 <4; iter2++)
		{
			//Operator < (actually <=) compares vertices coordinate-wise.
			if ((!(*Vertices[iter1] < *Vertices[iter2])) && (!(*Vertices[iter2] < *Vertices[iter1])))
				{
				Output::print(vert0);
				Output::print(vert1);
				Output::print(vert2);
				Output::print(vert3);
				ErrThrow("Two vertices are the same");
				}
		}
	}

	//Compute cross point (Solve LGS)
	double s,t;

	double x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3;
	vert0->GetCoords(x0,y0,z0);
	vert1->GetCoords(x1,y1,z1);
	vert2->GetCoords(x2,y2,z2);
	vert3->GetCoords(x3,y3,z3);

	double a = (x0-x1);
	double b = (x2-x3);
	double x = (x3-x1);
	double c = (y0-y1);
	double d = (y2-y3);
	double y = (y3-y1);
	double e = (z0-z1);
	double f = (z2-z3);
	double z = (z3-z1);

	double det = a*d-c*b;
	if (det!=0)
	{
		s=(d*x-b*y)/det;
		t=(a*y-c*x)/det;

		if (std::abs(e*s+f*t-z)>1e-10)
			return false; //Edges are parallel
	}
	else
	{
		det = c*f-e*d;
			if (det!=0)
			{
				s=(f*y-d*z)/det;
				t=(c*z-e*y)/det;

				if (std::abs(a*s+b*t-x)>1e-10)
					return false; //Edges are parallel
			}
			else
				return false; //Edges are parallel or Vertices are the same...

	}
	if (0<=s && s<=1 && 0<=t && t<=1)
		{
		Output::print("s ",s, " t ", t);
		return true;
		}

	return false;
}



bool TBaseCell::check_orientation() const
{
	#ifdef __3D__
	if (this->GetShapeDesc()->GetType()==Tetrahedron)
	{

		//calculate the triple product (volume of the hexahedron)
		double product =triple_product(
				this->GetVertex(0),this->GetVertex(1),
				this->GetVertex(2),this->GetVertex(3));
		return (product>0);
	}
	else if(this->GetShapeDesc()->GetType()==Hexahedron || this->GetShapeDesc()->GetType()==Brick)
	{
		//calculate the triple product (volume of the hexahedron)
		double product =triple_product(
				this->GetVertex(0),this->GetVertex(1),
				this->GetVertex(2),this->GetVertex(4));
		return (product>0);
	}
	else
	{
		Output::print("cell shape: ", this->GetShapeDesc()->GetType());
		ErrThrow("check_orientation: This cell shape is not implemented!");
		return false;
	}
    #endif

}

bool TBaseCell::check_shape() const
{
	if (this->GetShapeDesc()->GetType()==Tetrahedron)
	{
		//calculate the triple product (volume of the hexahedron)
		//if it is not equal 0 then the vertices are not in a plane
		double product =triple_product(
				this->GetVertex(0),this->GetVertex(1),
				this->GetVertex(2),this->GetVertex(3));
		return (std::abs(product)>1e-10);
	}
	else if(this->GetShapeDesc()->GetType()==Hexahedron|| this->GetShapeDesc()->GetType()==Brick)
	//TODO: Bricks have still more properties to be checked... use CheckHexa()
	{
		const int *vert_per_face;
		const int *n_vert_per_face;
	    int max_n;
		this->GetShapeDesc()->GetFaceVertex(vert_per_face, n_vert_per_face, max_n);

		for(int face_id=0; face_id<6; face_id++)
		{
			int vert0_id =vert_per_face[max_n*face_id +0];
			int vert1_id =vert_per_face[max_n*face_id +1];
			int vert2_id =vert_per_face[max_n*face_id +2];
			int vert3_id =vert_per_face[max_n*face_id +3];

			double product =triple_product(
				this->GetVertex(vert0_id),this->GetVertex(vert1_id),
				this->GetVertex(vert2_id),this->GetVertex(vert3_id));
			if(std::abs(product) >1e-10)
			{
				Output::print("Computing errors may cause that the cell is considered degenerated.");//TODO:
				return false;
			}

			//check if opposite edges cross
			if (find_cross_point(
					this->GetVertex(vert0_id),this->GetVertex(vert1_id),
					this->GetVertex(vert2_id),this->GetVertex(vert3_id)))
				ErrThrow("Opposite edges of a face cross.");

			if (find_cross_point(
					this->GetVertex(vert0_id),this->GetVertex(vert3_id),
					this->GetVertex(vert1_id),this->GetVertex(vert2_id)))
				ErrThrow("Opposite edges of a face cross.");

		} //end face
	}
	else
	{
		ErrThrow("check_shape: For this cell shape has not implemented a check method!");
		return false;
	}

	return true;
}

#endif //__3D__



