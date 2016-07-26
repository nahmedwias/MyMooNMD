#include <algorithm>
#include <iostream>

#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <TetGenMeshLoader.h>
#include <unordered_map>
// #include <TetGenMeshLoader_Test.h>

#include <Output3D.h>
#include <Joint.h>
#include <JointEqN.h>

/**
 *  Calculates the bary center of a cell
 */
void baryCentricCoords(const TBaseCell &gridCell, TVertex &baryCenter)
{
  double x=0., y=0., z=0.;
  int nPoints = gridCell.GetN_Vertices();

  for(int i=0; i<nPoints; i++)
  {
      x += gridCell.GetVertex(i)->GetX();
      y += gridCell.GetVertex(i)->GetY();
      z += gridCell.GetVertex(i)->GetZ();
  }
  baryCenter.SetCoords(x/nPoints,y/nPoints,z/nPoints);
}

/**
 *  get all vertices of a face 
 */
std::array<TVertex*, 3> getFaceVertices(TBaseCell &gridCell, int faceNo)
{
  const int *tempFaceVert, *tempLength;
  int maxLength;
  TShapeDesc *shapeDesc = gridCell.GetShapeDesc();

  shapeDesc->GetFaceVertex(tempFaceVert, tempLength, maxLength);
  std::array<TVertex*,3> vertices;

  vertices[0] = gridCell.GetVertex(tempFaceVert[faceNo*maxLength+0]);
  vertices[1] = gridCell.GetVertex(tempFaceVert[faceNo*maxLength+1]);
  vertices[2] = gridCell.GetVertex(tempFaceVert[faceNo*maxLength+2]);

  return vertices;
}

/**
 *  needed sets too 
 */
std::set<const TVertex*> getFaceVertices(TBaseCell* gridCell, int faceNo)
{
  const int *tempFaceVert, *tempLength;
  int maxLength;
  TShapeDesc *shapeDesc = gridCell->GetShapeDesc();

  shapeDesc->GetFaceVertex(tempFaceVert, tempLength, maxLength);

  std::array<TVertex*,3> vertices;

  vertices[0] = gridCell->GetVertex(tempFaceVert[faceNo*maxLength+0]);
  vertices[1] = gridCell->GetVertex(tempFaceVert[faceNo*maxLength+1]);
  vertices[2] = gridCell->GetVertex(tempFaceVert[faceNo*maxLength+2]);

  std::set<const TVertex*> ret = { vertices[0], vertices[1], vertices[2] };

  return ret;
}

/**
 *  get indices of all vertices of a face 
 */
std::array<int, 3> getFaceVertexIndices(TBaseCell &gridCell, int faceNo)
{
  const int *tempFaceVert, *tempLength;
  int maxLength;
  TShapeDesc *shapeDesc = gridCell.GetShapeDesc();

  shapeDesc->GetFaceVertex(tempFaceVert, tempLength, maxLength);

  std::array<int,3> indices;
  indices[0] = tempFaceVert[faceNo*maxLength+0];
  indices[1] = tempFaceVert[faceNo*maxLength+1];
  indices[2] = tempFaceVert[faceNo*maxLength+2];
  return indices;
}
 
/**
 *  constructs a cell of tetrahedron type
 */
TBaseCell* allocTetrahedron(TVertex* a, TVertex* b, TVertex* c, 
                            TVertex* d, TRefDesc *refDesc)
{
  TBaseCell *tetrahedron = new TGridCell(refDesc, 0);
  // set the vertices
  tetrahedron->SetVertex(0,a);
  tetrahedron->SetVertex(1,b);
  tetrahedron->SetVertex(2,c);
  tetrahedron->SetVertex(3,d);

  return tetrahedron;
}

/**
 *  sets joint and neighbours of a given face. Allocates joint if not exist. Face given
 *  in terms of three vertices.
 */ 
template <typename It>
void setNeighbours (It begin, It end, const TVertex* a, const TVertex* b, const TVertex* c)
{
  std::set<const TVertex*> face = { a, b, c };
  int n_neighs=0;
  TBaseCell* neighbs[2];
  int face_numbers[2];

  for(auto i=begin;i!=end;++i)
  {
    auto nFaces = (*i)->GetN_Faces();
    for(int j=0;j<nFaces;++j)
    {
      std::set<const TVertex*> current = getFaceVertices(*i, j);

      if(std::equal(std::begin(face), std::end(face), std::begin(current)) )
      {
        assert(n_neighs<2);

        neighbs[n_neighs] = *i;
        face_numbers[n_neighs] = j;
        ++n_neighs;
      }
    }
  }

  assert (n_neighs == 2);

  TJoint* joint = nullptr;
  if(neighbs[0]->GetJoint(face_numbers[0]) != nullptr )
  {
    joint = neighbs[0]->GetJoint(face_numbers[0]);
    neighbs[1]->SetJoint(face_numbers[1], joint);
    joint->SetMapType();
  }
  else if(neighbs[1]->GetJoint(face_numbers[1]) != nullptr )
  {
    joint = neighbs[1]->GetJoint(face_numbers[1]);
    neighbs[0]->SetJoint(face_numbers[0], joint);
    joint->SetMapType();
  }
  else
  {
    joint = new TJointEqN(neighbs[0], neighbs[1]);
    neighbs[0]->SetJoint(face_numbers[0], joint);
    neighbs[1]->SetJoint(face_numbers[1], joint);
    joint->SetMapType();
  }
}

/**
 *  Does bary centric refinement of one tetrahedron.
 *  Return newly created cell.
 */ 
std::array<TBaseCell*, 4> baryCentricRefine(TBaseCell &oldTetra)
{
  assert(oldTetra.GetType() == Tetrahedron);

  constexpr int nFaces = 4;
  TVertex* baryCenter = new TVertex(0,0,0);
  baryCentricCoords(oldTetra, *baryCenter);
  std::array<TBaseCell*, 4> newTetras;
  // create new tetrahedrons
  for(int i=0; i< nFaces; i++)
  {
    // allocating of tetrahedron
    auto faceVertices = getFaceVertices(oldTetra,i);
    newTetras[i] = allocTetrahedron(faceVertices[0], faceVertices[1],
                      faceVertices[2], baryCenter, oldTetra.GetRefDesc());

    auto joint = oldTetra.GetJoint(i);
    auto neighbour = joint->GetNeighbour(std::addressof(oldTetra));
    joint->SetNeighbour(0, neighbour);
    joint->SetNeighbour(1, newTetras[i]);
    // joint->SetMapType();
    // set the joint of new newTetras
    newTetras[i]->SetJoint(0, joint);
    newTetras[i]->SetJoint(1, nullptr);
    newTetras[i]->SetJoint(2, nullptr);
    newTetras[i]->SetJoint(3, nullptr);

    joint->SetMapType();
  }
  // set neighbours (joints)
  TVertex *a, *b;
  a = oldTetra.GetVertex(0);
  b = oldTetra.GetVertex(1);
  setNeighbours(std::begin(newTetras), std::end(newTetras), a, b, baryCenter);

  a = oldTetra.GetVertex(0);
  b = oldTetra.GetVertex(2);
  setNeighbours(std::begin(newTetras), std::end(newTetras), a, b, baryCenter);

  a = oldTetra.GetVertex(1);
  b = oldTetra.GetVertex(2);
  setNeighbours(std::begin(newTetras), std::end(newTetras), a, b, baryCenter);

  a = oldTetra.GetVertex(0);
  b = oldTetra.GetVertex(3);
  setNeighbours(std::begin(newTetras), std::end(newTetras), a, b, baryCenter);

  a = oldTetra.GetVertex(1);
  b = oldTetra.GetVertex(3);
  setNeighbours(std::begin(newTetras), std::end(newTetras), a, b, baryCenter);

  a = oldTetra.GetVertex(2);
  b = oldTetra.GetVertex(3);
  setNeighbours(std::begin(newTetras), std::end(newTetras), a, b, baryCenter);

  return newTetras;
}

/**
 *  predicate for mesh fix
 */ 
bool predAll (const TBaseCell& /*cell*/)
{
  return true;
}

/**
 *  predicate for mesh fix
 */ 
bool predTetgenFix (const TBaseCell& cell)
{
  int nBoundaryFaces=0;
  int nFaces = cell.GetN_Faces();
  for (int i=0;i<nFaces;++i)
  {
    auto jointType = cell.GetJoint(i)->GetType();

    if ( jointType == BoundaryFace || jointType == IsoBoundFace )
    {
      ++nBoundaryFaces;
    }
  }
  Output::print("Mesh Fixed: number of boundary faces: ", nBoundaryFaces);
  return nBoundaryFaces >= 3;
}

/**
 *  Does bary centric refinement of a cell if cell matches given predicate
 */
template <typename P>
void fixGrid (TDomain& domain, P predicate)
{
  TBaseCell** cellTree;
  int nCells, nFixed=0;
  
  domain.GetTreeInfo(cellTree, nCells);  

  std::vector<TBaseCell*> newCells;

  for(int i=0;i<nCells;++i)
  {
    if(predicate(*cellTree[i]) )
    {
      auto newTetras = baryCentricRefine(*cellTree[i]);
      
      newCells.push_back(newTetras[0]);
      newCells.push_back(newTetras[1]);
      newCells.push_back(newTetras[2]);
      newCells.push_back(newTetras[3]);
      
      ++nFixed;
    }
    else
    {
      newCells.push_back(cellTree[i]);
    }
  }

  delete [] cellTree;

  nCells = newCells.size();
  cellTree = new TBaseCell* [nCells];
  
  for(int i=0;i<nCells;++i)
  {
    cellTree[i] = newCells[i];
  }

  domain.SetTreeInfo(cellTree, nCells);
  
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_LE]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_Finest]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_Between]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_OCAF]->SetParam(std::addressof(domain));
  
  Output::print("Mesh Fix: ", nFixed, " Tetrahedron(s) refined");
}
// template <typename N>
void refineAll(TDomain& domain)
{
  TBaseCell** cellTree;
  int nCells, nFixed=0;

  domain.GetTreeInfo(cellTree, nCells);

  std::vector<TBaseCell*> newCells;

  for(int i=0;i<nCells;++i)
  {
    auto newTetras = baryCentricRefine(*cellTree[i]);

    newCells.push_back(newTetras[0]);
    newCells.push_back(newTetras[1]);
    newCells.push_back(newTetras[2]);
    newCells.push_back(newTetras[3]);

    ++nFixed;
  }
  delete [] cellTree;

  nCells = newCells.size();
  cellTree = new TBaseCell* [nCells];

  for(int i=0;i<nCells;++i)
  {
    cellTree[i] = newCells[i];
  }

  domain.SetTreeInfo(cellTree, nCells);
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_LE]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_Finest]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_Between]->SetParam(std::addressof(domain));
  TDatabase::IteratorDB[It_OCAF]->SetParam(std::addressof(domain));

  Output::print("Mesh Fix: ", nFixed, " Tetrahedron(s) refined");
}

void writeMesh(const char* fileName, TDomain& domain)
{
  TCollection* coll = domain.GetCollection(It_Finest,0);

  TOutput3D out(0,0,0,0, std::addressof(domain), coll);

  out.WriteVtk(fileName);
}

void usage (int /*argc*/, char** argv)
{
  Output::print("usage: \n");
  Output::print("\t", argv[0], " <in> <out>\n");
}

int main(int argc, char** argv)
{
  TDatabase Database;
  TFEDatabase3D feDatabase;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(ParameterDatabase::default_tetgen_database(), true);
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();

  TDomain domain(argv[1], parmoon_db);
  int nCells;
  // check for refined mesh
  // the initial mesh has one bad Tetrahedron
  /*writeMesh("initMesh.vtk", domain);
  nCells =domain.GetCollection(It_Finest,0)->GetN_Cells();
  Output::print("nCells Init: ", nCells,"\n");*/

  if(argc < 2 )
  {
    usage(argc, argv);
    return 1;
  }
  // refine that bad Tetrahedron using barycentric refinement
  // which replaces that Tetrahedron with four Tetrahedron
  /*fixGrid(domain, predTetgenFix);
  nCells =domain.GetCollection(It_Finest,0)->GetN_Cells();
  Output::print("nCells Fixed: ", nCells,"\n");
  writeMesh("vtk.vtk", domain);*/
  Output::print("refinement");
  // do refinement
  domain.RegRefineAll();
  nCells = domain.GetCollection(It_Finest,0)->GetN_Cells();
  Output::print("nCells Refined: ", nCells,"\n");
  // write the mesh to a file
  writeMesh("refinedmesh.vtk", domain);

  for(int i=0;i<2;i++)
  {
    nCells =domain.GetCollection(It_Finest,0)->GetN_Cells();
    Output::print("ref lev before: ",i , " nCells Refined fixed: ", nCells,"\n");
    refineAll(domain);
    nCells =domain.GetCollection(It_Finest,0)->GetN_Cells();
    Output::print("ref lev: ",i , " nCells Refined fixed: ", nCells,"\n");
  }
  writeMesh("ref.vtk", domain);
  exit(0);
  return 0;
}
