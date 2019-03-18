// =======================================================================
// @(#)FESpace.C        1.5 09/15/99
// 
// Class:       TFESpace
// Purpose:     general super class for all finite element spaces
//              special spaces are implemented in subclasses
//
// Author:      Gunar Matthies (04.11.97)
//
// History:     start of implementation 04.11.97 (Gunar Matthies)
//
//              split TFESpace into TFESpacexD (15.04.1998) Volker Behns
//
// =======================================================================

#include <Constants.h>
#include <FESpace.h>
#include <MooNMD_Io.h>

TFESpace::TFESpace(TCollection *coll,
                   const std::string& name,
                   const std::string& description)
:  Name(name), Description(description),
   Collection(coll)
{

  N_Cells=Collection->GetN_Cells();

  // The following values are dummies, they are re-set by
  // the ConstructSpace method of the daughter classes.
  N_DegreesOfFreedom = 0;
  GlobalNumbers = nullptr;
  BeginIndex = nullptr;
  N_UsedElements = 0;
  N_DiffBoundNodeTypes = 0;
  BoundaryNodeTypes = nullptr;
  N_Dirichlet = 0;
  N_BoundaryNodes = nullptr;
  N_Inner = 0;
  InnerBound = 0;
  BoundaryNodesBound = nullptr;
  DirichletBound = 0;
  ActiveBound = 0;
  is_discontinuous_galerkin_space = false;
}

/** destructor */
TFESpace::~TFESpace()
{
  if(GlobalNumbers)
    delete [] GlobalNumbers;
  if(BeginIndex)
    delete [] BeginIndex;
  if(BoundaryNodeTypes)
    delete [] BoundaryNodeTypes;
  if(N_BoundaryNodes)
    delete [] N_BoundaryNodes;
  if(BoundaryNodesBound)
    delete [] BoundaryNodesBound;
}

/** write info on fespace into file */
int TFESpace::Write(const std::string& filename)
{
  int header[4];
  int N_LocalDOF;

  std::ofstream dat(filename);
  if(!dat)
  {
    cerr << "cannot open file '" << filename << "' for output" << endl;
    return -1;
  }

  N_LocalDOF = BeginIndex[N_Cells];
  header[0] = N_Cells;
  header[1] = N_DegreesOfFreedom;
  header[2] = ActiveBound;
  header[3] = N_LocalDOF;

  dat.write((char *)header, sizeof(int)*4);
  
  dat.write((char *)BeginIndex, sizeof(int)*(N_Cells+1));
  dat.write((char *)GlobalNumbers, sizeof(int)*N_LocalDOF);
  
  dat.close();

  return 0;
}
