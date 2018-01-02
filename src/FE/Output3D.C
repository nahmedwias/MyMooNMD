#include <FEDatabase3D.h>
#include <Output3D.h>
#include <BaseCell.h>
#include <Joint.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <FEVectFunct3D.h>
#include <Domain.h>

#include <TetraAffin.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>

#include <sstream>
// #include <malloc.h>
#include <dirent.h> 
#include <unistd.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>

/** constructor maximum number of these things */
TOutput3D::TOutput3D(int maxn_fespaces, int maxn_scalar,
                     int maxn_vect, int maxn_parameters, TDomain *domain, TCollection *coll, const char *name)
{
  MaxN_FESpaces=maxn_fespaces;
  N_FESpaces=0;

  MaxN_ScalarVar=maxn_scalar;
  N_ScalarVar=0;

  MaxN_VectorVar=maxn_vect;
  N_VectorVar=0;

  MaxN_Parameters=maxn_parameters;
  N_Parameters=0;

  FESpaceArray=new const TFESpace3D*[MaxN_FESpaces];

  FEFunctionArray=new TFEFunction3D*[MaxN_ScalarVar];

  FEVectFunctArray=new TFEVectFunct3D*[MaxN_VectorVar];

  ParameterValues=new double[MaxN_Parameters];

  ParameterDescription=new const char*[MaxN_Parameters];

  Coll = coll;

  Domain = domain;

  Data = NULL;

  if (name) Name = strdup(name);
  else Name = strdup("none");
 
}

/** add a FESpace into this output object (internal use) */
int TOutput3D::AddFESpace(const TFESpace3D *fespace)
{
  const TFESpace3D **NewStorage;
  int i;

  if(N_FESpaces==0)
  { 
    // no fespace in this output object    
    Coll=fespace->GetCollection();
  }

  // check whether fespace is already stored
  for(i=0;i<N_FESpaces;i++)
    if(FESpaceArray[i]==fespace) return i+1;

  // this space is based on a different Collection
  if(fespace->GetCollection()!=Coll) return 0;

  if(MaxN_FESpaces<=N_FESpaces)
  {
    // enlarge storage
    NewStorage=new const TFESpace3D*[MaxN_FESpaces+5];
    memcpy(NewStorage, FESpaceArray, sizeof(const TFESpace3D *)*MaxN_FESpaces);
    MaxN_FESpaces +=5;
    delete FESpaceArray;
    FESpaceArray=NewStorage;
  }

  // store space on the next place
  FESpaceArray[N_FESpaces]=fespace;
  N_FESpaces++;

  return N_FESpaces;
}

/** add a FEFunction into this output object */
int TOutput3D::AddFEFunction(TFEFunction3D *fefunction)
{
  TFEFunction3D **NewStorage;

  if(!(AddFESpace(fefunction->GetFESpace3D()))) return 0;

  if(MaxN_ScalarVar<=N_ScalarVar)
  {
    // enlarge storage
    NewStorage=new TFEFunction3D*[MaxN_ScalarVar+5];
    memcpy(NewStorage, FEFunctionArray, MaxN_ScalarVar*sizeof(TFEFunction3D *));
    MaxN_ScalarVar +=5;
    delete FEFunctionArray;
    FEFunctionArray=NewStorage;
  }

  // store function on the next place
  FEFunctionArray[N_ScalarVar]=fefunction;
  N_ScalarVar++;

  return N_ScalarVar;
}

/** add a FEVectFunct into this output object */
int TOutput3D::AddFEVectFunct(TFEVectFunct3D *fevectfunct)
{
  TFEVectFunct3D **NewStorage;

  if(!(AddFESpace(fevectfunct->GetFESpace3D()))) return 0;

  if(MaxN_VectorVar<=N_VectorVar)
  {
    // enlarge storage
    NewStorage=new TFEVectFunct3D*[MaxN_VectorVar+5];
    memcpy(NewStorage, FEVectFunctArray, 
                MaxN_VectorVar*sizeof(TFEVectFunct3D *));
    MaxN_VectorVar +=5;
    delete FEVectFunctArray;
    FEVectFunctArray=NewStorage;
  }

  // store function on the next place
  FEVectFunctArray[N_VectorVar]=fevectfunct;
  N_VectorVar++;

  return N_VectorVar;
}

/** add parameter into this output object */
int TOutput3D::AddParameter(double value, const char *descr)
{
  double *NewValues;
  const char **NewDescr;

  if(MaxN_Parameters<=N_Parameters)
  {
    // enlarge storage
    NewValues=new double[MaxN_Parameters+5];
    NewDescr=new const char *[MaxN_Parameters+5];
    memcpy(NewValues, ParameterValues, sizeof(double)*MaxN_Parameters);
    memcpy(NewDescr, ParameterDescription, sizeof(const char *)*MaxN_Parameters);
    MaxN_Parameters +=5;
    delete ParameterValues;
    delete ParameterDescription;
    ParameterValues=NewValues;
    ParameterDescription=NewDescr;
  }

  // value and description on next place
  ParameterValues[N_Parameters]=value;
  ParameterDescription[N_Parameters]=strdup(descr);
  N_Parameters++;

  return N_Parameters;
}

static void Sort(TVertex **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, *rr, len;
  TVertex *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr=new int[2*len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete[] rr;

}

static int GetIndex(TVertex **Array, int Length, TVertex *Element)
{
  int l=0, r=Length, m=(r+l)/2;
  TVertex *Mid;

  Mid=Array[m];
  while(Mid!=Element)
  {
    if(Mid>Element)
    {
      l=m;
    }
    else
    {
      r=m;
    }
    m=(r+l)/2;
    Mid=Array[m];
  }
  return m;
}

/** write stored data into a VTK file */
/** start of implementation by Piotr Skrzypacz 24.03.04 */

int TOutput3D::WriteVtk(const char *name)
{
  int i,j,k,l,m,n;
  int N_Vertices, N_CellVertices, N_Comps;
  // *N_DOF, *N_Comp, *FESpaceNumber;
  int  N_, N_Elements, N_LocVertices;
  TVertex **Vertices, *Last, *Current;
  TBaseCell *cell;
  int *FESpaceNumber;
  const TFESpace3D *fespace;
  double xi=0, eta=0, zeta=0;
  double value;
  double *Coeffs;
  int N_LocDOF;
  int Length, N_Comp;
  double t;
  
  TBaseFunct3D *bf;
  FE3D FE_ID;
  double BFValues[3*MaxN_BaseFunctions3D]; // 3 for vector valued basis functions
  int *GlobalNumbers, *BeginIndex, *DOF;
  double *Coords;
  int *VertexNumbers;
  int *NumberVertex;
  int *IntArray;
  double *DoubleArray;

 // format (x_i, y_i, z_i)
  static double HexaCoords[] =
    { -1, -1, -1,
       1, -1, -1,
       1,  1, -1,
      -1,  1, -1,
      -1, -1,  1,  
       1, -1,  1,
       1,  1,  1,
      -1,  1,  1 
    };

  // format (x_i, y_i, z_i)
  static double TetraCoords[] =
    { 0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1 };


  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(4);

  FESpaceNumber = new int[N_ScalarVar+N_VectorVar];

  N_Comps = 0;
  for(i=0;i<N_ScalarVar;i++)
  {
    N_Comps++;
    fespace=FEFunctionArray[i]->GetFESpace3D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[i]=j;
  }

  k = N_ScalarVar;
  for(i=0;i<N_VectorVar;i++,k++)
  {
    N_Comps += FEVectFunctArray[i]->GetN_Components();
    fespace=FEVectFunctArray[i]->GetFESpace3D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[k]=j;
  }

  // determine data for vtk file

  N_Elements=Coll->GetN_Cells();
//   cout << "N_Elements: " <<  N_Elements << endl;
  N_LocVertices=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
  }
  Vertices=new TVertex*[N_LocVertices];
  N_=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Vertices[N_]=cell->GetVertex(j);
      N_++;
    }
  }
  
  // check for discontinuous scalar variables. In such a case write a new file
  // for this variable. ParaView will really display a discontinuous function
  // instead of projecting it onto P1/Q1 space.
  // However there are some drawbacks. 
  for(i = 0; i < N_ScalarVar; i++)
  {
    j = FEFunctionArray[i]->GetFESpace3D()->IsDGSpace();
    if(j==1)
    {
      // draw all discontinuous functions not just this i-th one.
      WriteVtkDiscontinuous(name,N_LocVertices,Vertices);
      break;
    }
  }
  
//   cout << "N_" << N_ << endl;
  Sort(Vertices, N_);
  //Sort(cell_types, N_);
  Last=NULL;
  N_Vertices=0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }
  //cout << "N_Vertices: " << N_Vertices << endl;
  Coords=new double[3*N_Vertices];
  VertexNumbers=new int[N_LocVertices];
  NumberVertex=new int[N_LocVertices];
  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[3*N_],Coords[3*N_+1], Coords[3*N_+2]);
      k++;
      N_++;
      Last=Current;
    }
    NumberVertex[i]=k;
  }
  
  m=0;
  for(i=0;i<N_Elements;i++)
    {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      //cout << "Vertex Number: " << VertexNumbers[m] << endl;
      m++;
    } // endfor j
  } //endfor i

 


  // one additional column for absolute values of velocity
  N_Comps++;

  // to check
  //cout << "MaxN_VerticesPerCell*N_Comps" << MaxN_VerticesPerCell*N_Comps << endl;
  //cout << "MaxN_VerticesPerCell" << MaxN_VerticesPerCell << endl;
  //cout << "N_Comps" << N_Comps << endl;
  
  dat << std::scientific;
  dat.precision(6);
  dat << "# vtk DataFile Version 4.0" << endl;
  dat << "file created by MooNMD" << endl;


  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl << endl;
  dat << "POINTS " << N_Vertices << " double" << endl;
  N_=0;
  //cout << "N_LocVertices: " << N_LocVertices << endl;
  for(i=0;i<N_Vertices;i++)
  {  
    dat << Coords[N_] << " " <<  Coords[N_+1] << " " << Coords[N_+2] << endl;
    N_ +=3;
  }
  dat << endl;
  dat << "CELLS " << N_Elements << " " <<  N_Elements+N_LocVertices << endl;
  l=0;
  for(i=0;i<N_Elements;i++)
  {
    N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
    dat <<  N_CellVertices << " ";
    for(j=0;j<N_CellVertices;j++)
    {
      dat << VertexNumbers[l] << " ";
      l++;
    }
    dat << endl;
  }
  dat << endl;
  dat << "CELL_TYPES " << N_Elements << endl;
  for(i=0;i<N_Elements;i++)
  {  
    N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
    switch(N_CellVertices)
    {
    case 4: dat << 10 << " ";
      break;
    case 8: dat << 12 << " ";
      break; 
    }
  }
  dat << endl << endl;
  dat << "POINT_DATA " << N_Vertices << endl;  

  // function values
  
  DoubleArray = new double[3*N_Vertices];
  IntArray = new int[N_Vertices];

   // write scalar variables into file
  for(k=0;k<N_ScalarVar;k++)
  {
    fespace = FEFunctionArray[k]->GetFESpace3D();
    Coeffs = FEFunctionArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    memset(DoubleArray, 0, SizeOfDouble*N_Vertices);
    memset(IntArray, 0, SizeOfInt*N_Vertices);
    m = 0;
    
    // will be set to true for vector valued basis functions (for example 
    // Raviart-Thomas or Brezzi-Douglas-Marini)
    bool VectOutput = false;
      
    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE3D(i, cell);
      bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      int BaseVectDim = bf->GetBaseVectDim();
      if(BaseVectDim == 3) 
        VectOutput = true;
      else if(BaseVectDim != 1)
        ErrMsg("unkown number of basis function components, assume 1.");
      for(j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
	  //  Tetrahedron
          case 4: 
	    xi   = TetraCoords[3*j];
            eta  = TetraCoords[3*j+1];
            zeta = TetraCoords[3*j+2];
	    break;
	  // Hexahedron
	  case 8: 
            xi   = HexaCoords[3*j];
            eta  = HexaCoords[3*j+1];
            zeta = HexaCoords[3*j+2];
	    break;
        }
        bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
        if(!VectOutput)
        {
          value = 0;
          for(l=0;l<N_LocDOF;l++)
            value += BFValues[l] * Coeffs[DOF[l]];
          DoubleArray[VertexNumbers[m]] += value;
        }
        else // VectOutput
        {
          // transform values using the Piola transform
          RefTrans3D RefTrans = TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FE_ID);
          TRefTrans3D *F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
          TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);
          double *BFValuesOrig = new double[3*N_LocDOF];
          switch(RefTrans)
          {
            case TetraAffin:
            case HexaAffin:
              F_K->PiolaMapOrigFromRef(N_LocDOF, BFValues, BFValuesOrig);
              break;
            case HexaTrilinear:
              ErrMsg("Piola transform for trilinear reference map not yet " << 
                     "implemented");
              break;
            default:
              ErrMsg("unknown reference transformation");
              exit(0);
              break;
          }
          
          double value_x = 0, value_y = 0, value_z = 0;
          for( l = 0; l < N_LocDOF; l++)
          {
            int face = TFEDatabase3D::GetFE3D(FE_ID)->GetFEDesc3D()
                ->GetJointOfThisDOF(l);
            int nsign = 1;
            if(face != -1)
              nsign = cell->GetNormalOrientation(face);
            value_x += BFValuesOrig[l             ] * Coeffs[DOF[l]]*nsign;
            value_y += BFValuesOrig[l +   N_LocDOF] * Coeffs[DOF[l]]*nsign;
            value_z += BFValuesOrig[l + 2*N_LocDOF] * Coeffs[DOF[l]]*nsign;
          }
          DoubleArray[BaseVectDim*VertexNumbers[m] + 0] += value_x;
          DoubleArray[BaseVectDim*VertexNumbers[m] + 1] += value_y;
          DoubleArray[BaseVectDim*VertexNumbers[m] + 2] += value_z;
          
          delete [] BFValuesOrig;
        }
        IntArray[VertexNumbers[m]]++;
        m++;
      } // endfor j
    } // endfor i

    if(!VectOutput)
    {
      // non conforming
      for(i=0;i<N_Vertices;i++)
        DoubleArray[i] /= IntArray[i];
    }
    else // VectOutput
    {
      for(i = 0; i < N_Vertices; i++)
      {
        DoubleArray[3*i    ] /= IntArray[i];
        DoubleArray[3*i + 1] /= IntArray[i];
        DoubleArray[3*i + 2] /= IntArray[i];
      }
    }

    if(!VectOutput)
    {
      dat << "SCALARS " << FEFunctionArray[k]->GetName();
      dat << " double"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(j=0;j<N_Vertices;j++)
        dat << DoubleArray[j] << endl;
      dat << endl;
      dat << endl;
    }
    else
    {
      // vector output, we don't write each component individually, 
      dat << "VECTORS " << FEFunctionArray[k]->GetName() << " double\n";
      for(i = 0; i < N_Vertices; i++)
      {
        for(j = 0; j < 3; j++)
        {
          dat << DoubleArray[3 * i + j] << " ";
        }
        dat << endl;
      }
      dat << endl;
      // reset
      VectOutput = false;
    }
  } // endfor k

 
  for(k=0;k<N_VectorVar;k++)
  {
    fespace = FEVectFunctArray[k]->GetFESpace3D();
    N_Comp = FEVectFunctArray[k]->GetN_Components();
    Length = FEVectFunctArray[k]->GetLength();
    Coeffs = FEVectFunctArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    
    memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
    memset(IntArray, 0, SizeOfInt*N_Vertices);
    m = 0;

    
    //for(k=0;k<FEVectFunctArray[i]->GetN_Components();k++)
    //{
      
    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE3D(i, cell);
      bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
          case 4: 
	    xi   = TetraCoords[3*j];
            eta  = TetraCoords[3*j+1];
            zeta = TetraCoords[3*j+2];
	    break;

	  case 8: 
	    xi   = HexaCoords[3*j];
            eta  = HexaCoords[3*j+1];
            zeta = HexaCoords[3*j+2];
	    break;
        }
        bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
       
	for(n=0;n<N_Comp;n++)
        {
	  value = 0;
	  for(l=0;l<N_LocDOF;l++)
	    value += BFValues[l] * Coeffs[DOF[l]+n*Length];
	  DoubleArray[N_Comp*VertexNumbers[m] + n] += value;
	} 
        IntArray[VertexNumbers[m]]++;
        m++;
      } // endfor j
    } // endfor i

    // midle value

    l = 0;
    for(i=0;i<N_Vertices;i++)
    {
      for(j=0;j<N_Comp;j++)
      {
        DoubleArray[l] /= IntArray[i];
	l++;
      }
    } // endfor l
    /*
    for(i=0;i<2*N_Vertices;i++)
    {
      cout << "Do[" << i << "]" << DoubleArray[i] << endl;
    }
    */
    
    for(j=0;j<N_Comp;j++)
    {
      dat << "SCALARS " << FEVectFunctArray[k]->GetName() << j;
      dat << " double"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(i=0;i<N_Vertices;i++)
      {
	dat << DoubleArray[i*N_Comp+j] << endl;
      }
      dat << endl << endl;
    }
    
    dat << "SCALARS " << "|" << FEVectFunctArray[k]->GetName() << "|";
    dat << " double"<< endl;
    dat << "LOOKUP_TABLE " << "default" << endl;
    l=0;
    for(i=0;i<N_Vertices;i++)
    {
      t=0;
      for(j=0;j<N_Comp;j++)
      {	
       t+=DoubleArray[l]*DoubleArray[l];
       l++;
      }
      dat << sqrt(t)<< endl;
    }
    dat << endl << endl;

    dat << "VECTORS " << FEVectFunctArray[k]->GetName();
    dat << " double"<< endl;
    
    l=0;
    for(i=0;i<N_Vertices;i++)
    {
      for(j=0;j<N_Comp;j++)
      {
	dat << DoubleArray[N_Comp*i+j] << " ";
      }
      dat << endl;
    }
    dat << endl;
  } // endfor k

  dat.close();
  
  delete [] IntArray;
  delete [] DoubleArray;
  delete [] NumberVertex;
  delete [] VertexNumbers;
  delete [] Vertices;
  delete [] Coords;
  delete [] FESpaceNumber;

  Output::print("wrote output into vtk file: ", name);
  return 0;
}


/*
  Ulrich Wilbrandt, December 2013.

  writes an extra vtk-file for (scalar) discontinuous functions. ParaView can
  then display discontinuous data. The "POINTS" in the resulting vtk-file are 
  the vertices of the mesh, but every vertex is put here as many times as there 
  are cells this vertex belongs to. That means if a vertex belongs to 4 cells 
  in the mesh it will appear 4 times in the list "POINTS" in the resulting 
  vtk-file. The input "TVertex **Vertices" should be in this pattern. It can be
  generated by

  N_Elements=Coll->GetN_Cells();
  N_LocVertices=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
  }
  Vertices=new TVertex*[N_LocVertices];
  N_=0;
  
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Vertices[N_]=cell->GetVertex(j);
      N_++;
    }
  }

  as is done in WriteVtk(cont char *name).
  WARNING: This destroys the topology of the mesh. Some filters in ParaView 
  might work incorrectly. Warp by scalar works though.
  The way this is done here is not very elegant, but I couldn't find a better
  solution.
  Also note that in each element the function is projected onto P1/Q1.
*/
void TOutput3D::WriteVtkDiscontinuous(const char *fileName,
int N_LocVertices, TVertex **Vertices)
{
  char Disc[80];             // the new file name
  strcpy(Disc,fileName);     // copy the file name ...
  strcat(Disc,"_disc.vtk");  // ... add a string to the output file name
  
  std::ofstream dat(Disc);
  if (!dat)
  {
    Error("cannot open file for output\n");
    exit(-1);
  }
  
  dat.setf(std::ios::fixed);
  dat << setprecision(9);

  dat << "# vtk DataFile Version 4.0" << endl;
  dat << "file created by MooNMD." << endl;

  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl;
  dat << "POINTS " << N_LocVertices << " float" << endl;

  for(int i = 0; i < N_LocVertices; i++)
  {
    double x,y,z;
    Vertices[i]->GetCoords(x,y,z);
    dat << x << " " <<  y << " " << z << endl;
  }
  dat << endl;
  int N_Elements = Coll->GetN_Cells();
  
  // writing which vertices belong to which cells, here it is ignored that
  // a vertex might belong to multiple cells.
  dat << "CELLS " << N_Elements << " " <<  N_Elements+N_LocVertices << endl;
  int l = 0;
  for(int i = 0; i < N_Elements; i++)
  {
    TBaseCell* current_cell = Coll->GetCell(i);
    int N_CellVertices = current_cell->GetN_Vertices();
    dat <<  N_CellVertices << " ";
    for(int j = 0; j < N_CellVertices; j++)
    {
      dat << l << " ";
      l++;
    }
    dat << endl;
  }
  dat << endl;
  
  // the cell types tell paraview if this is a tetrahedron or a hexahedron
  // (export of other types is not supported here)
  dat << "CELL_TYPES " << N_Elements << endl;
  for(int i = 0; i < N_Elements; i++)
  {
    int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
    switch(N_CellVertices)
    {
    case 4: dat << 10 << " ";
      break;
    case 8: dat << 12 << " ";
      break; 
    }
  }
  dat << endl << endl;
  
  // write the function values, only for scalar functions, which includes
  // vector valued basis functions (such as Raviart-Thomas), because these 
  // are handled like scalar basis functions
  dat << "POINT_DATA " << N_LocVertices << endl;
  for(int space_number = 0; space_number < N_ScalarVar; space_number++)
  {
    TFEFunction3D* fefunction = FEFunctionArray[space_number];
    if(!fefunction->GetFESpace3D()->IsDGSpace())
      continue;
    // this is a discontinuous space
    int BaseVectDim = TFEDatabase3D::GetFE3D(
      fefunction->GetFESpace3D()->GetFE3D(0,Coll->GetCell(0)))
      ->GetBaseFunct3D()->GetBaseVectDim();   // ugly, but we need to know this
    // scalar valued basis functions (normal case)
    if (BaseVectDim == 1)
    {
      dat << endl << endl;
      dat << "SCALARS " << fefunction->GetName();
      dat << " double" << endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      double *function_value = new double[1];
      for(int i = 0; i < N_Elements; i++)
      {
        TBaseCell* current_cell = Coll->GetCell(i);
        int N_CellVertices = current_cell->GetN_Vertices();
        for(int j = 0; j < N_CellVertices; j++)
        {
          double x,y,z;
          current_cell->GetVertex(j)->GetCoords(x,y,z);
          fefunction->FindValueLocal(current_cell,i,x,y,z,function_value);
          dat << function_value[0] << endl;
        }
      }
      delete [] function_value;
    }
    // vector valued basis functions (e.g. Raviart-Thomas)
    else if(BaseVectDim == 3)
    {
      // find values for all components
      double* function_value = new double[3];   // 3==BaseVectDim
      // store function values at all vertices (three components)
      double **allValues = new double*[3];       // 3==BaseVectDim
      for(int l = 0; l < BaseVectDim; l++)
        allValues[l] = new double[N_LocVertices];
      for(int i = 0, k = 0; i < N_Elements; i++)
      {
        TBaseCell* current_cell = Coll->GetCell(i);
        int N_CellVertices = current_cell->GetN_Vertices();
        for(int j = 0; j < N_CellVertices; j++)
        {
          double x,y,z;
          current_cell->GetVertex(j)->GetCoords(x, y, z);
          // FindValueLocal includes the necessary sign changes due to global
          // normals (for Raviart-Thomas elements)
          fefunction->FindValueLocal(current_cell,i,x,y,z,function_value);
          //for(int l = 0; l < BaseVectDim; l++)
          //  allValues[l][k] = function_value[l];
          allValues[0][k] = function_value[0];
          allValues[1][k] = function_value[1];
          allValues[2][k] = function_value[2];
          k++;
        }
      }
      delete [] function_value;
      // write the function values to the vtk-file
      dat << endl << endl;
      dat << "VECTORS " << fefunction->GetName();
      dat << " double"<< endl;
      for(int i = 0, k = 0; i < N_Elements; i++)
      {
        int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
        for(int j = 0; j < N_CellVertices; j++)
        {
          dat << allValues[0][k] << "\t" << allValues[1][k] 
              << "\t" << allValues[2][k] << endl;
          k++;
        }
      }
      
      for(int l = 0; l < BaseVectDim; l++)
        delete [] allValues[l];
      delete [] allValues;
    }
    else
      ErrMsg("TOutput3D::WriteVtkDiscontinuous: Basis functions of dimension "
             << BaseVectDim << " are not supported");
    
    
  }
  dat.close();
}



/** write stored PARALLEL data into a pvtu and vtu files (XML files for paraview) (Sashikumaar Ganesan) */

int TOutput3D::Write_ParVTK(
#ifdef _MPI
                                MPI_Comm comm,
#endif
                               int img, char *subID,
							   std::string directory,
							   std::string basename)
{
  int i, j, k,l,m,n, rank, size, N_, N_Elements, N_LocVertices;
  int N_Vertices, N_CellVertices, N_Comps;
  int *FESpaceNumber, N_LocDOF, Length, N_Comp, *GlobalNumbers, *BeginIndex, *DOF;
  int *VertexNumbers=nullptr, *NumberVertex=nullptr, begin, ID;

  double xi=0, eta=0, zeta=0, value, *Coeffs, *WArray=nullptr, *DoubleArray=nullptr;
  double BFValues[MaxN_BaseFunctions3D];
  double *Coords=nullptr;
  static double HexaCoords[] = { -1, -1, -1, 1, -1, -1, 1,  1, -1, -1,  1, -1,
                                 -1, -1,  1, 1, -1,  1, 1,  1,  1, -1,  1,  1  };
  static double TetraCoords[] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  char Dquot;

#ifdef _MPI
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
#else
  rank = 0;
  size =1;
#endif

  std::string outputdir(directory);
  std::string vtu("VTU/");
  std::string vtudir = outputdir + std::string("/") + vtu;

  //get the proper path for the VTU directory (which contains the parts)
//  if (rank == 0)
//  {
//    Output::print(vtudir);
//  }

  time_t rawtime;
  struct tm * timeinfo;

  TVertex **Vertices=nullptr, *Last, *Current;
  TBaseCell *cell;
  const TFESpace3D *fespace;
  TBaseFunct3D *bf;
  FE3D FE_ID;

  Dquot = 34; //  see ASCII Chart
  const char* VtkBaseName = basename.c_str();
  const char* output_directory = directory.c_str();

  if(rank==0)
   {
//     remove(vtudir);
    mkdir(vtudir.c_str(), 0777); // create the folder to store SubDomain vtu files
   }

#ifdef _MPI
  MPI_Barrier (TDatabase::ParamDB->Comm);
#endif

  std::ostringstream os;
  os << " ";

  time ( &rawtime );
  timeinfo = localtime (&rawtime );


// write the master pvtu file
  if(rank==0)
   {
    Output::info("Output3D","writing output into ", output_directory, "/", VtkBaseName
           , subID, "*.", img, " xml vtk file");
    os.seekp(std::ios::beg);
    os << output_directory << "/" << VtkBaseName << subID;
    if(img<10) os << ".0000"<<img<<".pvtu" << ends;
    else if(img<100) os << ".000"<<img<<".pvtu" << ends;
    else if(img<1000) os << ".00"<<img<<".pvtu" << ends;
    else if(img<10000) os << ".0"<<img<<".pvtu" << ends;
    else  os << "."<<img<<".pvtu" << ends;
    std::ofstream dat(os.str().c_str());
    if (!dat)
     {
      cerr << "cannot open file for output" << endl;
#ifdef _MPI
      MPI_Abort(MPI_COMM_WORLD, 0);
#else
     exit(0);
#endif
     }

    dat <<  "<?xml version="<<Dquot<<"1.0"<<Dquot<<"?>"<< endl;
    dat << endl;
    dat <<  "<!--" << endl;
    dat <<  "      Title: Master file for parallel vtk data" << endl;
    dat <<  "    Program: ParMooN " << endl;
    dat <<  "    Version: v1.0.0 " << endl;
    dat <<  "Date & Time: " <<asctime (timeinfo) << endl;
    dat <<  "Problem Current Time " <<TDatabase::TimeDB->CURRENTTIME << endl;
    dat <<  "  -->" << endl;

    dat << endl;


    dat <<  "<VTKFile type="<<Dquot<<"PUnstructuredGrid"<<Dquot<<" version="<<Dquot<<"0.1"
             <<Dquot<<" byte_order="<<Dquot<<"LittleEndian"<<Dquot<<">"<<endl;
    dat <<  "<PUnstructuredGrid GhostLevel="<<Dquot<<0<<Dquot<<">"<<endl;
    dat << endl;

    dat <<  " <PPoints>"<<endl;
    dat <<  "   <PDataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
        <<"Position"<<Dquot<<" NumberOfComponents="<<Dquot<<"3"<<Dquot<<"/>"<<endl;
    dat <<   "</PPoints>"<<endl;
    dat << endl;

    dat <<   "<PCells>"<<endl;
    dat <<   "  <PDataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"connectivity"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
    dat <<   "  <PDataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"offsets"<<Dquot
        <<"      NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
    dat <<   "  <PDataArray type="<<Dquot<<"UInt8"<<Dquot<<" Name="<<Dquot<<"types"<<Dquot
        <<"        NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
    dat <<   "</PCells>"<<endl;
    dat << endl;

    dat <<   "<PPointData Vectors="<<Dquot<<"Vectors"<<Dquot<<" "<<"Scalars="<<Dquot<<"Scalars"<<Dquot<<">"<<endl;
    for(i=0;i<N_VectorVar;i++)
    dat <<   "  <PDataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot<<FEVectFunctArray[i]->GetName()<<Dquot
        <<" NumberOfComponents="<<Dquot<<"3"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<"/>"<<endl;
      dat << endl;

    for(i=0;i<N_VectorVar;i++)
    {
    for(j=0;j<FEVectFunctArray[i]->GetN_Components();j++)
     {
      dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
          <<FEVectFunctArray[i]->GetName()<<j<<Dquot<<" NumberOfComponents="<<Dquot
          <<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<"/>"<<endl;
      dat << endl;
     }
    }

    for(i=0;i<N_ScalarVar;i++)
    dat <<   "  <PDataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot<<FEFunctionArray[i]->GetName()<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<"/>"<<endl;
    dat <<   "</PPointData>"<<endl;
    dat << endl;

    dat <<   "<PCellData Scalars="<<Dquot<<"SubDomainAndRegionID"<<Dquot<<">"<<endl;
#ifdef _MPI    
    dat <<   "  <PDataArray type="<<Dquot<<"Int32"<<Dquot<<"   Name="<<Dquot<<"SubDomain"<<Dquot
        <<"  NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
#endif
    dat <<   "  <PDataArray type="<<Dquot<<"Int32"<<Dquot<<"   Name="<<Dquot<<"RegionID"<<Dquot
        <<"  NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;

    dat <<   "</PCellData>"<<endl;
    dat << endl; 

 
    // root not take part in computation
//     begin = 1;
    begin = 0; // root take part in computation
    
    for(i=begin;i<size;i++)
     {
      ID = i; // root not take part in computation

      if(img<10)
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".000"<<ID<<".0000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".00" <<ID<<".0000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<".0000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<"."   <<ID<<".0000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
      else if(img<100)
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".000"<<ID<<".000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".00" <<ID<<".000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<".000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<"."   <<ID<<".000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
      else if(img<1000)
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".000"<<ID<<".00"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".00" <<ID<<".00"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<".00"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<"."   <<ID<<".00"<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
      else if(img<10000)
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".000"<<ID<<".0"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".00" <<ID<<".0"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<".0"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<"."   <<ID<<".0"<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
      else
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".000"<<ID<<"."<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".00" <<ID<<"."<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<"."<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<vtu<<VtkBaseName<<subID<<"."   <<ID<<"."<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
     }
    dat << endl;

    dat <<  "</PUnstructuredGrid>"<<endl;
    dat <<  "</VTKFile>"<<endl;
    dat.close();
   } // if(rank==0
   
  // root take part in computation 
//   else
   {

    // determine data for vtu file of each processor
    FESpaceNumber = new int[N_ScalarVar+N_VectorVar];
    //cout << "N_ScalarVar: " <<  N_ScalarVar << endl;
    N_Comps = 0;
    for(i=0;i<N_ScalarVar;i++)
    {
     N_Comps++;
     fespace=FEFunctionArray[i]->GetFESpace3D();
     j=0;
     while(FESpaceArray[j]!=fespace) j++;
     FESpaceNumber[i]=j;
    }

    k = N_ScalarVar;
    for(i=0;i<N_VectorVar;i++,k++)
    {
     N_Comps += FEVectFunctArray[i]->GetN_Components();
     fespace=FEVectFunctArray[i]->GetFESpace3D();
     j=0;
     while(FESpaceArray[j]!=fespace) j++;
     FESpaceNumber[k]=j;
    }


#ifdef _MPI
  N_Elements=Coll->GetN_OwnCells();
#else
  N_Elements=Coll->GetN_Cells();
#endif


  N_LocVertices=0;
  for(i=0;i<N_Elements;i++)
   {
    cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
   }

  if(N_LocVertices)
   Vertices=new TVertex*[N_LocVertices];
  N_=0;
  for(i=0;i<N_Elements;i++)
   {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
     {
      Vertices[N_]=cell->GetVertex(j);
      N_++;
     }
   }


  if(N_)
   Sort(Vertices, N_);

  Last=NULL;
  N_Vertices=0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }

  if(N_LocVertices)
   {
    Coords=new double[3*N_Vertices];
    VertexNumbers=new int[N_LocVertices];
    NumberVertex=new int[N_LocVertices];
   }
   Last=NULL;
   N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[3*N_],Coords[3*N_+1], Coords[3*N_+2]);
      k++;
      N_++;
      Last=Current;
    }
    NumberVertex[i]=k;
  }

  m=0;
  for(i=0;i<N_Elements;i++)
    {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      m++;
    } // endfor j
  } //endfor i


  // additional column for absolute values of velocity
  N_Comps++;

  // to check
  //cout << "MaxN_VerticesPerCell*N_Comps" << MaxN_VerticesPerCell*N_Comps << endl;
  //cout << "MaxN_VerticesPerCell" << MaxN_VerticesPerCell << endl;
  //cout << "N_Comps" << N_Comps << endl;

  if(N_LocVertices)
   {
    DoubleArray = new double[3*N_Vertices];
    WArray = new double[N_Vertices];
   }

  ID = rank;
    os.seekp(std::ios::beg);
      if(img<10)
       {
        if(ID<10)        os<<vtudir<<VtkBaseName<<subID<<".000"<<ID<<".0000"<<img<<".vtu" <<ends;
        else if(ID<100)  os<<vtudir<<VtkBaseName<<subID<<".00" <<ID<<".0000"<<img<<".vtu" <<ends;
        else if(ID<1000) os<<vtudir<<VtkBaseName<<subID<<".0"  <<ID<<".0000"<<img<<".vtu" <<ends;
        else            os<<vtudir<<VtkBaseName<<subID<<"."   <<ID<<".0000"<<img<<".vtu" <<ends;
       }
      else if(img<100)
       {
        if(ID<10)        os<<vtudir<<VtkBaseName<<subID<<".000"<<ID<<".000"<<img<<".vtu" <<ends;
        else if(ID<100)  os<<vtudir<<VtkBaseName<<subID<<".00" <<ID<<".000"<<img<<".vtu" <<ends;
        else if(ID<1000) os<<vtudir<<VtkBaseName<<subID<<".0"  <<ID<<".000"<<img<<".vtu" <<ends;
        else            os<<vtudir<<VtkBaseName<<subID<<"."   <<ID<<".000"<<img<<".vtu" <<ends;
       }
      else if(img<1000)
       {
        if(ID<10)        os<<vtudir<<VtkBaseName<<subID<<".000"<<ID<<".00"<<img<<".vtu" <<ends;
        else if(ID<100)  os<<vtudir<<VtkBaseName<<subID<<".00" <<ID<<".00"<<img<<".vtu" <<ends;
        else if(ID<1000) os<<vtudir<<VtkBaseName<<subID<<".0"  <<ID<<".00"<<img<<".vtu" <<ends;
        else            os<<vtudir<<VtkBaseName<<subID<<"."   <<ID<<".00"<<img<<".vtu" <<ends;
       }
      else if(img<10000)
       {
        if(ID<10)        os<<vtudir<<VtkBaseName<<subID<<".000"<<ID<<".0"<<img<<".vtu" <<ends;
        else if(ID<100)  os<<vtudir<<VtkBaseName<<subID<<".00" <<ID<<".0"<<img<<".vtu" <<ends;
        else if(ID<1000) os<<vtudir<<VtkBaseName<<subID<<".0"  <<ID<<".0"<<img<<".vtu" <<ends;
        else            os<<vtudir<<VtkBaseName<<subID<<"."   <<ID<<".0"<<img<<".vtu" <<ends;
       }
      else
       {
        if(ID<10)        os<<vtudir<<VtkBaseName<<subID<<".000"<<ID<<"."<<img<<".vtu" <<ends;
        else if(ID<100)  os<<vtudir<<VtkBaseName<<subID<<".00" <<ID<<"."<<img<<".vtu" <<ends;
        else if(ID<1000) os<<vtudir<<VtkBaseName<<subID<<".0"  <<ID<<"."<<img<<".vtu" <<ends;
        else            os<<vtudir<<VtkBaseName<<subID<<"."   <<ID<<"."<<img<<".vtu" <<ends;
       }

    std::ofstream dat(os.str().c_str());
//    Output::print(os.str().c_str());
    if (!dat)
     {
      cerr << "cannot open file for output" << endl;
#ifdef _MPI
      MPI_Abort(MPI_COMM_WORLD, 0);
#else
     exit(0);
#endif
     }
    dat << setprecision(8);

    dat <<  "<?xml version="<<Dquot<<"1.0"<<Dquot<<"?>" << endl;
    dat << endl;
    dat <<  "<!--" << endl;
    dat <<  "      Title: SubDomain data for master ptvu file" << endl;
    dat <<  "    Program: ParMooN " << endl;
    dat <<  "    Version: v1.0.0 " << endl;
    dat <<  "Date & Time: " <<asctime (timeinfo) << endl;
    dat <<  "  -->" << endl;
    dat << endl;

    dat <<  "<VTKFile type="<<Dquot<<"UnstructuredGrid"<<Dquot<<" version="<<Dquot<<"0.1"
             <<Dquot<<" byte_order="<<Dquot<<"LittleEndian"<<Dquot<<">"<<endl;
    dat <<  "<UnstructuredGrid>"<<endl;
    dat << endl;

    dat <<  "<Piece NumberOfPoints="<<Dquot<<N_Vertices<<Dquot<<" NumberOfCells="<<Dquot<<N_Elements<<Dquot<<">"<<endl;
    dat <<  "<Points>"<<endl;
    dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot<<"Position"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"3"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    N_=0;
    for(i=0;i<N_Vertices;i++)
    {
     dat <<  "   " << Coords[N_] << " " <<  Coords[N_+1] << " " << Coords[N_+2] << endl;
     N_ +=3;
    }
    dat <<  "  </DataArray>"<<endl;
    dat <<  "</Points>"<<endl;


    dat <<  "<Cells>"<<endl;
    dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"connectivity"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    l=0;
    for(i=0;i<N_Elements;i++)
     {
      N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
      dat << "       ";
      for(j=0;j<N_CellVertices;j++)
       {
        dat << VertexNumbers[l] << " ";
        l++;
       }
      dat << endl;
      }
    dat <<  "  </DataArray>"<<endl;
    dat << endl;

    dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"offsets"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    for(i=1;i<=N_Elements;i++)
      {
        N_CellVertices=Coll->GetCell(i-1)->GetN_Vertices();
        dat <<  i*N_CellVertices<<"  " ;
      }
    dat <<  "   </DataArray>"<<endl;
    dat << endl;

    dat <<  "  <DataArray type="<<Dquot<<"UInt8"<<Dquot<<"  Name="<<Dquot<<"types"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    for(i=0;i<N_Elements;i++)
     {
      N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
      switch(N_CellVertices)
       {
        case 4: dat << 10 << " ";
           break;
        case 8: dat << 12 << " ";
           break; 
    }
     }
    dat <<  "  </DataArray>"<<endl;
    dat <<  "</Cells>"<<endl;
    dat <<endl;
    dat <<  "<PointData Vectors="<<Dquot<<"Velocity"<<Dquot<<" "<<"Scalars="<<Dquot<<"Scalars"<<Dquot<<">"<<endl;
    dat << endl;


  // write vector variables into file
    if(N_LocVertices)
    for(k=0;k<N_VectorVar;k++)
     {
      fespace = FEVectFunctArray[k]->GetFESpace3D();
      N_Comp = FEVectFunctArray[k]->GetN_Components();
      Length = FEVectFunctArray[k]->GetLength();
      Coeffs = FEVectFunctArray[k]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();

      memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
      memset(WArray, 0, SizeOfDouble*N_Vertices);
      m = 0;

      for(i=0;i<N_Elements;i++)
       {
        cell = Coll->GetCell(i);
        N_ = cell->GetN_Vertices();

        // find FE data for this element
        FE_ID = fespace->GetFE3D(i, cell);
        bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
        DOF = GlobalNumbers+BeginIndex[i];
        N_LocDOF = bf->GetDimension();
        for(j=0;j<N_;j++)
         {
          switch(cell->GetN_Vertices())
           {
            case 4: 
	      xi   = TetraCoords[3*j];
              eta  = TetraCoords[3*j+1];
              zeta = TetraCoords[3*j+2];
            break;

	    case 8: 
	      xi   = HexaCoords[3*j];
              eta  = HexaCoords[3*j+1];
              zeta = HexaCoords[3*j+2];
	    break;
           }
          bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
          bf->ChangeBF(Coll, cell, BFValues);

	  for(n=0;n<N_Comp;n++)
           {
	    value = 0;
	    for(l=0;l<N_LocDOF;l++)
	      value += BFValues[l] * Coeffs[DOF[l]+n*Length];
	    DoubleArray[N_Comp*VertexNumbers[m] + n] += value;
	   } 
          WArray[VertexNumbers[m]] +=1.;
          m++;
         } // endfor j
        } // endfor i


       // midle value
       l = 0;
       for(i=0;i<N_Vertices;i++)
        {
         for(j=0;j<N_Comp;j++)
          {
           if(WArray[i]>1.)
            DoubleArray[l] /= WArray[i];
	   l++;
          }
         } // endfor l


//        for(i=0;i<2*N_Vertices;i++)
//         cout << "Do[" << i << "]" << DoubleArray[i] << endl;


      dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
          <<FEVectFunctArray[k]->GetName()<<Dquot<<" NumberOfComponents="<<Dquot
          <<"3"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
      l=0;
      l=0;
      for(i=0;i<N_Vertices;i++)
       {
        for(j=0;j<N_Comp;j++)
	 dat << DoubleArray[N_Comp*i+j] << " ";
       }
      dat << endl;

      dat <<  "  </DataArray>"<<endl;
      dat << endl;

      for(j=0;j<N_Comp;j++)
       {
        dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
          <<FEVectFunctArray[k]->GetName()<<j<<Dquot<<" NumberOfComponents="<<Dquot
          <<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
        for(i=0;i<N_Vertices;i++)
	 dat << DoubleArray[i*N_Comp+j] << " ";

        dat << endl;
        dat <<  "  </DataArray>"<<endl;
        dat << endl;
      }

     } // for(k=0;k<N_Vec

    // write scalar variables into file
    if(N_LocVertices)
    for(k=0;k<N_ScalarVar;k++)
     {
      fespace = FEFunctionArray[k]->GetFESpace3D();
      Coeffs = FEFunctionArray[k]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();

      memset(DoubleArray, 0, SizeOfDouble*N_Vertices);
      memset(WArray, 0, SizeOfDouble*N_Vertices);
      m = 0;

      for(i=0;i<N_Elements;i++)
       {
        cell = Coll->GetCell(i);
        N_ = cell->GetN_Vertices();

        // find FE data for this element
        FE_ID = fespace->GetFE3D(i, cell);
        bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
        DOF = GlobalNumbers+BeginIndex[i];
        N_LocDOF = bf->GetDimension();
        for(j=0;j<N_;j++)
         {
          switch(cell->GetN_Vertices())
           {
	    //  Tetrahedron
            case 4: 
	      xi   = TetraCoords[3*j];
              eta  = TetraCoords[3*j+1];
              zeta = TetraCoords[3*j+2];
	    break;
	    // Hexahedron
	    case 8: 
              xi   = HexaCoords[3*j];
              eta  = HexaCoords[3*j+1];
              zeta = HexaCoords[3*j+2];
	    break;
           }
          bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
          bf->ChangeBF(Coll, cell, BFValues);
          value = 0;
          for(l=0;l<N_LocDOF;l++)
            value += BFValues[l] * Coeffs[DOF[l]];
          DoubleArray[VertexNumbers[m]] += value;
          WArray[VertexNumbers[m]] +=1.;
          m++;
        } // endfor j
       } // endfor i

      // non conforming
      for(i=0;i<N_Vertices;i++)
       {
        if(WArray[i]>1.)
         DoubleArray[i] /= WArray[i];
       }

      dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
          <<FEFunctionArray[k]->GetName()<<Dquot<<" NumberOfComponents="<<Dquot
          <<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
      for(j=0;j<N_Vertices;j++)
        dat << DoubleArray[j]<< " ";
       dat << endl;
      dat <<  "  </DataArray>"<<endl;

       dat << endl;

     }// for(k=0;k<N_ScalarVar

     
//     dat <<  "</PointData>"<<endl;
//     dat << endl;
// 
//     dat <<  "<CellData Scalars="<<Dquot<<"SubDomain"<<Dquot<<">"<<endl;
//     dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"Region"<<Dquot
//         <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
//     for(i=0;i<N_Elements;i++)
//       dat << (Coll->GetCell(i))->GetRegionID()   << " ";
//     dat <<  "  </DataArray>"<<endl;
//     dat <<  "</CellData>"<<endl;
// 
//  
     
    dat <<  "</PointData>"<<endl;
    dat << endl;

    dat <<  "<CellData Scalars="<<Dquot<<"SubDomain"<<Dquot<<">"<<endl;
#ifdef _MPI    
    dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"SubDomain"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    for(i=0;i<N_Elements;i++)     
      dat << (Coll->GetCell(i))->GetSubDomainNo() << " ";
    dat <<  "  </DataArray>"<<endl;
#endif    
    dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"RegionID"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
     for(i=0;i<N_Elements;i++)
      dat << (Coll->GetCell(i))->GetRegionID() << " ";       
    dat <<  "  </DataArray>"<<endl;    

    dat <<  "</CellData>"<<endl;

    dat <<  "</Piece>"<<endl;
    dat <<  "</UnstructuredGrid>"<<endl;
    dat <<  "</VTKFile>"<<endl;

    dat.close();



    delete [] FESpaceNumber;

  if(N_LocVertices)
   {
    delete [] NumberVertex;
    delete [] VertexNumbers;
    delete [] Vertices;
    delete [] DoubleArray;
    delete [] WArray;
    delete [] Coords;
   }
  } //else root
//  if(rank==1)


  return 0;
} //TOutput3D::Write_ParVTK

TOutput3D::~TOutput3D()
{
  int i;

  delete [] FESpaceArray;
  delete [] FEFunctionArray;
  delete [] FEVectFunctArray;
  delete [] ParameterValues;

  for(i=0;i<N_Parameters;i++)
    free((void*) ParameterDescription[i]);

  delete [] ParameterDescription;

  if (Data) delete Data;

  free(Name);
}

TOutput3D::TOutputData::~TOutputData()
{
  if (Nodes) delete [] Nodes;
  if (ConList) delete [] ConList;
  if (FEFuncValues) 
  {
    delete [] FEFuncValues[0];
    delete [] FEFuncValues;
  }
}

void TOutput3D::ComputeOutputData()
{
  TBaseCell *Cell;
  TVertex *Vert;
  int N_Cells, N_=0, counter;
  
  if (Data) delete Data;
  Data = new TOutputData();

  // reset clipboard
  N_Cells = Coll->GetN_Cells();
  for(int i=0; i<N_Cells;i++)
  {
    Cell = Coll->GetCell(i);
    N_ = Cell->GetN_Vertices();

    switch(N_)
    {
      case 4:
	Data->Type = TOutputData::TETRAHEDRON;
	break;
      case 8:
	Data->Type = TOutputData::BRICK;
	break;
    }

    for(int j=0;j<N_;j++)
    {
      Cell->GetVertex(j)->SetClipBoard(-1);
    }
  }

  // count vertices
  counter = 0;
  for(int i=0; i<N_Cells;i++)
  {
    Cell = Coll->GetCell(i);
    N_ = Cell->GetN_Vertices();

    for(int j=0;j<N_;j++)
    {
      if( (Vert = Cell->GetVertex(j))->GetClipBoard() == -1 )
      {
	Vert->SetClipBoard(counter++);
      } 
    }
  }
  Data->N_Nodes = counter;

  // find all vertices and put into an array, create conn list
  Data->Nodes = new TVertex* [Data->N_Nodes];
  Data->ConList = new int [N_*N_Cells];
  for(int i=0;i<N_Cells;i++)
  {
    Cell = Coll->GetCell(i);
    N_ = Cell->GetN_Vertices();
    for(int j=0;j<N_;j++)
    {
      counter = (Vert = Cell->GetVertex(j))->GetClipBoard();
      Data->Nodes[counter] = Vert;
      Data->ConList[i*N_+j] = counter+1;
    }
  }
    
  ComputeFEValues();

}

void TOutput3D::ComputeFEValues()
{
  TBaseCell *Cell;
  const TFESpace3D *fespace;
  TBaseFunct3D *BaseFunc;
  FE3D FeID;
  double x, y, z, *xi, *eta, *zeta;
  double *FECoeffs;
  double BaseFuncValues[MaxN_BaseFunctions3D];
  int N_ , N_Cells, N_BF, counter;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int N_Comps, Length, index, N_VectorTotal;

  double xi_tetra[]   = {0, 1, 0, 0};
  double eta_tetra[]  = {0, 0, 1, 0};
  double zeta_tetra[] = {0, 0, 0, 1};
  double xi_brick[]   = {-1,  1,  1, -1, -1,  1,  1, -1};
  double eta_brick[]  = {-1, -1,  1,  1, -1, -1,  1,  1};
  double zeta_brick[] = {-1, -1, -1, -1,  1,  1,  1,  1};

  N_VectorTotal = 4*N_VectorVar;

  N_ = 3 + N_ScalarVar + N_VectorTotal;
  Data->N_Data = Data->N_Nodes * N_;

  N_Cells = Coll->GetN_Cells();

  switch(Data->Type)
  {
    case TOutputData::TETRAHEDRON:
      xi = xi_tetra;
      eta = eta_tetra;
      zeta = zeta_tetra;
    case TOutputData::BRICK:
      xi = xi_brick;
      eta = eta_brick;
      zeta = zeta_brick;
      break;
  }

  double *tmp = new double [Data->N_Data];
  Data->FEFuncValues = new double* [N_];
  for(int i=0;i<N_;i++)
    Data->FEFuncValues[i] = tmp + i*Data->N_Nodes;

  // coords
  for(int i=0;i<Data->N_Nodes;i++)
  {
    Data->Nodes[i]->GetCoords(x,y,z);

    Data->FEFuncValues[0][i] = x;
    Data->FEFuncValues[1][i] = y;
    Data->FEFuncValues[2][i] = z;
  }

  // find fe values (scalar)
  for(int i=0;i<N_ScalarVar;i++)
  {
    counter = 0;
    fespace = FEFunctionArray[i]->GetFESpace3D();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    FECoeffs = FEFunctionArray[i]->GetValues();

    for(int j=0;j<N_Cells;j++)
    {
      Cell = Coll->GetCell(j);
      FeID = fespace->GetFE3D(j, Cell);
      N_ = Cell->GetN_Vertices();
      DOF = GlobalNumbers + BeginIndex[j];
      BaseFunc = TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID);
      N_BF = BaseFunc->GetDimension();

      for(int k=0;k<N_;k++)
      {
	double val = 0;
	
	if ( Cell->GetVertex(k)->GetClipBoard() == counter )
	{
	  BaseFunc->GetDerivatives(D000, xi[k], eta[k], zeta[k],
				   BaseFuncValues);
	  for(int l=0;l<N_BF;l++)
	  {
	    val += BaseFuncValues[l] * FECoeffs[DOF[l]];
	  }
	  Data->FEFuncValues[i+3][counter++] = val;
	}   
      }
    } // end for j (Cells)
  } // end for i (ScalarVar)

  // find fe values (vector)
  index = 3 + N_ScalarVar;
  for(int i=0;i<N_VectorVar;i++)
  {
    counter = 0;
    fespace = FEVectFunctArray[i]->GetFESpace3D();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    N_Comps = FEVectFunctArray[i]->GetN_Components();
    Length = FEVectFunctArray[i]->GetLength();
    FECoeffs = FEVectFunctArray[i]->GetValues();

    for(int j=0;j<N_Cells;j++)
    {	
      Cell = Coll->GetCell(j);
      FeID = fespace->GetFE3D(j, Cell);
      N_ = Cell->GetN_Vertices();
      DOF = GlobalNumbers + BeginIndex[j];
      BaseFunc = TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID);
      N_BF = BaseFunc->GetDimension();

      for(int k=0;k<N_;k++)
      {
	if ( Cell->GetVertex(k)->GetClipBoard() == counter )
	{
	  double abs = 0;

	  BaseFunc->GetDerivatives(D000, xi[k], eta[k], zeta[k],
				   BaseFuncValues);

	  for(int p=0;p<N_Comps;p++)
	  {
	    double val = 0;
	    for(int l=0;l<N_BF;l++)
	    {
	      val += BaseFuncValues[l] * FECoeffs[p*Length+DOF[l]];
	    }
	    abs += val*val;
	    Data->FEFuncValues[index+p][counter] = val;
	  } // end for p (Components)
	  Data->FEFuncValues[index+N_Comps][counter] = sqrt(abs);
	  counter++;
	}
	
      } // end for k (Vertices)
    } // end for j (Cells)
    index += N_Comps + 1;
  } // end for i (VectorVar)
}
