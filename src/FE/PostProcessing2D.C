#include <Database.h>
#include <MainUtilities.h>
#include <FEDatabase2D.h>
#include <PostProcessing2D.h>
#include <QuadAffin.h>
#include <TriaAffin.h>
#include <QuadBilinear.h>

using namespace std;

PostProcessing2D::PostProcessing2D(){
  //  init(param_db);
};

// Possible variant for handling more output objects in a single problem
// e.g. Stokes-Darcy
//void PostProcessing2D::PostProcessing2D(string basename){
//  init();
//  testcaseName = basename;
//};


void PostProcessing2D::init(const ParameterDatabase& param_db)
{
  writeVTK = param_db["write_vtk"];
  writeCASE = TDatabase::ParamDB->WRITE_CASE;

  string tfile = param_db["base_name"];
  testcaseName=tfile;
  string tdir = param_db["output_directory"];
  testcaseDir=tdir;

  Coll=NULL;
    
  // for time dependent problems
  dt = param_db["time_step_length"];
  t0 = param_db["time_start"];
  period = TDatabase::TimeDB->STEPS_PER_IMAGE;
  timeValues.clear();
  
}



/** add a FEFunction into this output object */
void PostProcessing2D::addFEFunction(const TFEFunction2D* fefunction)
{
  FEFunctionArray.push_back(fefunction);

  // check that FE functions have the same collection
  if (Coll==NULL) {
    Coll = fefunction->GetFESpace2D()->GetCollection();
  } else {
    if (Coll != fefunction->GetFESpace2D()->GetCollection()) {
      OutPut(" ** WARNING ** in file : " << __FILE__ 
	     << ", line " << __LINE__ 
	     << " new FE function might have a different collection " << endl);
    }
  }

}
/** add a FEVectFunct into this output object */
void PostProcessing2D::addFEVectFunct(const TFEVectFunct2D* fevectfunction)
{
  FEVectFunctArray.push_back(fevectfunction);

  // check that FE functions have the same collection
  if (Coll==NULL) {
    Coll = fevectfunction->GetFESpace2D()->GetCollection();
  } else {
    if (Coll != fevectfunction->GetFESpace2D()->GetCollection()) {
      OutPut(" ** WARNING ** in file : " << __FILE__ 
	     << ", line " << __LINE__ 
	     << " new FE function might have a different collection " << endl);
    }
  }
}


void PostProcessing2D::write(const char *name, int i, double _current_time)
{
  write((string)name, i, _current_time);
}


void PostProcessing2D::write(int i,double _current_time)
{
  std::ostringstream os;

  if(writeVTK) {
    os.seekp(std::ios::beg);
    if (i>=0) { 
      os << testcaseDir << "/" << testcaseName << i << ".vtk"<< ends;
    } else {
      os << testcaseDir << "/" << testcaseName << ".vtk"<< ends;
    }
    cout << " PostProcessing2D:: writing " << os.str() << endl;
    writeVtk(os.str().c_str());
  }

  if(TDatabase::ParamDB->WRITE_GNU)
  {
    os.seekp(std::ios::beg);
    os << testcaseDir << "/" << testcaseName << i << ".gnu"<<ends;
    cout << " Output2D:: writing " << os.str() << endl;
    writeGnuplot(os.str().c_str());
  }

  if(writeCASE) {
    // note: i<0 is used to avoid suffix in vtk output.
    // it shall be disregarded for case output
    if (i<0) i=0;
    // store new time step value
    timeValues.push_back(_current_time);
    if (i==0) {
      // write geometry only in the first iteration
      writeCaseGeo();
    }
    writeCaseVars(i);
    writeCaseFile();

  }

}


void PostProcessing2D::write(string basename, int i,double _current_time)
{
  std::ostringstream os;

  if(writeVTK) {
    os.seekp(std::ios::beg);
    os << basename << i << ".vtk"<< ends;
    cout << " PostProcessing2D:: writing " << os.str() << endl;
    writeVtk(os.str().c_str());
  }

  if(TDatabase::ParamDB->WRITE_GNU)
  {
    os.seekp(std::ios::beg);
    os << basename << i << ".gnu"<<ends;
    cout << " Output2D:: writing " << os.str() << endl;
    writeGnuplot(os.str().c_str());
  }

  if(writeCASE) {
    cout << " --- WARNING: the basename of output files is set during the initialization of Postprocessing2D --- " 
	 << endl;
    // store new time step value
    timeValues.push_back(_current_time);
    if (i==0) {
      // write geometry only in the first iteration
      writeCaseGeo();
    }
    writeCaseVars(i);
    writeCaseFile();

  }

}



  


/** write stored data into a gunplot file */
void PostProcessing2D::writeGnuplot(const char *name)
{
  cout << " void PostProcessing2D::writeGnuplot(const char *name) not yet implemented " << endl;
  cout << " in File: " << __FILE__ << " line: " << __LINE__ << endl;
  exit(1);

}


/** write stored data into a VTK file */
///@todo this functions can be simplified using the new features of Collection class
void PostProcessing2D::writeVtk(const char *name)
{
  
  int l,m;
  int edge, nsign;
  int N_Vertices, N_CellVertices, MaxN_VerticesPerCell;
  int  N_, N_Elements, N_LocVertices;
  int N_LocDOF, Length, N_Comp;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int *VertexNumbers, *NumberVertex;

  double xi, eta, value, t;
  const double *Coeffs;
  double value_y;
  double BFValues[MaxN_BaseFunctions2D];
  double *Coords;
  //this is needed to hanlde vector fields approximated
   //  with vector FE (scalar unknown postprocessed as vectors)
   
  double BFValuesOrig[MaxN_BaseFunctions2D];
  bool VectOutput = false;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;

  double QuadCoords[] = { -1, -1, 1, -1, 1, 1, -1, 1};
  double TriaCoords[] = { 0, 0, 1, 0,  0, 1};

  const TFESpace2D *fespace;
  TBaseFunct2D *bf;
  FE2D FE_ID;
  int BaseVectDim;

  MaxN_VerticesPerCell = 4;                       // 2D case

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << " *** ERROR, file: " << __FILE__ << ", line: " << __LINE__ << ": cannot open file for output" << endl;
    exit(1);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(9);


  // determine data for vtk file
  N_Elements=Coll->GetN_Cells();
  //
  N_LocVertices=0;
  for(int i=0;i<N_Elements;i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
  }
  TVertex **Vertices =new TVertex*[N_LocVertices];
  N_=0;

  for(int i=0;i<N_Elements;i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    for(int j=0;j< cell->GetN_Vertices(); j++)
    {
      Vertices[N_]=cell->GetVertex(j);
      N_++;
    }
  }

  // check for discontinuous scalar variables. In such a case write a new file
  // for this variable. ParaView will really display a discontinuous function
  // instead of projecting it onto P1/Q1 space. However there are some
  // drawbacks.
  for(unsigned int i=0;i<FEFunctionArray.size();i++)
  {
    if ( FEFunctionArray[i]->GetFESpace2D()->IsDGSpace() ) {
      writeVtkDiscontinuous(name,N_LocVertices,Vertices);
      break;
    }
  }

  sort(Vertices, N_);
  TVertex *Last, *Current;

  Last=NULL;
  N_Vertices=0;
  for(int i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
  {
    N_Vertices++;
    Last=Current;
  }

  Coords=new double[2*N_Vertices];
  VertexNumbers=new int[N_LocVertices];
  NumberVertex=new int[N_LocVertices];
  Last=NULL;
  N_=0; 
  int k1=-1;
#ifdef __3D__
  double z = 0.;
#endif
  for(int i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
#ifdef __3D__
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1], z);
#else
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1]);
#endif
       k1++;
      N_ += 2;
      Last=Current;
    }
    NumberVertex[i]=k1;
  }

  m=0;
  for(int i=0;i<N_Elements;i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    for(int j=0;j<cell->GetN_Vertices(); j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      l=getIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      m++;
    }                                             // endfor j
  }                                               //endfor i



  dat << "# vtk DataFile Version 4.2" << endl;
  dat << "file created by ParMooN"
      << " Time < " << TDatabase::TimeDB->CURRENTTIME <<" >" << " L < " << TDatabase::ParamDB->REACTOR_P29 <<" >" <<endl;

  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl;
  dat << "POINTS " << N_Vertices << " float" << endl;
  N_=0;
  for(int i=0;i<N_Vertices;i++)
  {
    dat << Coords[N_] << " " <<  Coords[N_+1] << " ";
    dat << double(0) << endl;
    N_ +=2;
  }
  dat << endl;
  dat << "CELLS " << N_Elements << " " <<  N_Elements+N_LocVertices << endl;
  l=0;
  for(int i=0;i<N_Elements;i++)
  {
    N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
    dat <<  N_CellVertices << " ";
    for(int j=0;j<N_CellVertices;j++)
    {
      dat << VertexNumbers[l] << " ";
      l++;
    }
    dat << endl;
  }
  dat << endl;
  dat << "CELL_TYPES " << N_Elements << endl;
  for(int i=0;i<N_Elements;i++)
  {
    N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
    switch(N_CellVertices)
    {
      case 4: dat << 9 << " ";
      break;
      case 3: dat << 5 << " ";
      break;
    }
  }
  dat << endl << endl;
  dat << "POINT_DATA " << N_Vertices << endl;

  // function values
  double *WArray, *DoubleArray;
  DoubleArray = new double[2*N_Vertices];
  WArray = new double[N_Vertices];

  // write scalar variables into file
  for(unsigned int k=0;k<FEFunctionArray.size();k++)
  {
    fespace = FEFunctionArray[k]->GetFESpace2D();
    Coeffs = FEFunctionArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    // get dimension of basis functions
    TBaseCell *cell = Coll->GetCell(0);
    FE_ID = fespace->GetFE2D(0, cell);
    bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
    BaseVectDim = bf->GetBaseVectDim();
    N_Comp = BaseVectDim;

    memset(DoubleArray, 0, SizeOfDouble*N_Comp*N_Vertices);
    memset(WArray, 0, SizeOfDouble*N_Vertices);
    m = 0;

    // set to TRUE if basis functions are vectors
    if (BaseVectDim>1)  VectOutput = true;
    else VectOutput = false;

    for(int i=0;i<N_Elements;i++)
    {
      TBaseCell *cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FE_ID);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FE_ID);
      F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
      switch(RefTrans)
      {
        case TriaAffin:
          ((TTriaAffin*)F_K)->SetCell(cell);
          break;
        case QuadAffin:
          ((TQuadAffin*)F_K)->SetCell(cell);
          break;
        case QuadBilinear:
          ((TQuadBilinear*)F_K)->SetCell(cell);
          break;
        default: 
          cout << "Output2D():: no such reference transformation allowed" 
               << endl;
          break;
      }                                           // endswitch

      for(int j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
          case 3:
            xi = TriaCoords[2*j];
            eta = TriaCoords[2*j+1];
            break;

          case 4:
            xi = QuadCoords[2*j];
            eta = QuadCoords[2*j+1];
            break;
        }
        bf->GetDerivatives(D00, xi, eta, BFValues);

        value = 0;
        value_y = 0;
        if (VectOutput)
        {
          // apply Piola transform 
          switch(RefTrans)
          {
            case TriaAffin:
            case QuadAffin:
              F_K->PiolaMapOrigFromRef(N_LocDOF,BFValues,BFValuesOrig);
              break;
            case QuadBilinear: 
              // for non affine reference transformations one needs to know the 
              // point, because the determinant is not constant in this case.
              ((TQuadBilinear*)F_K)->PiolaMapOrigFromRefNotAffine(
                                        N_LocDOF,BFValues,BFValuesOrig,xi,eta);
              break;
            default:
              cout << "Output2D():: no such reference transformation allowed" 
               << endl;
          }
          for(l=0;l<N_LocDOF;l++)
          {
            // change sign of basis functions according to global normal
            edge=TFEDatabase2D::GetFE2D(FE_ID)->GetFEDesc2D()->GetJointOfThisDOF(l);
            if (edge != -1)
            {
              nsign = cell->GetNormalOrientation(edge);
            }
            else
              nsign=1;

            value += BFValuesOrig[l] * Coeffs[DOF[l]]*nsign;
            value_y += BFValuesOrig[N_LocDOF+l] * Coeffs[DOF[l]]*nsign;
          }
          DoubleArray[N_Comp*VertexNumbers[m] + 0] += value;
          DoubleArray[N_Comp*VertexNumbers[m] + 1] += value_y;
        } else
        {                                         // standard
          for(l=0;l<N_LocDOF;l++)
          {
            value += BFValues[l] * Coeffs[DOF[l]];
          }
          DoubleArray[VertexNumbers[m]] += value;
        }
        WArray[VertexNumbers[m]] +=1.;
        m++;
      }                                           // endfor j
    }                                             // endfor i

    // non conforming
    if (VectOutput)
    {
      l = 0;
      for(int i=0;i<N_Vertices;i++)
      {
        for(int j=0;j<N_Comp;j++)
        {
          if(WArray[i]!=0.)
            DoubleArray[l] /= WArray[i];
          l++;
        }
      }                                           // endfor i
    } else
    {
      for(int i=0;i<N_Vertices;i++)
      {
        if(WArray[i]!=0.)
        {
          DoubleArray[i] /= WArray[i];
        }
      }
    }

    if (!VectOutput)
    {
      // standard output writing
      dat << "SCALARS " << FEFunctionArray[k]->GetName();
      dat << " float"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(int j=0;j<N_Vertices;j++)
      {
        dat << DoubleArray[j] << endl;
      }
      dat << endl;
      dat << endl;
    }

    // write output as vector
    if (VectOutput)
    {
      // scalar components
      for(int j=0;j<N_Comp;j++)
      {
        dat << "SCALARS " << FEFunctionArray[k]->GetName() << j;
        dat << " float"<< endl;
        dat << "LOOKUP_TABLE " << "default" << endl;
        for(int i=0;i<N_Vertices;i++)
        {
          dat << DoubleArray[i*N_Comp+j] << endl;
        }
        dat << endl << endl;
      }

      // absolute value
      dat << "SCALARS " << "|" << FEFunctionArray[k]->GetName() << "|";
      dat << " float"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      l=0;
      for(int i=0;i<N_Vertices;i++)
      {
        t=0;
        for(int j=0;j<N_Comp;j++)
        {
          t+=DoubleArray[l]*DoubleArray[l];
          l++;
        }
        dat << sqrt(t)<< endl;
      }
      dat << endl << endl;

      dat << "VECTORS " << FEFunctionArray[k]->GetName();
      dat << " float"<< endl;

      l=0;
      for(int i=0;i<N_Vertices;i++)
      {
        for(int j=0;j<N_Comp;j++)
        {
          dat << DoubleArray[N_Comp*i+j] << " ";
        }
        dat << double(0) << " " << endl;
      }
      dat << endl;

      VectOutput = false;
    }

  }                                               // endfor k

  for(unsigned int k=0;k<FEVectFunctArray.size();k++) {
    fespace = FEVectFunctArray[k]->GetFESpace2D();
    N_Comp = FEVectFunctArray[k]->GetN_Components();
    Length = FEVectFunctArray[k]->GetLength();
    Coeffs = FEVectFunctArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    // cout << "N_Comp  " << N_Comp << endl;
    memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
    memset(WArray, 0, SizeOfDouble*N_Vertices);
    m = 0;

    for(int i=0;i<N_Elements;i++) {
      TBaseCell *cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(int j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices()) {
	case 3:
	  xi = TriaCoords[2*j];
	  eta = TriaCoords[2*j+1];
	  break;
	  
	case 4:
	  xi = QuadCoords[2*j];
	  eta = QuadCoords[2*j+1];
	  break;
        }
        bf->GetDerivatives(D00, xi, eta, BFValues);

        for(int n=0;n<N_Comp;n++) {
          value = 0;
          for(l=0;l<N_LocDOF;l++)
            value += BFValues[l] * Coeffs[DOF[l]+n*Length];
          DoubleArray[N_Comp*VertexNumbers[m] + n] += value;
        }
        WArray[VertexNumbers[m]] +=1.;
        m++;
      }                                           // endfor j
    }                                             // endfor i

    // mean value

    l = 0;
    for(int i=0;i<N_Vertices;i++) {
      for(int j=0;j<N_Comp;j++)
	{
	  if(WArray[i]!=0.)
	    DoubleArray[l] /= WArray[i];
	  l++;
	}
    }                                             // endfor l

    for(int j=0;j<N_Comp;j++)
    {
      dat << "SCALARS " << FEVectFunctArray[k]->GetName() << j;
      dat << " float"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(int i=0;i<N_Vertices;i++)
      {
        dat << DoubleArray[i*N_Comp+j] << endl;
      }
      dat << endl << endl;
    }

    // ***************
    // absolute value of vector variables
    // ***************
    dat << "SCALARS " << "|" << FEVectFunctArray[k]->GetName() << "|";
    dat << " float"<< endl;
    dat << "LOOKUP_TABLE " << "default" << endl;
    l=0;
    for(int i=0;i<N_Vertices;i++) {
      t=0;
      for(int j=0;j<N_Comp;j++) {
        t+=DoubleArray[l]*DoubleArray[l];
        l++;
      }
      dat << sqrt(t)<< endl;
    }
    dat << endl << endl;

    // ***************
    // VECTORS
    // ***************
    dat << "VECTORS " << FEVectFunctArray[k]->GetName();
    dat << " float"<< endl;

    l=0;
    for(int i=0;i<N_Vertices;i++) {
      for(int j=0;j<N_Comp;j++) {
        dat << DoubleArray[N_Comp*i+j] << " ";
      }
      dat << double(0) << " " << endl;
    }
    dat << endl;
  }                                               // endfor k

  dat << endl;

  delete [] NumberVertex;
  delete [] VertexNumbers;
  delete [] Vertices;
  delete [] DoubleArray;
  delete [] WArray;
  delete [] Coords;

  dat.close();
  //   cout << endl;
  if( TDatabase::ParamDB->SC_VERBOSE > 0 )
    OutPut("wrote output into vtk file: " << name << endl);


}




/** @brief
    writes an extra vtk-file for (scalar) discontinuous functions. ParaView can
    then display discontinuous data. The "POINTS" in the resulting vtk-file are the
    vertices of the mesh, but every vertex is put here as many times as there are
    cells this vertex belongs to. That means if a vertex belongs to 4 cells in the
    mesh it will appear 4 times in the list "POINTS" in the resulting vtk-file.
    The input "TVertex **Vertices" should be in this pattern. It can be generated by
    N_Elements=Coll->GetN_Cells();
    //
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
    
    as is done in writeVtk(cont char *name).
    Also note that in each element the function is projected onto P1/Q1.
    @warning This destroys the topology of the mesh. 
    @warning Some filters in ParaView might work incorrectly. 
    Warp by scalar works though. 
    The way this is done here is not very elegant, but I couldn't find a better
    solution.
    
*/
void PostProcessing2D::writeVtkDiscontinuous(const char *fileName,
int N_LocVertices, TVertex **Vertices)
{
  double x,y;                // coordinates of a vertex
  int N_Elements;            // number of all elements in this mesh
  int N_CellVertices;        // number of all vertices in this cell
  int i,j,k,l;               // loop variables
  const TFEFunction2D *fefunction; // this function, which is discontinuous
  TBaseCell *current_cell;   
  double *function_value;    // value of function at a particular vertex
  double **allValues;        // in case of vector valued basis functions, store
                             // all values of all components
  char Disc[80];             // the new file name
  strcpy(Disc,fileName);     // copy the file name ...
  strcat(Disc,"_disc.vtk");  // ... add a string to the output file name

  N_Elements = Coll->GetN_Cells();

  std::ofstream dat(Disc);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    exit(-1);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(9);

  dat << "# vtk DataFile Version 4.2" << endl;
  dat << "file created by ParMooN." << endl;

  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl;
  dat << "POINTS " << N_LocVertices << " float" << endl;

  for(i=0;i<N_LocVertices;i++)
  {
#ifdef __3D__
    //Vertices[i]->GetCoords(x,y,z);
    OutPut("writeVtkDiscontinuous has to be checked for 3D"<<endl);
    exit(4711);
#else
    Vertices[i]->GetCoords(x,y);
#endif
    dat << x << " " <<  y << " " << double(0) << endl;
    // in 2D: set last entry 0
  }
  dat << endl;
  dat << "CELLS " << N_Elements << " " <<  N_Elements+N_LocVertices << endl;
  l=0;
  for(i=0;i<N_Elements;i++)
  {
    current_cell = Coll->GetCell(i);
    N_CellVertices = current_cell->GetN_Vertices();
    dat <<  N_CellVertices << " ";
    for(j=0;j<N_CellVertices;j++)
    {
      dat << l << " ";
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
      case 4: dat << 9 << " ";
      break;
      case 3: dat << 5 << " ";
      break;
    }
  }
  dat << endl << endl;
  dat << "POINT_DATA " << N_LocVertices << endl;
  for(unsigned int space_number=0; space_number<FEFunctionArray.size(); space_number++)
  {
    fefunction = FEFunctionArray[space_number];
    if(fefunction->GetFESpace2D()->IsDGSpace() != 1)
      continue;
    // this is a discontinuous space

    int BaseVectDim =
      TFEDatabase2D::GetFE2D(
      fefunction->GetFESpace2D()->GetFE2D(0,Coll->GetCell(0)))
      ->GetBaseFunct2D()->GetBaseVectDim();   // ugly, but we need to know this
    if (BaseVectDim==1)
    {
      dat << endl << endl;
      dat << "SCALARS " << fefunction->GetName();
      dat << " float" << endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      function_value = new double[1];  // only need to allocate this first entry
      for(i=0;i<N_Elements;i++)
      {
        current_cell = Coll->GetCell(i);
        N_CellVertices=current_cell->GetN_Vertices();
        for(j=0;j<N_CellVertices;j++)
        {
#ifdef __3D__
          ;
#else
          current_cell->GetVertex(j)->GetCoords(x,y);
#endif
          fefunction->FindValueLocal(current_cell,i,x,y,function_value);
          dat << function_value[0] << endl;
        }
      }
    }
    else if(BaseVectDim==2)
    {
      // find values for all components
      function_value = new double[2];             // 2==BaseVectDim
      allValues = new double*[BaseVectDim];       // 2==BaseVectDim
      for(l=0; l<BaseVectDim; l++)
        allValues[l] = new double[N_LocVertices];
      k=0;
      for(i=0;i<N_Elements;i++)
      {
        current_cell = Coll->GetCell(i);
        N_CellVertices=current_cell->GetN_Vertices();
        for(j=0;j<N_CellVertices;j++)
        {
#ifdef __3D__
          ;
#else
          current_cell->GetVertex(j)->GetCoords(x,y);
#endif
          // FindValueLocal includes the necessary sign changes due to global
          // normals (for Raviart-Thomas elements)
          fefunction->FindValueLocal(current_cell,i,x,y,function_value);
          for(l=0; l<BaseVectDim; l++)
            allValues[l][k] = function_value[l];
          k++;
        }
      }
      // write the function values to the vtk-file
      for(l=0; l<BaseVectDim; l++)
      {
        dat << endl << endl;
        dat << "SCALARS " << fefunction->GetName() << l;
        dat << " float" << endl;
        dat << "LOOKUP_TABLE " << "default" << endl;
        k=0;
        for(i=0;i<N_Elements;i++)
        {
          N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
          for(j=0;j<N_CellVertices;j++)
          {
            dat << allValues[l][k] << endl;
            k++;
          }
        }
      }
    }
    else
      OutPut("TOutput2D::writeVtkDiscontinuous: Basis functions of dimension "
        << BaseVectDim << " are not supported." << endl);
  }

}




/** @brief Output in .case format

   1 file .case (main file)
   1 file .geo (mesh)
   scalar and vector output files at each time step
   Note that currently we only write the P1 projection (node values) of the results
   However, the case format allows also element-wise and higher order
   visualizations (but not with paraview) 
   @todo implement visualization for different element type
*/
void PostProcessing2D::writeCaseFile()
{
  ofstream casf;
  string filename = testcaseDir+"/"+testcaseName;
  string casename= filename + ".case";
  cout << " ** PostProcessing2D::writeCaseFile - write " << casename << endl;
  casf.open(casename.c_str());
  
  string geoname = testcaseName + ".geo";
  
  casf << "FORMAT\n";
  casf << "type: ensight\n";
  casf << "GEOMETRY\n";
  casf << "model: 1 " << geoname << endl;
  casf << "VARIABLE\n";
  for (unsigned int j=0; j<FEFunctionArray.size(); j++) {
    string sclname =  "scalar";
    casf << "scalar per node: 1 " << FEFunctionArray[j]->GetName()
	 << " " << testcaseName << "_"
	 << FEFunctionArray[j]->GetName() << ".****.scl" << endl;
  }
  for (unsigned int j=0; j<FEVectFunctArray.size(); j++){
    casf << "vector per node: 1 "<< FEVectFunctArray[j]->GetName() << " " 
	 << testcaseName << "_" 
	 << FEVectFunctArray[j]->GetName() << ".****.vct" << endl;  
  }
  // time values
  unsigned int N = timeValues.size();
  casf << "TIME\n";  
  casf << "time set: 1\n";
  casf << "number of steps: " << N << endl;
  casf << "filename start number: 0\n";
  casf << "filename increment: " << period << "\n";  
  casf << "time values:\n";
  for(unsigned int i=0;i<N;i++){
    casf.precision(5);	
    casf.width(12);    
    casf << timeValues[i];
    if( (i && !(i%5)) || (i==N) )
       casf << endl;
    else casf << " ";
  }
  casf.close();
}

// Geometry
void PostProcessing2D::writeCaseGeo()
{
  int dimension=2;
  
  Coll->createElementLists();
  unsigned int N_Vertices = Coll->NodesReferences.size();
  
  // write geo file
  ofstream f;
  string filename = testcaseDir+"/"+testcaseName;
  string geoname = filename + ".geo";
  f.open(geoname.c_str());

  f << testcaseName << " geometry\n";
  f << "t = 0.0000" << "\n";
  f << "node id assign" << "\n";
  f << "element id assign" << "\n";
  f << "coordinates" << "\n";
  f.width(8);
  // write points
  f << N_Vertices << endl;
  for(unsigned int i=0;i<N_Vertices;i++) {
    f.setf(ios_base::scientific);
    f.precision(5);	
    f.width(12);
    f <<  Coll->NodesCoords[i*dimension];
    f.setf(ios_base::scientific);
    f.precision(5);	
    f.width(12);      
    f  << Coll->NodesCoords[i*dimension+1]; 
    f.setf(ios_base::scientific);
    f.precision(5);	
    f.width(12);
    if (dimension==2) {
      f << 0.;
    } else {
      f  << Coll->NodesCoords[i*dimension+2]; 
    }
    f << endl;
  }

  // write elements
  int nVE = Coll->GetCell(0)->GetN_Vertices();
  string ensight_type;

  /// @warning this works now only for a single (bulk) domain
  unsigned int nParts = 1; // number of (bulk) subdomains
  string partName = "inner"; // name of collection subdomain
  
  for(unsigned int ifig=1;ifig<=nParts;ifig++){ 
    f << "part";
    f.width(8);
    f << ifig << endl;
    f << partName << endl; // name of the figure
    switch(dimension){
    case 2:
      switch(nVE){
      case 3:
	ensight_type="tria3";
	break;
      case 4:
	ensight_type="quad4";
	break;
      }
      break;
    case 3:
      cout << " ERROR file " << __FILE__ << ", line " 
	   << __LINE__ << " writeCaseGeo currently only for 2D output " << endl;
      break;
    default:
      ensight_type = ""; 
    }
    f << ensight_type << endl;
    
    f.width(8);
    f << Coll->GetN_Cells() << endl;
    for(int k=0;k<Coll->GetN_Cells();k++){
      for(int i=0;i<nVE;i++){
	f.width(8);
	f << Coll->ElementNodes[nVE*k+i];
      }
      f << endl;
    }// k
  }//ifig
  
  // write boundary elements
  unsigned int nBdParts = 1; // number of subdomains
  ///@warning this works now only for a single boundary domain
  partName = "boundary"; // name of collection subdomain
  int nVertexPerFace = 2;
  for(unsigned int ifig=1;ifig<=nBdParts;ifig++){ 
    f << "part";
    f.width(8);
    f << nParts+ifig << endl;
    f << partName << endl; // name of the figure
    ///@warning only for P1 in 1D
    ensight_type="bar2";
    f << ensight_type << endl;
    
    f.width(8);
    f << Coll->BdFacesReferences.size() << endl;
    for(unsigned int k=0;k<Coll->BdFacesReferences.size();k++){
      for(int i=0;i<nVertexPerFace;i++){
    	f.width(8);
   	f << Coll->BdFacesNodes[nVertexPerFace*k+i];
     }
     f << endl;
    }// k
  }//ifig

}


void PostProcessing2D::writeCaseVars(int iter)
{
  string filename = testcaseDir+"/"+testcaseName;
  char numstr[4];
  sprintf(numstr,"%d",iter);
  string number = "000"+string(numstr);
  if (iter>9)
    number = "00" + string(numstr);
  if (iter>99)
    number = "0" + string(numstr);
  if (iter>999)
    number = string(numstr);

  string ensight_type;
  int nParts = 1;
  int nVE = Coll->GetCell(0)->GetN_Vertices();
  int n;
  ///@todo do not recreate lists if they exist already
  Coll->createElementLists();

  // scalars
  for (unsigned int i=0;i<FEFunctionArray.size(); i++){
   vector<double> uP1;
   FEFunctionArray[i]-> computeNodeValues(uP1);

    // write scalars
    ofstream sclf;
    ostringstream sclname;
    sclname << filename << "_" << FEFunctionArray[i]->GetName()
				     <<  "." << number << ".scl";
    string fname = sclname.str();
    cout << " ** PostProcessing2D::write File - write " << sclname.str() << endl;
    sclf.open(fname.c_str());
    sclf << FEFunctionArray[i]->GetName() << " step = " << iter << endl; 
  
    for(int ifig=1;ifig<=nParts;ifig++){
     n = 0;
     for(unsigned int k=0;k< Coll->NodesReferences.size(); k++){
       double x = uP1[k];
       sclf.setf(ios_base::scientific);	
       sclf.precision(5);	
       sclf.width(12);
       if (fabs(x)<1e-18) x = 0.;
       sclf << x;
       ++n;
       if ( n == 6 ) {
	 sclf << std::endl;
	 n=0;
       }	
     }// k
     
     sclf << endl;
    }//ifig
    
  }

  // vectors
  for(unsigned int i=0;i<FEVectFunctArray.size();i++) {
    int nComp = FEVectFunctArray[i]->GetN_Components();
    vector< vector<double> > uP1vect;
    uP1vect.resize(nComp);
    for(int j=0;j<nComp;j++) {
      FEVectFunctArray[i]->GetComponent(j)->computeNodeValues(uP1vect[j]);
    }
    
    ofstream vctf;
    ostringstream vctname;
    vctname << filename << "_" << FEVectFunctArray[i]->GetName()
				     <<  "." << number << ".vct";
    string fname = vctname.str();
    cout << " ** PostProcessing2D::write File - write " << vctname.str() << endl;
    vctf.open(fname.c_str());
    vctf <<  FEVectFunctArray[i]->GetName() << " step = " << iter << endl; 
  
    for(int ifig=1;ifig<=nParts;ifig++){
      n = 0;
      for(int k=0;k<Coll->NodesReferences.size();k++){
	double x,y,z;
	x = uP1vect[0][k];
	vctf.setf(ios_base::scientific);	
	vctf.precision(5);	
	vctf.width(12);
	if (fabs(x)<1e-18) x = 0.;
	vctf << x;
	++n;
	if ( n == 6 ) {
	  vctf << std::endl;
	  n=0;
	}	
	
	y = 0.;
	if (nComp>1) {
	  y = uP1vect[1][k];
	}
	vctf.setf(ios_base::scientific);	
	vctf.precision(5);	
	vctf.width(12);
	if (fabs(y)<1e-18) y = 0.;
	vctf << y;
	++n;
	if ( n == 6 ) {
	  vctf << std::endl;
	  n=0;
	}	
	
	z = 0.;
	if (nComp>2) {
	  double z = uP1vect[2][k];
	}
	vctf.setf(ios_base::scientific);	
	vctf.precision(5);	
	vctf.width(12);
	if (fabs(z)<1e-18) z = 0.;
	vctf << z;
	++n;
	if ( n == 6 ) {
	  vctf << std::endl;
	  n=0;
	}	
	
	
      }// k
      
      vctf << endl;
    }//ifig
   
  }

}



/*
void PostProcessing2D::writeVtkNew(const char *name)
{
  
  int l,m;
  int edge, nsign;
  int N_Vertices, N_CellVertices, MaxN_VerticesPerCell;
  int  N_, N_Elements, N_LocVertices;
  int N_LocDOF, Length, N_Comp;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int *VertexNumbers, *NumberVertex;

  double xi, eta, value, t;
  const double *Coeffs;
  double value_y;
  double BFValues[MaxN_BaseFunctions2D];
  double *Coords;
  //this is needed to hanlde vector fields approximated
   //  with vector FE (scalar unknown postprocessed as vectors)
   
  double BFValuesOrig[MaxN_BaseFunctions2D];
  bool VectOutput = false;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;

  double QuadCoords[] = { -1, -1, 1, -1, 1, 1, -1, 1};
  double TriaCoords[] = { 0, 0, 1, 0,  0, 1};

  const TFESpace2D *fespace;
  TBaseFunct2D *bf;
  FE2D FE_ID;
  int BaseVectDim;

  MaxN_VerticesPerCell = 4;                       // 2D case

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << " *** ERROR, file: " << __FILE__ << ", line: " << __LINE__ << ": cannot open file for output" << endl;
    exit(1);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(9);


  // determine data for vtk file
  N_Elements=Coll->GetN_Cells();
  //
  N_LocVertices=0;
  for(int i=0;i<N_Elements;i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
  }
  TVertex **Vertices =new TVertex*[N_LocVertices];
  N_=0;

  for(int i=0;i<N_Elements;i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    for(int j=0;j< cell->GetN_Vertices(); j++)
    {
      Vertices[N_]=cell->GetVertex(j);
      N_++;
    }
  }

  // check for discontinuous scalar variables. In such a case write a new file
  // for this variable. ParaView will really display a discontinuous function
  // instead of projecting it onto P1/Q1 space. However there are some
  // drawbacks.
  for(unsigned int i=0;i<FEFunctionArray.size();i++)
  {
    if ( FEFunctionArray[i]->GetFESpace2D()->IsDGSpace() ) {
      writeVtkDiscontinuous(name,N_LocVertices,Vertices);
      break;
    }
  }

  sort(Vertices, N_);
  TVertex *Last, *Current;

  Last=NULL;
  N_Vertices=0;
  for(int i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
  {
    N_Vertices++;
    Last=Current;
  }

  Coords=new double[2*N_Vertices];
  VertexNumbers=new int[N_LocVertices];
  NumberVertex=new int[N_LocVertices];
  Last=NULL;
  N_=0; 
  int k1=-1;
#ifdef __3D__
  double z = 0.;
#endif
  for(int i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
#ifdef __3D__
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1], z);
#else
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1]);
#endif
       k1++;
      N_ += 2;
      Last=Current;
    }
    NumberVertex[i]=k1;
  }

  m=0;
  for(int i=0;i<N_Elements;i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    for(int j=0;j<cell->GetN_Vertices(); j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      l=getIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      m++;
    }                                             // endfor j
  }                                               //endfor i


  // ------ start writing file ---------------
  Coll->createElementLists();
  unsigned int nNodes = Coll->NodesReferences.size();
  unsigned int nElements = Coll->ElementReferences.size();
  dat << "# vtk DataFile Version 4.2" << endl;
  dat << "file created by ParMooN"
      << " Time < " << TDatabase::TimeDB->CURRENTTIME <<" >" << " L < " << TDatabase::ParamDB->REACTOR_P29 <<" >" <<endl;
  ///@todo what is Database::ParamDB->REACTOR_P29?
  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl;

  dat << "POINTS " << nNodes << " float" << endl;
  for(unsigned int i=0;i<nNodes;i++)
  {
    ///@attention this line now depend on dimension (dim=2)
    dat << Coll->NodesCoords[2*i] << " " <<  Coords[2*i+1] << " ";
    // note: we write a 3D file anyway
    dat << double(0) << endl;
  }
  dat << endl;
  
  dat << "CELLS " << nElements << " " <<  nElements+N_LocVertices << endl;
  l=0;
  for(unsigned int i=0;i<nElements;i++)
  {
    N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
    dat <<  N_CellVertices << " ";
    for(unsigned int j=0;j<N_CellVertices;j++)
    {
      dat << VertexNumbers[l] << " ";
      l++;
    }
    dat << endl;
  }
  dat << endl;
  dat << "CELL_TYPES " << N_Elements << endl;
  for(int i=0;i<N_Elements;i++)
  {
    N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
    switch(N_CellVertices)
    {
      case 4: dat << 9 << " ";
      break;
      case 3: dat << 5 << " ";
      break;
    }
  }
  dat << endl << endl;
  dat << "POINT_DATA " << N_Vertices << endl;

  // function values
  double *WArray, *DoubleArray;
  DoubleArray = new double[2*N_Vertices];
  WArray = new double[N_Vertices];

  // write scalar variables into file
  for(unsigned int k=0;k<FEFunctionArray.size();k++)
  {
    fespace = FEFunctionArray[k]->GetFESpace2D();
    Coeffs = FEFunctionArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    // get dimension of basis functions
    TBaseCell *cell = Coll->GetCell(0);
    FE_ID = fespace->GetFE2D(0, cell);
    bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
    BaseVectDim = bf->GetBaseVectDim();
    N_Comp = BaseVectDim;

    memset(DoubleArray, 0, SizeOfDouble*N_Comp*N_Vertices);
    memset(WArray, 0, SizeOfDouble*N_Vertices);
    m = 0;

    // set to TRUE if basis functions are vectors
    if (BaseVectDim>1)  VectOutput = true;
    else VectOutput = false;

    for(int i=0;i<N_Elements;i++)
    {
      TBaseCell *cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FE_ID);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FE_ID);
      F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
      switch(RefTrans)
      {
        case TriaAffin:
          ((TTriaAffin*)F_K)->SetCell(cell);
          break;
        case QuadAffin:
          ((TQuadAffin*)F_K)->SetCell(cell);
          break;
        case QuadBilinear:
          ((TQuadBilinear*)F_K)->SetCell(cell);
          break;
        default: 
          cout << "Output2D():: no such reference transformation allowed" 
               << endl;
          break;
      }                                           // endswitch

      for(int j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
          case 3:
            xi = TriaCoords[2*j];
            eta = TriaCoords[2*j+1];
            break;

          case 4:
            xi = QuadCoords[2*j];
            eta = QuadCoords[2*j+1];
            break;
        }
        bf->GetDerivatives(D00, xi, eta, BFValues);

        value = 0;
        value_y = 0;
        if (VectOutput)
        {
          // apply Piola transform 
          switch(RefTrans)
          {
            case TriaAffin:
            case QuadAffin:
              F_K->PiolaMapOrigFromRef(N_LocDOF,BFValues,BFValuesOrig);
              break;
            case QuadBilinear: 
              // for non affine reference transformations one needs to know the 
              // point, because the determinant is not constant in this case.
              ((TQuadBilinear*)F_K)->PiolaMapOrigFromRefNotAffine(
                                        N_LocDOF,BFValues,BFValuesOrig,xi,eta);
              break;
            default:
              cout << "Output2D():: no such reference transformation allowed" 
               << endl;
          }
          for(l=0;l<N_LocDOF;l++)
          {
            // change sign of basis functions according to global normal
            edge=TFEDatabase2D::GetFE2D(FE_ID)->GetFEDesc2D()->GetJointOfThisDOF(l);
            if (edge != -1)
            {
              nsign = cell->GetNormalOrientation(edge);
            }
            else
              nsign=1;

            value += BFValuesOrig[l] * Coeffs[DOF[l]]*nsign;
            value_y += BFValuesOrig[N_LocDOF+l] * Coeffs[DOF[l]]*nsign;
          }
          DoubleArray[N_Comp*VertexNumbers[m] + 0] += value;
          DoubleArray[N_Comp*VertexNumbers[m] + 1] += value_y;
        } else
        {                                         // standard
          for(l=0;l<N_LocDOF;l++)
          {
            value += BFValues[l] * Coeffs[DOF[l]];
          }
          DoubleArray[VertexNumbers[m]] += value;
        }
        WArray[VertexNumbers[m]] +=1.;
        m++;
      }                                           // endfor j
    }                                             // endfor i

    // non conforming
    if (VectOutput)
    {
      l = 0;
      for(int i=0;i<N_Vertices;i++)
      {
        for(int j=0;j<N_Comp;j++)
        {
          if(WArray[i]!=0.)
            DoubleArray[l] /= WArray[i];
          l++;
        }
      }                                           // endfor i
    } else
    {
      for(int i=0;i<N_Vertices;i++)
      {
        if(WArray[i]!=0.)
        {
          DoubleArray[i] /= WArray[i];
        }
      }
    }

    if (!VectOutput)
    {
      // standard output writing
      dat << "SCALARS " << FEFunctionArray[k]->GetName();
      dat << " float"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(int j=0;j<N_Vertices;j++)
      {
        dat << DoubleArray[j] << endl;
      }
      dat << endl;
      dat << endl;
    }

    // write output as vector
    if (VectOutput)
    {
      // scalar components
      for(int j=0;j<N_Comp;j++)
      {
        dat << "SCALARS " << FEFunctionArray[k]->GetName() << j;
        dat << " float"<< endl;
        dat << "LOOKUP_TABLE " << "default" << endl;
        for(int i=0;i<N_Vertices;i++)
        {
          dat << DoubleArray[i*N_Comp+j] << endl;
        }
        dat << endl << endl;
      }

      // absolute value
      dat << "SCALARS " << "|" << FEFunctionArray[k]->GetName() << "|";
      dat << " float"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      l=0;
      for(int i=0;i<N_Vertices;i++)
      {
        t=0;
        for(int j=0;j<N_Comp;j++)
        {
          t+=DoubleArray[l]*DoubleArray[l];
          l++;
        }
        dat << sqrt(t)<< endl;
      }
      dat << endl << endl;

      dat << "VECTORS " << FEFunctionArray[k]->GetName();
      dat << " float"<< endl;

      l=0;
      for(int i=0;i<N_Vertices;i++)
      {
        for(int j=0;j<N_Comp;j++)
        {
          dat << DoubleArray[N_Comp*i+j] << " ";
        }
        dat << double(0) << " " << endl;
      }
      dat << endl;

      VectOutput = false;
    }

  }                                               // endfor k

  for(unsigned int k=0;k<FEVectFunctArray.size();k++) {
    fespace = FEVectFunctArray[k]->GetFESpace2D();
    N_Comp = FEVectFunctArray[k]->GetN_Components();
    Length = FEVectFunctArray[k]->GetLength();
    Coeffs = FEVectFunctArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    // cout << "N_Comp  " << N_Comp << endl;
    memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
    memset(WArray, 0, SizeOfDouble*N_Vertices);
    m = 0;

    for(int i=0;i<N_Elements;i++) {
      TBaseCell *cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(int j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices()) {
	case 3:
	  xi = TriaCoords[2*j];
	  eta = TriaCoords[2*j+1];
	  break;
	  
	case 4:
	  xi = QuadCoords[2*j];
	  eta = QuadCoords[2*j+1];
	  break;
        }
        bf->GetDerivatives(D00, xi, eta, BFValues);

        for(int n=0;n<N_Comp;n++) {
          value = 0;
          for(l=0;l<N_LocDOF;l++)
            value += BFValues[l] * Coeffs[DOF[l]+n*Length];
          DoubleArray[N_Comp*VertexNumbers[m] + n] += value;
        }
        WArray[VertexNumbers[m]] +=1.;
        m++;
      }                                           // endfor j
    }                                             // endfor i

    // mean value

    l = 0;
    for(int i=0;i<N_Vertices;i++) {
      for(int j=0;j<N_Comp;j++)
	{
	  if(WArray[i]!=0.)
	    DoubleArray[l] /= WArray[i];
	  l++;
	}
    }                                             // endfor l

    for(int j=0;j<N_Comp;j++)
    {
      dat << "SCALARS " << FEVectFunctArray[k]->GetName() << j;
      dat << " float"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(int i=0;i<N_Vertices;i++)
      {
        dat << DoubleArray[i*N_Comp+j] << endl;
      }
      dat << endl << endl;
    }

    // ***************
    // absolute value of vector variables
    // ***************
    dat << "SCALARS " << "|" << FEVectFunctArray[k]->GetName() << "|";
    dat << " float"<< endl;
    dat << "LOOKUP_TABLE " << "default" << endl;
    l=0;
    for(int i=0;i<N_Vertices;i++) {
      t=0;
      for(int j=0;j<N_Comp;j++) {
        t+=DoubleArray[l]*DoubleArray[l];
        l++;
      }
      dat << sqrt(t)<< endl;
    }
    dat << endl << endl;

    // ***************
    // VECTORS
    // ***************
    dat << "VECTORS " << FEVectFunctArray[k]->GetName();
    dat << " float"<< endl;

    l=0;
    for(int i=0;i<N_Vertices;i++) {
      for(int j=0;j<N_Comp;j++) {
        dat << DoubleArray[N_Comp*i+j] << " ";
      }
      dat << double(0) << " " << endl;
    }
    dat << endl;
  }                                               // endfor k

  dat << endl;

  delete [] NumberVertex;
  delete [] VertexNumbers;
  delete [] Vertices;
  delete [] DoubleArray;
  delete [] WArray;
  delete [] Coords;

  dat.close();

}

*/

void PostProcessing2D::sort(TVertex **Array, int length)
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

  delete [] rr;
}


int PostProcessing2D::getIndex(TVertex **Array, int Length, TVertex *Element)
{
  int l,r,m;
  l=0;
  r=Length;
  TVertex *Mid;

  m=Length/2;
  Mid=Array[m];
  while(Mid!=Element)
  {
    if(Mid>Element) {
      l=m;
    } else {
      r=m;
    }
    m=(r+l)/2;
    Mid=Array[m];
  }
  return m;
}
