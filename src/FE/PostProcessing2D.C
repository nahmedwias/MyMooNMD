#include <Database.h>
#include <FEDatabase2D.h>
#include <PostProcessing2D.h>
#include <QuadAffin.h>
#include <TriaAffin.h>
#include <QuadBilinear.h>

using namespace std;

PostProcessing2D::PostProcessing2D(const ParameterDatabase& param_db)
 : testcaseDir(), testcaseName(), period(1), FEFunctionArray(),
   FEVectFunctArray(), Coll(nullptr), timeValues()
{
  // set the variables depending on input parameters
  ParameterDatabase db = ParameterDatabase::default_output_database();
  db.merge(param_db,false);
  
  writeVTK = db["output_write_vtk"];
  writeCASE = db["output_write_case"];
  if(writeCASE)
  {
    Output::warn("PostProcessing2D", "case output is not working with older ",
                 "versions of ParaView on Linux. Try (at least) version 5.0");
  }

  testcaseName = db["output_basename"].get<std::string>();
  testcaseDir = db["output_vtk_directory"].get<std::string>();
  //period = db["steps_per_output"];
};


void PostProcessing2D::add_fe_function(const TFEFunction2D* fefunction)
{
  // check that FE functions have the same collection
  if(this->Coll == nullptr)
  {
    this->Coll = fefunction->GetFESpace2D()->GetCollection();
  }
  else
  {
    if(this->Coll != fefunction->GetFESpace2D()->GetCollection())
    {
      // we could also just refuse to add this fe function and return without an
      // exception.
      ErrThrow("new FE function has a different collection");
    }
  }
  FEFunctionArray.push_back(fefunction);
}

void PostProcessing2D::add_fe_vector_function(
  const TFEVectFunct2D* fevectfunction)
{
  // check that FE functions have the same collection
  if(this->Coll == nullptr)
  {
    this->Coll = fevectfunction->GetFESpace2D()->GetCollection();
  }
  else
  {
    if(this->Coll != fevectfunction->GetFESpace2D()->GetCollection())
    {
      // we could also just refuse to add this fe function and return without an
      // exception.
      ErrThrow("new FE vector function has a different collection");
    }
  }
  FEVectFunctArray.push_back(fevectfunction);
}

void PostProcessing2D::write(double current_time)
{
  timeValues.push_back(current_time);
  if(writeVTK)
  {
    std::string name;
    name += testcaseDir + "/" + testcaseName;
    name += std::to_string(timeValues.size());
    name += ".vtk";
    Output::print<2>(" PostProcessing2D:: writing ", name);
    writeVtk(name);
  }
  if(writeCASE)
  {
    if(timeValues.size() == 1) // first call to this method
    {
      // write geometry only in the first iteration
      writeCaseGeo();
    }
    writeCaseVars(timeValues.size()-1);
    writeCaseFile();
  }
}

void PostProcessing2D::write(int i)
{
  if(!timeValues.empty())
    Output::warn<1>("PostProcessing2D::write(int i)", "it seems you have "
                    "called PostProcessing2D::write(double) before, which is "
                    "used for time dependent problems.\nThis method is "
                    "typically for stationary problems. I am not sure if this "
                    "works");
  if(writeVTK)
  {
    std::string name;
    name += testcaseDir + "/" + testcaseName;
    if(i>=0)
    {
      name += std::to_string(i);
    }
    name += ".vtk";
    Output::print<2>(" PostProcessing2D:: writing ", name);
    writeVtk(name);
  }

  if(writeCASE)
  {
    // note: i<0 is used to avoid suffix in vtk output in steady problems
    // it shall be disregarded for case output
    if(i < 0) 
      i = 0;
    if(i == 0)
    {
      // write geometry only in the first iteration
      writeCaseGeo();
    }
    writeCaseVars(i);
    writeCaseFile();
  }
}


/**
   @brief write stored data into a VTK file (old  version)
   This function writes the mesh and the solution into a vtk file
   Note: it uses an old-fashion way to compute the P1 solution
   The new implementation should be used eventually.
*/
void PostProcessing2D::writeVtk(std::string name)
{
  double BFValues[MaxN_BaseFunctions2D];
  double BFValuesOrig[MaxN_BaseFunctions2D];
  double QuadCoords[] = { -1, -1, 1, -1, 1, 1, -1, 1};
  double TriaCoords[] = { 0, 0, 1, 0,  0, 1};

  std::ofstream dat(name);
  if (!dat)
  {
    ErrThrow("cannot open file for output. ", name);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(12);


  // determine data for vtk file
  int N_Elements = Coll->GetN_Cells();
  // number of vertices, counting each vertex n times where n is the number of 
  // cells this vertex belongs to.
  int N_LocVertices = 0;
  for(int i=0;i<N_Elements;i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
  }
  TVertex **Vertices =new TVertex*[N_LocVertices];
  int N_=0;
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
  int N_Vertices=0;
  for(int i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
  {
    N_Vertices++;
    Last=Current;
  }

  std::vector<double> Coords(2*N_Vertices);
  std::vector<int> VertexNumbers(N_LocVertices);
  std::vector<int> NumberVertex(N_LocVertices);
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

  for(int i = 0, m = 0; i < N_Elements; i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    for(int j=0;j<cell->GetN_Vertices(); j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      int l=getIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      m++;
    }                                             // endfor j
  }                                               //endfor i



  dat << "# vtk DataFile Version 4.0" << endl;
  dat << "file created by ParMooN"
      << " Time < " << TDatabase::TimeDB->CURRENTTIME <<" >" << endl;

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
  for(int i = 0, l = 0; i < N_Elements; i++)
  {
    int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
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
    int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
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
  std::vector<double> DoubleArray(2*N_Vertices);
  std::vector<double> WArray(N_Vertices);

  // write scalar variables into file
  for(unsigned int k=0;k<FEFunctionArray.size();k++)
  {
    const TFESpace2D *fespace = FEFunctionArray[k]->GetFESpace2D();
    const double *Coeffs = FEFunctionArray[k]->GetValues();
    const int *GlobalNumbers = fespace->GetGlobalNumbers();
    const int *BeginIndex = fespace->GetBeginIndex();

    // get dimension of basis functions
    TBaseCell *cell = Coll->GetCell(0);
    FE2D FE_ID = fespace->GetFE2D(0, cell);
    TBaseFunct2D *bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
    //this is needed to handle vector fields approximated with vector FE 
    // (scalar unknown postprocessed as vectors)
    int BaseVectDim = bf->GetBaseVectDim();
    int N_Comp = BaseVectDim;

    std::fill(DoubleArray.begin(), DoubleArray.end(), 0.0);
    std::fill(WArray.begin(), WArray.end(), 0.0);

    // set to TRUE if basis functions are vectors
    bool VectOutput = false;
    if (BaseVectDim>1) VectOutput = true;
    
    for(int i = 0, m = 0; i < N_Elements; i++)
    {
      TBaseCell *cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE2D FE_ID = fespace->GetFE2D(i, cell);
      TBaseFunct2D *bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      const int *DOF = GlobalNumbers+BeginIndex[i];
      int N_LocDOF = bf->GetDimension();
      RefTrans2D RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FE_ID);
      TRefTrans2D *F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
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
          ErrThrow("no such reference transformation allowed ", RefTrans);
          break;
      }                                           // endswitch

      for(int j=0;j<N_;j++)
      {
        double xi, eta;
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

        double value = 0;
        double value_y = 0;
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
          for(int l=0;l<N_LocDOF;l++)
          {
            // change sign of basis functions according to global normal
            int edge = 
             TFEDatabase2D::GetFE2D(FE_ID)->GetFEDesc2D()->GetJointOfThisDOF(l);
            int nsign = 1;
            if(edge != -1)
            {
              nsign = cell->GetNormalOrientation(edge);
            }
            value += BFValuesOrig[l] * Coeffs[DOF[l]]*nsign;
            value_y += BFValuesOrig[N_LocDOF+l] * Coeffs[DOF[l]]*nsign;
          }
          DoubleArray[N_Comp*VertexNumbers[m] + 0] += value;
          DoubleArray[N_Comp*VertexNumbers[m] + 1] += value_y;
        }
        else
        {                                         // standard
          for(int l=0;l<N_LocDOF;l++)
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
      for(int i = 0, l = 0; i < N_Vertices; i++)
      {
        for(int j=0;j<N_Comp;j++)
        {
          if(WArray[i]!=0.)
            DoubleArray[l] /= WArray[i];
          l++;
        }
      }                                           // endfor i
    }
    else
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
      for(int i = 0, l = 0; i < N_Vertices; i++)
      {
        double t=0;
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

      for(int i = 0; i < N_Vertices; i++)
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

  for(unsigned int k=0;k<FEVectFunctArray.size();k++)
  {
    const TFESpace2D *fespace = FEVectFunctArray[k]->GetFESpace2D();
    int N_Comp = FEVectFunctArray[k]->GetN_Components();
    int Length = FEVectFunctArray[k]->GetLength();
    const double * Coeffs = FEVectFunctArray[k]->GetValues();
    const int * GlobalNumbers = fespace->GetGlobalNumbers();
    const int * BeginIndex = fespace->GetBeginIndex();

    // cout << "N_Comp  " << N_Comp << endl;
    std::fill(DoubleArray.begin(), DoubleArray.end(), 0.0);
    std::fill(WArray.begin(), WArray.end(), 0.0);
    
    for(int i = 0, m = 0; i < N_Elements; i++)
    {
      TBaseCell *cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE2D FE_ID = fespace->GetFE2D(i, cell);
      TBaseFunct2D *bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      const int * DOF = GlobalNumbers+BeginIndex[i];
      int N_LocDOF = bf->GetDimension();
      for(int j=0;j<N_;j++)
      {
        double xi, eta;
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

        for(int n=0;n<N_Comp;n++)
        {
          double value = 0;
          for(int l = 0; l < N_LocDOF; l++)
            value += BFValues[l] * Coeffs[DOF[l]+n*Length];
          DoubleArray[N_Comp*VertexNumbers[m] + n] += value;
        }
        WArray[VertexNumbers[m]] +=1.;
        m++;
      }                                           // endfor j
    }                                             // endfor i

    // mean value
    for(int i = 0, l = 0; i < N_Vertices; i++)
    {
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
    for(int i = 0, l = 0; i < N_Vertices; i++)
    {
      double t=0;
      for(int j=0;j<N_Comp;j++)
      {
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
    for(int i = 0; i < N_Vertices; i++)
    {
      for(int j=0;j<N_Comp;j++)
      {
        dat << DoubleArray[N_Comp*i+j] << " ";
      }
      dat << double(0) << " " << endl;
    }
    dat << endl;
  }                                               // endfor k
  dat << endl;
  
  delete [] Vertices;
  dat.close();
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
    
    as is done in writeVtk(std::string name).
    Also note that in each element the function is projected onto P1/Q1.
    @warning This destroys the topology of the mesh. 
    @warning Some filters in ParaView might work incorrectly. 
    Warp by scalar works though. 
    The way this is done here is not very elegant, but I couldn't find a better
    solution.
    
*/
void PostProcessing2D::writeVtkDiscontinuous(std::string fileName,
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
  // copy the file name add a string to the output file name
  std::string disc = fileName + std::string("_disc.vtk");

  N_Elements = Coll->GetN_Cells();

  std::ofstream dat(disc);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    exit(-1);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(9);

  dat << "# vtk DataFile Version 4.0" << endl;
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
  string filename = testcaseDir + "/" + testcaseName;
  string casename = filename + ".case";
  Output::print<2>(" ** PostProcessing2D::writeCaseFile - write ", casename);
  casf.open(casename);
  
  string geoname = testcaseName + ".geo";
  
  casf << "FORMAT\n";
  casf << "type: ensight\n";
  casf << "GEOMETRY\n";
  casf << "model: 1 " << geoname << endl;
  casf << "VARIABLE\n";
  for (unsigned int j=0; j<FEFunctionArray.size(); j++)
  {
    string sclname =  "scalar";
    casf << "scalar per node: 1 " << FEFunctionArray[j]->GetName() << " "
         << testcaseName << "_" << FEFunctionArray[j]->GetName() 
         << ".****.scl\n";
  }
  for (unsigned int j=0; j<FEVectFunctArray.size(); j++)
  {
    casf << "vector per node: 1 "<< FEVectFunctArray[j]->GetName() << " "
         << testcaseName << "_" << FEVectFunctArray[j]->GetName() 
         << ".****.vct\n";
  }
  // time values
  unsigned int n_time_steps = timeValues.size();
  casf << "TIME\n";
  casf << "time set: 1\n";
  casf << "number of steps: " << std::max(1U, n_time_steps) << endl;
  casf << "filename start number: 0\n";
  casf << "filename increment: " << period << "\n";  
  casf << "time values:\n";
  for(unsigned int i = 0; i < n_time_steps; i++)
  {
    casf.precision(5);
    casf.width(12);
    casf << timeValues[i];
    if( (i && !(i%5)) || (i==n_time_steps) )
       casf << endl;
    else casf << " ";
  }
  if(n_time_steps == 0) // this is the case for stationary problems
    casf << 0.0;
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
  f.open(geoname);

  f << testcaseName << " geometry\n";
  f << "t = 0.0000" << "\n";
  f << "node id assign" << "\n";
  f << "element id assign" << "\n";
  f << "coordinates" << "\n";
  f.width(8);
  // write points
  f << N_Vertices << endl;
  for(unsigned int i=0;i<N_Vertices;i++)
  {
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
    if (dimension==2)
    {
      f << 0.;
    }
    else
    {
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
  
  for(unsigned int ifig=1;ifig<=nParts;ifig++)
  { 
    f << "part";
    f.width(8);
    f << ifig << endl;
    f << partName << endl; // name of the figure
    switch(dimension)
    {
      case 2:
        switch(nVE)
        {
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
    for(int k=0;k<Coll->GetN_Cells();k++)
    {
      // @warning we cannot write yet .case output with mixed meshes.
      if(Coll->GetCell(k)->GetN_Vertices() != nVE)
      {
        ErrThrow(" **ERROR: CASE output with mixed (tria and quad) meshes not "
                 "supported yet. Use VTK instead.");
      }
      
      for(int i=0;i<nVE;i++)
      {
        f.width(8);
        f << Coll->ElementNodes[k][i];
      }
      f << endl;
    }// k
  }//ifig
  
  // write boundary elements
  unsigned int nBdParts = 1; // number of subdomains
  ///@warning this works now only for a single boundary domain
  partName = "boundary"; // name of collection subdomain
  int nVertexPerFace = 2;
  for(unsigned int ifig=1;ifig<=nBdParts;ifig++)
  { 
    f << "part";
    f.width(8);
    f << nParts+ifig << endl;
    f << partName << endl; // name of the figure
    ///@warning only for P1 in 1D
    ensight_type="bar2";
    f << ensight_type << endl;
    
    f.width(8);
    f << Coll->BdFacesReferences.size() << endl;
    for(unsigned int k=0;k<Coll->BdFacesReferences.size();k++)
    {
      for(int i=0;i<nVertexPerFace;i++)
      {
        f.width(8);
        f << Coll->BdFacesNodes[nVertexPerFace*k+i];
      }
      f << endl;
    }// k
  }//ifig

}


void PostProcessing2D::writeCaseVars(int iter)
{
  std::string filename = testcaseDir+"/"+testcaseName;
  char numstr[4];
  sprintf(numstr,"%d",iter);
  std::string number = "000"+string(numstr);
  if (iter>9)
    number = "00" + string(numstr);
  if (iter>99)
    number = "0" + string(numstr);
  if (iter>999)
    number = string(numstr);

  std::string ensight_type;
  int nParts = 1;
  int n;
  ///@todo do not recreate lists if they exist already
  this->Coll->createElementLists();

  // scalars
  for (unsigned int i=0;i<FEFunctionArray.size(); i++)
  {
   std::vector<double> uP1;
   FEFunctionArray[i]->computeNodeValues(uP1);

    // write scalars
    std::ofstream sclf;
    std::ostringstream sclname;
    sclname << filename << "_" << FEFunctionArray[i]->GetName() <<  "."
            << number << ".scl";
    std::string fname = sclname.str();
    cout << " ** PostProcessing2D::write File - write " << sclname.str() << endl;
    sclf.open(fname);
    sclf << FEFunctionArray[i]->GetName() << " step = " << iter << endl; 
  
    for(int ifig=1;ifig<=nParts;ifig++)
    {
      n = 0;
      for(unsigned int k=0;k< Coll->NodesReferences.size(); k++)
      {
        double x = uP1[k];
        sclf.setf(ios_base::scientific);	
        sclf.precision(5);	
        sclf.width(12);
        if(fabs(x)<1e-18)
          x = 0.;
        sclf << x;
        ++n;
        if( n == 6 )
        {
          sclf << std::endl;
          n=0;
        }
      }// k
     
     sclf << endl;
    }//ifig
  }

  // vectors
  for(unsigned int i=0;i<FEVectFunctArray.size();i++)
  {
    int nComp = FEVectFunctArray[i]->GetN_Components();
    std::vector< std::vector<double> > uP1vect(nComp);
    for(int j=0;j<nComp;j++)
    {
      FEVectFunctArray[i]->GetComponent(j)->computeNodeValues(uP1vect[j]);
    }
    
    std::ofstream vctf;
    std::ostringstream vctname;
    vctname << filename << "_" << FEVectFunctArray[i]->GetName()
            <<  "." << number << ".vct";
    std::string fname = vctname.str();
    Output::print<2>(" ** PostProcessing2D::write File - write ",
                     vctname.str());
    vctf.open(fname.c_str());
    vctf <<  FEVectFunctArray[i]->GetName() << " step = " << iter << endl; 
  
    for(int ifig=1;ifig<=nParts;ifig++)
    {
      n = 0;
      for(int k=0;k<Coll->NodesReferences.size();k++)
      {
        double x,y,z;
        x = uP1vect[0][k];
        vctf.setf(ios_base::scientific);	
        vctf.precision(5);	
        vctf.width(12);
        if(fabs(x)<1e-18)
          x = 0.;
        vctf << x;
        ++n;
        if( n == 6 )
        {
          vctf << std::endl;
          n=0;
        }
        
        y = 0.;
        if (nComp>1)
        {
          y = uP1vect[1][k];
        }
        vctf.setf(ios_base::scientific);	
        vctf.precision(5);	
        vctf.width(12);
        if(fabs(y)<1e-18)
          y = 0.;
        vctf << y;
        ++n;
        if( n == 6 )
        {
          vctf << std::endl;
          n=0;
        }
        
        z = 0.;
        if(nComp>2)
        {
          z = uP1vect[2][k];
        }
        vctf.setf(ios_base::scientific);	
        vctf.precision(5);	
        vctf.width(12);
        if(fabs(z) < 1e-18)
          z = 0.;
        vctf << z;
        ++n;
        if( n == 6 )
        {
          vctf << std::endl;
          n=0;
        }
      }// k
      vctf << endl;
    }//ifig
  }
}


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
