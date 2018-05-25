#include <Database.h>
#include <FEDatabase2D.h>
#include <PostProcessing2D.h>

#include <sys/stat.h>
#include <algorithm>
#include <type_traits>

using namespace std;

PostProcessing2D::PostProcessing2D(const ParameterDatabase& param_db)
 : testcaseDir(), testcaseName(), n_steps_per_output(1), 
   n_steps_until_next_output(0), FEFunctionArray(), FEVectFunctArray(),
   Coll(nullptr), timeValues()
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
  n_steps_per_output = std::max<size_t>(1, db["steps_per_output"]); // avoid 0
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
  bool already_known = std::any_of(FEFunctionArray.begin(),
                                   FEFunctionArray.end(),
                                   [fefunction](const TFEFunction2D* f)
                                   {return f == fefunction;});
  if(!already_known )
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
  bool already_known = std::any_of(FEVectFunctArray.begin(),
                                   FEVectFunctArray.end(),
                                   [fevectfunction](const TFEVectFunct2D* f)
                                   {return f == fevectfunction;});
  if(!already_known )
    FEVectFunctArray.push_back(fevectfunction);
}

void PostProcessing2D::write(double current_time)
{
  if(n_steps_until_next_output != 0)
  {
    n_steps_until_next_output--;
  }
  else
  {
    n_steps_until_next_output = n_steps_per_output-1;
    timeValues.push_back(current_time);
    // index of current time step in 'timeValues'
    auto index = timeValues.size() - 1;
    if(writeVTK)
    {
      mkdir(testcaseDir.c_str(), 0777);
      std::string name;
      name += testcaseDir + "/" + testcaseName;
      name += std::to_string(index*n_steps_per_output);
      name += ".vtk";
      Output::print<2>(" PostProcessing2D:: writing ", name);
      writeVtk(name);
    }
    if(writeCASE)
    {
      if(index == 0) // first call to this method
      {
        // write geometry only in the first iteration
        writeCaseGeo();
      }
      writeCaseVars(index*n_steps_per_output);
      writeCaseFile();
    }
  }
}

void PostProcessing2D::write()
{
  if(writeVTK)
  {
    mkdir(testcaseDir.c_str(),0777);
    std::string name = testcaseDir + "/" + testcaseName + ".vtk";
    Output::print<2>(" PostProcessing2D:: writing ", name);
    writeVtk(name);
  }
  if(writeCASE)
  {
    writeCaseGeo();
    writeCaseVars(0);
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
  std::ofstream dat(name);
  if (!dat)
  {
    ErrThrow("cannot open file for output. ", name);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(12);
  
  // determine data for vtk file
  unsigned int N_Elements = Coll->GetN_Cells();

  // check for discontinuous scalar variables. In such a case write a new file
  // for this variable. ParaView will really display a discontinuous function
  // instead of projecting it onto P1/Q1 space. However there are some
  // drawbacks.
  std::vector <TVertex *> Vertices;
  for(unsigned int ncell = 0; ncell<N_Elements; ncell++)
  {
      TBaseCell* cell = Coll->GetCell(ncell);
      unsigned int N_Loc = cell->GetN_Vertices(); 
      for(unsigned int nvert =0; nvert < N_Loc; nvert++)
          Vertices.push_back(cell->GetVertex(nvert));
  }
  
  for(unsigned int i=0;i<FEFunctionArray.size();i++)
  {
    if ( FEFunctionArray[i]->GetFESpace2D()->IsDGSpace() ) 
   {
    writeVtkDiscontinuous(name,Vertices.size(),Vertices);
     break;
   }
  }

  dat << "# vtk DataFile Version 4.0" << endl;
  dat << "file created by ParMooN"
      << " Time < " << TDatabase::TimeDB->CURRENTTIME <<" >" << endl;

  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl;
  dat << "POINTS " << Coll->GetN_Vertices() << " double" << endl;

  writeCoord(dat, 2);

  unsigned int N_Vertices=Coll->GetN_Vertices();
  unsigned int N_LocVertices=Coll->GetNLocVertices();
  
  dat << "CELLS " << N_Elements << " " <<  N_Elements+ N_LocVertices<< endl;
  for(unsigned int i = 0; i < N_Elements; i++)
  {
    int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
    dat <<  N_CellVertices << " ";
    for(int j=0;j<N_CellVertices;j++)
    {
      dat << this->Coll->GetGlobalVerNo(i,j) << " "; 
    }
    dat << endl;
  }
  dat << endl;

  dat << "CELL_TYPES " << N_Elements << endl;
  for(unsigned int i=0;i<N_Elements;i++)
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

  unsigned int dimension;
  std::vector<double> solutionAtNode;
  
  for(unsigned int nfunction =0; nfunction< FEFunctionArray.size(); nfunction++)
  {   
    computeNodeValues<TFEFunction2D>(FEFunctionArray.at(nfunction), solutionAtNode, dimension);
    
    std::string name =FEFunctionArray.at(nfunction)->GetName();
   
    printVectCompwise(dat, name, N_Vertices, dimension, solutionAtNode); 
      
    if (dimension > 1)
    {     
      printVectAbsValue(dat, name, N_Vertices, dimension, solutionAtNode);
      printVectPointwise(dat, name, N_Vertices, dimension, solutionAtNode);      
    }
  }
  
  for(unsigned int nfunction =0; nfunction< FEVectFunctArray.size(); nfunction++)
  {    
    computeNodeValues<TFEVectFunct2D>(FEVectFunctArray.at(nfunction), solutionAtNode, dimension);
    
    std::string name =FEVectFunctArray.at(nfunction)->GetName();
    
    printVectCompwise(dat, name, N_Vertices, dimension, solutionAtNode);
    printVectAbsValue(dat, name, N_Vertices, dimension, solutionAtNode);
    printVectPointwise(dat, name, N_Vertices, dimension, solutionAtNode);
  } 
  
  dat << endl;
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
int N_LocVertices, std::vector<TVertex *> Vertices)
{
  double x,y;                // coordinates of a vertex
  int N_Elements;            // number of all elements in this mesh
  int N_CellVertices;        // number of all vertices in this cell
  int i,j,k,l;               // loop variables
  const TFEFunction2D *fefunction; // this function, which is discontinuous
  TBaseCell *current_cell;
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
  dat << "POINTS " << N_LocVertices << " double" << endl;

  for(i=0;i<N_LocVertices;i++)
  {
#ifdef __3D__
    double z;
    Vertices[i]->GetCoords(x,y,z);
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
    if(!fefunction->GetFESpace2D()->IsDGSpace())
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
      dat << " double" << endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      double function_value;    // value of function at a particular vertex
      for(i=0;i<N_Elements;i++)
      {
        current_cell = Coll->GetCell(i);
        N_CellVertices=current_cell->GetN_Vertices();
        for(j=0;j<N_CellVertices;j++)
        {
#ifdef __3D__
          double z;
          current_cell->GetVertex(j)->GetCoords(x, y, z);
#else
          current_cell->GetVertex(j)->GetCoords(x, y);
#endif
          fefunction->FindValueLocal(current_cell,i,x,y, &function_value);
          dat << function_value << endl;
        }
      }
      
    }
    else if(BaseVectDim==2)
    {
      // find values for all components
      double function_value[BaseVectDim]; // function values at a vertex
      // store all values of all components
      std::vector<std::vector<double>> 
          allValues(BaseVectDim, std::vector<double>(N_LocVertices));
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
        dat << " double" << endl;
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
      Output::print("TOutput2D::writeVtkDiscontinuous: Basis functions of "
                    "dimension ", BaseVectDim, " are not supported.");
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
    unsigned int dimension = TFEDatabase2D::GetFE2D(
      FEFunctionArray[j]->GetFESpace2D()->GetFE2D(0,Coll->GetCell(0)))
      ->GetBaseFunct2D()->GetBaseVectDim();
    
    if (dimension==1)
    {
      casf << "scalar per node: 1 " << FEFunctionArray[j]->GetName() << " "
           << testcaseName << "_" << FEFunctionArray[j]->GetName()
	   << ".****.scl\n";
    }
    else
    {
      casf << "vector per node: 1 "<< FEFunctionArray[j]->GetName() << " "
	   << testcaseName << "_" << FEFunctionArray[j]->GetName() 
	   << ".****.vct\n";
    }
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
  casf << "filename increment: " << n_steps_per_output << "\n";  
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
  
  unsigned int N_Vertices = Coll->GetN_Vertices();
  
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
  
  writeCoord(f, dimension);

  // write elements
  int nVE = Coll->GetCell(0)->GetN_Vertices();
  string ensight_type;

  /// @warning this works now only for a single (bulk) domain
  unsigned int nParts = 1; // number of (bulk) subdomains
  string partName = "inner"; // name of collection subdomain
  
  unsigned int VERTEX_OFFSET =1; //index change for the output
  
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
        f << Coll->GetGlobalVerNo(k,i)+VERTEX_OFFSET;
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
    f << Coll->GetN_BdFaces() << endl;
    for(unsigned int k=0;k<Coll->GetN_BdFaces();k++)
    {
      for(int i=0;i<nVertexPerFace;i++)
      {
        f.width(8);
        f << Coll->GetBdFacesNode(nVertexPerFace*k+i);
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
  
  unsigned int N_Vertices = Coll->GetN_Vertices();
  unsigned int dimension;

  // scalars
  for (unsigned int i=0;i<FEFunctionArray.size(); i++)
  {
    std::vector<double> solutionAtNode;
    computeNodeValues(FEFunctionArray[i], solutionAtNode, dimension);

    // write scalars
    std::ofstream sclf;
    std::ostringstream sclname;
    
    if(dimension==1)
    {    
      sclname << filename << "_" << FEFunctionArray[i]->GetName() <<  "."
	      << number << ".scl";    
    }
    else
    {
      sclname << filename << "_" << FEFunctionArray[i]->GetName() <<  "."
	      << number << ".vct"; 
    }
	      
    std::string fname = sclname.str();

    Output::print<3>(" ** PostProcessing2D::write File - write ", fname);
    sclf.open(fname);
    sclf << FEFunctionArray[i]->GetName() << " step = " << iter << endl; 
  
    writeVectCase(sclf, N_Vertices, dimension, solutionAtNode);
  }

  // vectors
  for(unsigned int i=0;i<FEVectFunctArray.size();i++)
  {
    //int nComp = FEVectFunctArray[i]->GetN_Components();
    std::vector<double > solutionAtNode;
    computeNodeValues(FEVectFunctArray[i], solutionAtNode, dimension); 
    
    std::ofstream vctf;
    std::ostringstream vctname;
    vctname << filename << "_" << FEVectFunctArray[i]->GetName()
            <<  "." << number << ".vct";
    std::string fname = vctname.str();
    Output::print<3>(" ** PostProcessing2D::write File - write ", fname);
    vctf.open(fname.c_str());
    vctf <<  FEVectFunctArray[i]->GetName() << " step = " << iter << endl; 
  
    writeVectCase(vctf, N_Vertices, dimension, solutionAtNode);
  }
}

//Check, if PostProcessing object should be generalized
void PostProcessing2D::writeCoord(ofstream & f, int dimension)
{
  unsigned int N_Vertices = Coll->GetN_Vertices();
  for(unsigned int i=0;i<N_Vertices;i++)
  {
    //f.setf(ios_base::scientific);
    f.precision(12);
    f.width(12);
    f <<  Coll->GetCoord(i*dimension);
    f << " ";
    //f.setf(ios_base::scientific);
    f.precision(12);
    f.width(12);
    f << Coll->GetCoord(i*dimension+1);
    f << " ";
    //f.setf(ios_base::scientific);
    f.precision(12);
    f.width(12);
    if (dimension==2)
    {
      f << 0.;
    }
    else
    {
      f  << Coll->GetCoord(i*dimension+2); 
    }
    f << endl;
  }
  f << endl;
}

template <class T>
void PostProcessing2D::computeNodeValues(const T* function, 
					 std::vector< double >& solutionAtNode,
					 unsigned int & dimension)
{
  auto FESpace2D = function->GetFESpace2D();
  TCollection* coll = FESpace2D->GetCollection();
  int nPoints = coll->GetN_Vertices();

  // compute FE type on first cell
  TBaseCell *cell = coll->GetCell(0);
  FE2D FE_ID = FESpace2D->GetFE2D(0, cell);
  TBaseFunct2D *bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
  
  // get function type of T
  bool IsVect = std::is_same< T, TFEVectFunct2D >::value;  
  dimension=1;
  if (IsVect || bf->GetBaseVectDim()>1)
    dimension = 2;

  solutionAtNode.assign(dimension*nPoints, 0.0);
  std::vector<double> WArray(nPoints, 0.);
  
  for(int ncell=0;ncell<coll->GetN_Cells();ncell++) 
  {
    cell = coll->GetCell(ncell);
    unsigned int nLocalVertices = cell->GetN_Vertices();

    bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();

    for(unsigned int nvert=0; nvert<nLocalVertices; nvert++) 
    {   
      unsigned int globalVert_index = coll->GetGlobalVerNo(ncell,nvert);
      
      double x,y;
#ifdef __3D__
          ;
#else
      cell->GetVertex(nvert)->GetCoords(x,y);
#endif
      double function_value[dimension];
    
      if (dimension>1) //computations for "vector valued scalar functions"
      {
	function->FindValueLocal(cell,ncell,x,y,function_value);	
	    
	solutionAtNode[dimension*globalVert_index + 0 ] += function_value[0];
	solutionAtNode[dimension*globalVert_index + 1 ] += function_value[1];
      }  
      else	// compute nodal value of FE function
      {
	function->FindValueLocal(cell, ncell, x,y, function_value); 
	solutionAtNode.at(dimension *globalVert_index)+=function_value[0];
      }
      WArray.at(globalVert_index)+=1.;

    } // endfor nvert
  }  // endfor ncell

  for(int nvert=0;nvert<nPoints;nvert++) 
  {
    if(WArray[nvert]!=0.) 
    {
      solutionAtNode[dimension*nvert] /= WArray[nvert];
      if (dimension>1) //dimension==2
      {
	solutionAtNode[dimension*nvert+1]/= WArray[nvert];
      }
    }
  }
}

// ****************
// scalar components of vector variables
// ****************
void PostProcessing2D::printVectCompwise(ofstream& dat,  std::string name,
					 unsigned int N_Vertices, unsigned int N_Comp,
					 vector< double > solutionAtNode)
{
  for(unsigned int ncomp=0;ncomp<N_Comp;ncomp++)
    {
      if(N_Comp==1)
	dat << "SCALARS " << name;
      else
	dat << "SCALARS " << name << ncomp;
      dat << " double"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(unsigned int nvert=0;nvert<N_Vertices;nvert++)
      {
        double value =solutionAtNode[nvert*N_Comp+ncomp];
	if(fabs(value)<1e-18)
	  value = 0.;
	dat <<  value << endl;
      }
      dat << endl << endl;
    }
}

// ***************
// absolute value of vector variables
// ***************
void PostProcessing2D::printVectAbsValue(ofstream& dat, std::string name,
		        unsigned int N_Vertices, unsigned int N_Comp,
					 vector< double > solutionAtNode)
{
  dat << "SCALARS " << "|" << name << "|";
  dat << " double"<< endl;
  dat << "LOOKUP_TABLE " << "default" << endl;
  for(unsigned int nvert = 0, l = 0; nvert < N_Vertices; nvert++)
  {
    double t=0;
    for(unsigned int ncomp=0;ncomp<N_Comp;ncomp++)
    {
      t+=solutionAtNode[l]*solutionAtNode[l];
      l++;
    }
    dat << sqrt(t)<< endl;
  }
  dat << endl << endl;
}

// ***************
// VECTORS
// ***************
void PostProcessing2D::printVectPointwise(ofstream& dat, std::string name,
		        unsigned int N_Vertices, unsigned int N_Comp,
					  vector< double > solutionAtNode)
{
  dat << "VECTORS " << name;
  dat << " double"<< endl;
  for(unsigned int nvert = 0; nvert < N_Vertices; nvert++)
  {
    for(unsigned int ncomp=0;ncomp<N_Comp;ncomp++)
    {
        double value =solutionAtNode[nvert*N_Comp+ncomp];
	if(fabs(value)<1e-18)
	  value = 0.;
	dat <<  value << " ";
    }
    dat << double(0) << " " << endl;
  }
  dat << endl;
}

//auxiliary function for writeVectCase
void PostProcessing2D::printEntry(ofstream & dat, double value, int counter)
{
    dat.setf(ios_base::scientific);	
    dat.precision(5);	
    dat.width(12);
    if(fabs(value)<1e-18)
      value = 0.;
    dat << value;

  // One line contains most six entries.    
    if( (counter+1) % 6 == 0 )
      dat << std::endl;
}

void PostProcessing2D::writeVectCase(ofstream& dat, 
				     unsigned int N_Vertices, unsigned int N_Comp, 
				     vector< double > solutionAtNode)
{
  for (unsigned int nvert =0; nvert < N_Vertices; nvert++)
  {
    unsigned int entry = N_Comp*nvert;
    
    double value_x = solutionAtNode.at(entry);
    
    if (N_Comp==1)
      printEntry(dat, value_x, entry);
    else
    {
      int counter = 3*nvert;
      
      printEntry(dat,value_x,counter);
      
      double value_y = solutionAtNode.at(entry+1);
      printEntry(dat,value_y,counter+1);
    
      double value_z = 0.;
      if (N_Comp>2)
	value_z = solutionAtNode.at(entry+2);  
      printEntry(dat,value_z,counter+2);    
    } 
  }
}