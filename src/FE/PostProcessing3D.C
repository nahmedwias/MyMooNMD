#include <Database.h>
#include <FEDatabase3D.h>
#include <PostProcessing3D.h>

#include <sys/stat.h>
#include <type_traits>
#include <algorithm>
#include <iostream>
#include <fstream>

#ifdef _MPI
#include "mpi.h"
#endif

using namespace std;

PostProcessing3D::PostProcessing3D(const ParameterDatabase& param_db)
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

void PostProcessing3D::add_fe_function(const TFEFunction3D* fefunction)
{
  // check that FE functions have the same collection
  if(this->Coll == nullptr)
  {
    this->Coll = fefunction->GetFESpace3D()->GetCollection();
  }
  else
  {
    if(this->Coll != fefunction->GetFESpace3D()->GetCollection())
    {
      // we could also just refuse to add this fe function and return without an
      // exception.
      ErrThrow("new FE function has a different collection");
    }
  }
  bool already_known = std::any_of(FEFunctionArray.begin(),
                                   FEFunctionArray.end(),
                                   [fefunction](const TFEFunction3D* f)
                                   {return f == fefunction;});
  if(!already_known )
    FEFunctionArray.push_back(fefunction);
}

void PostProcessing3D::add_fe_vector_function(const TFEVectFunct3D* fevectfunction)
{
  // check that FE functions have the same collection
  if(this->Coll == nullptr)
  {
    this->Coll = fevectfunction->GetFESpace3D()->GetCollection();
  }
  else
  {
    if(this->Coll != fevectfunction->GetFESpace3D()->GetCollection())
    {
      // we could also just refuse to add this fe function and return without an
      // exception.
      ErrThrow("new FE vector function has a different collection");
    }
  }
  bool already_known = std::any_of(FEVectFunctArray.begin(),
                                   FEVectFunctArray.end(),
                                   [fevectfunction](const TFEVectFunct3D* f)
                                   {return f == fevectfunction;});
  if(!already_known )
    FEVectFunctArray.push_back(fevectfunction);
}

void PostProcessing3D::write()
{
  if(writeVTK)
  {
#ifdef _MPI
     char SubID[] = "";
     int my_rank;
     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
     if(my_rank==0)
       mkdir(testcaseDir.c_str(), 0777);
     Write_ParVTK(MPI_COMM_WORLD, 0, SubID, testcaseDir, testcaseName);
#else    
    mkdir(testcaseDir.c_str(), 0777);
    std::string name = testcaseDir + "/" + testcaseName + ".vtk";
    Output::print<2>(" PostProcessing3D:: writing ", name);
    writeVtk(name);
#endif
    
  }
  if(writeCASE)
  {
#ifdef _MPI
    ErrThrow(" A parallel option for CASE-output is not implemented, yet.");
#endif
    writeCaseGeo();
    writeCaseVars(0);
    writeCaseFile();
  }
}

void PostProcessing3D::write(double current_time)
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
#ifdef _MPI
      char SubID[] = "";
      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      if(my_rank==0)
        mkdir(testcaseDir.c_str(), 0777);
      Write_ParVTK(MPI_COMM_WORLD, index*n_steps_per_output, SubID, testcaseDir, testcaseName);
#else      
      mkdir(testcaseDir.c_str(), 0777);
      std::string name;
      name += testcaseDir + "/" + testcaseName;
      name += std::to_string(index*n_steps_per_output);
      name += ".vtk";
      Output::print<2>(" PostProcessing3D:: writing ", name);
      writeVtk(name);
#endif
    }
    if(writeCASE)
    {
#ifdef _MPI
      ErrThrow(" A parallel option for CASE-output is not implemented, yet.");
#endif      
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

void PostProcessing3D::writeMesh(std::string name)
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
  
  dat << "# vtk DataFile Version 4.0" << endl;
  dat << "file created by ParMooN"
      << " Time < " << TDatabase::TimeDB->CURRENTTIME <<" >" << endl;

  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl;
  dat << "POINTS " << Coll->GetN_Vertices() << " double" << endl;

  writeCoord(dat, 3);

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
      case 4: dat << 10 << " ";
      break;
      case 8: dat << 12 << " ";
      break;
    }
  }
}

void PostProcessing3D::writeVtk(string name)
{  
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
    if ( FEFunctionArray[i]->GetFESpace3D()->IsDGSpace() ) 
   {
    writeVtkDiscontinuous(name,Vertices.size(),Vertices);
     break;
   }
  }
  
  writeMesh(name);
  ofstream dat(name, ios::out|ios::app);
  if (!dat)
  {
    ErrThrow("cannot open file for output. ", name);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(12);
  
  unsigned int N_Vertices=Coll->GetN_Vertices();

  dat << endl << endl;
  dat << "POINT_DATA " << N_Vertices << endl;

  unsigned int dimension;
  std::vector<double> solutionAtNode;
  
  for(unsigned int nfunction =0; nfunction< FEFunctionArray.size(); nfunction++)
  {   
    computeNodeValues<TFEFunction3D>(FEFunctionArray.at(nfunction), solutionAtNode, dimension);
    
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
    computeNodeValues<TFEVectFunct3D>(FEVectFunctArray.at(nfunction), solutionAtNode, dimension);
    
    std::string name =FEVectFunctArray.at(nfunction)->GetName();
    
    printVectCompwise(dat, name, N_Vertices, dimension, solutionAtNode);
    printVectAbsValue(dat, name, N_Vertices, dimension, solutionAtNode);
    printVectPointwise(dat, name, N_Vertices, dimension, solutionAtNode);
  } 
  
  dat << endl;
  dat.close();
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
void PostProcessing3D::writeVtkDiscontinuous(string fileName, 
					     int N_LocVertices, std::vector< TVertex* > Vertices)
{

  char Disc[80];             // the new file name
  strcpy(Disc,fileName.c_str());     // copy the file name ...
  strcat(Disc,"_disc.vtk");  // ... add a string to the output file name
  
  std::ofstream dat(Disc);
  if (!dat)
  {
    Error("cannot open file for output\n");
    exit(-1);
  }
  
  dat.setf(std::ios::fixed);
  dat << setprecision(9);

  dat << "# vtk DataFile Version 4.2" << endl;
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
  for(unsigned int space_number = 0; space_number < FEFunctionArray.size(); space_number++)
  {
    const TFEFunction3D* fefunction = FEFunctionArray[space_number];
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

//TODO: alter Write_ParVTK
/** write stored PARALLEL data into a pvtu and vtu files (XML files for paraview) 
 *(Sashikumaar Ganesan) */
int PostProcessing3D::Write_ParVTK(
#ifdef _MPI
                                MPI_Comm comm,
#endif
                               int img, char *subID,
							   std::string directory,
							   std::string basename)
{
  int i, j, k,l,m,n, rank, size, N_, N_Elements, N_LocVertices;
  int N_Vertices, N_CellVertices;
  int N_LocDOF, Length, N_Comp, *GlobalNumbers, *BeginIndex, *DOF;
  int *VertexNumbers=nullptr, *NumberVertex=nullptr, begin, ID;

  double xi=0, eta=0, zeta=0, value;
  const double *Coeffs;
  double *WArray=nullptr, *DoubleArray=nullptr;
  double BFValues[MaxN_BaseFunctions3D];
  double *Coords=nullptr;
  static double HexaCoords[] = { -1, -1, -1, 1, -1, -1, 1,  1, -1, -1,  1, -1,
                                 -1, -1,  1, 1, -1,  1, 1,  1,  1, -1,  1,  1  };
  static double TetraCoords[] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1 };
  
  int N_ScalarVar = FEFunctionArray.size();
  int N_VectorVar = FEVectFunctArray.size();

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

  std::vector<TVertex *> Vertices;
  TVertex *Last, *Current;
  TBaseCell *cell;
  const TFESpace3D *fespace;
  TBaseFunct3D *bf;
  FE3D FE_ID;

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

    dat <<  "<?xml version=\"1.0\"?>"<< endl;
    dat << endl;
    dat <<  "<!--" << endl;
    dat <<  "      Title: Master file for parallel vtk data" << endl;
    dat <<  "    Program: ParMooN " << endl;
    dat <<  "    Version: v1.0.0 " << endl;
    dat <<  "Date & Time: " <<asctime (timeinfo) << endl;
    dat <<  "Problem Current Time " <<TDatabase::TimeDB->CURRENTTIME << endl;
    dat <<  "  -->" << endl;

    dat << endl;


    dat <<  "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1"
             <<"\" byte_order=\"LittleEndian\">"<<endl;
    dat <<  "<PUnstructuredGrid GhostLevel=\""<<0<<"\">"<<endl;
    dat << endl;

    dat <<  " <PPoints>"<<endl;
    dat <<  "   <PDataArray type=\"Float32\" Name=\""
        <<"Position\" NumberOfComponents=\"3\"/>"<<endl;
    dat <<   "</PPoints>"<<endl;
    dat << endl;

    dat <<   "<PCells>"<<endl;
    dat <<   "  <PDataArray type=\"Int32\" Name=\"connectivity\""
        <<" NumberOfComponents=\"1\"/>"<<endl;
    dat <<   "  <PDataArray type=\"Int32\" Name=\"offsets\""
        <<"      NumberOfComponents=\"1\"/>"<<endl;
    dat <<   "  <PDataArray type=\"UInt8\" Name=\"types\""
        <<"        NumberOfComponents=\"1\"/>"<<endl;
    dat <<   "</PCells>"<<endl;
    dat << endl;

    dat <<   "<PPointData Vectors=\"Vectors\" Scalars=\"Scalars\">"<<endl;
    for(i=0;i<N_VectorVar;i++)
    dat <<   "  <PDataArray type=\"Float32\" Name=\""<<FEVectFunctArray[i]->GetName()<<"\""
        <<" NumberOfComponents=\"3\" format=\"ascii\"/>"<<endl;
      dat << endl;

    for(i=0;i<N_VectorVar;i++)
    {
    for(j=0;j<FEVectFunctArray[i]->GetN_Components();j++)
     {
      dat <<  "  <DataArray type=\"Float32\" Name=\""
          <<FEVectFunctArray[i]->GetName()<<j<<"\" NumberOfComponents=\""
          <<"1\" format=\"ascii\"/>"<<endl;
      dat << endl;
     }
    }

    for(i=0;i<N_ScalarVar;i++)
    dat <<   "  <PDataArray type=\"Float32\" Name=\""<<FEFunctionArray[i]->GetName()<<"\""
        <<" NumberOfComponents=\"1\" format=\"ascii\"/>"<<endl;
    dat <<   "</PPointData>"<<endl;
    dat << endl;

    dat <<   "<PCellData Scalars=\"SubDomainAndRegionID\">"<<endl;
#ifdef _MPI    
    dat <<   "  <PDataArray type=\"Int32\"   Name=\"SubDomain\""
        <<"  NumberOfComponents=\"1\"/>"<<endl;
#endif
    dat <<   "  <PDataArray type=\"Int32\"   Name=\"RegionID\""
        <<"  NumberOfComponents=\"1\"/>"<<endl;

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
        if(i<10)        dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".000"<<ID<<".0000"<<img<<".vtu\"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".00" <<ID<<".0000"<<img<<".vtu\"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<".0000"<<img<<".vtu\"/>"<<endl;
        else            dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<"."   <<ID<<".0000"<<img<<".vtu\"/>"<<endl;
       }
      else if(img<100)
       {
        if(i<10)        dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".000"<<ID<<".000"<<img<<".vtu\"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".00" <<ID<<".000"<<img<<".vtu\"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<".000"<<img<<".vtu\"/>"<<endl;
        else            dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<"."   <<ID<<".000"<<img<<".vtu\"/>"<<endl;
       }
      else if(img<1000)
       {
        if(i<10)        dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".000"<<ID<<".00"<<img<<".vtu\"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".00" <<ID<<".00"<<img<<".vtu\"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<".00"<<img<<".vtu\"/>"<<endl;
        else            dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<"."   <<ID<<".00"<<img<<".vtu\"/>"<<endl;
       }
      else if(img<10000)
       {
        if(i<10)        dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".000"<<ID<<".0"<<img<<".vtu\"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".00" <<ID<<".0"<<img<<".vtu\"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<".0"<<img<<".vtu\"/>"<<endl;
        else            dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<"."   <<ID<<".0"<<img<<".vtu\"/>"<<endl;
       }
      else
       {
        if(i<10)        dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".000"<<ID<<"."<<img<<".vtu\"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".00" <<ID<<"."<<img<<".vtu\"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<".0"  <<ID<<"."<<img<<".vtu\"/>"<<endl;
        else            dat <<   "  <Piece Source=\""<<vtu<<VtkBaseName<<subID<<"."   <<ID<<"."<<img<<".vtu\"/>"<<endl;
       }
     }
    dat << endl;

    dat <<  "</PUnstructuredGrid>"<<endl;
    dat <<  "</VTKFile>"<<endl;
    dat.close();
   } // if(rank==0)
   
  // root take part in computation 
//   else
   {

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
   Vertices.resize(N_LocVertices);
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
   std::sort(Vertices.begin(), Vertices.end());

  Last=nullptr;
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
   Last=nullptr;
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
      l=std::distance(Vertices.begin(), std::lower_bound(Vertices.begin(), Vertices.end(), Current));
      VertexNumbers[m]=NumberVertex[l];
      m++;
    } // endfor j
  } //endfor i

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

    dat <<  "<?xml version=\"1.0\"?>" << endl;
    dat << endl;
    dat <<  "<!--" << endl;
    dat <<  "      Title: SubDomain data for master ptvu file" << endl;
    dat <<  "    Program: ParMooN " << endl;
    dat <<  "    Version: v1.0.0 " << endl;
    dat <<  "Date & Time: " <<asctime (timeinfo) << endl;
    dat <<  "  -->" << endl;
    dat << endl;

    dat <<  "<VTKFile type=\"UnstructuredGrid\" version=\"0.1"
             <<"\" byte_order=\"LittleEndian\">"<<endl;
    dat <<  "<UnstructuredGrid>"<<endl;
    dat << endl;

    dat <<  "<Piece NumberOfPoints=\""<<N_Vertices<<"\" NumberOfCells=\""<<N_Elements<<"\">"<<endl;
    dat <<  "<Points>"<<endl;
    dat <<  "  <DataArray type=\"Float32\" Name=\"Position\""
        <<" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
    N_=0;
    for(i=0;i<N_Vertices;i++)
    {
     dat <<  "   " << Coords[N_] << " " <<  Coords[N_+1] << " " << Coords[N_+2] << endl;
     N_ +=3;
    }
    dat <<  "  </DataArray>"<<endl;
    dat <<  "</Points>"<<endl;


    dat <<  "<Cells>"<<endl;
    dat <<  "  <DataArray type=\"Int32\" Name=\"connectivity\""
        <<" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;
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

    dat <<  "  <DataArray type=\"Int32\" Name=\"offsets\""
        <<" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;
    for(i=1;i<=N_Elements;i++)
      {
        N_CellVertices=Coll->GetCell(i-1)->GetN_Vertices();
        dat <<  i*N_CellVertices<<"  " ;
      }
    dat <<  "   </DataArray>"<<endl;
    dat << endl;

    dat <<  "  <DataArray type=\"UInt8\"  Name=\"types\""
        <<" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;
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
    dat <<  "<PointData Vectors=\"Velocity\" Scalars=\"Scalars\">"<<endl;
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


      dat <<  "  <DataArray type=\"Float32\" Name=\""
          <<FEVectFunctArray[k]->GetName()<<"\" NumberOfComponents=\""
          <<"3\" format=\"ascii\">"<<endl;
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
        dat <<  "  <DataArray type=\"Float32\" Name=\""
          <<FEVectFunctArray[k]->GetName()<<j<<"\" NumberOfComponents=\""
          <<"1\" format=\"ascii\">"<<endl;
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

      dat <<  "  <DataArray type=\"Float32\" Name=\""
          <<FEFunctionArray[k]->GetName()<<"\" NumberOfComponents=\""
          <<"1\" format=\"ascii\">"<<endl;
      for(j=0;j<N_Vertices;j++)
        dat << DoubleArray[j]<< " ";
       dat << endl;
      dat <<  "  </DataArray>"<<endl;

       dat << endl;

     }// for(k=0;k<N_ScalarVar

     
//     dat <<  "</PointData>"<<endl;
//     dat << endl;
// 
//     dat <<  "<CellData Scalars=\"SubDomain\">"<<endl;
//     dat <<  "  <DataArray type=\"Int32\" Name=\"Region\""
//         <<" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;
//     for(i=0;i<N_Elements;i++)
//       dat << (Coll->GetCell(i))->GetRegionID()   << " ";
//     dat <<  "  </DataArray>"<<endl;
//     dat <<  "</CellData>"<<endl;
// 
//  
     
    dat <<  "</PointData>"<<endl;
    dat << endl;

    dat <<  "<CellData Scalars=\"SubDomain\">"<<endl;
#ifdef _MPI    
    dat <<  "  <DataArray type=\"Int32\" Name=\"SubDomain\""
        <<" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;
    for(i=0;i<N_Elements;i++)     
      dat << (Coll->GetCell(i))->GetSubDomainNo() << " ";
    dat <<  "  </DataArray>"<<endl;
#endif    
    dat <<  "  <DataArray type=\"Int32\" Name=\"RegionID\""
        <<" NumberOfComponents=\"1\" format=\"ascii\">"<<endl;
     for(i=0;i<N_Elements;i++)
      dat << (Coll->GetCell(i))->GetRegionID() << " ";       
    dat <<  "  </DataArray>"<<endl;    

    dat <<  "</CellData>"<<endl;

    dat <<  "</Piece>"<<endl;
    dat <<  "</UnstructuredGrid>"<<endl;
    dat <<  "</VTKFile>"<<endl;

    dat.close();

  if(N_LocVertices)
   {
    delete [] NumberVertex;
    delete [] VertexNumbers;
    delete [] DoubleArray;
    delete [] WArray;
    delete [] Coords;
   }
  } //else root
//  if(rank==1)


  return 0;
} //TOutput3D::Write_ParVTK


/** @brief Output in .case format

   1 file .case (main file)
   1 file .geo (mesh)
   scalar and vector output files at each time step
   Note that currently we only write the P1 projection (node values) of the results
   However, the case format allows also element-wise and higher order
   visualizations (but not with paraview) 
   @todo implement visualization for different element type
*/
void PostProcessing3D::writeCaseFile()
{
  ofstream casf;
  string filename = testcaseDir + "/" + testcaseName;
  string casename = filename + ".case";
  Output::print<2>(" ** PostProcessing3D::writeCaseFile - write ", casename);
  casf.open(casename);
  
  string geoname = testcaseName + ".geo";
  
  casf << "FORMAT\n";
  casf << "type: ensight\n";
  casf << "GEOMETRY\n";
  casf << "model: 1 " << geoname << endl;
  casf << "VARIABLE\n";
  for (unsigned int j=0; j<FEFunctionArray.size(); j++)
  {
    unsigned int dimension = TFEDatabase3D::GetFE3D(
      FEFunctionArray[j]->GetFESpace3D()->GetFE3D(0,Coll->GetCell(0)))
      ->GetBaseFunct3D()->GetBaseVectDim();
    
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
void PostProcessing3D::writeCaseGeo()
{
  int dimension=3;
  
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
        switch(nVE)
	{
	  case 4:
	    ensight_type ="tetra4";
	    break;
	  case 8:
	    ensight_type ="hexa8"; 
	    break;
	}
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
        ErrThrow(" **ERROR: CASE output with mixed (tets and hexas) meshes not "
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

void PostProcessing3D::writeCaseVars(int iter)
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

    Output::print<3>(" ** PostProcessing3D::write File - write ", fname);
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
    Output::print<3>(" ** PostProcessing3D::write File - write ", fname);
    vctf.open(fname.c_str());
    vctf <<  FEVectFunctArray[i]->GetName() << " step = " << iter << endl; 
  
    writeVectCase(vctf, N_Vertices, dimension, solutionAtNode);
  }
}


void PostProcessing3D::writeCoord(ofstream& f, int dimension)
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
void PostProcessing3D::computeNodeValues(const T* function, vector< double >& solutionAtNode, 
					 unsigned int& dimension)
{
  auto FESpace3D = function->GetFESpace3D();
  TCollection* coll = FESpace3D->GetCollection();
  int nPoints = coll->GetN_Vertices();

  // compute FE type on first cell
  TBaseCell *cell = coll->GetCell(0);
  FE3D FE_ID = FESpace3D->GetFE3D(0, cell);
  TBaseFunct3D *bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
  
  // get function type of T
  bool IsVect = std::is_same< T, TFEVectFunct3D >::value;  
  dimension=1;
  if (IsVect || bf->GetBaseVectDim()>1)
    dimension = 3;

  solutionAtNode.assign(dimension*nPoints, 0.0);
  std::vector<double> WArray(nPoints, 0.);
  
  for(int ncell=0;ncell<coll->GetN_Cells();ncell++) 
  {
    cell = coll->GetCell(ncell);
    unsigned int nLocalVertices = cell->GetN_Vertices();

    bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();

    for(unsigned int nvert=0; nvert<nLocalVertices; nvert++) 
    {   
      unsigned int globalVert_index = coll->GetGlobalVerNo(ncell,nvert);
      
      double x, y, z;
#ifdef __2D__
          ;
#else
      cell->GetVertex(nvert)->GetCoords(x,y,z);
#endif
      double function_value[dimension];
    
      if (dimension>1) //computations for "vector valued scalar functions"
      {
	function->FindValueLocal(cell,ncell,x,y,z,function_value);	
	    
	solutionAtNode[dimension*globalVert_index + 0 ] += function_value[0];
	solutionAtNode[dimension*globalVert_index + 1 ] += function_value[1];
	solutionAtNode[dimension*globalVert_index + 2 ] += function_value[2];
      }  
      else	// compute nodal value of FE function
      {
	function->FindValueLocal(cell, ncell, x,y,z, function_value); 
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
      if (dimension>1) //dimension==3
      {
	solutionAtNode[dimension*nvert+1]/= WArray[nvert];
	solutionAtNode[dimension*nvert+2]/=WArray[nvert];
      }
    }
  }
}

// ****************
// scalar components of vector variables
// ****************
void PostProcessing3D::printVectCompwise(ofstream& dat,  std::string name,
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
	//if(fabs(value)<1e-18)
	//  value = 0.;
	dat << scientific;
	dat <<  value << endl;
      }
      dat << endl << endl;
    }
}

// ***************
// absolute value of vector variables
// ***************
void PostProcessing3D::printVectAbsValue(ofstream& dat, std::string name,
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
void PostProcessing3D::printVectPointwise(ofstream& dat, std::string name,
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
	//if(fabs(value)<1e-18)
	  //value = 0.;
	dat << scientific;
	dat <<  value << " ";
    }
    dat << endl;
  }
  dat << endl;
}

//auxiliary function for writeVectCase
void printEntry(ofstream & dat, double value, int counter)
{
    dat.setf(ios_base::scientific);	
    dat.precision(5);	
    dat.width(12);
    dat << value;

  // One line contains most six entries.    
    if( (counter+1) % 6 == 0 )
      dat << std::endl;
}

void PostProcessing3D::writeVectCase(ofstream& dat, 
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