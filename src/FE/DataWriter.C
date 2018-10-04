#include <DataWriter.h>
#ifdef __2D__
#include <FEVectFunct2D.h>
#else
#include <FEVectFunct3D.h>
#endif
#include <Vertex.h>

#include <algorithm>
#include <sys/stat.h>
#include <type_traits>

#ifdef _MPI
#include "mpi.h"
#endif

template <int d>
DataWriter<d>::DataWriter(const ParameterDatabase& param_db)
  : testcaseDir(), testcaseName(), n_steps_per_output(1),
    n_steps_until_next_output(0), FEFunctionArray(), FEVectFunctArray(),
    Coll(nullptr), timeValues()
{
  // set the variables depending on input parameters
  ParameterDatabase db = ParameterDatabase::default_output_database();
  db.merge(param_db, false);

  writeVTK = db["output_write_vtk"];
  writeVTU = db["output_write_vtu"];
  writeCASE = db["output_write_case"];
  if(writeCASE)
  {
    Output::warn("DataWriter", "case output is not working with older ",
                 "versions of ParaView on Linux. Try (at least) version 5.0");
  }

  testcaseName = db["output_basename"].get<std::string>();
  testcaseDir = db["output_vtk_directory"].get<std::string>();
  n_steps_per_output = std::max<size_t>(1, db["steps_per_output"]); // avoid 0
};

template <int d>
void DataWriter<d>::add_fe_function(const FEFunction* fefunction)
{
  // check that FE functions have the same collection
  if(this->Coll == nullptr)
  {
    this->Coll = fefunction->GetFESpace()->GetCollection();
  }
  else
  {
    if(this->Coll != fefunction->GetFESpace()->GetCollection())
    {
      // we could also just refuse to add this fe function and return without an
      // exception.
      ErrThrow("new FE function has a different collection");
    }
  }
  bool already_known = std::any_of(
      FEFunctionArray.begin(), FEFunctionArray.end(),
      [fefunction](const FEFunction* f) { return f == fefunction; });
  if(!already_known)
    FEFunctionArray.push_back(fefunction);
}

template <int d>
void DataWriter<d>::add_fe_vector_function(const FEVectFunct* fevectfunction)
{
  // check that FE functions have the same collection
  if(this->Coll == nullptr)
  {
    this->Coll = fevectfunction->GetFESpace()->GetCollection();
  }
  else
  {
    if(this->Coll != fevectfunction->GetFESpace()->GetCollection())
    {
      // we could also just refuse to add this fe function and return without an
      // exception.
      ErrThrow("new FE vector function has a different collection");
    }
  }
  bool already_known = std::any_of(
      FEVectFunctArray.begin(), FEVectFunctArray.end(),
      [fevectfunction](const FEVectFunct* f) { return f == fevectfunction; });
  if(!already_known)
    FEVectFunctArray.push_back(fevectfunction);
}

template <int d>
void DataWriter<d>::write()
{
  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  if(my_rank == 0 && (writeVTK || writeVTU || writeCASE))
  {
    mkdir(testcaseDir.c_str(), 0777);
  }
  std::string name = testcaseDir + "/" + testcaseName;
  
  if(writeVTK)
  {
#ifdef _MPI
    char SubID[] = "";
    Write_ParVTK(MPI_COMM_WORLD, 0, SubID, testcaseDir, testcaseName);
#else
    Output::print<2>(" DataWriter:: writing ", name + ".vtk");
    writeVtk(name + ".vtk");
#endif
  }
  if(writeVTU)
  {
    Output::print<2>(" DataWriter:: writing ", name + ".vtu");
    writeVtu(name + ".vtu");
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

template <int d>
void DataWriter<d>::write(double current_time)
{
  if(n_steps_until_next_output != 0)
  {
    n_steps_until_next_output--;
  }
  else
  {
    n_steps_until_next_output = n_steps_per_output - 1;
    timeValues.push_back(current_time);
    // index of current time step in 'timeValues'
    auto index = timeValues.size() - 1;
    int my_rank = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    if(my_rank == 0 && (writeVTK || writeVTU || writeCASE))
    {
      mkdir(testcaseDir.c_str(), 0777);
    }
    std::string name;
    name += testcaseDir + "/" + testcaseName;
    name += std::to_string(index * n_steps_per_output);  
    if(writeVTK)
    {
#ifdef _MPI
      char SubID[] = "";
      Write_ParVTK(MPI_COMM_WORLD, index * n_steps_per_output, SubID,
                   testcaseDir, testcaseName);
#else
      Output::print<2>(" DataWriter:: writing ", name + ".vtk");
      writeVtk(name + ".vtk");
#endif
    }
    if(writeVTU)
    {
      writeVtu(name + ".vtu");
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
      writeCaseVars(index * n_steps_per_output);
      writeCaseFile();
    }
  }
}

template <int d>
int n_local_vertices_to_type(int n_loc_vert)
{
  switch(d)
  {
    case 2:
      switch(n_loc_vert)
      {
        case 4: return 9; break;
        case 3: return 5; break;
      }
      break;
    case 3:
      switch(n_loc_vert)
      {
        case 4: return 10; break;
        case 8: return 12; break;
      }
      break;
  }
  ErrThrow("a ", d, "D cell with ", n_loc_vert, " vertices is not supported");
}

template <int d>
void DataWriter<d>::writeMesh(std::string name)
{
  std::ofstream dat(name);
  if(!dat)
  {
    ErrThrow("cannot open file for output. ", name);
  }
  dat << setprecision(12);

  dat << "# vtk DataFile Version 4.0\n";
  dat << "file created by ParMooN"
      << " Time < " << (timeValues.empty() ? 0. : timeValues[0]) << " >\n";
  dat << "ASCII\n";
  dat << "DATASET UNSTRUCTURED_GRID\n";
  dat << "POINTS " << Coll->GetN_Vertices() << " double\n";

  writeCoord(dat);

  unsigned int N_LocVertices = Coll->GetNLocVertices();
  unsigned int N_Elements = Coll->GetN_Cells();

  dat << "CELLS " << N_Elements << " " << N_Elements + N_LocVertices << "\n";
  for(unsigned int i = 0; i < N_Elements; i++)
  {
    int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
    dat << N_CellVertices << " ";
    for(int j = 0; j < N_CellVertices; j++)
    {
      dat << this->Coll->GetGlobalVerNo(i, j) << " ";
    }
    dat << "\n";
  }
  dat << "\n";

  dat << "CELL_TYPES " << N_Elements << "\n";
  for(unsigned int i = 0; i < N_Elements; i++)
  {
    int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
    dat << n_local_vertices_to_type<d>(N_CellVertices) << " ";
  }
}

template <int d>
void DataWriter<d>::writeVtk(std::string name)
{
  // determine data for vtk file
  unsigned int N_Elements = Coll->GetN_Cells();

  // check for discontinuous scalar variables. In such a case write a new file
  // for this variable. ParaView will really display a discontinuous function
  // instead of projecting it onto P1/Q1 space. However there are some
  // drawbacks.
  std::vector<const TVertex*> Vertices;
  for(unsigned int ncell = 0; ncell < N_Elements; ncell++)
  {
    const TBaseCell* cell = Coll->GetCell(ncell);
    unsigned int N_Loc = cell->GetN_Vertices();
    for(unsigned int nvert = 0; nvert < N_Loc; nvert++)
      Vertices.push_back(cell->GetVertex(nvert));
  }

  for(unsigned int i = 0; i < FEFunctionArray.size(); i++)
  {
    bool isDiscont = FEFunctionArray[i]->GetFESpace()->IsDGSpace();
    if(isDiscont)
    {
      writeVtkDiscontinuous(name, Vertices.size(), Vertices);
      break;
    }
  }

  writeMesh(name);
  std::ofstream dat(name, std::ios::out | std::ios::app);
  if(!dat)
  {
    ErrThrow("cannot open file for output. ", name);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(12);

  unsigned int N_Vertices = Coll->GetN_Vertices();
  unsigned int dimension;

  dat << "\n\n";
  dat << "POINT_DATA " << N_Vertices << "\n";

  std::vector<double> solutionAtNode;
  for(unsigned int i = 0; i < FEFunctionArray.size(); ++i)
  {
    computeNodeValues(FEFunctionArray.at(i), solutionAtNode, dimension);
    std::string name = FEFunctionArray.at(i)->GetName();
    printVectCompwise(dat, name, N_Vertices, dimension, solutionAtNode);
    if(dimension > 1)
    {
      printVectAbsValue(dat, name, N_Vertices, dimension, solutionAtNode);
      printVectPointwise(dat, name, N_Vertices, dimension, solutionAtNode);
    }
  }
  for(unsigned int i = 0; i < FEVectFunctArray.size(); ++i)
  {
    computeNodeValues(FEVectFunctArray.at(i), solutionAtNode, dimension);
    std::string name = FEVectFunctArray.at(i)->GetName();
    printVectCompwise(dat, name, N_Vertices, dimension, solutionAtNode);
    printVectAbsValue(dat, name, N_Vertices, dimension, solutionAtNode);
    printVectPointwise(dat, name, N_Vertices, dimension, solutionAtNode);
  }
  dat << "\n";
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
template <int d>
void DataWriter<d>::writeVtkDiscontinuous(std::string fileName,
                                          int N_LocVertices,
                                          std::vector<const TVertex*> Vertices)
{
  // copy the file name add a string to the output file name
  std::string disc = fileName + std::string("_disc.vtk");

  std::ofstream dat(disc);
  if(!dat)
  {
    ErrThrow("cannot open file for output ", disc);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(12);

  dat << "# vtk DataFile Version 4.0\n";
  dat << "file created by ParMooN.\n";
  dat << "ASCII\n";
  dat << "DATASET UNSTRUCTURED_GRID\n";
  dat << "POINTS " << N_LocVertices << " double\n";

  for(int i = 0; i < N_LocVertices; i++)
  {
    double x, y, z;
    Vertices[i]->GetCoords(x, y, z);
    dat << x << " " << y << " " << z << "\n";
  }
  dat << "\n";
  int N_Elements = Coll->GetN_Cells();

  // writing which vertices belong to which cells, here it is ignored that
  // a vertex might belong to multiple cells.
  dat << "CELLS " << N_Elements << " " << N_Elements + N_LocVertices << "\n";
  int l = 0;
  for(int i = 0; i < N_Elements; i++)
  {
    auto current_cell = Coll->GetCell(i);
    int N_CellVertices = current_cell->GetN_Vertices();
    dat << N_CellVertices << " ";
    for(int j = 0; j < N_CellVertices; j++)
    {
      dat << l << " ";
      l++;
    }
    dat << "\n";
  }
  dat << "\n";

  // the cell types tell paraview if this is a tetrahedron or a hexahedron
  // (export of other types is not supported here)
  dat << "CELL_TYPES " << N_Elements << "\n";
  for(int i = 0; i < N_Elements; i++)
  {
    int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
    dat << n_local_vertices_to_type<d>(N_CellVertices) << " ";
  }
  dat << "\n\n";

  // write the function values, only for scalar functions, which includes
  // vector valued basis functions (such as Raviart-Thomas), because these
  // are handled like scalar basis functions
  dat << "POINT_DATA " << N_LocVertices << "\n";
  for(unsigned int i = 0; i < FEFunctionArray.size(); ++i)
  {
    const auto* fefunction = FEFunctionArray[i];
    const int BaseVectDim = fefunction->GetFESpace()->GetBaseVectDim();

    // scalar valued basis functions (normal case)
    if(BaseVectDim == 1)
    {
      dat << "\n\n";
      dat << "SCALARS " << fefunction->GetName();
      dat << " double\n";
      dat << "LOOKUP_TABLE default\n";
      double function_value;
      for(int i = 0; i < N_Elements; i++)
      {
        TBaseCell* current_cell = Coll->GetCell(i);
        int N_CellVertices = current_cell->GetN_Vertices();
        for(int j = 0; j < N_CellVertices; j++)
        {
          double x, y, z;
          current_cell->GetVertex(j)->GetCoords(x, y, z);
#ifdef __2D__
          fefunction->FindValueLocal(current_cell, i, x, y, &function_value);
#else
          fefunction->FindValueLocal(current_cell, i, x, y, z, &function_value);
#endif
          dat << function_value << "\n";
        }
      }
    }
    // vector valued basis functions (e.g. Raviart-Thomas)
    else if(BaseVectDim == 2 || BaseVectDim == 3)
    {
      // find values for all components
      double function_value[BaseVectDim];
      dat << "\n\n";
      dat << "VECTORS " << fefunction->GetName();
      dat << " double\n";
      for(int i = 0, k = 0; i < N_Elements; i++)
      {
        TBaseCell* current_cell = Coll->GetCell(i);
        int N_CellVertices = current_cell->GetN_Vertices();
        for(int j = 0; j < N_CellVertices; j++)
        {
          double x, y, z;
          current_cell->GetVertex(j)->GetCoords(x, y, z);
#ifdef __3D__
          fefunction->FindValueLocal(current_cell, i, x, y, z, function_value);
#else
          fefunction->FindValueLocal(current_cell, i, x, y, function_value);
#endif
          for(int l = 0; l < BaseVectDim; ++l)
            dat << function_value[l] << "\t";
          if(BaseVectDim == 2)
            dat << double(0.);
          dat << "\n";
          k++;
        }
      }
    }
    else
    {
      ErrThrow("DataWriter::WriteVtkDiscontinuous: Basis functions of "
               "dimension ",
               BaseVectDim, " are not supported");
    }
  }
  dat.close();
}

template<int d>
void DataWriter<d>::writeVtu(std::string name) const
{
  std::ofstream f(name);
  
  const std::string byte_order = "LittleEndian"; // "BigEndian"
  std::string ascii_or_binary = "ascii"; // "binary"
  bool compression = false;
  if(!f)
  {
    ErrThrow("cannot open file for output. ", name);
  }
  f << std::fixed << std::scientific << std::setprecision(16);
  f << "<?xml version=\"1.0\" ?>"             << "\n"
    << "<!--"                                 << "\n"
    << "# This file was generated by ParMooN" << "\n"
    << "-->"                                  << "\n"
    << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
  if(compression)
    f << " compressor=\"vtkZLibDataCompressor\"";
  f << " byte_order=\""  << byte_order << "\"";
  f << ">"                        << "\n"
    << "  <UnstructuredGrid>"     << "\n";
  unsigned int N_LocVertices = Coll->GetNLocVertices();
  unsigned int n_cells = Coll->GetN_Cells();
  f << "  <Piece NumberOfPoints=\"" << N_LocVertices
    << "\" NumberOfCells=\""      << n_cells << "\" >" << "\n";
  
  
  // write the vertices of the grid
  f << "  <Points>\n";
  f << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\""
    << ascii_or_binary << "\">\n";
  
  std::vector<const TVertex*> Vertices;
  for(unsigned int icell = 0; icell < n_cells; icell++)
  {
    const TBaseCell* cell = Coll->GetCell(icell);
    unsigned int n_Loc_vert = cell->GetN_Vertices();
    for(unsigned int ivert = 0; ivert < n_Loc_vert; ++ivert)
    {
      const TVertex * v = cell->GetVertex(ivert);
      Vertices.push_back(v);
      double x, y, z;
      v->GetCoords(x, y, z);
      f << setw(16) << x << " " << setw(16) << y << " " << setw(16) << z 
        << "\n";
    }
  }
  f << "    </DataArray>\n";
  f << "  </Points>\n";
  
  // Write cell connectivity
  f << "  <Cells>" << "\n";
  f << "    <DataArray  type=\"Int32\"  Name=\"connectivity\"  format=\""
    << ascii_or_binary << "\">" << "\n";
  for(unsigned int icell = 0, loc_vert = 0; icell < n_cells; icell++)
  {
    const TBaseCell* cell = Coll->GetCell(icell);
    unsigned int n_Loc_vert = cell->GetN_Vertices();
    f << loc_vert;
    for(unsigned int ivert = 1; ivert < n_Loc_vert; ++ivert)
    {
      f << " " << loc_vert + ivert;
    }
    loc_vert += n_Loc_vert;
    f << "\n";
  }
  f << "    </DataArray>" << "\n";
  // Write offset into connectivity array for the end of each cell
  f << "    <DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"" 
    << ascii_or_binary << "\">" << "\n";
  for(unsigned int icell = 0, loc_vert = 0; icell < n_cells; icell++)
  {
    const TBaseCell* cell = Coll->GetCell(icell);
    unsigned int n_Loc_vert = cell->GetN_Vertices();
    loc_vert += n_Loc_vert;
    f << loc_vert << ((icell !=0 && icell%10 == 0) ? "\n" : " ");
  }
  f << "\n" << "    </DataArray>" << "\n";
  // Write cell type
  f << "    <DataArray  type=\"UInt32\"  Name=\"types\"  format=\"" 
    << ascii_or_binary << "\">" << "\n";
  for(unsigned int icell = 0; icell < n_cells; icell++)
  {
    const TBaseCell* cell = Coll->GetCell(icell);
    unsigned int n_Loc_vert = cell->GetN_Vertices();
    f << n_local_vertices_to_type<d>(n_Loc_vert);
    f << ((icell !=0 && icell%10 == 0) ? "\n" : " ");
  }
  f << "\n" << "    </DataArray>" << "\n";
  f << "  </Cells>" << "\n";
  
  // writing data
  std::string name_default; // paraview will show this data at first
  if(!FEFunctionArray.empty())
    name_default = FEFunctionArray[0]->GetName();
  else if(!FEVectFunctArray.empty())
    name_default = FEVectFunctArray[0]->GetName();
  
  f << "  <PointData Scalars=\"" << name_default << "\">" << "\n";
  for(unsigned int i = 0; i < FEFunctionArray.size(); ++i)
  {
    const auto* fefunction = FEFunctionArray[i];
    const int BaseVectDim = fefunction->GetFESpace()->GetBaseVectDim();
    if(BaseVectDim == 1)
    {
      f << "    <DataArray  type=\"Float64\"  Name=\"" << fefunction->GetName()
        << "\" format=\"" << ascii_or_binary << "\">" << "\n";
      for(unsigned int icell = 0; icell < n_cells; icell++)
      {
        const TBaseCell* cell = Coll->GetCell(icell);
        unsigned int n_Loc_vert = cell->GetN_Vertices();
        double function_value;
        for(unsigned int j = 0; j < n_Loc_vert; j++)
        {
          double x, y, z;
          cell->GetVertex(j)->GetCoords(x, y, z);
#ifdef __2D__
          fefunction->FindValueLocal(cell, icell, x, y, &function_value);
#else
          fefunction->FindValueLocal(cell, icell, x, y, z, &function_value);
#endif
          f << function_value << "\n";
        }
      }
      f << "    </DataArray> "     << "\n";
    }
    // vector valued basis functions (e.g. Raviart-Thomas)
    else if(BaseVectDim == d)
    {
      f << "    <DataArray  type=\"Float64\"  Name=\"" << fefunction->GetName()
        << "\" NumberOfComponents=\"" << d << "\" format=\"" 
        << ascii_or_binary << "\">" << "\n";
      // find values for all components
      double function_value[BaseVectDim];
      for(unsigned int icell = 0; icell < n_cells; icell++)
      {
        const TBaseCell* cell = Coll->GetCell(icell);
        unsigned int n_Loc_vert = cell->GetN_Vertices();
        for(unsigned int j = 0; j < n_Loc_vert; j++)
        {
          double x, y, z;
          cell->GetVertex(j)->GetCoords(x, y, z);
#ifdef __3D__
          fefunction->FindValueLocal(cell, icell, x, y, z, function_value);
#else
          fefunction->FindValueLocal(cell, icell, x, y, function_value);
#endif
          for(int l = 0; l < BaseVectDim; l++)
          {
            f << function_value[l] << " ";
          }
          f << "\n";
        }
      }
      f << "    </DataArray> "     << "\n";
    }
    else
    {
      ErrThrow("DataWriter::WriteVtkDiscontinuous: Basis functions of "
               "dimension ",
               BaseVectDim, " are not supported");
    }
  }
  
  for(unsigned int i = 0; i < FEVectFunctArray.size(); i++)
  {
    const auto* fe_vect = FEVectFunctArray[i];
    f << "    <DataArray  type=\"Float64\"  Name=\"" << fe_vect->GetName()
      << "\" NumberOfComponents=\""<< d << "\" format=\"" 
      << ascii_or_binary << "\">" << "\n";
    for(unsigned int icell = 0; icell < n_cells; icell++)
    {
      const TBaseCell* cell = Coll->GetCell(icell);
      unsigned int n_Loc_vert = cell->GetN_Vertices();
      for(unsigned int j = 0; j < n_Loc_vert; j++)
      {
        double x, y, z;
        cell->GetVertex(j)->GetCoords(x, y, z);
        double function_value[d];
#ifdef __2D__
        fe_vect->FindValueLocal(cell, icell, x, y, function_value);
#else
        fe_vect->FindValueLocal(cell, icell, x, y, z, function_value);
#endif
        for(unsigned int dim = 0; dim < d; ++dim)
        {
          f << function_value[dim] << " ";
        }
        f << "\n";
      }
    }
    f << "    </DataArray> "     << "\n";
  }
  
  f << "  </PointData>" << "\n";
  f << "  </Piece>"            << "\n";
  f << "  </UnstructuredGrid>" << "\n";
  f << "</VTKFile>\n";
}

// TODO: alter Write_ParVTK
/** write stored PARALLEL data into a pvtu and vtu files (XML files for
 *paraview)
 *(Sashikumaar Ganesan) */
#ifdef _MPI
template <int d>
void DataWriter<d>::Write_ParVTK(MPI_Comm comm, int img, char* subID,
                                 std::string directory, std::string basename)
{
  int i, j, k, l, m, n, rank, size, N_, N_Elements, N_LocVertices;
  int N_Vertices, N_CellVertices;
  int N_LocDOF, Length, N_Comp, *GlobalNumbers, *BeginIndex, *DOF;
  int *VertexNumbers = nullptr, *NumberVertex = nullptr, begin, ID;

  double xi = 0, eta = 0, zeta = 0, value;
  const double* Coeffs;
  double *WArray = nullptr, *DoubleArray = nullptr;
  double BFValues[MaxN_BaseFunctions3D];
  double* Coords = nullptr;
  static double HexaCoords[] = {-1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1,
                                -1, -1, 1,  1, -1, 1,  1, 1, 1,  -1, 1, 1};
  static double TetraCoords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};

  int N_ScalarVar = FEFunctionArray.size();
  int N_VectorVar = FEVectFunctArray.size();

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  std::string outputdir(directory);
  std::string vtu("VTU/");
  std::string vtudir = outputdir + std::string("/") + vtu;

  // get the proper path for the VTU directory (which contains the parts)
  //  if (rank == 0)
  //  {
  //    Output::print(vtudir);
  //  }

  time_t rawtime;
  struct tm* timeinfo;

  std::vector<TVertex*> Vertices;
  TVertex *Last, *Current;
  TBaseCell* cell;
  const TFESpace3D* fespace;
  TBaseFunct3D* bf;

  const char* VtkBaseName = basename.c_str();
  const char* output_directory = directory.c_str();

  if(rank == 0)
  {
    //     remove(vtudir);
    mkdir(vtudir.c_str(),
          0777); // create the folder to store SubDomain vtu files
  }

  MPI_Barrier(comm);

  std::ostringstream os;
  os << " ";

  time(&rawtime);
  timeinfo = localtime(&rawtime);


  // write the master pvtu file
  if(rank == 0)
  {
    Output::info("Output3D", "writing output into ", output_directory, "/",
                 VtkBaseName, subID, "*.", img, " xml vtk file");
    os.seekp(std::ios::beg);
    os << output_directory << "/" << VtkBaseName << subID;
    if(img < 10)
      os << ".0000" << img << ".pvtu" << ends;
    else if(img < 100)
      os << ".000" << img << ".pvtu" << ends;
    else if(img < 1000)
      os << ".00" << img << ".pvtu" << ends;
    else if(img < 10000)
      os << ".0" << img << ".pvtu" << ends;
    else
      os << "." << img << ".pvtu" << ends;
    std::ofstream dat(os.str().c_str());
    if(!dat)
    {
      cerr << "cannot open file for output\n";
      MPI_Abort(MPI_COMM_WORLD, 0);
    }

    dat << "<?xml version=\"1.0\"?>\n";
    dat << "\n";
    dat << "<!--\n";
    dat << "      Title: Master file for parallel vtk data\n";
    dat << "    Program: ParMooN\n";
    dat << "    Version: v1.0.0\n";
    dat << "Date & Time: " << asctime(timeinfo) << "\n";
    dat << "Problem Current Time " << (timeValues.empty() ? 0. : timeValues[0])
        << "\n";
    dat << "  -->\n";
    dat << "\n";


    dat << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1"
        << "\" byte_order=\"LittleEndian\">\n";
    dat << "<PUnstructuredGrid GhostLevel=\"0\">\n";
    dat << "\n";

    dat << " <PPoints>\n";
    dat << "   <PDataArray type=\"Float32\" Name=\""
        << "Position\" NumberOfComponents=\"3\"/>\n";
    dat << "</PPoints>\n";
    dat << "\n";

    dat << "<PCells>\n";
    dat << "  <PDataArray type=\"Int32\" Name=\"connectivity\""
        << " NumberOfComponents=\"1\"/>\n";
    dat << "  <PDataArray type=\"Int32\" Name=\"offsets\""
        << "      NumberOfComponents=\"1\"/>\n";
    dat << "  <PDataArray type=\"UInt8\" Name=\"types\""
        << "        NumberOfComponents=\"1\"/>\n";
    dat << "</PCells>\n";
    dat << "\n";

    dat << "<PPointData Vectors=\"Vectors\" Scalars=\"Scalars\">\n";
    for(i = 0; i < N_VectorVar; i++)
      dat << "  <PDataArray type=\"Float32\" Name=\""
          << FEVectFunctArray[i]->GetName() << "\""
          << " NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    dat << "\n";

    for(i = 0; i < N_VectorVar; i++)
    {
      for(j = 0; j < FEVectFunctArray[i]->GetN_Components(); j++)
      {
        dat << "  <DataArray type=\"Float32\" Name=\""
            << FEVectFunctArray[i]->GetName() << j << "\" NumberOfComponents=\""
            << "1\" format=\"ascii\"/>\n";
        dat << "\n";
      }
    }

    for(i = 0; i < N_ScalarVar; i++)
      dat << "  <PDataArray type=\"Float32\" Name=\""
          << FEFunctionArray[i]->GetName() << "\""
          << " NumberOfComponents=\"1\" format=\"ascii\"/>\n";
    dat << "</PPointData>\n";
    dat << "\n";

    dat << "<PCellData Scalars=\"SubDomainAndRegionID\">\n";
    dat << "  <PDataArray type=\"Int32\"   Name=\"SubDomain\""
        << "  NumberOfComponents=\"1\"/>\n";
    dat << "  <PDataArray type=\"Int32\"   Name=\"RegionID\""
        << "  NumberOfComponents=\"1\"/>\n";

    dat << "</PCellData>\n";
    dat << "\n";


    // root not take part in computation
    //     begin = 1;
    begin = 0; // root take part in computation

    for(i = begin; i < size; i++)
    {
      ID = i; // root not take part in computation

      if(img < 10)
      {
        if(i < 10)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".000"
              << ID << ".0000" << img << ".vtu\"/>\n";
        else if(i < 100)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".00"
              << ID << ".0000" << img << ".vtu\"/>\n";
        else if(i < 1000)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".0"
              << ID << ".0000" << img << ".vtu\"/>\n";
        else
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << "."
              << ID << ".0000" << img << ".vtu\"/>\n";
      }
      else if(img < 100)
      {
        if(i < 10)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".000"
              << ID << ".000" << img << ".vtu\"/>\n";
        else if(i < 100)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".00"
              << ID << ".000" << img << ".vtu\"/>\n";
        else if(i < 1000)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".0"
              << ID << ".000" << img << ".vtu\"/>\n";
        else
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << "."
              << ID << ".000" << img << ".vtu\"/>\n";
      }
      else if(img < 1000)
      {
        if(i < 10)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".000"
              << ID << ".00" << img << ".vtu\"/>\n";
        else if(i < 100)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".00"
              << ID << ".00" << img << ".vtu\"/>\n";
        else if(i < 1000)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".0"
              << ID << ".00" << img << ".vtu\"/>\n";
        else
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << "."
              << ID << ".00" << img << ".vtu\"/>\n";
      }
      else if(img < 10000)
      {
        if(i < 10)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".000"
              << ID << ".0" << img << ".vtu\"/>\n";
        else if(i < 100)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".00"
              << ID << ".0" << img << ".vtu\"/>\n";
        else if(i < 1000)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".0"
              << ID << ".0" << img << ".vtu\"/>\n";
        else
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << "."
              << ID << ".0" << img << ".vtu\"/>\n";
      }
      else
      {
        if(i < 10)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".000"
              << ID << "." << img << ".vtu\"/>\n";
        else if(i < 100)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".00"
              << ID << "." << img << ".vtu\"/>\n";
        else if(i < 1000)
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << ".0"
              << ID << "." << img << ".vtu\"/>\n";
        else
          dat << "  <Piece Source=\"" << vtu << VtkBaseName << subID << "."
              << ID << "." << img << ".vtu\"/>\n";
      }
    }
    dat << "\n";

    dat << "</PUnstructuredGrid>\n";
    dat << "</VTKFile>\n";
    dat.close();
  } // if(rank==0)

  N_Elements = Coll->GetN_OwnCells();
  N_LocVertices = 0;
  for(i = 0; i < N_Elements; i++)
  {
    cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
  }

  if(N_LocVertices)
    Vertices.resize(N_LocVertices);
  N_ = 0;
  for(i = 0; i < N_Elements; i++)
  {
    cell = Coll->GetCell(i);
    k = cell->GetN_Vertices();
    for(j = 0; j < k; j++)
    {
      Vertices[N_] = cell->GetVertex(j);
      N_++;
    }
  }


  if(N_)
    std::sort(Vertices.begin(), Vertices.end());

  Last = nullptr;
  N_Vertices = 0;
  for(i = 0; i < N_LocVertices; i++)
  {
    if((Current = Vertices[i]) != Last)
    {
      N_Vertices++;
      Last = Current;
    }
  }

  if(N_LocVertices)
  {
    Coords = new double[3 * N_Vertices];
    VertexNumbers = new int[N_LocVertices];
    NumberVertex = new int[N_LocVertices];
  }
  Last = nullptr;
  N_ = 0;
  k = -1;
  for(i = 0; i < N_LocVertices; i++)
  {
    if((Current = Vertices[i]) != Last)
    {
      Vertices[i]->GetCoords(Coords[3 * N_], Coords[3 * N_ + 1],
                             Coords[3 * N_ + 2]);
      k++;
      N_++;
      Last = Current;
    }
    NumberVertex[i] = k;
  }

  m = 0;
  for(i = 0; i < N_Elements; i++)
  {
    cell = Coll->GetCell(i);
    k = cell->GetN_Vertices();
    for(j = 0; j < k; j++)
    {
      Current = cell->GetVertex(j);
      l = std::distance(
          Vertices.begin(),
          std::lower_bound(Vertices.begin(), Vertices.end(), Current));
      VertexNumbers[m] = NumberVertex[l];
      m++;
    } // endfor j
  }   // endfor i

  // to check
  // cout << "MaxN_VerticesPerCell*N_Comps" << MaxN_VerticesPerCell*N_Comps <<
  // "\n";
  // cout << "MaxN_VerticesPerCell" << MaxN_VerticesPerCell << "\n";
  // cout << "N_Comps" << N_Comps << "\n";

  if(N_LocVertices)
  {
    DoubleArray = new double[3 * N_Vertices];
    WArray = new double[N_Vertices];
  }

  ID = rank;
  os.seekp(std::ios::beg);
  if(img < 10)
  {
    if(ID < 10)
      os << vtudir << VtkBaseName << subID << ".000" << ID << ".0000" << img
         << ".vtu" << ends;
    else if(ID < 100)
      os << vtudir << VtkBaseName << subID << ".00" << ID << ".0000" << img
         << ".vtu" << ends;
    else if(ID < 1000)
      os << vtudir << VtkBaseName << subID << ".0" << ID << ".0000" << img
         << ".vtu" << ends;
    else
      os << vtudir << VtkBaseName << subID << "." << ID << ".0000" << img
         << ".vtu" << ends;
  }
  else if(img < 100)
  {
    if(ID < 10)
      os << vtudir << VtkBaseName << subID << ".000" << ID << ".000" << img
         << ".vtu" << ends;
    else if(ID < 100)
      os << vtudir << VtkBaseName << subID << ".00" << ID << ".000" << img
         << ".vtu" << ends;
    else if(ID < 1000)
      os << vtudir << VtkBaseName << subID << ".0" << ID << ".000" << img
         << ".vtu" << ends;
    else
      os << vtudir << VtkBaseName << subID << "." << ID << ".000" << img
         << ".vtu" << ends;
  }
  else if(img < 1000)
  {
    if(ID < 10)
      os << vtudir << VtkBaseName << subID << ".000" << ID << ".00" << img
         << ".vtu" << ends;
    else if(ID < 100)
      os << vtudir << VtkBaseName << subID << ".00" << ID << ".00" << img
         << ".vtu" << ends;
    else if(ID < 1000)
      os << vtudir << VtkBaseName << subID << ".0" << ID << ".00" << img
         << ".vtu" << ends;
    else
      os << vtudir << VtkBaseName << subID << "." << ID << ".00" << img
         << ".vtu" << ends;
  }
  else if(img < 10000)
  {
    if(ID < 10)
      os << vtudir << VtkBaseName << subID << ".000" << ID << ".0" << img
         << ".vtu" << ends;
    else if(ID < 100)
      os << vtudir << VtkBaseName << subID << ".00" << ID << ".0" << img
         << ".vtu" << ends;
    else if(ID < 1000)
      os << vtudir << VtkBaseName << subID << ".0" << ID << ".0" << img
         << ".vtu" << ends;
    else
      os << vtudir << VtkBaseName << subID << "." << ID << ".0" << img << ".vtu"
         << ends;
  }
  else
  {
    if(ID < 10)
      os << vtudir << VtkBaseName << subID << ".000" << ID << "." << img
         << ".vtu" << ends;
    else if(ID < 100)
      os << vtudir << VtkBaseName << subID << ".00" << ID << "." << img
         << ".vtu" << ends;
    else if(ID < 1000)
      os << vtudir << VtkBaseName << subID << ".0" << ID << "." << img << ".vtu"
         << ends;
    else
      os << vtudir << VtkBaseName << subID << "." << ID << "." << img << ".vtu"
         << ends;
  }

  std::ofstream dat(os.str().c_str());
  if(!dat)
  {
    Error("cannot open file for output");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  dat << setprecision(8);

  dat << "<?xml version=\"1.0\"?>\n";
  dat << "\n";
  dat << "<!--\n";
  dat << "      Title: SubDomain data for master ptvu file\n";
  dat << "    Program: ParMooN\n";
  dat << "    Version: v1.0.0\n";
  dat << "Date & Time: " << asctime(timeinfo) << "\n";
  dat << "  -->\n";
  dat << "\n";

  dat << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1"
      << "\" byte_order=\"LittleEndian\">\n";
  dat << "<UnstructuredGrid>\n";
  dat << "\n";

  dat << "<Piece NumberOfPoints=\"" << N_Vertices << "\" NumberOfCells=\""
      << N_Elements << "\">\n";
  dat << "<Points>\n";
  dat << "  <DataArray type=\"Float32\" Name=\"Position\""
      << " NumberOfComponents=\"3\" format=\"ascii\">\n";
  N_ = 0;
  for(i = 0; i < N_Vertices; i++)
  {
    dat << "   " << Coords[N_] << " " << Coords[N_ + 1] << " " << Coords[N_ + 2]
        << "\n";
    N_ += 3;
  }
  dat << "  </DataArray>\n";
  dat << "</Points>\n";


  dat << "<Cells>\n";
  dat << "  <DataArray type=\"Int32\" Name=\"connectivity\""
      << " NumberOfComponents=\"1\" format=\"ascii\">\n";
  l = 0;
  for(i = 0; i < N_Elements; i++)
  {
    N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
    dat << "       ";
    for(j = 0; j < N_CellVertices; j++)
    {
      dat << VertexNumbers[l] << " ";
      l++;
    }
    dat << "\n";
  }
  dat << "  </DataArray>\n";
  dat << "\n";

  dat << "  <DataArray type=\"Int32\" Name=\"offsets\""
      << " NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(i = 1; i <= N_Elements; i++)
  {
    N_CellVertices = Coll->GetCell(i - 1)->GetN_Vertices();
    dat << i * N_CellVertices << "  ";
  }
  dat << "   </DataArray>\n";
  dat << "\n";

  dat << "  <DataArray type=\"UInt8\"  Name=\"types\""
      << " NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(i = 0; i < N_Elements; i++)
  {
    N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
    dat << n_local_vertices_to_type<d>(N_CellVertices) << " ";
  }
  dat << "  </DataArray>\n";
  dat << "</Cells>\n";
  dat << "\n";
  dat << "<PointData Vectors=\"Velocity\" Scalars=\"Scalars\">\n";
  dat << "\n";

  // write vector variables into file
  if(N_LocVertices)
    for(k = 0; k < N_VectorVar; k++)
    {
      fespace = FEVectFunctArray[k]->GetFESpace3D();
      N_Comp = FEVectFunctArray[k]->GetN_Components();
      Length = FEVectFunctArray[k]->GetLength();
      Coeffs = FEVectFunctArray[k]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();

      memset(DoubleArray, 0, SizeOfDouble * N_Vertices * N_Comp);
      memset(WArray, 0, SizeOfDouble * N_Vertices);
      m = 0;

      for(i = 0; i < N_Elements; i++)
      {
        cell = Coll->GetCell(i);
        N_ = cell->GetN_Vertices();

        // find FE data for this element
        bf = fespace->get_fe(i).GetBaseFunct3D();
        DOF = GlobalNumbers + BeginIndex[i];
        N_LocDOF = bf->GetDimension();
        for(j = 0; j < N_; j++)
        {
          switch(cell->GetN_Vertices())
          {
            case 4:
              xi = TetraCoords[3 * j];
              eta = TetraCoords[3 * j + 1];
              zeta = TetraCoords[3 * j + 2];
              break;

            case 8:
              xi = HexaCoords[3 * j];
              eta = HexaCoords[3 * j + 1];
              zeta = HexaCoords[3 * j + 2];
              break;
          }
          bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
          bf->ChangeBF(Coll, cell, BFValues);

          for(n = 0; n < N_Comp; n++)
          {
            value = 0;
            for(l = 0; l < N_LocDOF; l++)
              value += BFValues[l] * Coeffs[DOF[l] + n * Length];
            DoubleArray[N_Comp * VertexNumbers[m] + n] += value;
          }
          WArray[VertexNumbers[m]] += 1.;
          m++;
        } // endfor j
      }   // endfor i

      // midle value
      l = 0;
      for(i = 0; i < N_Vertices; i++)
      {
        for(j = 0; j < N_Comp; j++)
        {
          if(WArray[i] > 1.)
            DoubleArray[l] /= WArray[i];
          l++;
        }
      } // endfor l

      //      for(i=0;i<2*N_Vertices;i++)
      //       cout << "Do[" << i << "]" << DoubleArray[i] << "\n";

      dat << "  <DataArray type=\"Float32\" Name=\""
          << FEVectFunctArray[k]->GetName() << "\" NumberOfComponents=\""
          << "3\" format=\"ascii\">\n";
      l = 0;
      l = 0;
      for(i = 0; i < N_Vertices; i++)
      {
        for(j = 0; j < N_Comp; j++)
          dat << DoubleArray[N_Comp * i + j] << " ";
      }
      dat << "\n";

      dat << "  </DataArray>\n";
      dat << "\n";

      for(j = 0; j < N_Comp; j++)
      {
        dat << "  <DataArray type=\"Float32\" Name=\""
            << FEVectFunctArray[k]->GetName() << j << "\" NumberOfComponents=\""
            << "1\" format=\"ascii\">\n";
        for(i = 0; i < N_Vertices; i++)
          dat << DoubleArray[i * N_Comp + j] << " ";
        dat << "\n";
        dat << "  </DataArray>\n";
        dat << "\n";
      }
    } // for(k=0;k<N_Vec

  // write scalar variables into file
  if(N_LocVertices)
    for(k = 0; k < N_ScalarVar; k++)
    {
      fespace = FEFunctionArray[k]->GetFESpace3D();
      Coeffs = FEFunctionArray[k]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();

      memset(DoubleArray, 0, SizeOfDouble * N_Vertices);
      memset(WArray, 0, SizeOfDouble * N_Vertices);
      m = 0;

      for(i = 0; i < N_Elements; i++)
      {
        cell = Coll->GetCell(i);
        N_ = cell->GetN_Vertices();

        // find FE data for this element
        bf = fespace->get_fe(i).GetBaseFunct3D();
        DOF = GlobalNumbers + BeginIndex[i];
        N_LocDOF = bf->GetDimension();
        for(j = 0; j < N_; j++)
        {
          switch(cell->GetN_Vertices())
          {
            case 4: // Tetrahedron
              xi = TetraCoords[3 * j];
              eta = TetraCoords[3 * j + 1];
              zeta = TetraCoords[3 * j + 2];
              break;
            case 8: // Hexahedron
              xi = HexaCoords[3 * j];
              eta = HexaCoords[3 * j + 1];
              zeta = HexaCoords[3 * j + 2];
              break;
          }
          bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
          bf->ChangeBF(Coll, cell, BFValues);
          value = 0;
          for(l = 0; l < N_LocDOF; l++)
            value += BFValues[l] * Coeffs[DOF[l]];
          DoubleArray[VertexNumbers[m]] += value;
          WArray[VertexNumbers[m]] += 1.;
          m++;
        } // endfor j
      }   // endfor i

      // non conforming
      for(i = 0; i < N_Vertices; i++)
      {
        if(WArray[i] > 1.)
          DoubleArray[i] /= WArray[i];
      }

      dat << "  <DataArray type=\"Float32\" Name=\""
          << FEFunctionArray[k]->GetName() << "\" NumberOfComponents=\""
          << "1\" format=\"ascii\">\n";
      for(j = 0; j < N_Vertices; j++)
        dat << DoubleArray[j] << " ";
      dat << "\n";
      dat << "  </DataArray>\n";
      dat << "\n";
    } // for(k=0;k<N_ScalarVar

  //     dat <<  "</PointData>\n";
  //     dat << "\n";
  //
  //     dat <<  "<CellData Scalars=\"SubDomain\">\n";
  //     dat <<  "  <DataArray type=\"Int32\" Name=\"Region\""
  //         <<" NumberOfComponents=\"1\" format=\"ascii\">\n";
  //     for(i=0;i<N_Elements;i++)
  //       dat << (Coll->GetCell(i))->GetRegionID()   << " ";
  //     dat <<  "  </DataArray>\n";
  //     dat <<  "</CellData>\n";
  //
  //

  dat << "</PointData>\n";
  dat << "\n";

  dat << "<CellData Scalars=\"SubDomain\">\n";
  dat << "  <DataArray type=\"Int32\" Name=\"SubDomain\""
      << " NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(i = 0; i < N_Elements; i++)
    dat << (Coll->GetCell(i))->GetSubDomainNo() << " ";
  dat << "  </DataArray>\n";
  dat << "  <DataArray type=\"Int32\" Name=\"RegionID\""
      << " NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(i = 0; i < N_Elements; i++)
    dat << (Coll->GetCell(i))->GetRegionID() << " ";
  dat << "  </DataArray>\n";

  dat << "</CellData>\n";

  dat << "</Piece>\n";
  dat << "</UnstructuredGrid>\n";
  dat << "</VTKFile>\n";

  dat.close();

  if(N_LocVertices)
  {
    delete[] NumberVertex;
    delete[] VertexNumbers;
    delete[] DoubleArray;
    delete[] WArray;
    delete[] Coords;
  }
} // DataWriter::Write_ParVTK
#endif // _MPI

/** @brief Output in .case format

   1 file .case (main file)
   1 file .geo (mesh)
   scalar and vector output files at each time step
   Note that currently we only write the P1 projection (node values) of the
   results
   However, the case format allows also element-wise and higher order
   visualizations (but not with paraview)
   @todo implement visualization for different element type
*/
template <int d>
void DataWriter<d>::writeCaseFile()
{
  std::string casename = testcaseDir + "/" + testcaseName + ".case";
  Output::print<2>(" ** DataWriter::writeCaseFile - write ", casename);
  std::ofstream casf(casename);

  std::string geoname = testcaseName + ".geo";

  casf << "FORMAT\n";
  casf << "type: ensight\n";
  casf << "GEOMETRY\n";
  casf << "model: 1 " << geoname << "\n";
  casf << "VARIABLE\n";
  for(unsigned int j = 0; j < FEFunctionArray.size(); j++)
  {
    unsigned int dimension = FEFunctionArray[j]->GetFESpace()->GetBaseVectDim();
    if(dimension == 1)
    {
      casf << "scalar per node: 1 " << FEFunctionArray[j]->GetName() << " "
           << testcaseName << "_" << FEFunctionArray[j]->GetName()
           << ".****.scl\n";
    }
    else
    {
      casf << "vector per node: 1 " << FEFunctionArray[j]->GetName() << " "
           << testcaseName << "_" << FEFunctionArray[j]->GetName()
           << ".****.vct\n";
    }
  }
  for(unsigned int j = 0; j < FEVectFunctArray.size(); j++)
  {
    casf << "vector per node: 1 " << FEVectFunctArray[j]->GetName() << " "
         << testcaseName << "_" << FEVectFunctArray[j]->GetName()
         << ".****.vct\n";
  }
  // time values
  unsigned int n_time_steps = timeValues.size();
  casf << "TIME\n";
  casf << "time set: 1\n";
  casf << "number of steps: " << std::max(1U, n_time_steps) << "\n";
  casf << "filename start number: 0\n";
  casf << "filename increment: " << n_steps_per_output << "\n";
  casf << "time values:\n";
  for(unsigned int i = 0; i < n_time_steps; i++)
  {
    casf.precision(5);
    casf.width(12);
    casf << timeValues[i];
    if((i && !(i % 5)) || (i == n_time_steps))
      casf << "\n";
    else
      casf << " ";
  }
  if(n_time_steps == 0) // this is the case for stationary problems
    casf << 0.0;
  casf << "\n";
  casf.close();
}

// Geometry
template <int d>
void DataWriter<d>::writeCaseGeo()
{
  // write geo file
  std::string geoname = testcaseDir + "/" + testcaseName + ".geo";
  std::ofstream f(geoname);

  f << testcaseName << " geometry\n";
  f << "t = 0.0000"
    << "\n";
  f << "node id assign"
    << "\n";
  f << "element id assign"
    << "\n";
  f << "coordinates"
    << "\n";
  f.width(8);

  // write points
  f << Coll->GetN_Vertices() << "\n";

  writeCoord(f);

  // write elements
  int nVE = Coll->GetCell(0)->GetN_Vertices();
  std::string ensight_type;

  /// @warning this works now only for a single (bulk) domain
  unsigned int nParts = 1;        // number of (bulk) subdomains
  std::string partName = "inner"; // name of collection subdomain

  unsigned int VERTEX_OFFSET = 1; // index change for the output

  for(unsigned int ifig = 1; ifig <= nParts; ifig++)
  {
    f << "part";
    f.width(8);
    f << ifig << "\n";
    f << partName << "\n"; // name of the figure
    switch(d)
    {
      case 2:
        switch(nVE)
        {
          case 3:
            ensight_type = "tria3";
            break;
          case 4:
            ensight_type = "quad4";
            break;
        }
        break;
      case 3:
        switch(nVE)
        {
          case 4:
            ensight_type = "tetra4";
            break;
          case 8:
            ensight_type = "hexa8";
            break;
        }
        break;
      default:
        ErrThrow("unnkown number of vertices in ", d, "D: ", nVE);
    }
    f << ensight_type << "\n";

    f.width(8);
    f << Coll->GetN_Cells() << "\n";
    for(int k = 0; k < Coll->GetN_Cells(); k++)
    {
      // @warning we cannot write yet .case output with mixed meshes.
      if(Coll->GetCell(k)->GetN_Vertices() != nVE)
      {
        ErrThrow(" **ERROR: CASE output with mixed (tets and hexas) meshes not "
                 "supported yet. Use VTK instead.");
      }

      for(int i = 0; i < nVE; i++)
      {
        f.width(8);
        f << Coll->GetGlobalVerNo(k, i) + VERTEX_OFFSET;
      }
      f << "\n";
    } // k
  }   // ifig

  // write boundary elements
  unsigned int nBdParts = 1; // number of subdomains
  if(d == 3)                 // this does not work in 3D
    nBdParts = 0;
  ///@warning this works now only for a single boundary domain
  partName = "boundary"; // name of collection subdomain
  int nVertexPerFace = 2;
  for(unsigned int ifig = 1; ifig <= nBdParts; ifig++)
  {
    f << "part";
    f.width(8);
    f << nParts + ifig << "\n";
    f << partName << "\n"; // name of the figure
    ///@warning only for P1 in 1D
    ensight_type = "bar2";
    f << ensight_type << "\n";

    f.width(8);
    f << Coll->GetN_BdFaces() << "\n";
    for(unsigned int k = 0; k < Coll->GetN_BdFaces(); k++)
    {
      for(int i = 0; i < nVertexPerFace; i++)
      {
        f.width(8);
        f << Coll->GetBdFacesNode(nVertexPerFace * k + i);
      }
      f << "\n";
    } // k
  }   // ifig
}

template <int d>
void DataWriter<d>::writeCaseVars(int iter)
{
  std::string filename = testcaseDir + "/" + testcaseName;
  char numstr[4];
  sprintf(numstr, "%d", iter);
  std::string number = "000" + std::string(numstr);
  if(iter > 9)
    number = "00" + std::string(numstr);
  if(iter > 99)
    number = "0" + std::string(numstr);
  if(iter > 999)
    number = std::string(numstr);

  std::string ensight_type;

  unsigned int N_Vertices = Coll->GetN_Vertices();
  unsigned int dimension;

  // scalars
  for(unsigned int i = 0; i < FEFunctionArray.size(); i++)
  {
    std::vector<double> solutionAtNode;
    computeNodeValues(FEFunctionArray[i], solutionAtNode, dimension);

    // write scalars
    std::ofstream sclf;
    std::ostringstream sclname;

    if(dimension == 1)
    {
      sclname << filename << "_" << FEFunctionArray[i]->GetName() << "."
              << number << ".scl";
    }
    else
    {
      sclname << filename << "_" << FEFunctionArray[i]->GetName() << "."
              << number << ".vct";
    }
    std::string fname = sclname.str();

    Output::print<3>(" ** DataWriter::write File - write ", fname);
    sclf.open(fname);
    sclf << FEFunctionArray[i]->GetName() << " step = " << iter << "\n";

    writeVectCase(sclf, N_Vertices, dimension, solutionAtNode);
  }

  // vectors
  for(unsigned int i = 0; i < FEVectFunctArray.size(); i++)
  {
    // int nComp = FEVectFunctArray[i]->GetN_Components();
    std::vector<double> solutionAtNode;
    computeNodeValues(FEVectFunctArray[i], solutionAtNode, dimension);

    std::ofstream vctf;
    std::ostringstream vctname;
    vctname << filename << "_" << FEVectFunctArray[i]->GetName() << "."
            << number << ".vct";
    std::string fname = vctname.str();
    Output::print<3>(" ** DataWriter::write File - write ", fname);
    vctf.open(fname.c_str());
    vctf << FEVectFunctArray[i]->GetName() << " step = " << iter << "\n";

    writeVectCase(vctf, N_Vertices, dimension, solutionAtNode);
  }
}

template <int d>
void DataWriter<d>::writeCoord(std::ofstream& f)
{
  f.precision(8);   
  unsigned int N_Vertices = Coll->GetN_Vertices();
  for(unsigned int i = 0; i < N_Vertices; i++)
  {
    ///@attention case output works only with this format
    f << Coll->GetCoord(i * d);
    f << " ";
     f << Coll->GetCoord(i * d + 1);
    f << " ";
    if(d == 2)
    {
      f << 0.;
    }
    else
    {
      f << Coll->GetCoord(i * d + 2);
    }
    f << "\n";
  }
  f << "\n";
}

template <int d>
template <class T>
void DataWriter<d>::computeNodeValues(const T* function,
                                      std::vector<double>& solutionAtNode,
                                      unsigned int& dimension)
{
  auto FESpace = function->GetFESpace();
  TCollection* coll = FESpace->GetCollection();
  int nPoints = coll->GetN_Vertices();

  // get function type of T
  bool IsVect = std::is_same<T, TFEVectFunct3D>::value
                || std::is_same<T, TFEVectFunct2D>::value;
  dimension = 1;
  if(IsVect || FESpace->GetBaseVectDim() > 1)
  {
    dimension = d;
  }

  solutionAtNode.assign(dimension * nPoints, 0.0);
  std::vector<int> WArray(nPoints, 0);

  for(int ncell = 0; ncell < coll->GetN_Cells(); ncell++)
  {
    const TBaseCell* cell = coll->GetCell(ncell);
    unsigned int nLocalVertices = cell->GetN_Vertices();

    for(unsigned int nvert = 0; nvert < nLocalVertices; nvert++)
    {
      unsigned int globalVert_index = coll->GetGlobalVerNo(ncell, nvert);

      double function_value[dimension];
      double x, y, z;
      cell->GetVertex(nvert)->GetCoords(x, y, z);
#ifdef __2D__
      function->FindValueLocal(cell, ncell, x, y, function_value);
#else
      function->FindValueLocal(cell, ncell, x, y, z, function_value);
#endif
      for(unsigned int dim = 0; dim < dimension; dim++)
      {
        solutionAtNode[dimension * globalVert_index + dim]
            += function_value[dim];
      }
      WArray.at(globalVert_index) += 1;
    } // endfor nvert
  }   // endfor ncell

  for(int nvert = 0; nvert < nPoints; nvert++)
  {
    if(WArray[nvert] != 0)
    {
      for(unsigned int dim = 0; dim < dimension; dim++)
      {
        solutionAtNode[dimension * nvert + dim] /= WArray[nvert];
      }
    }
  }
}

// ****************
// scalar components of vector variables
// ****************
template <int d>
void DataWriter<d>::printVectCompwise(std::ofstream& dat, std::string name,
                                      unsigned int N_Vertices,
                                      unsigned int N_Comp,
                                      std::vector<double> solutionAtNode)
{
  for(unsigned int ncomp = 0; ncomp < N_Comp; ncomp++)
  {
    dat << "SCALARS " << name;
    if(N_Comp != 1)
      dat << ncomp;
    dat << " double\n";
    dat << "LOOKUP_TABLE default\n";
    for(unsigned int nvert = 0; nvert < N_Vertices; nvert++)
    {
      double value = solutionAtNode[nvert * N_Comp + ncomp];
      dat << std::scientific << value << "\n";
    }
    dat << "\n\n";
  }
}

// ***************
// absolute value of vector variables
// ***************
template <int d>
void DataWriter<d>::printVectAbsValue(std::ofstream& dat, std::string name,
                                      unsigned int N_Vertices,
                                      unsigned int N_Comp,
                                      std::vector<double> solutionAtNode)
{
  dat << "SCALARS |" << name << "|";
  dat << " double\n";
  dat << "LOOKUP_TABLE default\n";
  for(unsigned int nvert = 0, l = 0; nvert < N_Vertices; nvert++)
  {
    double t = 0;
    for(unsigned int ncomp = 0; ncomp < N_Comp; ncomp++)
    {
      t += solutionAtNode[l] * solutionAtNode[l];
      l++;
    }
    dat << std::sqrt(t) << "\n";
  }
  dat << "\n\n";
}

// ***************
// VECTORS
// ***************
template <int d>
void DataWriter<d>::printVectPointwise(std::ofstream& dat, std::string name,
                                       unsigned int N_Vertices,
                                       unsigned int N_Comp,
                                       std::vector<double> solutionAtNode)
{
  dat << "VECTORS " << name;
  dat << " double\n";
  for(unsigned int nvert = 0; nvert < N_Vertices; nvert++)
  {
    for(unsigned int ncomp = 0; ncomp < N_Comp; ncomp++)
    {
      double value = solutionAtNode[nvert * N_Comp + ncomp];
      dat << std::scientific << value << " ";
    }
    if(N_Comp == 2)
      dat << double(0);
    dat << "\n";
  }
  dat << "\n";
}

// auxiliary function for writeVectCase
void printEntry(std::ofstream& dat, double value, int counter)
{
  dat.setf(std::ios_base::scientific);
  dat.precision(5);
  dat.width(12);
  dat << value;

  // One line contains at most six entries.
  if((counter + 1) % 6 == 0)
    dat << "\n";
}

template <int d>
void DataWriter<d>::writeVectCase(std::ofstream& dat, unsigned int N_Vertices,
                                  unsigned int N_Comp,
                                  std::vector<double> solutionAtNode)
{
  for(unsigned int nvert = 0; nvert < N_Vertices; nvert++)
  {
    unsigned int entry = N_Comp * nvert;
    double value_x = solutionAtNode.at(entry);
    if(N_Comp == 1)
      printEntry(dat, value_x, entry);
    else
    {
      int counter = 3 * nvert;
      printEntry(dat, value_x, counter);

      double value_y = solutionAtNode.at(entry + 1);
      printEntry(dat, value_y, counter + 1);

      double value_z = 0.;
      if(N_Comp > 2)
        value_z = solutionAtNode.at(entry + 2);
      printEntry(dat, value_z, counter + 2);
    }
  }
}

#ifdef __3D__
template class DataWriter<3>;
#else
template class DataWriter<2>;
#endif
