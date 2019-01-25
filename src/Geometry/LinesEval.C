/******************************************************************************* 
*
* @class LineEval
* @date  19.12.18
* @brief Some tools for function evaluation over line
* @author Baptiste Moreau
* @History:
*
*******************************************************************************/

#include <LinesEval.h>

#include <algorithm>
#include <Domain.h>
#include <BaseCell.h>
#include <FEDatabase3D.h>
#ifdef __2D__
#include "FEFunction2D.h"
#else
#include "FEFunction3D.h"
#endif
#ifdef _MPI
  #include <numeric> // std::accumulate
#endif // _MPI
#include <limits>
#include <sys/stat.h>

// see https://stackoverflow.com/a/8016853
template <int d>
constexpr char LinesEval<d>::required_database_name[];

/* ************************************************************************** */
bool operator < (const LineEval::point_on_line& P,
                 const LineEval::point_on_line& Q)
{
  if( P.lmin_cell < Q.lmin_cell)
  {
    return true;
  }
  else
  {
    return false;
  }
}

/* ************************************************************************** */
bool operator == (const LineEval::point_on_line& P,
                  const LineEval::point_on_line& Q)
{
  if( fabs(P.lmin_cell - Q.lmin_cell) < 1e-16
  &&  fabs(P.lmax_cell - Q.lmax_cell) < 1e-16 )
  {
    return true;
  }
  else
  {
    return false;
  }
}

/* ************************************************************************** */
template <int d>
ParameterDatabase LinesEval<d>::default_lineseval_parameters()
{
  ParameterDatabase db(LinesEval<d>::required_database_name);

  db.add("line_position", std::vector<double>(1, std::nan("")),
         "List of points coordinates throuh which the lines to be defined pass."
         " In 3D the number of coordinates must be divisible by three.");

  db.add("line_direction", 0, "Line direction, only one direction possible for "
         "all lines (possible direction are paralell to axis x=0 y=1 z=2).",
         {0, 1, 2});

  db.add("line_refinement", 0, "Line refinement.\n0 no refinement\n"
         "i additional point: between two points for line defined by file or\n"
         "                    in cells for line defined by nested database",
         0, 1000);

  db.add("position_file", "", "Name of the file containing the points "
         "coordinates for the line evaluation.");
  
  db.add("directory_name", "line_evaluations", "Name of the directory where "
         "the files of LinesEval::write_fe_values are written to.");

  return db;
}

/* ************************************************************************** */
template <int d>
LinesEval<d>::LinesEval():db(default_lineseval_parameters())
{  
}

/* ************************************************************************** */
template <int d>
LinesEval<d>::LinesEval(const TDomain&           domain,
                     const ParameterDatabase& param_db):
                     db(default_lineseval_parameters())
{
  // Try to find a nested database with the correct name
  try
  {
    this->db.merge(param_db.get_nested_database(LinesEval
                                                ::required_database_name));
  }
  catch(...)
  {
    if(param_db.get_name() == std::string(LinesEval::required_database_name))
    {
      this->db.merge(param_db, false);
    }
    else
    {
      ErrThrow("LineEval: ", "If the values should be evaluated on a line, "
               "you must give a (nested) database named 'LineEval Database'.");
    }
  }

  int refi = this->db["line_refinement"];

  // create a directory for the data to be saved
  bool i_am_root = true;
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == 0);
#endif
  if( i_am_root )
  {
    std::string directory_name = db["directory_name"];
    mkdir(directory_name.c_str(), 0777);
  }
  
  std::string filename = this->db["position_file"];

  // define lines from file
  if( ! filename.empty() )
  {
    std::vector<int>    n_coord; 
    std::vector<int>    direction;
    std::vector<double> position;
    std::vector<double> coord;

    int offset = 0;

    Output::info<3>("LineEval:", "Coordinates for line evaluation are read "
                "from the file ", filename);

    read_position(filename, direction, position, coord, n_coord);
    check_position(filename, direction, position, n_coord);

    // loop over all lines to be defined
    for( unsigned int i=0 ; i< n_coord.size() ; i++ )
    {
      double P[3];
      std::copy(position.begin()+3*i,
                position.begin()+(3*i+3),
                P);

      std::vector<double> coordinates;
      std::copy(coord.begin()+offset,
                coord.begin()+(offset+n_coord[i]),
                back_inserter(coordinates));

      LineEval line(domain, direction[i], P, coordinates, refi);    
      lines_for_postprocess.emplace_back(line);

      offset += n_coord[i];
    }
  }

  // define lines from nested database
  std::vector<double> line_position = LinesEval::db["line_position"];

  // skip definition of lines from nested database if no positions are given
  if( (line_position.size()==1) && (std::isnan(line_position.at(0))) )
  {
    return;
  }
  // Check parameters compatibility and define N_line
  if( (line_position.size())%3 != 0 )
  {
    ErrThrow("LineEval: ", "Some coordinates in parameter: line_position "
                  "are missing.");
  }

  double P[3];
  int N_line = line_position.size()/3;

  // Loop over all lines to be defined
  for( int i=0 ; i<N_line ; i++ )
  {
    std::copy(line_position.begin()+3*i,
              line_position.begin()+(3*i+3),
              P);

    LineEval line(domain, LinesEval::db["line_direction"], P, refi);

    lines_for_postprocess.emplace_back(line);
  }
}

/* ************************************************************************** */
template <int d>
void LinesEval<d>::read_position(const std::string    filename,
                              std::vector<int>&    direction,
                              std::vector<double>& position,
                              std::vector<double>& coord,
                              std::vector<int>&    n_coord)
{
  Output::info<3>("LineEval:",
                  "Reading points coordinates from file ", filename);

  std::string line;
  std::ifstream file(filename.c_str());
  if( !file )
  {
    ErrThrow("LineEval: File ", filename, " could not be opened.");
  }

  int line_i = 1;

  while( std::getline(file, line) )
  {
    if(line_i == 1)
    {
      int dir;
      std::istringstream iss(line);
      while (iss >> dir)
      {
        direction.push_back(dir);
      }
      line_i++;
      continue;
    }
    else if( line_i == 2 )
    {
      double p;
      std::istringstream iss(line);
      while (iss >> p)
      {  
        position.push_back(p);
      }
      line_i++;
      continue;
    }
    else
    {
      double c;
      int n = 0;
      std::istringstream iss(line);
      while (iss >> c)
      {  
        coord.push_back(c);
        n++;
      }
      n_coord.push_back(n);
    }
  }
  file.close();
}

/* ************************************************************************** */
template <int d>
void LinesEval<d>::check_position(const std::string         filename,
                               const std::vector<int>    direction,
                               const std::vector<double> position,
                               std::vector<int>&         n_coord)
{
  unsigned int N_dir = direction.size();
  unsigned int N_pos = position.size()/3;
  unsigned int N_coo = n_coord.size();
  
#ifdef __3D__  
  if( (position.size())%3 != 0 )
  {
    ErrThrow("LineEval :", "Some position coordinates in file: ", filename,
             " are missing.");
  }
#endif

  if( N_dir != N_pos )
  {
    ErrThrow("LineEval :", "Numbers of directions and position coordinates "
             "in file: ", filename, " are mismaching.");
  }

  if( N_dir != N_coo )
  {
    ErrThrow("LineEval :", "Numbers of directions and lines in file: ",
              filename, " are mismaching.");
    n_coord.resize(std::min(N_dir, N_coo));
  }

  if( N_pos != N_coo )
  {
    ErrThrow("LineEval :", "Numbers of position coordinates and lines "
             "in file: ", filename, " are mismaching.");
    N_coo = n_coord.size();
    n_coord.resize(std::min(N_pos, N_coo));
  }
}

/* ************************************************************************** */
template <int d>
void LinesEval<d>::write_fe_values(std::vector<const FEFunction*> fe_functions,
                                   double      time_step,
                                   std::string file_prefix) const
{
  bool i_am_root = true;

#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == 0);
#endif

  double default_val    = -1.e32;
  auto   n_fe_functions = fe_functions.size();

  if( n_fe_functions==0 )
  {
    return;
  }

  // loop over all lines
  for(int l_i=0 ; l_i < this->GetLength() ; l_i++)
  {
    int j_cell;

    auto line      = this->GetLine(l_i);
    int  size_li   = line.GetNbPoints();
    int  direction = line.GetDirection();

    std::vector<std::vector<double>> values(
      n_fe_functions+d, std::vector<double>(size_li, default_val));

    std::array<double, 3> P = line.GetBasePoint();

    // loop over all points of line l_i to get the values
    for( int j=0 ; j<size_li ; j++ )
    {
      j_cell = line.GetCellIdx(j);
      P[direction] = line.GetPosition(j);

      for(int i = 0; i < d; ++i)
      {
        values[i][j] = P[i];
      }
      
      if( j_cell==-1 ) //only append in case of many processes
      {
        continue;
      }

      const TBaseCell* cell = line.GetCell(j_cell);

      for(int i = 0; i < n_fe_functions; ++i)
      {
#ifdef __3D__
        fe_functions[i]->FindValueLocal(cell, j_cell, P[0], P[1], P[2],
                                        &values[i+d][j]);
#else
        fe_functions[i]->FindValueLocal(cell, j_cell, P[0], P[1], 
                                        &values[i+d][j]);
#endif
      }
    }

#ifdef _MPI
    //communicate the values to root
    double* rbuf;

    if( i_am_root )
    {
      for(int i = 0; i < n_fe_functions; ++i)
      {
        MPI_Reduce(MPI_IN_PLACE, &values[i+d][0], size_li, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);
      }
    }
    else
    {
      for(int i = 0; i < n_fe_functions; ++i)
      {
        MPI_Reduce(&values[i+d][0], &rbuf, size_li, MPI_DOUBLE, MPI_MAX, 0,
                   MPI_COMM_WORLD);
      }
    }
#endif /** #ifdef _MPI */

    if( i_am_root )
    {
      std::stringstream s_out;

      bool        flag_print  = false;
      std::string dir_name    = db["directory_name"];

      std::ostringstream l_i_format;
      l_i_format << std::setfill('0') << std::setw(2) << std::to_string(l_i);

      std::string file_name   = dir_name + "/" + file_prefix + "_l"
                              + l_i_format.str()
                              + "_t" + std::to_string(time_step) + ".txt";

#ifdef _MPI
      s_out << "#" << setw(10) << "x0" << setw(15) << "x1" << setw(15) << "x2";
#else
      s_out << "#" << setw(13) << "x0" << setw(15) << "x1" << setw(15) << "x2";
#endif

      for(int i = 0; i < n_fe_functions; ++i)
      {
        s_out << setw(15) << fe_functions[i]->GetName();
      }
      s_out << "\n";
      
      for( int j=0 ; j < size_li ; j++ )
      {
        if(values[d][j] != default_val)
        {
          flag_print = true;
          for(int i = 0; i < d + n_fe_functions; ++i)
          {
            s_out << setprecision(7) << setw(14) << values[i][j] << " ";
          }
          s_out << "\n";
        }
      }

      if( flag_print )
      {
        Output::redirect(file_name);
        Output::print<1>(s_out.str());
        Output::resetOutfile();
      }
    }
  }
}

/* ************************************************************************** */
LineEval::LineEval(const TDomain&      domain,
                   const int           direction,
                   const double        P[3],
                   std::vector<double> coord,
                   const int           refi) :
                   coll(domain.get_grid_collections().front()),
                   direction(direction),
                   base_point{P[0], P[1],P[2]},
                   n_refine(refi),
                   fromDB(false)
{  
  int        my_rank;
  TBaseCell* cell;

  int n_coord = coord.size();
  int N_Cells = coll->GetN_Cells();

#ifdef _MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
#else
  my_rank = 0;
#endif

  if( my_rank == 0 )
  {
    Output::info<3>("LineEval:",
                    "Constructing a list of cells for line post-processing.");
  }

  // add points to coord if refine > 0
  if( n_refine>0 )
  {
    std::vector<double> coordinate;

    for( int i_pnt=0 ; i_pnt<n_coord-1 ; i_pnt++ )
    {
      double p0 = coord[i_pnt];
      double p1 = coord[i_pnt+1];

      for( int i=0 ; i<n_refine+1 ; i++ )
      {
        double pi = p0 + i*(p1-p0)/(n_refine+1);
        coordinate.push_back(pi);
      }
    }
    // add the last point
    coordinate.push_back(coord[n_coord-1]);

    coord.clear();
    std::copy(coordinate.begin(),
              coordinate.end(),
              back_inserter(coord));
    n_coord = coord.size();
  }

  // add points to line
  for( int i_pnt=0 ; i_pnt<n_coord ; i_pnt++ )
  {
    int j_cell = -1;
    std::array<double,3> pnt_eval = this->base_point;
    pnt_eval[direction] = coord[i_pnt];

    // loop to get the extremities of the line
    for( int i_cell=0 ; i_cell<N_Cells ; i_cell++ )
    {
      cell = coll->GetCell(i_cell);

#ifdef _MPI
      //only perform the following calculations on OwnCells
      if( cell->IsHaloCell() )
      {
        continue;
      }
#endif

      if( cell->PointInCell(pnt_eval[0],pnt_eval[1],pnt_eval[2]) )
      {
        j_cell = i_cell;
        break;
      }
    }
    // to keep the point ordering in case of many processes
    line_for_postprocess.emplace_back(j_cell, coord[i_pnt], 0.);
  }
}

/* ************************************************************************** */
LineEval::LineEval(const TDomain& domain,
                   const int      direc,
                   const double   P[3],
                   const int      refi) :
                   coll(domain.get_grid_collections().front()),
                   direction(direc),
                   base_point{P[0], P[1],P[2]},
                   n_refine(refi),
                   fromDB(true)
{  
  int        my_rank;
  double     lmin_cell;
  double     lmax_cell;
  TBaseCell* cell;

  int N_Cells = coll->GetN_Cells();

#ifdef _MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
#else
  my_rank = 0;
#endif

  if( my_rank == 0 )
  {
    Output::info<3>("LineEval:",
                    "Constructing a list of cells for line post-processing.");
  }

  // loop over all cells to get the intersected cells
  for( int i_cell=0 ; i_cell<N_Cells ; i_cell++ )
  {
    cell = coll->GetCell(i_cell);

    // Error message if the cell is not a thetrahedron
    if( cell->GetType() != Tetrahedron )
    {
      ErrThrow("LineEval:", "Only for tetrahedron.");
    }

#ifdef _MPI
    //only perform the following calculations on OwnCells
    if( cell->IsHaloCell() )
    {
      continue;
    }
#endif

    if( cell->IsLineCutingCell(direction, base_point, lmin_cell, lmax_cell) )
    {
      // only keep cells with more than one intersection point
      if( (lmax_cell - lmin_cell) > 1e-16)
      {
        line_for_postprocess.emplace_back(i_cell, lmin_cell, lmax_cell);
      }
    }
  }

#ifdef _MPI
  // add points from other processes in order to keep the point ordering
  int mpi_size;
  
  int totlength = 0;
  int size_l    = line_for_postprocess.size();

  std::vector<double> slmin(size_l,0.);
  std::vector<double> slmax(size_l,0.);

  for( int i=0 ; i<size_l ; i++ )
  {
    slmin[i] = line_for_postprocess.at(i).lmin_cell;
    slmax[i] = line_for_postprocess.at(i).lmax_cell;
  }
  
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
  std::vector<int> displs(mpi_size, 0);
  std::vector<int> rlengths(mpi_size, 0);

  // gather lengths of line on every process
  MPI_Allgather(&size_l, 1, MPI_INT, &rlengths[0], 1, MPI_INT, MPI_COMM_WORLD);

  for( int i_mpi=0 ; i_mpi<mpi_size ; i_mpi++ )
  {
    displs[i_mpi] = totlength;
    totlength    += rlengths[i_mpi];
  }

  std::vector<double> rlmin(totlength, 0.);
  std::vector<double> rlmax(totlength, 0.);

  MPI_Allgatherv(&slmin[0], size_l, MPI_DOUBLE,
                 &rlmin[0], &rlengths[0], &displs[0], MPI_DOUBLE,
                 MPI_COMM_WORLD);
  MPI_Allgatherv(&slmax[0], size_l, MPI_DOUBLE,
                 &rlmax[0], &rlengths[0], &displs[0], MPI_DOUBLE,
                 MPI_COMM_WORLD);

  // only add new points by skipping points already in process according to
  // ordering from MPI_Allgatherv (i.e. according to rank order)
  int start_pos = std::accumulate(rlengths.begin(),
                                  rlengths.begin()+my_rank,
                                  0);

  for( int i=0 ; i<totlength ; i++ )
  {
    if( (i>=start_pos) && (i<start_pos+rlengths[my_rank]) )
    {
      continue;
    }
    line_for_postprocess.emplace_back(-1, rlmin[i], rlmax[i]);
  }
#endif /** #ifdef _MPI */

  // sort points by lmin_cell in ascending order.
  std::sort(line_for_postprocess.begin(),line_for_postprocess.end());
  
  // and keep only one cell if line on joint of many cells
  auto it = std::unique(line_for_postprocess.begin(),
                        line_for_postprocess.end());
  line_for_postprocess.erase(it, line_for_postprocess.end());
}

/* ************************************************************************** */
std::array<double, 3> LineEval::GetBasePoint() const
{
  return base_point;
}

/* ************************************************************************** */
int LineEval::GetDirection() const
{
  return direction;
}

/* ************************************************************************** */
int LineEval::GetNbPoints() const
{
  int size = this->line_for_postprocess.size();

  if( this->fromDB )
  {
    return (size + this->n_refine*size + 1);
  }
  else
  {
    return size;
  }
}

/* ************************************************************************** */
int LineEval::GetCellIdx(int i) const
{
  if( this->fromDB )
  {
    int j = (i==this->GetNbPoints()-1) ? this->line_for_postprocess.size()-1
                                       : i/(this->n_refine + 1);
    return this->line_for_postprocess.at(j).cell_index;
  }
  else
  {
    return this->line_for_postprocess.at(i).cell_index;
  }
}

/* ************************************************************************** */
const TBaseCell* LineEval::GetCell(int cell_idx) const
{
  if( cell_idx>=0 )
  {
    return this->coll->GetCell(cell_idx);
  }
  else
  {
    return nullptr; 
  }
}

/* ************************************************************************** */
double LineEval::GetPosition(int i) const
{
  if( this->fromDB )
  {
    if( i==this->GetNbPoints()-1 )
    {
      return this->line_for_postprocess.at(this->line_for_postprocess.size()-1)
                                       .lmax_cell;
    }
    else
    {
      int j = i/(this->n_refine + 1);
      int k = i%(this->n_refine + 1);
      double p0 = this->line_for_postprocess.at(j).lmin_cell;
      double p1 = this->line_for_postprocess.at(j).lmax_cell;
      return (p0 + k*(p1-p0)/(this->n_refine+1));
    }
  }
  else
  {
    return this->line_for_postprocess.at(i).lmin_cell;
  }
}

/* ************************************************************************** */
template <int d>
int LinesEval<d>::GetLength() const
{
  return this->lines_for_postprocess.size();
}

/* ************************************************************************** */
template <int d>
LineEval LinesEval<d>::GetLine(int i) const
{
  if(i >= (int)this->lines_for_postprocess.size() || i < 0)
  {
    ErrThrow("LinesEval::GetLine: index ", i, " out of range. There are only ",
             this->lines_for_postprocess.size(), " lines available.");
  }
  return this->lines_for_postprocess.at(i);
}

/* ************************************************************************** */
bool LineEval::IsFromDB() const
{
  return this->fromDB;
}

/* ************************************************************************** */
double LineEval::mean_value(const TFEFunction3D& f) const
{
  double val_p0;
  double val_p1;
  double dlength;
  std::array<double,3> p0;
  std::array<double,3> p1;

  double val    = 0.;
  double length = 0.;

  p0 = base_point;
  p1 = base_point;
  
  Output::info<3>("Post-processing function", "Space averaging using trapezoid "
                  "rule (holes in the domain are not consider).");

#ifdef _MPI
  Output::warn("Post-processing function", "Space averaging is not yet "
               "implemented for MPI.");
  return val;
#endif

  for( int i=0 ; i<GetNbPoints()-1 ; i++ )
  {
    int i_cell0 = GetCellIdx(i);
    int i_cell1 = GetCellIdx(i+1);

    if( i_cell0==-1 || i_cell1==-1 )
    {
      continue;
    }

    p0[direction] = GetPosition(i);
    p1[direction] = GetPosition(i+1);

    dlength = p1[direction] - p0[direction];

    f.FindValueLocal(GetCell(i_cell0), i_cell0, p0[0], p0[1], p0[2], &val_p0);
    f.FindValueLocal(GetCell(i_cell1), i_cell1, p1[0], p1[1], p1[2], &val_p1);

    val += (val_p0 + val_p1)/2. * dlength;
    length += dlength;
  }

  val /= (length==0.) ? 1. : length;

  return val;
}

/* ************************************************************************** */
double LineEval::space_average_value(const TFEFunction3D& f) const
{
  double val = 0.;

  if( ! IsFromDB() )
  {
    val = mean_value(f);
  }
  else
  {
    int                  l;
    int                  i_cell;
    int                  N_LinePoints;
    double               p0;
    double               p1;
    double               val_tmp;
    double               hE;
    std::array<double,3> pnt_eval;
    const double*        zeta;
    const double*        LineWeights;
    const TFESpace3D*    FSpace;
    TBaseCell*           cell;
    FE3D                 CurrentElement;
    QuadFormula1D        LineQuadFormula;
    TQuadFormula1D*      qf1;

    FSpace = f.GetFESpace3D();

    pnt_eval = base_point;

    Output::info<3>("Post-processing function", "Space averaging using "
                     "quadrature formula.");

    // loop over all cells cuting the line
    for( unsigned int i=0 ; i<line_for_postprocess.size() ; i++ )
    {
      i_cell = line_for_postprocess.at(i).cell_index;

      if( i_cell<0 ) //only append in case of many processes
      {
        continue;
      }

      cell = coll->GetCell(i_cell);

      p0 = line_for_postprocess.at(i).lmin_cell;
      p1 = line_for_postprocess.at(i).lmax_cell;

      // get polynomial degree of finite element in current cell
      CurrentElement = FSpace->GetFE3D(i_cell, cell);
      l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(CurrentElement);

      // define quadrature formula for triangular face
      LineQuadFormula = TFEDatabase3D::GetQFLineFromDegree(l);

      // get quadrature formula for 1D (line)
      qf1 = TFEDatabase3D::GetQuadFormula1D(LineQuadFormula);

      // get number of quad points, weights and positions
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

      // half of the length of the segment
      hE = 0.5 * (p1 - p0);

      // sum over all quad points
      for( int j=0 ; j<N_LinePoints ; j++ )
      {
        // get quadrature point in original element
        pnt_eval[direction] = p0 + 0.5*(p1-p0)*(zeta[j]+1);

        // get value at this quadrature point (in s)
        f.FindValueLocal(cell, i_cell,
                        pnt_eval[0], pnt_eval[1], pnt_eval[2],
                        &val_tmp);

        // multiply value with weight from quadrature formula and determinant
        // from integral transformation to the unit edge (-1,1)
        val += val_tmp * LineWeights[j] * hE;
      }
    }
  }

#ifdef _MPI
  MPI_Allreduce(MPI_IN_PLACE, &val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  return val;
}

#ifdef __3D__
template class LinesEval<3>;
#else
template class LinesEval<2>;
#endif
