/**
 * @file Implementation of namespace and methods in CoiledPipe.h.
 * @author Clemens Bartsch
 * @date 2016/04/28
 */

#include <CoiledPipe.h>
#include <Collection.h>

namespace CoiledPipe
{

namespace GeoConsts
{
// directly input consts
size_t n_twists;
size_t n_segments_per_twist;

double l_inflow;
size_t n_segments_inflow;

double l_outflow;
size_t n_segments_outflow;

double tube_radius;
double twist_radius;
double space_between_twists;

// dependent consts
double l_tube;
double l_twisted_part;
std::vector<double> segment_marks;
double h;
}

void set_up_geoconsts(
    size_t in_n_twists,
    size_t in_n_segments_per_twist,
    double in_l_inflow,
    size_t in_n_segments_inflow,
    double in_l_outflow,
    size_t in_n_segments_outflow,
    double in_tube_radius,
    double in_twist_radius,
    double in_space_between_twists)
{
  using namespace GeoConsts;

  // 1) just set them input parameters
  n_twists = in_n_twists;
  n_segments_per_twist = in_n_segments_per_twist;
  l_inflow = in_l_inflow;
  n_segments_inflow = in_n_segments_inflow;
  l_outflow = in_l_outflow;
  n_segments_outflow = in_n_segments_outflow;
  tube_radius = in_tube_radius;
  twist_radius = in_twist_radius;
  space_between_twists = in_space_between_twists;


  // 2) compute the dependent quantities
  double PI = 3.14159265358979323846;
  h =  2*tube_radius +  space_between_twists; //ganghoehe der rohrmitte-helix
  double k = h / (2*PI*twist_radius); //steigung der rohrmitte-helix

  double s_1 = 2*PI*twist_radius * sqrt(1 + k*k); //Bogenlaenge einer vollen Rundung

  l_twisted_part = n_twists * s_1;

  //length of the tube in total (will form DRIFT_Z)
  l_tube = l_inflow + l_twisted_part + l_outflow;

  // set up the vector of segment marks
  size_t n_segment_marks =
      n_segments_inflow
      + n_twists * n_segments_per_twist
      + n_segments_outflow + 1;

  segment_marks = std::vector<double>(n_segment_marks , 0.0);

  double inflow_fraction = l_inflow / l_tube;
  double twisted_fraction = l_twisted_part / l_tube;
  double outflow_fraction = l_outflow / l_tube;

  for (size_t i =0 ; i < n_segment_marks ; ++i)
  {
    if(i < n_segments_inflow)
    {//inflow part
      segment_marks.at(i) = (double) i / n_segments_inflow * inflow_fraction;
    }
    else if (i < n_segments_inflow + n_twists * n_segments_per_twist )
    {//twisted part
      size_t i_twist = i - n_segments_inflow;
      segment_marks.at(i) = inflow_fraction + (double) i_twist / (n_twists * n_segments_per_twist) * twisted_fraction;
    }
    else
    {//outflow part
      size_t i_out = i - n_segments_inflow - n_twists * n_segments_per_twist;
      segment_marks.at(i) = inflow_fraction + twisted_fraction
          + (double) i_out / n_segments_outflow * outflow_fraction;
    }
  }
}

/// Determine in which part of the pipe a certain point lies.
/// Note that this is very coarse - it checks whether (x,y,z) is contained
/// in a closed rectangular bounding box of either in- or outflow and if not so
/// assume it must lie in the coiled part.
PIPE_PIECE in_which_piece_is_this_point(double x, double y, double z)
{
  using namespace GeoConsts;

  double tol = 1e-10; //TODO hope that tolerance is good enough to distinguish parts!

  if(  -l_inflow <= x     && x <= 0
      && -tube_radius <= y  && y <= tube_radius
      && -tube_radius <= z  && z <= tube_radius) //check inflow bounding box
  {

    if(fabs(- l_inflow - x ) < tol
        && sqrt(y*y + z*z) < tube_radius - tol) //check inflow face
    {
      return PIPE_PIECE::INFLOW_FACE;
    }
    else //inflow piece, but not on inflow face
      return PIPE_PIECE::INFLOW_PIECE;
  }

  double height = GeoConsts::h*n_twists;
  if(  0 <= x     && x <= l_outflow
      && -tube_radius <= y  && y <= tube_radius
      && -tube_radius + height <= z  && z <= tube_radius + height) //check outflow bounding box
  {
    if(fabs(l_outflow - x ) < tol
        && sqrt(y*y + (z-height)*(z-height)) < tube_radius - tol) //check outflow face
    {
      return PIPE_PIECE::OUTFLOW_FACE;
    }
    else
      return PIPE_PIECE::OUTFLOW_PIECE;
  }
  //we assume that no garbage input goes in here - so everything which is neither
  //inflow nor outflow is assumed to lie in the coil
  return PIPE_PIECE::COIL;
}

void swap_x_and_z_coordinates(TCollection *coll)
{

  int N_Cells = coll->GetN_Cells();

  // initialise all Vertex ClipBoards to 0
  for(int i=0 ; i<N_Cells ; i++)
  {
    TBaseCell* cell = coll->GetCell(i);
    int n_v = cell->GetN_Vertices();

    for (int j=0 ; j<n_v ; j++)
    {
      TVertex* vertex = cell->GetVertex(j);
      vertex->SetClipBoard(0);
    }
  }

  // swap x and z coordinates if vertex->GetClipBoard() = 0
  for(int i=0 ; i<N_Cells ; i++)
  {
    TBaseCell* cell = coll->GetCell(i);
    int n_v = cell->GetN_Vertices();

    for (int j=0 ; j<n_v ; j++)
    {
      TVertex* vertex = cell->GetVertex(j);
      if (!vertex->GetClipBoard())
      {
        double x, y, z;
        vertex->GetCoords(x, y, z);
        //swithc coordinates and place the start of the coiled tube at (0,0,0)
        vertex->SetCoords(z - GeoConsts::l_inflow, y, x);
        vertex->SetClipBoard(1);
      }
    }
  }

}

/**
 * This method calculates the new position of a given point in a
 * pipe coiled according to the parameters given in namespace TwistedPipeConstants.
 *
 * @param[in] x Old x value (in straight setup).
 * @param[in] y Old y value (in straight setup).
 * @param[in] z Old z value (in straight setup).
 * @param[out] x_coiled New x value (in coiled setup).
 * @param[out] y_coiled New y value (in coiled setup).
 * @param[out] z_coiled New z value (in coiled setup).
 */
void compute_position_in_coiled_pipe(
    double x, double y,double z,
    double& x_coiled, double& y_coiled, double& z_coiled)
{
  using namespace GeoConsts;

  // three cases - inflow piece, coil, outflow piece
  if( x  < 0)
  {// inflow piece - let the point untouched
    x_coiled = x;
    y_coiled = y;
    z_coiled = z;
  }
  else if ( x < l_twisted_part)
  {//coil
    double t_of_x = x / (l_twisted_part / n_twists);
    double r_of_yz = twist_radius - y;
    double c_of_yz = z;

    double PI = 3.14159265358979323846; //...constant

    x_coiled = r_of_yz * cos(2*PI*t_of_x - PI/2);
    y_coiled = r_of_yz * sin(2*PI*t_of_x - PI/2);
    z_coiled = GeoConsts::h * t_of_x + c_of_yz;

    //move the center of the coil in accordance with the fixed inflow position
    y_coiled += twist_radius;
  }
  else
  {// outflow piece - attach the outflow to the coil
    x_coiled = x - l_twisted_part;
    y_coiled = y;
    z_coiled = z + GeoConsts::h * n_twists;
  }
}

/**
 * This is the second method to call from the program -
 * it takes care of the coiling of the pipe.
 */
void coil_pipe_helically(TCollection *coll)
{
  int N_Cells = coll->GetN_Cells();

  // initialise ClipBoard
  for(int i=0 ; i<N_Cells ; i++)
  {
    TBaseCell* cell = coll->GetCell(i);
    int n_verts = cell->GetN_Vertices();

    for (int j=0 ; j<n_verts ; j++)
    {
      TVertex* vertex = cell->GetVertex(j);
      vertex->SetClipBoard(0);
    }
  }

  // visit each vertex, reset its coordinates to its new position
  for(int i=0 ; i < N_Cells ; i++)
  {
    TBaseCell* cell = coll->GetCell(i);
    int n_verts = cell->GetN_Vertices();

    for (int j=0 ; j < n_verts ; j++)
    {
      TVertex* vertex = cell->GetVertex(j);
      if (!vertex->GetClipBoard())
      {
        double x, y, z;
        double x_coiled, y_coiled, z_coiled;

        vertex->GetCoords(x, y, z);
        compute_position_in_coiled_pipe(x, y, z, x_coiled, y_coiled, z_coiled);
        vertex->SetCoords(x_coiled, y_coiled, z_coiled);

        //mark this vertex as treated
        vertex->SetClipBoard(1);
      }
    }
  }
}

}
