/**
 * @file This example file contains an upward
 * directed (stationary) flow through a wound up pipe.
 * A large portion of the file deals with methods which are used to
 * transform a cylindrical sandwich geometry to a helically coiled tube
 * geometry. These methods have to be calles onto the TCollection
 * before the rest of the problem (FESpaces, FEFUnctions,
 * Matrices, Parallel Communicators etc.) is set up.
 *
 * @date 2016/04/26
 * @author Various, imported to ParMooN by Clemens Bartsch
 */
#ifndef TWISTED_PIPE_FLOW_
#define TWISTED_PIPE_FLOW_

namespace TwistedPipeConstants
{
size_t n_twists = 5;
size_t n_segments_per_twist = 50;

double l_inflow = 10;   // (cm)length of the straight inflow part
size_t n_segments_inflow = 2;

double l_outflow = 20;  // (cm) length of the straight outflow part
size_t n_segments_outflow = 4;

double tube_radius = 0.3;  // (cm) must fit the .PRM and .GEO!!!!!!
double twist_radius = 5.9; //5.9; // (cm) from center of tube to central axis
double space_between_twists = 0.3; //(cm) vertical distance between two coils

}

namespace FluidProperties
{
double eta = 1.19e-3;       // the dynamic viscosity ( kg / (m*s) ) here: of Kalialaun, as reported by V.Wiedmeyer
double rho = 1100;          // the density           ( kg /  m^3) ) here: of a Kalialaun solution, as reported by V.Wiedmeyer

double u_infty = 0.01;      // the characteristic velocity of the fluid (was UREA_u_infty before)
double l_infty = 0.01;      // the characteristic length scale of the tube (was UREA_l_infty before)

// note: in the coefficients function the diffusion coefficient will
// be calculated as:
//      eps = (eta/rho) / (u_infty*l_infty);
}

class derived_properties
{
  public:
    static double l_tube;

    static double l_twisted_part;

    static std::vector<double> segment_marks;

    static void set_statics()
    {
      using namespace TwistedPipeConstants;

      double PI = 3.14159265358979323846;
      double h =  2*tube_radius +  space_between_twists; //ganghoehe der rohrmitte-helix
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
};

//init static members
double derived_properties::l_tube = 0;
double derived_properties::l_twisted_part = 0;
std::vector<double> derived_properties::segment_marks = {};

enum class PIPE_PIECE {INFLOW, COIL, OUTFLOW}; //TODO _PIECE and _FACE!

/// Determine in which part of the pipe a certain point lies.
/// Note that this is very coarse - it checks whether (x,y,z) is contained
/// in a closed rectangular bounding box of either in- or outflow and if not so
/// assume it must lie in the coiled part.
PIPE_PIECE in_which_piece_is_this_point(double x, double y, double z)
{
  using namespace TwistedPipeConstants;

  if(  -l_inflow <= x     && x <= 0
    && -tube_radius <= y  && y <= tube_radius
    && -tube_radius <= z  && z <= tube_radius) //check inflow bounding box
  {
    return PIPE_PIECE::INFLOW;
  }
  double height = (2*tube_radius +  space_between_twists)*n_twists;
  if(  0 <= x     && x <= l_outflow
      && -tube_radius <= y  && y <= tube_radius
      && -tube_radius + height <= z  && z <= tube_radius + height) //check outflow bounding box
  {
   return PIPE_PIECE::OUTFLOW;
  }
  //we assume that no garbage goes in here - so everything which is neither
  //inflow nor outflow is assumed to lie in the coil
  return PIPE_PIECE::COIL;
}

/// This sets some parameters which are relevant for the Domain and must thus
/// be called before Domain is initialized!
void ExampleFile()
{
  // set global parameters which are too deeply rooted in the
  // code to be easily removed at the moment
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1356; //FIXME Remove this from the entire code!
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;

  derived_properties::set_statics();

  Output::print(" > Example: twisted_pipe_flow.h.");

  using namespace FluidProperties;
  double eps = (eta/rho) / (u_infty*l_infty);
  Output::print(" > diffusion coefficient eps = ", eps);
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
//  // Poiseuille flow with R = 1, dp/dz = -0.1
//  // v = -1/(4 eta) (R^2-r^2) * nabla p
//  // ParamDB->UREA_INFLOW_SCALE  in ml/min !!!
//  // based on derivation 13/12/06
//  double R, R2, r2yz;
//  double l_infty = TDatabase::ParamDB->UREA_l_infty;
//
//  R = 0.3;
//  R2=R*R;
//
//  //if (y > 0)
//  //  values[0] =
//  //    (R2-(y-R_bent)*(y-R_bent)-z*z) * V_in /(3000 * Pi * R2*R2)/TDatabase::ParamDB->UREA_u_infty;
//  //else
//  //  values[0] =
//  //    -(R2-(y+R_bent)*(y+R_bent)-z*z) * V_in /(3000 * Pi * R2*R2)/TDatabase::ParamDB->UREA_u_infty;
//  double vFlux = 7.2;
//  r2yz = (y-5.65)*(y-5.65)+z*z;
//
//  // factor covering inflow from boundary, dimensionless representation, etc.
//  double uFactor = 2. * vFlux * TDatabase::ParamDB->UREA_INFLOW_SCALE
//      / (Pi * R2 * TDatabase::ParamDB->UREA_u_infty);
//  values[0] = uFactor * (1-r2yz/R2);
//  //if (TDatabase::TimeDB->CURRENTTIME < 1.0)
//  //  values[0] *= TDatabase::TimeDB->CURRENTTIME;
//
//  values[1] = 0;
//  values[2] = 0;
//  values[3] = 0;
//  values[4] = 0;
//
//  //OutPut(values[0] << " ");
//  // if (abs(y-5.65)>0.3)
//  //  exit(1);
//  //OutPut(" in " << x << " "  << y << " " << z << " " << values[0] << " : ");

  PIPE_PIECE piece = in_which_piece_is_this_point(x,y,z);
  switch(piece)
  {
    case PIPE_PIECE::INFLOW:
      values[0]=1; //TODO Hagen-P.!
      break;
    default:
      values[0]=0;
      break;
  }
  values[1] = 0; //TODO Hagen-P.!
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}



/// Initial solution in U1 direction - Hagen-Poiseuille in inflow piece.
void InitialU1(double x, double y, double z, double *values)
{
  // Poiseuille flow with R = 1, dp/dz = -0.1
  // v = -1/(4 eta) (R^2-r^2) * nabla p
  double val[5];

  ExactU1(x,y,z,val);
  //if ((x<50)&& (y > 0))
  if (x < 50)
  {
    values[0] = val[0];
    //OutPut(y << " " << val[0] << "::");
  }
  else
    values[0] = 0;

  values[0] = 0.0;
  //OutPut(values[0] << " ");
}

/// Initial solution in U2 direction. No flow, nowhere.
void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

/// Initial solution in U3 direction. No flow, nowhere.
void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
}

/// Initial pressure. No pressure, nowhere.
void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}


// ========================================================================
// boundary conditions
// ========================================================================

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
//  double eps = 1e-6, r2yz;
//
//  cond = DIRICHLET;
//
//  r2yz =  (y-5.65)*(y-5.65)+(z-5.3333333333333325)*(z-5.3333333333333325);
//
//  if (x>53.75-25)
//  {
//    cond = NEUMANN;
//    OutPut("neum " << x <<" " << y << " " << z << endl);
//  }
//  //if ((fabs(x-43.43137254901961)< eps))//&& ((0.09-r2yz) < 1e-3))
//  //if (x>53.75)
//  //  OutPut("bdry " <<  setprecision(16) << x <<" " << y << " " <<  setprecision(16) << z << " "  << r2yz << endl);
//  //if ((fabs(60-x)<eps))
//  //if (z>5.03333)
//  // {
//  // outflow
//  // cond = NEUMANN;
//  //OutPut("neum " << x <<" " << y << " " << z << endl);
//  // TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
//  // }

  PIPE_PIECE piece = in_which_piece_is_this_point(x,y,z);
  switch(piece)
  {
    case PIPE_PIECE::INFLOW:
      cond = DIRICHLET;
      break;
    default:
      cond = DIRICHLET;
      break;
  }
  //TODO Add Dirichlet at inflow and Neumann at outflow here!
}


// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
//  double eps = 1e-6, val[5];
//
//  value = 0.0;
//
//  // inflow
//  //if ((fabs(x)<eps)&&(y>0))
//  if ((fabs(x)<eps))
//  {
//    ExactU1(x,y,z,val);
//    value = val[0];
//  }
//  if (x>53.75-25)
//  {
//    ExactU1(0,y,z-5.3333333333333325,val);
//    value = val[0];
//    value = 0;
//    OutPut("out " << x << " "  << y << " " << z << " " << value << endl);
//  }
//  //OutPut(value << " ");

  PIPE_PIECE piece = in_which_piece_is_this_point(x,y,z);
  switch(piece)
  {
    case PIPE_PIECE::INFLOW:
      value = 1; //TODO Hagen-P.!
      break;
    case PIPE_PIECE::OUTFLOW:
      value = -1;
      break;
    default:
      value =0;
      break;
  }
}


// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  using namespace FluidProperties;

  double eps =  (eta / rho) / (l_infty * u_infty);

  for(int i=0; i < n_points; i++)
  {
    double* coeff = coeffs[i];

    coeff[0] = eps;       // momentum diffusion
    coeff[1] = 0;         // f1
    coeff[2] = 0;         // f2
    coeff[3] = 0;         // f3
    coeff[4] = 0;         // g
  }
}


/**
 * Traverse a collection of cells and in each of its vertices
 * swap the x and z coordinate.
 *
 * As these Cartesian coordinates are only stored in the vertices,
 * it is enough to do the swapping there in order to swap it all.
 *
 * THis is used to "lie down" the initially standing sandwich grid.
 *
 * @param coll The collection to treat in such a way.
 */
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
        vertex->SetCoords(z-TwistedPipeConstants::l_inflow, y, x);
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
  using namespace TwistedPipeConstants;

  // three cases - inflow piece, coil, outflow piece
  if( x  < 0)
  {// inflow piece - let the point untouched
    x_coiled = x;
    y_coiled = y;
    z_coiled = z;
  }
  else if ( x < derived_properties::l_twisted_part)
  {//coil
    double t_of_x = x / (derived_properties::l_twisted_part / n_twists);
    double r_of_yz = twist_radius - y;
    double h =  2*tube_radius +  space_between_twists; //constant...
    double c_of_yz = z;

    double PI = 3.14159265358979323846; //...constant

    x_coiled = r_of_yz * cos(2*PI*t_of_x - PI/2);
    y_coiled = r_of_yz * sin(2*PI*t_of_x - PI/2);
    z_coiled = h * t_of_x + c_of_yz;

    //move the center of the coil in accordance with the fixed inflow position
    y_coiled += twist_radius;
  }
  else
  {// outflow piece - attach the outflow to the coil
    double h =  2*tube_radius +  space_between_twists; //constant...
    x_coiled = x - derived_properties::l_twisted_part;
    y_coiled = y;
    z_coiled = z + h * n_twists;
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
#endif /*TWISTED_PIPE_FLOW_*/
