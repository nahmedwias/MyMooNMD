/**
 * @file Holds routines which are used to turn a cylindrical
 * sandwich geometry into a coiled pipe.
 *
 * Is used in the project "twisted_pipe" (modelling of the flow in a tube
 * crystallizer) and could be used as a model of how to control and
 * modify a ParMooN sandwich geometry.
 *
 * @author various, rewritten by Clemens Bartsch
 * @date 2016/04/28
 */

#ifndef COILED_PIPE_H_
#define COILED_PIPE_H_

#include <cstddef>
#include <vector>

//forwared declaration
class TCollection;

namespace CoiledPipe
{

/// Constants describing the geometry
namespace GeoConsts
{
//all following variables are implemented as input data
extern size_t n_twists;
extern size_t n_segments_per_twist;

extern double l_inflow;   // (cm)length of the straight inflow part
extern size_t n_segments_inflow;

extern double l_outflow;  // (cm) length of the straight outflow part
extern size_t n_segments_outflow;

extern double tube_radius; // (cm) must fit the .PRM and .GEO!!!!!!
extern double twist_radius; //5.9; // (cm) from center of tube to central axis
extern double space_between_twists; //(cm) vertical distance between two coils

//all following variables are implemented as dependently computed data
extern double l_tube;
extern double l_twisted_part;
extern std::vector<double> segment_marks;
extern double h;

}

/**
 *
 */
void set_up_geoconsts(
    size_t n_twists,
    size_t n_segments_per_twist,
    double l_inflow,
    size_t n_segments_inflow,
    double l_outflow,
    size_t n_segments_outflow,
    double tube_radius,
    double twist_radius,
    double space_between_twists
);

/// Enum class naming the sections of the pipe.
enum class PIPE_PIECE {INFLOW_FACE, INFLOW_PIECE, COIL,
    OUTFLOW_PIECE, OUTFLOW_FACE};

/// Determine in which part of the pipe a certain point lies.
/// Note that this is very coarse - it checks whether (x,y,z) is contained
/// in a closed rectangular bounding box of either in- or outflow and if not so
/// assume it must lie in the coiled part.
PIPE_PIECE in_which_piece_is_this_point(double x, double y, double z);


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
void swap_x_and_z_coordinates(TCollection *coll);

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
    double& x_coiled, double& y_coiled, double& z_coiled);

/**
 * This is the second method to call from the program -
 * it takes care of the coiling of the pipe.
 */
void coil_pipe_helically(TCollection *coll);

}

#endif
