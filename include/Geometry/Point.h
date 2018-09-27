#ifndef __POINT__
#define __POINT__

#include <vector>


/* ****************************************************************************
*
* @class      Point
* @brief      represent a point in space
*
* This is not a Vertex. A Vertex is part of a Mesh, while a Point describes an
* arbitrary position in the computational domain or a reference element.
*
* Usually a Point will be a quadrature point either in the computational domain
* or in a reference cell. This is where basis functions, coefficients, boundary
* conditions, boundary data, exact solutions and such are evaluated.
* 
* It is possible to construct a 1D Point, however this is not the usual case.
* 
* All constructors are explicit, to prevent accidentally creating Point objects.
*
* @author     Ulrich Wilbrandt
* @date       07.11.2014
*
***************************************************************************** */

class Point
{
 public:
  /** @brief construct a 1D Point x */
  explicit Point(double x);
  /** @brief construct a 2D Point (x,y) */
  explicit Point(double x, double y);
  /** @brief construct a 3D Point (x,y,z) */
  explicit Point(double x, double y, double z);
  
  /** @brief construct a point from a vector 
   * 
   * If p.size()==1, this is equivalent to Point(p[0]).
   * If p.size()==2, this is equivalent to Point(p[0], p[1]).
   * if p.size()>=3, this is equivalent to Point(p[0], p[1], p[2]).
   * If none of the above is true, an exception is thrown, as this is probably
   * not intended.
   */
  explicit Point(std::vector<double> p);
  
  /** @brief construct a Point with the given dimension (1,2, or 3) 
   * 
   * If the dimension is not 1,2, or 3, an exception is thrown.
   */
  explicit Point(unsigned int dimension);
  
  /** @brief copy constructor */
  Point(const Point &) = default;
  
  /** @brief move constructor */
  Point(Point &&) noexcept = default;
  
  /** @brief copy assignment */
  Point& operator=(Point&) noexcept = default;
  
  /** @brief move assignment */
  Point& operator=(Point&&) noexcept = default;
  
  /** @brief destructor */
  ~Point() noexcept  = default;
  
  
  /** @brief return a constant reference to the x-coordinate */
  const double& x() const
  {return X; }
  /** @brief return a reference to the x-coordinate */
  double& x()
  {return X; }
  
  /** @brief return a constant reference to the y-coordinate */
  const double& y() const
  {return Y; }
  /** @brief return a reference to the y-coordinate */
  double& y()
  {return Y; }
  
  /** @brief return a constant reference to the z-coordinate 
   * 
   * Note that this will return std::numeric_limits<double>::quiet_NaN() for a 
   * 2D Point. (Should we throw an exception in such a case?)
   */
  const double& z() const
  {return Z; }
  /** @brief return a reference to the z-coordinate 
   * 
   * Note that this will return std::numeric_limits<double>::quiet_NaN() for a 
   * 2D Point. (Should we throw an exception in such a case?)
   */
  double& z()
  {return Z; }
  
  /** @brief return the space dimension (either 1,2, or 3) this point is in 
   * 
   * Here we use that if a coordinate is set to nan, it is not used.
   */
  unsigned int dimension() const;
  
 private:
  /** @brief the x-coordinate of the point */
  double X;
  /** @brief the y-coordinate of the point 
   * 
   * In 1D this is set to be std::numeric_limits<double>::quiet_NaN().
   */
  double Y;
  /** @brief the z-coordinate of the point 
   * 
   * In 1D and 2D this is set to be std::numeric_limits<double>::quiet_NaN().
   */
  double Z;
};

#endif // __POINT__