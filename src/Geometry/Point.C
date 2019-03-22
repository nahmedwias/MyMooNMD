#include <limits> // std::numeric_limits<double>::quiet_NaN
#include <cmath>  // std::isnan
#include <MooNMD_Io.h>
#include <Point.h>

/* ************************************************************************** */
Point::Point(double x)
 : X(x), Y(std::numeric_limits<double>::quiet_NaN()), 
   Z(std::numeric_limits<double>::quiet_NaN())
{
  
}

/* ************************************************************************** */
Point::Point(double x, double y)
 : X(x), Y(y), Z(std::numeric_limits<double>::quiet_NaN())
{
  
}

/* ************************************************************************** */
Point::Point(double x, double y, double z)
 : X(x), Y(y), Z(z)
{
  
}

/* ************************************************************************** */
Point::Point(const std::vector<double>& p)
{
  const unsigned int p_size = p.size();
  if(p_size == 2)
  {
    this->X = p[0];
    this->Y = p[1];
    this->Z = std::numeric_limits<double>::quiet_NaN();
  }
  else if(p_size == 3)
  {
    this->X = p[0];
    this->Y = p[1];
    this->Z = p[2];
  }
  else if(p_size == 1)
  {
    this->X = p[0];
    this->Y = std::numeric_limits<double>::quiet_NaN();
    this->Z = std::numeric_limits<double>::quiet_NaN();
  }
  else
  {
    ErrThrow("unable to create a Point with a vector of size ", p_size);
  }
}

/* ************************************************************************** */
Point::Point(unsigned int dimension)
 : X(0.0), Y(0.0), Z(0.0)
{
  if(dimension == 0 || dimension > 3)
  {
    ErrThrow("cannot create a Point of dimension ", dimension);
  }
  if(dimension < 3)
  {
    this->Z = std::numeric_limits<double>::quiet_NaN();
    if(dimension < 2)
      this->Y = std::numeric_limits<double>::quiet_NaN();
  }
}

/* ************************************************************************** */
unsigned int Point::dimension() const
{
  if(!std::isnan(this->Z))
    return 3;
  else if(!std::isnan(this->Y))
    return 2;
  else
    return 1;
}

