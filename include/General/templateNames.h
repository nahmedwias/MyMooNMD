#ifndef INCLUDE_GENERAL_TEMPLATENAMES_H_
#define INCLUDE_GENERAL_TEMPLATENAMES_H_

class Example_Darcy2D;
class Example_Darcy3D;
class TFEFunction2D;
class TFEFunction3D;
class TFESpace2D;
class TFESpace3D;
class TFEVectFunct2D;
class TFEVectFunct3D;

/**
 * @brief this struct gives us the opportunity to write code which does not
 * depend on the space dimension.
 *
 * Basically, classes which are implemented separately for 2D and 3D can be 
 * used in the same code via 
 *     using FEFunction = typename Template_names<d>::FEFunction;
 * where `d` is either 2 or 3.
 * 
 */
template <int d>
struct Template_names
{
};
template <>
struct Template_names<2>
{
  typedef Example_Darcy2D Example_Darcy;
  typedef TFEFunction2D FEFunction;
  typedef TFESpace2D FESpace;
  typedef TFEVectFunct2D FEVectFunct;
};
template <>
struct Template_names<3>
{
  typedef Example_Darcy3D Example_Darcy;
  typedef TFEFunction3D FEFunction;
  typedef TFESpace3D FESpace;
  typedef TFEVectFunct3D FEVectFunct;
};

#endif // INCLUDE_GENERAL_TEMPLATENAMES_H_
