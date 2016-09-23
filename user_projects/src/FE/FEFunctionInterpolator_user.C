/**
 * Implements class FEFunctionInterpolator declared in FEFunctioninterpolator.h
 *
 * @author Clemens Bartsch
 * @date 2016/02/05
 */

#include <FEFunctionInterpolator_user.h>
#include <MooNMD_Io.h>

#include <FESpace1D.h>
#include <FESpace2D.h>
#include <FESpace3D.h>

#include <FEFunction2D_user.h>

#include <memory>

FEFunctionInterpolator::FEFunctionInterpolator(
    const TFESpace* fe_space)
{
  // try to downcast the input pointer to figure out the dimension of the space
  // (a nice alternative would be a virtual get_dimension() method in TFESpace
  //
  // TODO: syntax for shared pointers would be
  //  if(std::dynamic_pointer_cast<TFESpace1D>(fe_space))
  //  {
  //    dimension_ = 3;
  //  }
  //
  //
  if(dynamic_cast<const TFESpace1D*>(fe_space))
  {
    dimension_ = 1;
  }
  else if(dynamic_cast<const TFESpace2D*>(fe_space))
  {
    dimension_ = 2;
  }
#ifdef __3D__
  else if(dynamic_cast<const TFESpace3D*>(fe_space))
  {
    dimension_ = 3;
  }
#endif
  else
  {
    ErrThrow("The fe space to interpolate to must be of one of the"
        " derived types TFESpace1D, TFESpace2D or TFESpace3D!");
  }

  //store the shared_ptr to the space to interpolate to
  into_space_ = fe_space;

}

/** ************************************************************************ */

void FEFunctionInterpolator::check() const
{
  //call check for casting of the stored space
  check_for_shearing();
}

/** ************************************************************************ */

TFEFunction2D FEFunctionInterpolator::interpolate(
    const TFEFunction2D& original_funct,
    std::vector<double>& values_memory ) const
{
  //check if there is enough memory given


  double* values = &values_memory.front();
  int length = values_memory.size();
  const TFESpace2D* into_space_cast = dynamic_cast<const TFESpace2D*>(into_space_);

  TFEFunction2D interpolation(
      into_space_cast, (char*) "interpolated", (char*) "interpolated", values, length);

  // make use of method "Interpolate" for the moment
  // - TODO that code should be moved here entirely in time
  interpolation.Interpolate(&original_funct);

  return interpolation;
}

// ////////////////////// //    private method(s)   // ////////////////////// //

void FEFunctionInterpolator::check_for_shearing() const
{
  // switch over the dimension and check if the dynamic cast
  // to the correct derived FESpace can be done
  switch (dimension_)
  {
    case 1:
      if(!dynamic_cast<const TFESpace1D*>(into_space_))
      {
        ErrThrow("In dimension 1 the cast to TFESpace1D must be possible!");
      }
      else
      {
        Output::print<1>("Correct cast to TFESpace1D.");
      }
      break;
    case 2:
      if(!dynamic_cast<const TFESpace2D*>(into_space_))
      {
        ErrThrow("In dimension 2 the cast to TFESpace2D must be possible!");
      }
      else
      {
        Output::print<1>("Correct cast to TFESpace2D.");
      }
      break;
#ifdef __3D__
    case 3:
      if(!dynamic_cast<const TFESpace3D*>(into_space_))
      {
        ErrThrow("In dimension 3 the cast to TFESpace3D must be possible!");
      }
      else
      {
        Output::print<1>("Correct cast to TFESpace3D.");
      }
      break;
#endif
    default:
      ErrThrow("What is that? Neither dimension 1, 2 nor 3?"
          " Something's terribly wrong here.")
  }
}
