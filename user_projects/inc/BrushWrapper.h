/**
 * BrushWrapper.h
 *
 * This wraps up Brushs interface to ParMooN in order to use it natively in
 * ParMooN. Especially it should be able to transform the output of Brush into
 * ParMooN style FEFunctions, so that they can be
 *  - visualized with Paraview
 *  - used in the assembling process
 *  - ...
 *
 */

#ifndef USER_PROJECTS_INC_BRUSHWRAPPER_H_
#define USER_PROJECTS_INC_BRUSHWRAPPER_H_

#include <Domain.h>
#include <ParameterDatabase.h>
#include <PostProcessing2D.h>
#ifdef __2D__
#include <FEFunction2D.h>
#elif __3D__
#include <FEFunction3D.h>
#endif

#include <valarray>

namespace Brush
{
  class InterfacePM;
}

class BrushWrapper
{
  public:

    /// Constructor taking a domain and a database.
    BrushWrapper(TCollection* coll, const ParameterDatabase& db);

    /// Destructor.
    ~BrushWrapper();

    /// Get the k-moment of the particle distribution as piecewise constant
    /// fe function on the domain.
    const TFEFunction2D& get_moment_fe(int k) const;

    /// Give the values of velocity u to Brush.
    void set_velocity(const TFEVectFunct2D& u);
    void set_temperature(const TFEFunction2D& T);
    void set_concentrations(std::vector<const TFEFunction2D*> c);

    /// Run the particle solver.
    void solve(double t_start, double t_end);

    /// Write moments of the particle distribution to a .vtk-file, which was
    /// specified in the database given to the constructor.
    void output(double t);

  private:

    ParameterDatabase db_;

    // TODO Ownership for the object is taken, but there is some trouble with the
    // destructor.
    Brush::InterfacePM* interface_;

    PostProcessing2D output_writer_;

    TCollection* pd_moments_grid_; //NO OWNERSHIP TAKEN.

    // This is all adapted to 0th order elements.
    // The function values of the moments which come
    // from Bruhs are directly written into pd_moments_values_,
    // i.e. they can be directly used as fe values of 0th order moments.
    // Anything else would require an extrapolation.
    TFESpace2D pd_moments_space_;
    std::vector<TFEFunction2D*> pd_moments_;
    std::vector<std::vector<double>> pd_moments_values_;

    /// Those are the points in space at which
    /// Brush expects function values as input.
    std::vector<std::valarray<double>> output_sample_points_;


};



#endif /* USER_PROJECTS_INC_BRUSHWRAPPER_H_ */
