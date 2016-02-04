/** ************************************************************************
 *
 * @name         FEFunctionInterpolator
 *
 * @brief        Functions like a factory for FE functions. Construct it with
 * a certain FESpace, give it a reference to an FE Function and it will give
 * you an interpolation of the given function onto its stored FESpace.
 * (The details are still under construction).
 *
 * This replaces the (wicked!) former "Interpolate" method of the FEFunctions.
 *
 * TODO Implement first version.
 * TODO Unit test this thing
 * TODO Put it into the coupling of TNSE2D and TCD2D.
 *
 * @author       Clemens Bartsch
 * @date         2016/02/05
 **************************************************************************/

#ifndef USER_PROJECTS_INC_FEFUNCTIONINTERPOLATOR_H_
#define USER_PROJECTS_INC_FEFUNCTIONINTERPOLATOR_H_

#include <memory>

class FEFunctionInterpolator{
  public:
    /**
     * Construct and Interpolater object whith a given FESpace. That is the
     * space the interpolator is going to project into.
     * @param fe_space
     */
    FEFunctionInterpolator(std::shared_ptr<TFESpace> fe_space);

    //TFEFunction1D& interpolate(const TFEFunction1D) const;
    /**
     * Project a given TFEFunction onto our space and return a reference to
     * that object.
     * So far the only sensible way to call this is using the operator= of
     * TFEFunction2D. I hesistate calling it a copy assignment...
     * @param
     * @return
     */
    TFEFunction2D& interpolate(const TFEFunction2D) const;
    //TFEFunction3D& interpolate(const TFEFunction3D) const;


    /**
     * As for special member functions, since the class stores nothing than
     * a shared pointer we think for the moment, that the default ones
     * are created by all compilers and suffice.
     *
     * NEED BE TESTED - WHAT ABOUT DERIVED TFESPACExD classes?
     */

    //! Default copy constructor.
    FEFunctionInterpolator(const FEFunctionInterpolator&) = default;

    //! DelDefaultete move constructor.
    FEFunctionInterpolator(FEFunctionInterpolator&&) = default;

    //! Default copy assignment operator.
    FEFunctionInterpolator& operator=(const FEFunctionInterpolator&) = default;

    //! Default move assignment operator.
    FEFunctionInterpolator& operator=(FEFunctionInterpolator&&) = default;

    //! Default destructor.
    ~FEFunctionInterpolator() = default;

  private:
    /** The FESpace this produces interpolations to.
     * (Maybe later its more sensible to give this as a parameter every time..
     * we will see).
     */
    std::shared_ptr<TFESpace> projection_space_;

    //! The dimension of the stored space to project into.
    size_t dimension_;

    /**
     * Use this as a developmental tool to see, if we can maintain the 1D/2D/3D
     * properties of the
     */
    void check_for_shearing();


};



#endif /* USER_PROJECTS_INC_FEFUNCTIONINTERPOLATOR_H_ */
