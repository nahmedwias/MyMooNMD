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
#include <vector>
#include <FESpace.h>

// Experimental Macro to avoid double code.
#ifdef __2D__
class TFESpace2D;
class TFEFunction2D;
#define TFESpaceXD TFESpace2D
#elif defined(__3D__)
class TFESpace3D;
class TFEFunction3D;
#define TFESpaceXD TFESpace3D
#endif

class FEInterpolationCheatSheet{

  public:
    typedef std::vector<int> ContainingCells;

    /// Constructor. Needs both old and new space
    FEInterpolationCheatSheet(
        const TFESpaceXD* old_fe_space, const TFESpaceXD* new_fe_space);

    const ContainingCells& get_cheats(int new_cell, int new_point) const
    {
      return cheat_sheet.at(new_cell).at(new_point);
    }

    size_t get_n_points_to_cell(int cell) const
    {
      return cheat_sheet.at(cell).size();
    }

    size_t get_n_cells() const
    {
      return cheat_sheet.size();
    }

  private:
    /// This is the data.
    std::vector<std::vector<ContainingCells>> cheat_sheet;
};

class FEFunctionInterpolator{
  public:
    /**
     * Construct an Interpolater object whith a given FESpace. That is the
     * space the interpolator is going to project into.
     * @param fe_space
     */
    FEFunctionInterpolator(const TFESpace* into_space);

    /**
     * Project a given TFEFunction onto our space and return a reference to
     * that object.
     * Note that in the current setup a TFEFunction2D does not take care of
     * the value for its own memory. Usually the values are administered by
     * a BlockVector, which also takes care of the deletion.
     * For the new, interpolating fe function one has to provide memory and
     * must take care that it does not go out of scope before the fe function.
     * This is currently (2016/02/05) under investigation.
     *
     * So far the only sensible way to call this method is using the
     * (temporary) copy constructor of TFEFunction2D.
     *
     * That call would look like this:
     * \code
     * std::vector<double> entries_memory(length, 0.0);
     * TFEFunction2D new_function =
     * fe_function_interpolator.interpolate(old_function, entries_memory);
     * \endcode
     * Hereby length determines the expected size and old_function is an
     * TFEFunction2D to be interpolated.
     *
     * The method relies on a function Interpolate method in TFEFunction2D
     * - sooner or later that method should be reworked, tested and put here
     * entirely.
     *
     * @param original_funct The function to be interpolated.
     * @param values_memory Memory which will be used for the new fe function.
     *
     * @return A reference to a newly constructed fe function which interpolates
     * original_function on the stored fe space.
     */
#ifdef __2D__
    TFEFunction2D interpolate(
        const TFEFunction2D& original_funct,
        std::vector<double>& values_memory,
        bool use_cheat_sheet) const;
#elif defined(__3D__)
    TFEFunction3D& interpolate(
      const TFEFunction3D& original_funct,
      std::vector<double>& values_memory,
      bool use_cheat_sheet) const;
#endif

    //! Call (different?) check methods to test the state of the object.
    //! Will throw an error if something is wrong.
    void check() const;

    const TFESpace* get_from_space(){return from_space_;};

  private:
    /**
     * The FESpace this produces interpolations to.
     * (Maybe later its more sensible to give this as a parameter every time..
     * we will see).
     *
     * TODO Change to shared_ptr as soon as that concept is enforced in the
     * entire ParMooN.
     */
    //std::shared_ptr<TFESpace> into_space_;
    const TFESpace* into_space_;

    //! The dimension of the stored space to project into.
    size_t dimension_;


    mutable std::shared_ptr<FEInterpolationCheatSheet> cheat_sheet_;
    mutable const TFESpace* from_space_; //to check which space the interpolation cheat_sheet belongs to

    /**
     * Use this as a developmental tool to see, if we can maintain the 1D/2D/3D
     * properties of the
     */
    void check_for_shearing() const;

};



#endif /* USER_PROJECTS_INC_FEFUNCTIONINTERPOLATOR_H_ */
