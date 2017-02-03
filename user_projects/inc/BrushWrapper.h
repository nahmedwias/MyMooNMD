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
  class ParMooNData;
}

namespace Exmpl
{
  enum class SourceAndSinkTerms;
}

class BrushWrapper
{
  public:

    /// Constructor taking a domain and a database.
    BrushWrapper(TCollection* coll, const ParameterDatabase& db);

    /// Destructor.
    ~BrushWrapper();

    ///Implicitely deletes copying and moving.

    /// Let Brush compute the source and sink terms which will
    /// be added to the right hand side of the CDR equations.
    std::vector<TFEFunction2D*> sources_and_sinks();

    /**
     * Reset the fluid phase data in Brush, i.e., velocity, pressure, temperature
     * and all species concentrations.
     * @param[in] u The fluid velocity (all components).
     * @param[in] p The pressure.
     * @param[in] species All transported species, the first of which
     *            is supposed to be temperature. The following must be in that
     *            order in which Brush expects them.
     * TODO Determine, fix and comment on units for all these quantities.
     */
    void reset_fluid_phase(
        const TFEVectFunct2D& u,
        const TFEFunction2D& p,
        std::vector<const TFEFunction2D*> species);

    /// Run the particle solver.
    void solve(double t_start, double t_end);

    /// Write moments of the particle distribution to a .vtk-file, which was
    /// specified in the database given to the constructor.
    void output(double t);

  private:

    /// The database for control parameters.
    ParameterDatabase db_;

    /// The ParMooN interface of Brush which is wrapped up by 'this'.
    /// TODO Ownership for the object is taken, but there is some trouble with the
    /// destructor.
    Brush::InterfacePM* interface_;


    PostProcessing2D output_writer_;

    std::ofstream moment_stats_file_;

    /// All functions which come from Brush start their life in ParMooN
    /// as 0th order fe functions on the following grid and fe space.
    TCollection* from_brush_grid_; //NO OWNERSHIP TAKEN.
    TFESpace2D from_brush_space_;

    std::vector<Exmpl::SourceAndSinkTerms> source_and_sink_requests_;
    std::vector<TFEFunction2D*> source_and_sink_fcts_;
    std::vector<std::vector<double>> source_and_sink_fcts_values_;

    // The function values of the moments which come
    // from Bruhs are directly written into pd_moments_values_,
    // i.e. they can be directly used as fe values of 0th order moments.
    // Anything else would require an extrapolation.
    std::vector<TFEFunction2D*> pd_moments_;
    std::vector<std::vector<double>> pd_moments_values_;

    /// Those are the points in space at which
    /// Brush expects function values as input.
    std::vector<std::valarray<double>> output_sample_points_;


};



#endif /* USER_PROJECTS_INC_BRUSHWRAPPER_H_ */
