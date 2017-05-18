/**
 * BrushWrapper.h
 *
 * This wraps up Brushs interface to ParMooN in order to use it natively in
 * ParMooN. Especially it is able to transform the output of Brush into
 * ParMooN style FEFunctions, so that they can be, e.g., visualized with Paraview
 * and used in the assembling of other systems.
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
/**
 * This enum class is implemented in Brush. It holds names for all those
 * example specific source and sink terms that Brush is able to compute. Each
 * BrushWrapper example will request its own source and sink terms from Brush.
 */
  enum class SourceAndSinkTerms;
}

class BrushWrapper
{
  public:

    /// Constructor taking a TCollection and a database.
    BrushWrapper(TCollection* coll, const ParameterDatabase& db);

    /// Destructor.
    ~BrushWrapper();

    ///Implicitely deletes copying and moving.

    /// Let Brush compute the source and sink terms which will
    /// be added to the right hand side of the CDR equations.
    std::vector<const TFEFunction2D*> sources_and_sinks();

    /**
     * Reset the fluid phase data in Brush, i.e., velocity, pressure, temperature
     * and all species concentrations.
     * @param[in] u The fluid velocity (all components).
     * @param[in] p The pressure.
     * @param[in] species All transported species, the first of which
     *            is supposed to be temperature. The following must be in that
     *            order in which Brush expects them.
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
    // Fill example depending data. This is like an "init" function,
    // called only by the constructor.
    //  0 - Eder's ASA flow crystallizer example.
    //
    void pick_example(int exmpl_code);

    //call this in order to fill the reset_fluid_phase_cheats_ cache
    void cache_output_point_containing_cells(const TFESpace2D& one_space);

    /// The database for control parameters.
    ParameterDatabase db_;

    /// The ParMooN interface of Brush which is wrapped up by 'this'.
    /// TODO Ownership for the object is taken, but there is some trouble with the
    /// destructor.
    Brush::InterfacePM* interface_;

    /* ****** CONTROL STRUCTURE FOR PARAMETERS THAT GO TO BRUSH*/
    /// Spatial dimension - must equal number of velocity components.
    size_t parameter_spatial_dimension_;
    /// Number of primary species that go to Brush as parameters.
    /// ('primary': the species is represented as an fe function in ParMooN)
    size_t parameter_n_specs_primary_;
    /// Number of derived species that go to Brush as parameters.
    /// ('derived': the species is computed from primary species in ParMooN)
    size_t parameter_n_specs_derived_;
    /// The vector of the names of the parameters. Must be of size
    /// parameter_spatial_dimension_ + 1 + parameter_n_specs_primary_ +
    /// parameter_n_specs_derived_. Where "+1" is for the fluid pressure.
    std::vector<std::string> parameter_function_names_;
    /// The real valued functions which compute the derived parameters from
    /// the primary parameters. One for each derived parameter.
    std::vector<std::function<double(const std::vector<double>&)>>
        parameter_specs_derived_fcts_;

    /// Those are the points in space at which
    /// Brush expects function values as input.
    std::vector<std::valarray<double>> parameter_sample_points_;

    /// Resetting the fluid phase in Brush requires point evaluation of fe
    /// functions (FindGradient). The cells in which to find the respective cells
    /// are cached here, which saves a lot of runtime.
    /// reset_fluid_phase_cheats_[p] holds the numbers of all cells that
    /// contain output point nr p.
    /// Note: THIS WIL ONLY WORK IF THE GRID IS THE SAME FOR ALL PARMOON FUNCTIONS!
    typedef std::vector<int> ContainingCells;
    std::vector<ContainingCells> reset_fluid_phase_cheats_;

    /* ****** CONTROL AND DATA STRUCTURE FOR SOURCES AND SINKS (COME FROM BRUSH)*/
    /// The following members bundle all that must be known for a representation
    /// of Brush's return values in ParMooN.
    /// All functions which come from Brush start their life in ParMooN
    /// as 0th order fe functions on the following grid and fe space.
    TCollection* from_brush_grid_; //NO OWNERSHIP TAKEN.
    TFESpace2D from_brush_space_;
    std::vector<const TFEFunction2D*> source_and_sink_fcts_;
    std::vector<std::vector<double>> source_and_sink_fcts_values_;
    // the following two informations come from the example
    std::vector<std::string> source_and_sink_function_names_;
    std::vector<Exmpl::SourceAndSinkTerms> source_and_sink_requests_;

    /* ****** PROGRAM OUTPUT OBJECTS */
    // (their sole purpose is post-processing.)

    // Files into which Brush writes its output. TODO Can't they be handled in Brush?
    std::ofstream moment_stats_file_;
    std::ofstream outflow_particles_file_;
    std::ofstream inflow_particles_file_;
    // Moments of the PSD which can be calculated by Brush. Used only for visual-
    // ization with paraview, which is performed by the output_writer_.
    PostProcessing2D output_writer_;
    std::vector<TFEFunction2D*> pd_moments_;
    std::vector<std::vector<double>> pd_moments_values_;




};



#endif /* USER_PROJECTS_INC_BRUSHWRAPPER_H_ */
