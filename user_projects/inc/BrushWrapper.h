/**
 * BrushWrapper.h
 *
 * This wraps up Brushs interface to ParMooN in order to use it natively in
 * ParMooN. Especially it is able to transform the output of Brush into
 * ParMooN style FEFunctions, so that they can be, e.g., visualized with Paraview
 * and used in the assembling of other systems.
 *
 * TODO When debugging is complete, disable fetch_moments - this might speed things up somewhat.
 */
#ifndef USER_PROJECTS_INC_BRUSHWRAPPER_H_
#define USER_PROJECTS_INC_BRUSHWRAPPER_H_

#include <Domain.h>
#include <GridTransferTool.h>
#include <Output3D.h>
#include <ParameterDatabase.h>
#include <PostProcessing2D.h>
#ifdef __2D__
#include <FEFunction2D.h>
#elif defined(__3D__)
#include <FEFunction3D.h>
#endif

#ifdef __2D__
#define TFEFunctionXD TFEFunction2D
#define TFESpaceXD TFESpace2D
#define GetFESpaceXD() GetFESpace2D()
#elif defined(__3D__)
#define TFEFunctionXD TFEFunction3D
#define TFESpaceXD TFESpace3D
#define GetFESpaceXD() GetFESpace3D()
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
    BrushWrapper(TCollection* brush_grid,
                 TCollection* parmoon_grid,
                 const ParameterDatabase& db);

    /// Destructor.
    ~BrushWrapper();

    ///Implicitely deletes copying and moving.

    /// Let Brush compute the source and sink terms which will
    /// be added to the right hand side of the CDR equations.
    std::vector<TFEFunctionXD*> sources_and_sinks();

    /**
     * Reset the fluid phase data in Brush, i.e., velocity, pressure, temperature
     * and all species concentrations.
     * @param[in] u1 The fluid velocity (first component)
     * @param[in] u2 The fluid velocity (second component)
     * @param[in] p The pressure.
     * @param[in] species All transported species, the first of which
     *            is supposed to be temperature. The following must be in that
     *            order in which Brush expects them.
     */
    void reset_fluid_phase(
        const TFEFunctionXD& u1,
    	  const TFEFunctionXD& u2,
#ifdef __3D__
    	  const TFEFunctionXD& u3,
#endif
        const TFEFunctionXD& p,
        std::vector<TFEFunctionXD*> species);

    /// Run the particle solver.
    void solve(double t_start, double t_end);

    /// Write moments of the particle distribution to a .vtk-file, which was
    /// specified in the database given to the constructor.
    void output(int& image_nr, double t);

  private:
    // Fill example depending data. This is like an "init" function,
    // called only by the constructor.
    //  eder_crystallizer_axis - Eder's ASA flow crystallizer example (axisymmetric).
    //  eder_crystallizer - Eder's ASA flow crystallizer example (non-axisymmetric, not working).
    //  wiedmeyer_crystallizer - 2017 Wiedmeyer Batch Crystallizer, 3d example
    // viscosity is a out-parameter which then goes to Brush
    void pick_example(const std::string& exmpl_name, double& viscosity);

    /// The ParMooN interface of Brush which is wrapped up by 'this'.
    /// TODO Ownership for the object is taken, but there is some trouble with the
    /// destructor.
    Brush::InterfacePM* interface_;

    /// The database for control parameters.
    ParameterDatabase db_;

    // BRUSH FE OBJECTS
    //grid
    TCollection* brush_grid_;
    //fe space (Q0 or P0)
    TFESpaceXD br_grid_space_;
    //fe functions for source and sink data (as returned from Brush)
    std::vector<TFEFunctionXD*> br_grid_source_fcts_;
    std::vector<std::vector<double>> br_grid_source_fcts_values_;
    //fe functions for parameter data (as handed over to Brush)
    std::vector<TFEFunctionXD*> br_grid_param_fcts_;
    std::vector<std::vector<double>> br_grid_param_fcts_values_;
    //fe functions for moments of the psd data (as returned from Brush) (used for output only)
    std::vector<TFEFunctionXD*> br_grid_psdmom_fcts_;
    std::vector<std::vector<double>> br_grid_psdmom_fcts_values_;

    // PARMOON FE OBJECTS
    //grid
    TCollection* parmoon_grid_;
    //fe space on parmoon_grid_, for return values to parmoon
    TFESpaceXD pm_grid_space_;
    //fe functions for source and sink data (as handed over to other ParMooN parts)
    std::vector<TFEFunctionXD*> pm_grid_source_fcts_;
    std::vector<std::vector<double>> pm_grid_source_fcts_values_;



    // OBJECTS FOR TRANSFERRING DATA BETWEEN BOTH 'WORLDS'
    GridTransferTool pm_to_brush_tool_;
    GridTransferTool brush_to_pm_tool_;


    // CONTROL STRUCTURE FOR PARAMETERS THAT GO TO BRUSH
    /// The vector of the names of the parameters. Must be of size
    /// parameter_spatial_dimension_ + 1 + parameter_n_specs_primary_ +
    /// parameter_n_specs_derived_. Where "+1" is for the fluid pressure.
    std::vector<std::string> parameter_function_names_;
    /// Number of primary species that go to Brush as parameters.
    /// ('primary': the species is represented as an fe function in ParMooN)
    size_t parameter_n_specs_primary_;
    /// Spatial dimension - must equal number of velocity components.
     size_t parameter_spatial_dimension_;
    /// Number of derived species that go to Brush as parameters.
    /// ('derived': the species is computed from primary species in ParMooN)
    size_t parameter_n_specs_derived_;
    /// The real valued functions which compute the derived parameters from
    /// the primary parameters. One for each derived parameter.
    std::vector<std::function<double(const std::vector<double>&)>>
        parameter_specs_derived_fcts_;


    // CONTROL STRUCTURE FOR SOURCES AND SINKS THAT COME FROM BRUSH
    // the following two informations come from the example
    std::vector<std::string> source_and_sink_function_names_;
    std::vector<Exmpl::SourceAndSinkTerms> source_and_sink_requests_;


    // OBJECTS FOR PROGRAM OUTPUT (their sole purpose is post-processing.)
    // Files into which Brush writes its output.
    std::ofstream moment_stats_file_;
    std::ofstream inflow_particles_file_;
    std::ofstream snapshot_particles_file_;
    std::ofstream outflow_particles_file_;
    // Moments of the PSD which can be calculated by Brush. Used only for visual-
    // ization with paraview, which is performed by the output_writer_.
#ifdef __2D__
    PostProcessing2D output_writer_;
#elif defined(__3D__)
    TOutput3D output_writer_;
    // cache some values that are needed for discontinuous .vtk output
    int disc_output_n_loc_verts;
    TVertex** disc_output_loc_verts;
#endif


};



#endif /* USER_PROJECTS_INC_BRUSHWRAPPER_H_ */
