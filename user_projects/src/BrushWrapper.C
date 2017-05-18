/*
 * BrushWrapper.C
 *
 *  Created on: Dec 14, 2016
 *      Author: bartsch
 */

#include <BrushWrapper.h>
//Brush code for the particles
#include <parmoon_interface.h>
#include <parmoon_data.h>

#include <algorithm>

namespace ASA_crystallizer
{
  #include <ASA_crystallizer.h>
}

void DirichletBoundaryConditions(int BdComp, double t, BoundCond &cond)
{
      cond = DIRICHLET;
}
void ZeroBoundaryValues(int BdComp, double Param, double &value)
{
      value = 0;
}

//Get the barycenter of a ParMooN grid cell.
std::valarray<double> center_point(const TBaseCell& cell)
{
#ifdef __3D__
  ErrThrow("Does center point calculation work in 3D?");
#endif
  std::valarray<double> p(0.0,3);

  unsigned int n_verts = cell.GetN_Vertices();
  for(unsigned int v = 0; v < n_verts; v++)
  {
    cell.GetVertex(v)->GetX();
    p[0] += cell.GetVertex(v)->GetX();
    p[1] += cell.GetVertex(v)->GetY();
  }
  p[0] /= n_verts;
  p[1] /= n_verts;

  return p;
}

void BrushWrapper::pick_example(int exmpl_code)
{
  switch(exmpl_code)
  {
    case 0:
    {
      using namespace ASA_crystallizer;
      //parameters information ('to Brush')
      parameter_spatial_dimension_ = parameter_spatial_dimension;
      parameter_n_specs_primary_ = parameter_n_specs_primary;
      parameter_n_specs_derived_ = parameter_n_specs_derived;
      parameter_function_names_ = parameter_term_names;
      parameter_specs_derived_fcts_ = parameter_specs_derived_fcts;
      //source and sink information ('from Brush')
      source_and_sink_function_names_ = source_and_sink_term_names;
      source_and_sink_requests_ = source_and_sink_fct_requests;

      // this stuff must be done for the Eder example, which has 4 different
      // parameter sets. It is not nice, but it does the trick.
      if(db_["sweep_file_depends_on_velocity_code"].is(true) && db_.contains("velocity_code"))
      {
        size_t parameter_set = db_["velocity_code"];
        //pick a sweep file according to the requested parameter set
        db_["sweep_file"].impose(
            Parameter("sweep_file", db_["sweep_file"].value_as_string()
                      + ".param_set_"+ std::to_string(parameter_set), ""));
        Output::info("Sweep File", "Due to choice of parameter set ", parameter_set,
                     " ('velocity_code'), the sweep file was changed"
                     " to ", db_["sweep_file"], ".");
      }

      break;
    }
    default:
      ErrThrow("Not implemented: example ", exmpl_code);
  }
}

// use 0-order elements, because that is the 'language' Brush speaks
int fe_order = 0;
BrushWrapper::BrushWrapper(TCollection* coll, const ParameterDatabase& db)
: db_(db),
  from_brush_grid_(coll),
  from_brush_space_(from_brush_grid_, (char*)"psd-moms", (char*)"psd-moms",
                    DirichletBoundaryConditions, fe_order, nullptr),
  moment_stats_file_(db_["out_part_moments_file"].get<std::string>()),
  outflow_particles_file_(db_["out_part_lists_file"].get<std::string>()),
  inflow_particles_file_(db_["in_part_lists_file"].get<std::string>()),
  output_writer_(db_)
{
  // get and store example specific information
  pick_example(db_["example"]);

  // set up Brushs ParMooN interface
  interface_ = new Brush::InterfacePM(
      db_["geo_file"], db_["third_dim_stretch"],
      db_["sweep_file"], " ",
      db_["therm_file"], db_["chem_file"],
      db_["max_sp_per_cell"], db_["max_m0_per_cell"]
  );

  // load the initial particle solution
  interface_->set_initial_particles(db_["init_partsol_file"]);

  // get those points where Brush expects values from ParMooN
  parameter_sample_points_ = interface_->get_cell_centers();

  // Set up and initialize FE functions which receive their data from Brush
  int fe_length = from_brush_space_.GetN_DegreesOfFreedom();
  std::vector<double> dummy_fe_values(fe_length, 0.0);

  // tell the interface which sample points parMooN is interested in
  // (for a start: the cell centers of all cells)
  int n_cells = coll->GetN_Cells();
  std::vector<std::valarray<double>> parmoon_cell_centers(n_cells,{0,0,0}); //z value will stay 0 in 2D
  for(int c = 0; c < n_cells; ++c)
  {
    parmoon_cell_centers[c] = center_point(*coll->GetCell(c));
  }
  interface_->set_output_sample_points(parmoon_cell_centers); //set ParMooN points in Brush

  // *** The source-and-sink functions contain the back-coupling from Brush to ParMooN
  int n_s_and_s_terms = source_and_sink_function_names_.size();
  source_and_sink_fcts_values_ =
      std::vector<std::vector<double>>(n_s_and_s_terms, dummy_fe_values);
  source_and_sink_fcts_.resize(n_s_and_s_terms);
  for(int s=0; s < n_s_and_s_terms; ++s)
  {
    source_and_sink_fcts_.at(s) =
        new TFEFunction2D(&from_brush_space_, (char*) &source_and_sink_function_names_.at(s),
                          (char*) " ", &source_and_sink_fcts_values_[s].at(0), fe_length);
  }

  // *** The moments functions are only used for output and visualization.
  pd_moments_values_ = std::vector<std::vector<double>>(3, dummy_fe_values);

  pd_moments_.resize(3);
  pd_moments_.at(0) = new TFEFunction2D(&from_brush_space_,
        (char*)"pd-m0", (char*)"", &pd_moments_values_[0].at(0), fe_length);
  pd_moments_.at(1) = new TFEFunction2D(&from_brush_space_,
        (char*)"pd-m1", (char*)"", &pd_moments_values_[1].at(0), fe_length);
  pd_moments_.at(2) = new TFEFunction2D(&from_brush_space_,
        (char*)"pd-m2", (char*)"", &pd_moments_values_[2].at(0), fe_length);

  //add the moments functions to the output writer
  output_writer_.add_fe_function(pd_moments_[0]);
  output_writer_.add_fe_function(pd_moments_[1]);
  output_writer_.add_fe_function(pd_moments_[2]);

  //force an update of all ensemble stats
  interface_->update_stats();

  // store moments m0, m1 and m2.
  interface_->fetch_moment(0, &pd_moments_values_[0].at(0));
  interface_->fetch_moment(1, &pd_moments_values_[1].at(0));
  interface_->fetch_moment(2, &pd_moments_values_[2].at(0));

  //open the file streams
  if(moment_stats_file_.is_open())
    moment_stats_file_.close();
  moment_stats_file_.open(db_["out_part_moments_file"].get<std::string>(), std::ofstream::out);
  if(outflow_particles_file_.is_open())
    outflow_particles_file_.close();
  outflow_particles_file_.open(db_["out_part_lists_file"].get<std::string>(), std::ofstream::out);
  if(inflow_particles_file_.is_open())
    inflow_particles_file_.close();
  inflow_particles_file_.open(db_["in_part_lists_file"].get<std::string>(), std::ofstream::out);

  // Finally write a nice header into the particle stats file,
  // so that it can be filled with data from now on.
  interface_->write_headers(moment_stats_file_, outflow_particles_file_, inflow_particles_file_);
}

BrushWrapper::~BrushWrapper()
{
  //FIXME THERE IS A BUG (CONNECTED TO THE MOON-GEOMETRY) WHEN
  // CALLING THE FOLLOWING DESTRUCTOR!!!
  //delete interface_;

  moment_stats_file_.close();
  outflow_particles_file_.close();
  inflow_particles_file_.close();

  for (auto f : pd_moments_)
  {
    delete f;
  }
  for (auto f : source_and_sink_fcts_)
  {
    delete f;
  }
}

std::vector<const TFEFunction2D*> BrushWrapper::sources_and_sinks()
{

  // Let the interface recalculate the function values.
  interface_->fetch_sources_and_sinks(
      source_and_sink_requests_,
      source_and_sink_fcts_values_
        );

  return source_and_sink_fcts_;
}

void BrushWrapper::reset_fluid_phase(
    const TFEVectFunct2D& u,
    const TFEFunction2D& p,
    std::vector<const TFEFunction2D*> species
    )
{
  //check input
  if(species.size() != parameter_n_specs_primary_)
   ErrThrow("Incorrect number of species fe functions given.");
  if(u.GetN_Components() != (int) this->parameter_spatial_dimension_)
    ErrThrow("Incorrect number of velocity components given.");

  size_t n_points = parameter_sample_points_.size();

  Brush::DataPM pm_data(parameter_function_names_, n_points);

  //For the velocity we need two xtra fe functions
  TFEFunction2D* u0_fe = u.GetComponent(0);
  TFEFunction2D* u1_fe = u.GetComponent(1);

  //prepare a vector of pointers to FEFunctions
  size_t dim = parameter_spatial_dimension_;
  size_t n_specs_primary = parameter_n_specs_primary_;
  size_t n_specs_derived = parameter_n_specs_derived_;

  std::vector<const TFEFunction2D*> fe_functs(dim + 1 + n_specs_primary);
  fe_functs[0] = u0_fe;
  fe_functs[1] = u1_fe;
  fe_functs[2] = &p;
  std::copy( species.begin() , species.end(), &fe_functs[3]);

  // This is a good place to cache containing cells of the output points
  if(reset_fluid_phase_cheats_.empty())
  {
    cache_output_point_containing_cells(*p.GetFESpace2D());
  }

  //loop over all evaluation points (i.e., Brushs cell midpoints)
  for (size_t point = 0 ; point < n_points ; ++point)
  {
    //x and y value of the current output sample point
    double x = parameter_sample_points_[point][0];
    double y = parameter_sample_points_[point][1];

    size_t data_set_size = dim + 1 + n_specs_primary + n_specs_derived;
    std::vector<double> data_set(data_set_size , 0.0);

    // point evaluations of all fe functions
    for(size_t i = 0 ; i < dim + 1 + n_specs_primary ; ++i)
    {
      double eval[3]; // will include differentials, thus length is '3')
      fe_functs[i]->FindGradient(x, y, eval, reset_fluid_phase_cheats_.at(point));
      data_set[i] = eval[0];
    }
    // computation for all derived species concentrations
    for(size_t i = dim + 1 + n_specs_primary; i < data_set_size; ++i)
    {
      size_t index = i - (dim + 1 + n_specs_primary);
      double f = parameter_specs_derived_fcts_[index](data_set);
      data_set[i] = f;
    }

    //now the data set is finished and can be put to the data vector
    pm_data.add_data_set(data_set);
  }

  delete u0_fe;
  delete u1_fe;

  //give the data to Brush
  interface_->reset_fluid_phase(pm_data);

}


void BrushWrapper::solve(double t_start, double t_end)
{
  //now call the solver
  size_t n_steps = db_["n_solves_per_time_step"];

  int rs = db_["random_seed"];
  //db_["random_seed"] = rs + 1; //count up the rng for use in the next step

  interface_->run_particle_phase(t_start, t_end, n_steps, rs);

  // Updating stats and fetching moments is only relevant for
  // visualization and the output
  // TODO Make this runtime-controllable.
  interface_->update_stats();
  interface_->fetch_moment(0, &pd_moments_values_[0].at(0));
  interface_->fetch_moment(1, &pd_moments_values_[1].at(0));
  interface_->fetch_moment(2, &pd_moments_values_[2].at(0));

}


void BrushWrapper::output(double t)
{
  output_writer_.write(t);

  interface_->write_particle_stats(t, moment_stats_file_);

  interface_->write_outlet_particle_list(outflow_particles_file_);

  interface_->write_inlet_particle_list(inflow_particles_file_);

}


void BrushWrapper::cache_output_point_containing_cells(const TFESpace2D& one_space)
{

  for (auto p : parameter_sample_points_)
  {
    double x = p[0];
    double y = p[1];
    std::vector<int> cheat;
    const TCollection& coll = *one_space.GetCollection();
    int N_Cells = coll.GetN_Cells();
    for(int i=0;i<N_Cells;i++)
    {
      TBaseCell* cell = coll.GetCell(i);
      if(cell->PointInCell(x,y))
      {
        cheat.push_back(i);
      }
    }
    reset_fluid_phase_cheats_.push_back(cheat);
  }
}
