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

namespace Axisymmetric_ASA_crystallizer
{
  #include <Axisymmetric_ASA_Crystallizer.h>
}

//CB DEBUG
#include "fstream"
#include "cmath"
//END DEBUG

void DirichletBoundaryConditions(int BdComp, double t, BoundCond &cond)
{
      cond = DIRICHLET;
}
void ZeroBoundaryValues(int BdComp, double Param, double &value)
{
      value = 0;
}



void BrushWrapper::pick_example(int exmpl_code)
{
  switch(exmpl_code)
  {
    case 0:
    {//Axisymmetric ASA crystallizer
      using namespace Axisymmetric_ASA_crystallizer;
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
    case 1:
    {//non-axisymmetric (==false) ASA crystallizer
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

/**
 *
 * brush_grid is the grid where the brush simulation should run on it should
 * (but need not) stand in a refinement relation to the grid of the parmoon space
 *
 * parmoon_space is the space of the transported species which is used in parmoon
 * it gets re-constructed for this wrapper, because there is currently no copy
 * constructor in the fespace - THIS MEANS THAT ALL RELEVANT GLOBAL PARAMETERS
 * MUST NOT CHANGE
 */
BrushWrapper::BrushWrapper(TCollection* brush_grid,
                           TCollection* parmoon_grid,
                           const ParameterDatabase& db,
						   bool axisymmetric)
: //database
  db_(db),
  //brush fe
  brush_grid_(brush_grid),
  br_grid_space_(brush_grid_, (char*)"brush-space", (char*)"Space for the "
      "direct representation of Brush's 0 order fe functions in ParMooN.",
       DirichletBoundaryConditions, 0, nullptr),
  //parmoon fe
  parmoon_grid_(parmoon_grid),
  pm_grid_space_(parmoon_grid_,(char*)"parmoon-space", (char*)"Space for the "
                 "ParMooN-representation of the functions which Brush returns",
                 DirichletBoundaryConditions , 0,
                 nullptr),
  //grid trafo and control
  //pm_to_brush_tool_(&br_grid_space_, GridTransferType::Interpolation),
  //brush_to_pm_tool_(&pm_grid_space_, GridTransferType::Interpolation),
  pm_to_brush_tool_(&br_grid_space_, GridTransferType::MultiGrid),
  brush_to_pm_tool_(&pm_grid_space_, GridTransferType::MultiGrid),
  parameter_n_specs_primary_(0),   // gets reset by pick_example
  parameter_spatial_dimension_(0), // gets reset by pick_example
  parameter_n_specs_derived_(0),   // gets reset by pick_example
  //output files and vtk object
  moment_stats_file_(db_["out_part_moments_file"].get<std::string>()),
  outflow_particles_file_(db_["out_part_lists_file"].get<std::string>()),
  inflow_particles_file_(db_["in_part_lists_file"].get<std::string>()),
  output_writer_(db_)
{

  // get and store example specific information
  pick_example(db_["example"]);

  //write the brush grid to a file, which can be read by Brush
  brush_grid_->writeMesh("brush_mesh.mesh", 2);

  // set up Brushs ParMooN interface
  interface_ = new Brush::InterfacePM(
      "brush_mesh.mesh", db_["third_dim_stretch"],
      db_["sweep_file"], " ",
      db_["therm_file"], db_["chem_file"],
      db_["max_sp_per_cell"], db_["max_m0_per_cell"],
	  axisymmetric
  );

  // check if the numbering of cells in the brush grid is the same in Brush and ParMooN
  std::vector<std::valarray<double>> centers = interface_->get_cell_centers();
  for(size_t brush_cell = 0 ; brush_cell < centers.size()  ; ++brush_cell)
  {
    std::valarray<double> point = centers.at(brush_cell);
    double x = point[0];
    double y = point[1];
    //Output::print("Remote cell ", brush_cell, " midpoint (", point[0],",",point[1],")");
    std::vector<int> found_in;
    for(int loc_cell = 0 ; loc_cell < brush_grid_->GetN_Cells() ;++loc_cell)
    {
      TBaseCell* cell = brush_grid_->GetCell(loc_cell);

      if( cell->PointInCell(x,y) )
      {
        found_in.push_back(loc_cell);
        //Output::print("Midpoint of remote cell ", brush_cell, " found in local cell ", loc_cell );
      }
    }
    //check the vector found_in - point found in and only found in the right local cell?
    if(found_in.size() == 0)
      ErrThrow("Did not find mid point of Brush cell ", brush_cell, "in any local cell.");
    if(found_in.size() > 1)
      ErrThrow("Found mid point of Brush cell ", brush_cell, " in ",
               found_in.size(), " local cells , which is too much." );
    if(found_in.at(0) != (int) brush_cell)
      ErrThrow("Found mid point of Brush cell ", brush_cell, " in "
               "local cell ", found_in.at(0), " which is unexpected.");
    // if these checks passed, everything is fine.
  }

  // load the initial particle solution
  interface_->set_initial_particles(db_["init_partsol_file"]);

  // *** The source-and-sink functions contain the back-coupling from Brush to ParMooN
  int n_s_and_s_terms = source_and_sink_function_names_.size(); //their number
  //...here their ParMooN-form gets initialized
  int pm_space_fe_length = pm_grid_space_.GetN_DegreesOfFreedom();
  std::vector<double> dummy_fe_values(pm_space_fe_length, 0.0);
  pm_grid_source_fcts_values_ =
      std::vector<std::vector<double>>(n_s_and_s_terms, dummy_fe_values);
  pm_grid_source_fcts_.resize(n_s_and_s_terms);
  for(int s=0; s < n_s_and_s_terms; ++s)
  {
    pm_grid_source_fcts_.at(s) =
        new TFEFunction2D(&pm_grid_space_,source_and_sink_function_names_[s].c_str(),
                          (char*) "Brush function transferred to ParMooN grid.",
                          &pm_grid_source_fcts_values_[s].at(0), pm_space_fe_length);
  }
  //...and here is their Brush-form
  int br_space_fe_length = br_grid_space_.GetN_DegreesOfFreedom();
  dummy_fe_values.resize(br_space_fe_length, 0.0);
  br_grid_source_fcts_values_ =
      std::vector<std::vector<double>>(n_s_and_s_terms, dummy_fe_values);
  br_grid_source_fcts_.resize(n_s_and_s_terms);
  for(int s=0; s < n_s_and_s_terms; ++s)
  {
    br_grid_source_fcts_.at(s) =
        new TFEFunction2D(&br_grid_space_, (char*) &source_and_sink_function_names_.at(s),
                          (char*) "Brush function representation on Brush grid.",
                          &br_grid_source_fcts_values_[s].at(0), br_space_fe_length);
  }

  // *** The parameter functions contain the fluid information that will be given to Brush.
  int n_param_fcts = parameter_spatial_dimension_ + 1 + parameter_n_specs_primary_; // +1 is for the pressure
  br_grid_param_fcts_values_ =
      std::vector<std::vector<double>>(n_param_fcts, dummy_fe_values);
  br_grid_param_fcts_.resize(n_param_fcts);
  for(int p = 0; p < n_param_fcts; ++p)
  {
    br_grid_param_fcts_.at(p) =
        new TFEFunction2D(&br_grid_space_, (char*) &parameter_function_names_.at(p),
                          (char*) "ParMooN function representation on Brush grid.",
                          &br_grid_param_fcts_values_[p].at(0), br_space_fe_length);
  }

  // *** The moments functions on the brush grid are only used for output and visualization.
  br_grid_psdmom_fcts_values_ = std::vector<std::vector<double>>(3, dummy_fe_values);

  br_grid_psdmom_fcts_.resize(3);
  br_grid_psdmom_fcts_.at(0) = new TFEFunction2D(&br_grid_space_,
        (char*)"pd-m0", (char*)"", &br_grid_psdmom_fcts_values_[0].at(0), br_space_fe_length);
  br_grid_psdmom_fcts_.at(1) = new TFEFunction2D(&br_grid_space_,
        (char*)"pd-m1", (char*)"", &br_grid_psdmom_fcts_values_[1].at(0), br_space_fe_length);
  br_grid_psdmom_fcts_.at(2) = new TFEFunction2D(&br_grid_space_,
        (char*)"pd-m2", (char*)"", &br_grid_psdmom_fcts_values_[2].at(0), br_space_fe_length);

  if(db_["output_write_vtk"].is(true))
  {
	  //add the moments functions to the output writer
	  output_writer_.add_fe_function(br_grid_psdmom_fcts_[0]);
	  output_writer_.add_fe_function(br_grid_psdmom_fcts_[1]);
	  output_writer_.add_fe_function(br_grid_psdmom_fcts_[2]);
	  // store moments m0, m1 and m2.
	  interface_->update_stats();
	  interface_->fetch_moment(0, &br_grid_psdmom_fcts_values_[0].at(0));
	  interface_->fetch_moment(1, &br_grid_psdmom_fcts_values_[1].at(0));
	  interface_->fetch_moment(2, &br_grid_psdmom_fcts_values_[2].at(0));
  }

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

  //CB DEBUG
  remove("mass_balance.csv"); //file from old run, clean away
  std::ofstream myfile;
  myfile.open ("mass_balance.csv",std::ios_base::app);
  myfile << "t" << ",";
  myfile << "cell #" << ",";
  myfile << "z" << ",";
  myfile << "r" << ",";
  myfile << "ASA crys [kg/m^3]" << ",";
  myfile << "ASA diss [kg/m^3]" << "\n";
  myfile.close();
  //END DEBUG
}

BrushWrapper::~BrushWrapper()
{
  //FIXME THERE IS A BUG (CONNECTED TO THE MOON-GEOMETRY) WHEN
  // CALLING THE FOLLOWING DESTRUCTOR!!!
  //delete interface_;

  moment_stats_file_.close();
  outflow_particles_file_.close();
  inflow_particles_file_.close();

  for (auto f : br_grid_source_fcts_)
    delete f;

  for (auto f : br_grid_param_fcts_)
    delete f;

  for (auto f : br_grid_psdmom_fcts_)
  {
	  if(f != NULL)
		  delete f;
  }

  for (auto f : pm_grid_source_fcts_)
    delete f;

}

std::valarray<double> center_point_calc(const TBaseCell& cell);

std::vector<TFEFunction2D*> BrushWrapper::sources_and_sinks()
{

  // The sources and sinks from Brush are picked and stored as function values
  // on the Brush grid.
  // In order for this to work, two class invariants must be fulfilled:
  //  1) Brush grid and the cells in Brush must be in the same order
  //  2) The Brush grid source function must be of order 0
  interface_->fetch_sources_and_sinks(
      source_and_sink_requests_,
      br_grid_source_fcts_values_
        );

  int n_source_fcts = source_and_sink_function_names_.size();
  for(int f = 0 ; f < n_source_fcts ; ++f)
  {
    brush_to_pm_tool_.transfer(*br_grid_source_fcts_.at(f),
                               *pm_grid_source_fcts_.at(f),
                               pm_grid_source_fcts_values_.at(f));
  }


//  //CB DEBUG
//  double out [3] = {0,0,0};
//  double diff=0;
//  br_grid_source_fcts_[1]->GetMassAndMean(out, true, 'x');
//  Output::print("Sinks on Brush grid: ", out[0], " ", out[1], " ", out[2]);
//  diff = out[2];
//  pm_grid_source_fcts_[1]->GetMassAndMean(out, true, 'x');
//  Output::print("Sinks on ParMooN grid: ", out[0], " ", out[1], " ", out[2]);
//  diff -= out[2];
//  //END DEBUG



 return pm_grid_source_fcts_;
}

//CB DEBUG
std::valarray<double> center_point_calc(const TBaseCell& cell)
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
//END DEBUG

void BrushWrapper::reset_fluid_phase(
		const TFEFunction2D& u1,
		const TFEFunction2D& u2,
		const TFEFunction2D& p,
		std::vector<TFEFunction2D*> species
)
{
  //check input
  if(species.size() != parameter_n_specs_primary_)
   ErrThrow("Incorrect number of species fe functions given.");

  size_t n_data_sets = brush_grid_->GetN_Cells();

  Brush::DataPM pm_data(parameter_function_names_, n_data_sets);

  //prepare a vector of pointers to FEFunctions
  size_t dim = parameter_spatial_dimension_;
  size_t n_specs_primary = parameter_n_specs_primary_;
  size_t n_specs_derived = parameter_n_specs_derived_;

  std::vector<const TFEFunction2D*> fe_fcts_original(dim + 1 + n_specs_primary);
  fe_fcts_original[0] = &u1;
  fe_fcts_original[1] = &u2;
  //if(dim==3) fe_fcts_original[2] = &u3; //in 3D case...
  fe_fcts_original[dim] = &p;
  std::copy( species.begin() , species.end(), &fe_fcts_original[dim + 1]);

  // transfer all 'original' (ParMooN) functions to the Brush grid
  for(size_t f = 0; f<fe_fcts_original.size() ;++f)
  {
    pm_to_brush_tool_.transfer(
        *fe_fcts_original[f],*br_grid_param_fcts_[f], br_grid_param_fcts_values_[f]);
  }

//  //CB DEBUG
//  double out [3] = {0,0,0};
//  fe_fcts_original[4]->GetMassAndMean(out, true, 'x');
//  Output::print("Values on ParMooN grid: ", out[0], " ", out[1], " ", out[2]);
//  br_grid_param_fcts_[4]->GetMassAndMean(out, true, 'x');
//  Output::print("Values on Brush grid: ", out[0], " ", out[1], " ", out[2]);
//  //END DEBUG

  //loop over all evaluation points (i.e., Brushs cell midpoints)
  for (size_t d = 0 ; d < n_data_sets ; ++d)
  {

    size_t data_set_size = dim + 1 + n_specs_primary + n_specs_derived;
    std::vector<double> data_set(data_set_size , 0.0);

    // point evaluations of all fe functions
    for(size_t i = 0 ; i < dim + 1 + n_specs_primary ; ++i)
    {
      data_set[i] = br_grid_param_fcts_values_[i][d];
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

  //give the data to Brush
  interface_->reset_fluid_phase(pm_data);

  //CB DEBUG
  double del = fabs(TDatabase::TimeDB->CURRENTTIME - std::round(TDatabase::TimeDB->CURRENTTIME));
  if(del < 1e-3)
  {
  interface_->update_stats();
  interface_->fetch_moment(1, &br_grid_psdmom_fcts_values_[1].at(0));
  std::ofstream myfile;
  myfile.open ("mass_balance.csv",std::ios_base::app);

  for(int c=0; c < brush_grid_->GetN_Cells(); ++c)
  {

	  auto cntr = center_point_calc(*brush_grid_->GetCell(c));
	  myfile << TDatabase::TimeDB->CURRENTTIME << ",";
	  myfile << c << ",";
	  myfile << cntr[0] << ",";
	  myfile << cntr[1] << ",";
	  //Brush ASA values
	  auto bgn_ind= br_grid_space_.GetBeginIndex();
	  auto dofs = br_grid_space_.GetGlobalNumbers();
	  auto dof = dofs[bgn_ind[c]];
	  myfile << br_grid_psdmom_fcts_values_[1].at(dof) << ",";
	  //ParMooN ASA values
	  bgn_ind= br_grid_param_fcts_[4]->GetFESpace2D()->GetBeginIndex();
	  dofs = br_grid_param_fcts_[4]->GetFESpace2D()->GetGlobalNumbers();
	  dof = dofs[bgn_ind[c]];
	  double val_in_kg = br_grid_param_fcts_[4]->GetValues()[dof] * 0.18016;
	  myfile << val_in_kg << "\n";
  }
  myfile.close();
  }
  //END DEBUG

}


void BrushWrapper::solve(double t_start, double t_end)
{
  //now call the solver
  size_t n_steps = db_["n_solves_per_time_step"];

  int rs = db_["random_seed"];
  db_["random_seed"] = rs + 1; //count up the rng for use in the next step

  interface_->run_particle_phase(t_start, t_end, n_steps, rs);

}



void BrushWrapper::output(double t)
{
	// Updating stats and fetching moments is only relevant for
	// visualization and the output
	if(db_["output_write_vtk"].is(true))
	{
		interface_->update_stats();
		interface_->fetch_moment(0, &br_grid_psdmom_fcts_values_[0].at(0));
		interface_->fetch_moment(1, &br_grid_psdmom_fcts_values_[1].at(0));
		interface_->fetch_moment(2, &br_grid_psdmom_fcts_values_[2].at(0));
		output_writer_.write(t);
	}

  interface_->write_particle_stats(t, moment_stats_file_);

  interface_->write_outlet_particle_list(outflow_particles_file_);

  interface_->write_inlet_particle_list(inflow_particles_file_);

}
