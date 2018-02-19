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
#include <cmath>

namespace ASA_crystallizer
{
  #include <ASA_crystallizer.h>
}

namespace Axisymmetric_ASA_crystallizer
{
  #include <Axisymmetric_ASA_Crystallizer.h>
}

namespace Wiedmeyer_Batch_crystallizer
{
#include "TNSE_3D/WiedmeyerBatchCrystallizer.h"
}

// Hashing function for strings, allows switch over string.
// https://stackoverflow.com/questions/16388510/evaluate-a-string-with-a-switch-in-c
constexpr unsigned int string_hash(const char* str, int h = 0)
{
    return !str[h] ? 5381 : (string_hash(str, h+1) * 33) ^ str[h];
}

void DirichletBoundaryConditions(int BdComp, double t, BoundCond &cond)
{
      cond = DIRICHLET;
}
void DirichletBoundaryConditions_3D(double x, double y, double z, BoundCond &cond)
{
      cond = DIRICHLET;
}
void ZeroBoundaryValues_3D(double x, double y, double z, double &value)
{
      value = 0;
}

std::valarray<double> simplex_barycenter(const TBaseCell& cell);

void BrushWrapper::pick_example(const std::string& exmpl_name,
                                double& viscosity)
{
  switch(string_hash(exmpl_name.c_str()))
  {
    case string_hash("eder_crystallizer_axis"):
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
      viscosity = Physics::mu;

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
    case string_hash("eder_crystallizer"):
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

      viscosity = 0.01; //dummy

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
    case string_hash("wiedmeyer_crystallizer"):
    {
      using namespace Wiedmeyer_Batch_crystallizer;

      parameter_spatial_dimension_ = BrushInfo::parameter_spatial_dimension;
      parameter_n_specs_primary_ = BrushInfo::parameter_n_specs_primary;
      parameter_n_specs_derived_ = BrushInfo::parameter_n_specs_derived;
      parameter_function_names_ = BrushInfo::parameter_term_names;
      parameter_specs_derived_fcts_ = BrushInfo::parameter_specs_derived_fcts;
      //source and sink information ('from Brush')
      source_and_sink_function_names_ = BrushInfo::source_and_sink_term_names;
      source_and_sink_requests_ = BrushInfo::source_and_sink_fct_requests;

      //TODO this might change in case of de-dimensionalization
      viscosity = FluidProperties::eta/FluidProperties::rho;

      break;
    }
    default:
      ErrThrow("Not implemented: example ", exmpl_name);
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
                           const ParameterDatabase& db)
: //database
  db_(db),
  //brush fe
  brush_grid_(brush_grid),
#ifdef __2D__
  br_grid_space_(brush_grid_, "brush-space", "Space for the "
      "direct representation of Brush's 0 order fe functions in ParMooN.",
       DirichletBoundaryConditions, 0, nullptr),
#elif defined(__3D__)
  br_grid_space_(brush_grid_, "brush-space", "Space for the "
      "direct representation of Brush's 0 order fe functions in ParMooN.",
      DirichletBoundaryConditions_3D, 0),
#endif
  //parmoon fe
  parmoon_grid_(parmoon_grid),
#ifdef __2D__
  pm_grid_space_(parmoon_grid_,"parmoon-space", "Space for the "
                 "ParMooN-representation of the functions which Brush returns",
                 DirichletBoundaryConditions , 0,
                 nullptr),
#elif defined(__3D__)
  pm_grid_space_(parmoon_grid_,"parmoon-space", "Space for the "
                 "ParMooN-representation of the functions which Brush returns",
                 DirichletBoundaryConditions_3D, 0),
#endif
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
  inflow_particles_file_(db_["in_part_lists_file"].get<std::string>())
#ifdef __2D__
  ,output_writer_(db_)
#elif defined(__3D__)
  , output_writer_(1,3,0,0,brush_grid) //...bunch of magic numbers for control of the (almost deprecated) 3DOutput object
#endif
{

  // get and store example specific information
  double viscosity=0;
  pick_example(db_["example"], viscosity);

  //write the brush grid to a file, which can be read by Brush
  std::string out_dir(db_["output_vtk_directory"].value_as_string());
  std::size_t pos = out_dir.find("VTK");// Lil hack - the directory above 'VTK',
  out_dir = out_dir.substr(0,pos); 		// ...should be the general output directory.
#ifdef __2D
  int dim = 2;
#elif defined(__3D__)
  int dim = 3;
#endif
  std::string brush_grid_file = "./brush_grid.mesh";
  brush_grid_->writeMesh(brush_grid_file.c_str(), dim);
#ifdef __2D__
  double third_dim_stretch = 0;
  if(db_.contains("third_dim_stretch"))
	  third_dim_stretch=db_["third_dim_stretch"];
  bool axisymmetric = (db_["example"] == eder_crystallizer_axis);
#endif

  std::vector<std::valarray<double>> cell_centers = {};
// it seems this is unnecessary, tetgen preserves the input order of the medit file.
//  for(int c=0; c < brush_grid_->GetN_Cells(); ++c)
//  {
//    cell_centers.push_back(simplex_barycenter(*brush_grid_->GetCell(c)));
//    Output::print("Cell ", c, " center ", cell_centers[c][0], " ", cell_centers[c][1], " ", cell_centers[c][2]);
//  }

  // set up Brush's ParMooN interface
  interface_ = new Brush::InterfacePM(
      brush_grid_file, cell_centers,
      db_["sweep_file"], " ",
      db_["therm_file"], db_["chem_file"],
      db_["max_sp_per_cell"], db_["max_m0_per_cell"],
      db_["example"],
      viscosity
      //,db_["coagulation_parameter"]
  );

  Output::print("Setting up Brush::InterfacePM SUCCESS.");
//  Output::print("Performing cell checks");
//
//  // check if the numbering of cells in the brush grid is the same in Brush and ParMooN
//  // This scales quadratically in the number of cells,
//  // and should therefore only run in test scenarios.
//  // BEGIN CHECK
//  std::vector<std::valarray<double>> centers = interface_->get_cell_centers();
//  for(size_t brush_cell = 0 ; brush_cell < centers.size()  ; ++brush_cell)
//  {
//    std::valarray<double> point = centers.at(brush_cell);
//    double x = point[0];
//    double y = point[1];
//    double z;
//    if(dim == 3)
//     z = point[2];
//    //Output::print("Remote cell ", brush_cell, " midpoint (", point[0],",",point[1],")");
//    std::vector<int> found_in;
//    for(int loc_cell = 0 ; loc_cell < brush_grid_->GetN_Cells() ;++loc_cell)
//    {
//      TBaseCell* cell = brush_grid_->GetCell(loc_cell);
//
//      if(dim == 2)
//      {
//        if(cell->PointInCell(x,y))
//          found_in.push_back(loc_cell);
//      }
//      else if(dim == 3)
//      {
//        if( cell->PointInCell(x,y,z) )
//          found_in.push_back(loc_cell);
//      }
//    }
//    //check the vector found_in - point found in and only found in the right local cell?
//    if(found_in.size() == 0)
//      ErrThrow("Did not find mid point of Brush cell ", brush_cell, "in any local cell.");
//    if(found_in.size() > 1)
//      ErrThrow("Found mid point of Brush cell ", brush_cell, " in ",
//               found_in.size(), " local cells , which is too much." );
//    if(found_in.at(0) != (int) brush_cell)
//      ErrThrow("Found mid point of Brush cell ", brush_cell, " in "
//               "local cell ", found_in.at(0), " which is unexpected.");
//    // if these checks passed, everything is fine.
//  }
//  Output::print("Performing cell checks SUCCESS.");
//  // END CHECK

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
        new TFEFunctionXD(&pm_grid_space_,source_and_sink_function_names_[s].c_str(),
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
        new TFEFunctionXD(&br_grid_space_, (char*) &source_and_sink_function_names_.at(s),
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
        new TFEFunctionXD(&br_grid_space_, (char*) &parameter_function_names_.at(p),
                          (char*) "ParMooN function representation on Brush grid.",
                          &br_grid_param_fcts_values_[p].at(0), br_space_fe_length);
  }

  // *** The moments functions on the brush grid are only used for output and visualization.
  br_grid_psdmom_fcts_values_ = std::vector<std::vector<double>>(3, dummy_fe_values);

  br_grid_psdmom_fcts_.resize(3);
  br_grid_psdmom_fcts_.at(0) = new TFEFunctionXD(&br_grid_space_,
        (char*)"pd-m0", (char*)"", &br_grid_psdmom_fcts_values_[0].at(0), br_space_fe_length);
  br_grid_psdmom_fcts_.at(1) = new TFEFunctionXD(&br_grid_space_,
        (char*)"pd-m1", (char*)"", &br_grid_psdmom_fcts_values_[1].at(0), br_space_fe_length);
  br_grid_psdmom_fcts_.at(2) = new TFEFunctionXD(&br_grid_space_,
        (char*)"pd-m2", (char*)"", &br_grid_psdmom_fcts_values_[2].at(0), br_space_fe_length);

  if(db_["output_write_vtk"].is(true))
  {
#ifdef __2D__
	  //add the moments functions to the output writer
	  output_writer_.add_fe_function(br_grid_psdmom_fcts_[0]);
	  output_writer_.add_fe_function(br_grid_psdmom_fcts_[1]);
	  output_writer_.add_fe_function(br_grid_psdmom_fcts_[2]);
#elif defined(__3D__)
	  output_writer_.AddFEFunction(br_grid_psdmom_fcts_[0]);
    output_writer_.AddFEFunction(br_grid_psdmom_fcts_[1]);
    output_writer_.AddFEFunction(br_grid_psdmom_fcts_[2]);

    //CB EXPERIMENTAL - this is for discontinuous .vtk output
    disc_output_n_loc_verts=0;
    for(int i=0;i<brush_grid_->GetN_Cells();i++)
    {
      TBaseCell* cell = brush_grid_->GetCell(i);
      disc_output_n_loc_verts += cell->GetN_Vertices();
    }
    disc_output_loc_verts =new TVertex*[disc_output_n_loc_verts];
    int N_=0;

    for(int i=0;i<brush_grid_->GetN_Cells();i++)
    {
      TBaseCell* cell = brush_grid_->GetCell(i);
      int k=cell->GetN_Vertices();
      for(int j=0;j<k;j++)
      {
        disc_output_loc_verts[N_]=cell->GetVertex(j);
        N_++;
      }
    }
    //END EXPERIMENTAL

#endif
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

  //CB DEBUG - Total mass control file, write header.
  std::string mass_bal_file = out_dir + "mass_balance.csv";
  remove(mass_bal_file.c_str()); //file from old run, clean away
  std::ofstream myfile;
  myfile.open (mass_bal_file.c_str(),std::ios_base::app);
  myfile << "t" << ",";
  myfile << "cell #" << ",";
  myfile << "x" << ",";
  myfile << "y" << ",";
  myfile << "z" << ",";
  myfile << "Volume [m^3]" << ",";
  myfile << "POTASHALUM crys [kg/m^3]" << ",";
  //myfile << "ASA diss [kg/m^3]" << ",";
  myfile << "x velo [m/s] " << ",";
  myfile << "y velo [m/s] " << ",";
  myfile << "z velo [m/s] " << "\n";
  myfile.close();
  //END DEBUG
}

BrushWrapper::~BrushWrapper()
{
  delete interface_;

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

std::vector<TFEFunctionXD*> BrushWrapper::sources_and_sinks()
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

 return pm_grid_source_fcts_;
}

std::valarray<double> simplex_barycenter(const TBaseCell& cell)
{
  std::valarray<double> p(0.0,3);

  unsigned int n_verts = cell.GetN_Vertices();
#ifdef __2D__
  if(n_verts != 3)
    Output::warn("Calling simplex_barycenter in 2d on non-triangle!");
#elif defined (__3D__)
  if(n_verts != 4)
    Output::warn("Calling simplex_barycenter in 3d on non-tetrahedron!");
#endif

  for(unsigned int v = 0; v < n_verts; v++)
  {
    cell.GetVertex(v)->GetX();
    p[0] += cell.GetVertex(v)->GetX();
    p[1] += cell.GetVertex(v)->GetY();
#ifdef __3D__
    p[2] += cell.GetVertex(v)->GetZ();
#endif
  }
  p[0] /= n_verts;
  p[1] /= n_verts;
  p[2] /= n_verts;

  return p;
}

void BrushWrapper::reset_fluid_phase(
		const TFEFunctionXD& u1,
		const TFEFunctionXD& u2,
#ifdef __3D__
    const TFEFunctionXD& u3,
#endif
		const TFEFunctionXD& p,
		std::vector<TFEFunctionXD*> species
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

  std::vector<const TFEFunctionXD*> fe_fcts_original(dim + 1 + n_specs_primary);
  fe_fcts_original[0] = &u1;
  fe_fcts_original[1] = &u2;
#ifdef __3D__
  fe_fcts_original[2] = &u3;
#endif
  fe_fcts_original[dim] = &p;
  std::copy( species.begin() , species.end(), &fe_fcts_original[dim + 1]);

  // transfer all 'original' (ParMooN) functions to the Brush grid
  for(size_t f = 0; f<fe_fcts_original.size() ;++f)
  {
    pm_to_brush_tool_.transfer(
        *fe_fcts_original[f],*br_grid_param_fcts_[f], br_grid_param_fcts_values_[f]);
  }

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

}


void BrushWrapper::solve(double t_start, double t_end)
{
  //now call the solver
  size_t n_steps = db_["n_solves_per_time_step"];

  int rs = db_["random_seed"];
  db_["random_seed"] = rs + 1; //count up the rng for use in the next step

  interface_->run_particle_phase(t_start, t_end, n_steps, rs);

}



void BrushWrapper::output(int &image, double t)
{
	// Updating stats and fetching moments is only relevant for
	// visualization and the output
	if(db_["output_write_vtk"].is(true))
	{
		interface_->update_stats();
		interface_->fetch_moment(0, &br_grid_psdmom_fcts_values_[0].at(0));
		interface_->fetch_moment(1, &br_grid_psdmom_fcts_values_[1].at(0));
		interface_->fetch_moment(2, &br_grid_psdmom_fcts_values_[2].at(0));
#ifdef __2D__
		output_writer_.write(t);
#elif defined(__3D__)
    std::string filename = db_["output_vtk_directory"];
    filename += "/" + db_["output_basename"].value_as_string();

    if(image<10) filename += ".0000";
    else if(image<100) filename += ".000";
    else if(image<1000) filename += ".00";
    else if(image<10000) filename += ".0";
    else filename += ".";
    filename += std::to_string(image) + ".vtk";

    //CB EXPERIMENTAL
    output_writer_.WriteVtk(filename.c_str());
    //output_writer_.WriteVtkDiscontinuous(filename.c_str(), disc_output_n_loc_verts ,disc_output_loc_verts);
    //END EXPERIMENTAL

    //CB DEBUG - Total mass control file, output every 1s
    //double del = fabs(t - std::round(t));
    //if(del < 1e-3)
    {
      std::ofstream myfile;
      std::string out_dir(db_["output_vtk_directory"].value_as_string());
      std::size_t pos = out_dir.find("VTK");// Lil hack - the directory above 'VTK',
      out_dir = out_dir.substr(0,pos);    // ...should be the general output directory.
      std::string mass_bal_file = out_dir + "mass_balance.csv";
      myfile.open (mass_bal_file.c_str(),std::ios_base::app);

      for(int c=0; c < brush_grid_->GetN_Cells(); ++c)
      {

        // Brush function representation
        auto bgn_ind= br_grid_space_.GetBeginIndex();
        auto dofs = br_grid_space_.GetGlobalNumbers();
        auto dof = dofs[bgn_ind[c]];
        if (br_grid_psdmom_fcts_values_[1].at(dof)!= 0.0) //get active only if there are particles in this cell
        {
          auto cntr = simplex_barycenter(*brush_grid_->GetCell(c));
          myfile << t << ",";
          myfile << c << ",";
          myfile << cntr[0] << ",";
          myfile << cntr[1] << ",";
          myfile << cntr[2] << ",";
          //Cell volume
          myfile << brush_grid_->GetCell(c)->GetMeasure() << ",";
          //Brush POTASHALUM mass (kg/m^3)
          myfile << br_grid_psdmom_fcts_values_[1].at(dof) << ",";
//          //ParMooN dissolved POTASHALUM values
//          bgn_ind= br_grid_param_fcts_[4]->GetFESpaceXD()->GetBeginIndex();
//          dofs = br_grid_param_fcts_[4]->GetFESpaceXD()->GetGlobalNumbers();
//          dof = dofs[bgn_ind[c]];
//          double val_in_kg = br_grid_param_fcts_[4]->GetValues()[dof] * density;
//          myfile << val_in_kg << ",";
          //velocity in the current ambient
          bgn_ind= br_grid_param_fcts_[0]->GetFESpaceXD()->GetBeginIndex();
          dofs = br_grid_param_fcts_[0]->GetFESpaceXD()->GetGlobalNumbers();
          dof = dofs[bgn_ind[c]];
          //x velo
          double velo = br_grid_param_fcts_[0]->GetValues()[dof];
          myfile << velo << ",";
          //y velo
          velo = br_grid_param_fcts_[1]->GetValues()[dof];
          myfile << velo << ",";
          //z velo
          velo = br_grid_param_fcts_[2]->GetValues()[dof];
          myfile << velo << "\n";
        }
      }
      myfile.close();
    }
    //END DEBUG

    image++;
#endif
	}

  interface_->write_particle_stats(t, moment_stats_file_);

  interface_->write_inlet_particle_list(inflow_particles_file_);

  double it = 0;
  if(std::abs(std::modf(t,&it)) < 1e-10 || std::abs(std::modf(t,&it) -1 ) < 1e-10) //this is a very inelegant way of checking if t is close to an integer
  {
    Output::info("OUTPUT","Taking particle snapshot at time ", it ," s");
    interface_->particle_population_snapshot(outflow_particles_file_,t);
  }

  // use this for the Eder flow crystallizer example
  //interface_->write_outlet_particle_list(outflow_particles_file_);

}
