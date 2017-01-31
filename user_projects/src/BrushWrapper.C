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

#include <ASA_crystallizer.h> //TODO This should be handled at runtime instead!

#include <algorithm>

// TODO Before we think about proper boundary conditions we go for Dirichlet zero.
void DirichletBoundaryConditions(int BdComp, double t, BoundCond &cond)
{
      cond = DIRICHLET;
}
void ZeroBoundaryValues(int BdComp, double Param, double &value)
{
      value = 0;
}

//TODO Use an std::function instead
typedef double my_funct_type (const std::vector<double>&);

// With given values ux, uy, p, T and ASA conc in a certain point,
// will evaluate the molar concentration of EtOH, which is a derived quantity.
//TODO Move to an example!
double derived_concentration_EtOH(const std::vector<double>& data)
{
  if (data.size() != 7)
    throw std::runtime_error("derived_concentration_EtOH: expected 7 data points."
        " ux, uy, p, T, ASA, 0(CH3CH2OH) , 0(ASASUP)");

  return 1; //TODO CODE!
}
//TODO Move to an example!
// With given values ux, uy, p, T and ASA conc in a certain point,
// will evaluate the supersaturation concentration of ASA, which is a
// derived quantity.
// Temperature must be in K and ASA concentration in mol/m^3. Output
// will also be in mol/m^3.
double derived_concentration_ASASUP(const std::vector<double>& data)
{
  if (data.size() != 7)
    throw std::runtime_error("derived_concentration_EtOH: expected 7 data points."
        " ux, uy, p, T, ASA, CH3CH2OH , 0(ASASUP)");

    double T = data[3];     // grab temperature
    double c_asa = data[4]; // grab ASA concentration

    //theoretical supersaturation in terms of ASA mole fraction, as given by Eder et al.
    double chi_sat = pow(10, 27.769 - (2500.906/T) - 8.323 * log10(T));

    //now this has to be transformed to supersaturation in
    // terms of molar concentration, which is a bit cumbersome
    double M_S = chi_sat*Physics::M_ASA + (1-chi_sat)*Physics::M_Ethanol; //molar mass of solution
    double w_Ethanol = 1 - chi_sat*(Physics::M_ASA/M_S); //Ethanol mass fraction
    //TODO ideal solution assumption is questionable here!
    double rho_S = Physics::rho_E / w_Ethanol; //density of solution calculated under "ideal solution" asumption

    //now here comes supersaturation in terms of molar concentration
    double c_sat = chi_sat * rho_S / M_S;
    double c_supsat = std::max(0.0, c_asa - c_sat);

//    if (c_supsat > 0)
//      Output::print("c_asa: ", c_asa, ", c_sat: ", c_sat , ", supsat: ", c_supsat, " mol/m^3");

    return c_supsat;
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

int fe_order = 0;

BrushWrapper::BrushWrapper(TCollection* coll, const ParameterDatabase& db)
: db_(db), output_writer_(db_), moment_stats_file_(db_["out_part_moments_file"].get<std::string>()),
  pd_moments_grid_(coll),
  pd_moments_space_(pd_moments_grid_, (char*)"psd-moms", (char*)"psd-moms",
                    DirichletBoundaryConditions, fe_order, nullptr)
{

  // set up Brushs ParMooN interface
  interface_ = new Brush::InterfacePM(
      db_["geo_file"], db_["third_dim_stretch"],
      db_["sweep_file"], " ",
      db_["therm_file"], db_["chem_file"],
      db_["max_sp_per_cell"], db_["max_m0_per_cell"]
  );

  // load the initial particle solution
  interface_->set_initial_particles(db_["init_partsol_file"]);

  // set up FE functions and their values

  int fe_length = pd_moments_space_.GetN_DegreesOfFreedom();

  std::vector<double> dummy_fe_values(fe_length, 0.0);
  pd_moments_values_ = std::vector<std::vector<double>>(3, dummy_fe_values);

  pd_moments_.resize(3);
  pd_moments_.at(0) = new TFEFunction2D(&pd_moments_space_,
        (char*)"pd-m0", (char*)"", &pd_moments_values_[0].at(0), fe_length);
  pd_moments_.at(1) = new TFEFunction2D(&pd_moments_space_,
        (char*)"pd-m1", (char*)"", &pd_moments_values_[1].at(0), fe_length);
  pd_moments_.at(2) = new TFEFunction2D(&pd_moments_space_,
        (char*)"pd-m2", (char*)"", &pd_moments_values_[2].at(0), fe_length);

  //add the fe function to the output writer for later use
  output_writer_.add_fe_function(pd_moments_[0]);
  output_writer_.add_fe_function(pd_moments_[1]);
  output_writer_.add_fe_function(pd_moments_[2]);

  //write a nice header into the particle stats file
  interface_->write_header(moment_stats_file_);

  // fill them and add them to the output and print them (which is a temporary test)
  // tell the interface which sample points parMooN is interested in
  // (for a start: the cell centers of all cells)
  int n_cells = coll->GetN_Cells();
  std::vector<std::valarray<double>> parmoon_cell_centers(n_cells,{0,0,0}); //z value will stay 0
  for(int c = 0; c < n_cells; ++c)
  {
    parmoon_cell_centers[c] = center_point(*coll->GetCell(c));
  }

  interface_->set_output_sample_points(parmoon_cell_centers);

  output_sample_points_ = interface_->get_cell_centers();

  //force an update of all ensemble stats
  interface_->update_stats();

  // store moments m0, m1 and m2.
  interface_->fetch_moment(0, &pd_moments_values_[0].at(0));
  interface_->fetch_moment(1, &pd_moments_values_[1].at(0));
  interface_->fetch_moment(2, &pd_moments_values_[2].at(0));

}

BrushWrapper::~BrushWrapper()
{
  //FIXME THERE IS A BUG (CONNECTED TO THE MOON-GEOMETRY) WHEN
  // CALLING THE FOLLOWING DESTRUCTOR!!!
  //delete interface_;
  for (auto f : pd_moments_)
  {
    delete f;
  }
}

void BrushWrapper::reset_fluid_phase(
    const TFEVectFunct2D& u,
    const TFEFunction2D& p,
    std::vector<const TFEFunction2D*> species
    )
{
  size_t n_points = output_sample_points_.size();
  //TODO "decoder-vektor" koennte im Beispiel gespeichert sein
  Brush::DataPM pm_data({"ux","uy","p","T","ASA","CH3CH2OH","ASASUP"}, n_points);

  //For the velocity we need two xtra fe functions
  TFEFunction2D* u0_fe = u.GetComponent(0);
  TFEFunction2D* u1_fe = u.GetComponent(1);

  //prepare a vector of pointers to FEFunctions
  size_t dim = 2;
  size_t n_specs = species.size();
  size_t n_specs_derived = 2; //TODO EtOH and ASASUP - but do not hard code!

  std::vector<const TFEFunction2D*> fe_functs(dim + 1 + n_specs);
  fe_functs[0] = u0_fe;
  fe_functs[1] = u1_fe;
  fe_functs[2] = &p;
  std::copy( species.begin() , species.end(), &fe_functs[3]);

  std::vector<my_funct_type*> derived_concentrations(n_specs_derived);
  derived_concentrations[0] = derived_concentration_EtOH;   //TODO do not hard code
  derived_concentrations[1] = derived_concentration_ASASUP; //TODO do not hard code

  //loop over all evaluation points (i.e., Brushs cell midpoints)
  for (size_t point = 0 ; point < n_points ; ++point)
  {
    //x and y value of the current output sample point
    double x = output_sample_points_[point][0];
    double y = output_sample_points_[point][1];

    size_t data_set_size = dim + 1 + n_specs + n_specs_derived;
    std::vector<double> data_set(data_set_size , 0.0);

    // point evaluations of all fe functions
    for(size_t i = 0 ; i < dim + 1 + n_specs ; ++i)
    {
      double eval[3]; // will include differentials, thus length is '3')
      fe_functs[i]->FindGradient(x,y,eval);
      data_set[i] = eval[0];
    }

    for(size_t i = dim + 1 + n_specs; i < data_set_size; ++i)
    {
      size_t index = i - (dim + 1 + n_specs);
      double f = derived_concentrations[index](data_set); //TODO evaluate the 'index'th derived concentration function
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
  interface_->run_particle_phase(t_start, t_end);

  // Updating stats and fetching moments is only relevant for
  // visualization and the output!
  // TODO Get this (Brush to ParMooN) right after ParMooN to Brush
  interface_->update_stats();
  interface_->fetch_moment(0, &pd_moments_values_[0].at(0));
  interface_->fetch_moment(1, &pd_moments_values_[1].at(0));
  interface_->fetch_moment(2, &pd_moments_values_[2].at(0));

}


void BrushWrapper::output(double t)
{
  output_writer_.write(t);

  //TODO enable this for only certain steps
  interface_->write_particle_stats(t, moment_stats_file_);

}



