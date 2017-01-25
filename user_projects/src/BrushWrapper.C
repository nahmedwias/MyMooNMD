/*
 * BrushWrapper.C
 *
 *  Created on: Dec 14, 2016
 *      Author: bartsch
 */

#include <BrushWrapper.h>
//Brush code for the particles
#include <parmoon_interface.h>


// TODO Before we think about proper boundary conditions we go for Dirichlet zero.
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
      db_["max_sp_per_cell"], db_["max_m0_per_cell"],
      1 // number of species in the continuous phase TODO hard coded so far
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
  // CALLING THE FOLLOWING DESTRUCTOR EXPLICITELY!!!
  //delete interface_;
  for (auto f : pd_moments_)
  {
    delete f;
  }
}

void BrushWrapper::set_velocity(const TFEVectFunct2D& u)
{ // NOTE: This will only set the velocity at the center points
  // of Brush's cells, as Brush expects. I am aware that this point
  // evaluation means loss of information - this whole process of exchanging
  // information might be reworked later.
  size_t n_points = output_sample_points_.size();
  std::vector<double> u0(n_points,0);
  std::vector<double> u1(n_points,0);

  TFEFunction2D* u0_fe = u.GetComponent(0);
  TFEFunction2D* u1_fe = u.GetComponent(1);

  for (size_t p = 0 ; p < n_points ; ++p)
  {
    double vals0[5];
    double vals1[5];
    double x = output_sample_points_[p][0];
    double y = output_sample_points_[p][1];
    u0_fe->FindGradient(x,y,vals0);
    u1_fe->FindGradient(x,y,vals1);
    u0.at(p)=vals0[0];
    u1.at(p)=vals1[0];
  }

  interface_->fill_velocity_buffer(u0,0);
  interface_->fill_velocity_buffer(u1,1);

  delete u0_fe;
  delete u1_fe;
}

void BrushWrapper::set_temperature(const TFEFunction2D& T)
{
  size_t n_points = output_sample_points_.size();
  std::vector<double> T_values(n_points,0);

  for (size_t p = 0 ; p < n_points ; ++p)
   {
     double vals[5];
     double x = output_sample_points_[p][0];
     double y = output_sample_points_[p][1];
     T.FindGradient(x,y,vals);
     T_values.at(p)=vals[0];
   }

   interface_->fill_temperature_buffer(T_values);
}

void BrushWrapper::set_concentrations(std::vector<const TFEFunction2D*> c)
{
  // FIXME Note that the first element is supposed to be temperature and is disregarded.
  for( size_t i=1 ; i< c.size() ; ++i )
  {
    size_t n_points = output_sample_points_.size();
    std::vector<double> ci_values(n_points,0);

    for (size_t p = 0 ; p < n_points ; ++p)
    {
      double vals[5];
      double x = output_sample_points_[p][0];
      double y = output_sample_points_[p][1];
      c[i]->FindGradient(x,y,vals);
      ci_values.at(p)=vals[0];
    }

    interface_->fill_concentration_buffer(ci_values,i-1);
  }
}

void BrushWrapper::solve(double t_start, double t_end)
{
  //now call the solver
  interface_->run_particle_phase(t_start, t_end);

  //and get the ensemble statistics updated
  interface_->update_stats();

  // store the updated moments
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




