/*
 * BrushWrapper.C
 *
 *  Created on: Dec 14, 2016
 *      Author: bartsch
 */

#include <BrushWrapper.h>


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
  //TODO 3D?!
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

  //std::cout << p.size() << std::endl;
  return p;
}

int fe_order = 0;

BrushWrapper::BrushWrapper(TCollection* coll, const ParameterDatabase& db)
: db_(db), output_writer_(db_),
  pd_moments_grid_(coll),
  pd_moments_space_(pd_moments_grid_, (char*)"psd-moms", (char*)"psd-moms",
                    DirichletBoundaryConditions, fe_order, nullptr)
{

  // set up Brushs ParMooN interface
  interface_ = new Brush::InterfacePM(
      db_["geo_file"], db_["sweep_file"], " ",
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

  // fill them and add them to the output and print them (which is a temporary test)
  // tell the interface which sample points parMooN is interested in
  // (for a start: the cell centers of all cells)
  int n_cells = coll->GetN_Cells();
  std::vector<std::valarray<double>> sample_points(n_cells,{0,0,0}); //z value will stay 0
  for(int c = 0; c < n_cells; ++c)
  {
    sample_points[c] = center_point(*coll->GetCell(c));
  }

  interface_->set_output_sample_points(sample_points);
  input_sample_points_ = sample_points;

  output_sample_points_ = interface_->get_input_sample_points();

  //force an update of all ensemble stats
  interface_->update_stats();

  interface_->write_moment(0, &pd_moments_values_[0].at(0));
  interface_->write_moment(1, &pd_moments_values_[1].at(0));
  interface_->write_moment(2, &pd_moments_values_[2].at(0));

  output_writer_.add_fe_function(pd_moments_[0]);
  output_writer_.add_fe_function(pd_moments_[1]);
  output_writer_.add_fe_function(pd_moments_[2]);

  output_writer_.write();

}

BrushWrapper::~BrushWrapper()
{
  //delete interface_; //FIXME THERE IS A BUG (CONNECTED TO THE MOON-GEOMETRY) WHEN CALLING THIS DESTRUCTOR EXPLICITELY!!!
  for (auto f : pd_moments_)
  {
    delete f;
  }
}


void BrushWrapper::write_vtk()
{
  ;
}
