#include "CD3D_Temperature.hpp"
#include "Multigrid.h"
#include "LocalAssembling.h"
#include "Assemble3D.h"




/* **************************************************************************** */



#ifdef _MPI
CD3D_Temperature::CD3D_Temperature(TDomain& domain,
                                   const ParameterDatabase& param_db,
                                   const Example_TimeCD3D& example, 
                                   int maxSubDomainPerDof)
 : Time_CD3D(domain.refine_and_get_hierarchy_of_collections(param_db, maxSubDomainPerDof), param_db, example, maxSubDomainPerDof)
{
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("assembling with a finite element function as convection is not "
             "implemented for multigrid");
  }
}
#else
CD3D_Temperature::CD3D_Temperature( TDomain& domain,
                                   const ParameterDatabase& param_db,
                                   const Example_TimeCD3D& example)
:Time_CD3D(domain.refine_and_get_hierarchy_of_collections(param_db), param_db, Example_TimeCD3D(param_db))
{
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("assembling with a finite element function as convection is not "
             "implemented for multigrid");
  }
}
#endif


#ifdef _MPI
CD3D_Temperature::CD3D_Temperature( TDomain &domain, const ParameterDatabase& param_db, int maxSubDomainPerDof) 
:Time_CD3D(domain.refine_and_get_hierarchy_of_collections(param_db, maxSubDomainPerDof), param_db, Example_TimeCD3D(param_db), maxSubDomainPerDof)
}
#else
CD3D_Temperature::CD3D_Temperature(TDomain &domain, const ParameterDatabase& param_db) 
:Time_CD3D(domain.refine_and_get_hierarchy_of_collections(param_db), param_db, Example_TimeCD3D(param_db))
{
}
#endif



/*
CD3D_Temperature::CD3D_Temperature(std::list< TCollection* > collections, //const TDomain& domain,
                                   const ParameterDatabase& param_db,
                                   const Example_TimeCD3D& example_)
 : Time_CD3D(collections, //domain, 
param_db, example_)
{
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("assembling with a finite element function as convection is not "
             "implemented for multigrid");
  }
}
*/

void mapping_local_parameters(double *in, double *out)
{
  // coordinates:  x at in[0], y at in[1], z at in[3]
  out[0] = in[2];
  out[1] = in[3];
  out[2] = in[4];
}

void temperature_coefficients(int n_points, double *x, double *y, double *z,
                              double **parameters, double **coeffs,
                              double distance, double nu)
{
  for(int i = 0; i < n_points; ++i)
  {
    //another approx. for domain [0, 10] x [0, 6]
    double a = 0.05;
    double T_in = 50;
    std::array<double, 2> center_source = {{5.0 - distance/2., 3}};
    double x_distance_to_source = std::pow(std::abs(x[i]-center_source[0]), 2);
    double y_distance_to_source = std::pow(std::abs(y[i]-center_source[1]), 2);
    bool at_source = x_distance_to_source + y_distance_to_source < a*a;
    coeffs[i][0] = nu; // diffusion
    coeffs[i][1] = parameters[i][0]; // convection, x-direction
    coeffs[i][2] = parameters[i][1]; // convection, y-direction
    coeffs[i][3] = 0.; // reaction
    coeffs[i][4] = 0.; // right-hand side
    if(at_source)
    {
      //double magnitude = cos(Pi*x_distance_to_source/a) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/a) + 1;
      //magnitude /= 4.*a*a;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      double penalty_factor = 1000.;
      coeffs[i][3] = 1 * penalty_factor; // reaction
      coeffs[i][4] = T_in * penalty_factor; // right-hand side
    }
  }
}


void CD3D_Temperature::assemble(const TFEVectFunct3D& convection,
                                const double * x, double nu)
{
  if(db["space_discretization_type"].is("supg"))
  {
    ErrThrow("SUPG is not yet supported");
  }
  double distance = x[0];
  auto u1 = convection.GetComponent(0);
  auto u2 = convection.GetComponent(1);
  auto u3 = convection.GetComponent(2);  
  std::array<TFEFunction3D*, 3> fe_functions_pointers{{u1, u2, u3}};
  SystemPerGrid& s = this->systems_.front();
  
  LocalAssembling3D la(this->db, LocalAssembling_type::TCD3D,
                       fe_functions_pointers.data(),
                       this->example_.get_coeffs(), 0);
  
  using namespace std::placeholders;
  la.ResetCoeffFct(std::bind(temperature_coefficients,_1, _2, _3, _4, _5, _6,
                             distance, nu));
  
  la.setBeginParameter({0});
  la.SetN_Parameters(3);
  la.setN_ParamFct(1);
  la.setParameterFct({mapping_local_parameters});
  la.setN_FeValues(3);
  la.setFeValueFctIndex({0, 1, 2});
  la.setFeValueMultiIndex({D000, D000, D000});
  
  auto fe_space = &this->get_space();
  auto blocks = s.stiffMatrix_.get_blocks_uniquely();
  auto matrix = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
  double * rhs_entries = s.rhs_.get_entries();
  auto * boundary_conditions = fe_space->get_boundary_condition();
  BoundValueFunct3D * non_const_bound_value[1] {example_.get_bd()[0]};
  
  s.rhs_.reset();
  s.stiffMatrix_.reset();
  Assemble3D(1, &fe_space, 1, &matrix, 0, nullptr, 1, &rhs_entries,
             &fe_space, &boundary_conditions, non_const_bound_value, la);
  delete u1;
  delete u2;
  delete u3;

  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  if(TDatabase::TimeDB->THETA4)
  {
    // scale by time step length and theta4 (only active dofs)
    s.rhs_.scaleActive(tau*TDatabase::TimeDB->THETA4);    
    // add old right hand side scaled by time step length and theta3 (only 
    // active dofs)
    if (TDatabase::TimeDB->THETA3 != 0.)
      s.rhs_.addScaledActive((this->old_rhs), tau * TDatabase::TimeDB->THETA3);	

    // save old right hand side (only if THETA3 != 0)
    if(TDatabase::TimeDB->THETA3)
    {
      this->old_rhs.addScaledActive(s.rhs_, -1./(tau*TDatabase::TimeDB->THETA3));
      this->old_rhs.scaleActive(-TDatabase::TimeDB->THETA3/TDatabase::TimeDB->THETA4);
    }
  }
  else
  {
    if(TDatabase::TimeDB->TIME_DISC == 0)
    {
      ErrThrow("Forward Euler method is not supported. "
          "Choose TDatabase::TimeDB->TIME_DISC as 1 (bw Euler)"
          " or 2 (Crank-Nicoloson)");
    }
  }
  
  // rhs += M*uold
  s.massMatrix_.apply_scaled_add_actives(s.solution_, s.rhs_, 1.0);
  // rhs -= tau*theta2*A_old*uold
  s.rhs_.addScaledActive(s.old_Au, -tau*TDatabase::TimeDB->THETA2);
  
  // preparing the left hand side, i.e., the system matrix
  // stiffness matrix is scaled by tau*THETA1, after solving 
  // the matrix needs to be descaled if the coeffs does not depends
  // on time

  //scale  stiffness matrix...
  const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
  s.stiffMatrix_.scale_blocks_actives(tau*TDatabase::TimeDB->THETA1, cell_positions);
  // ...and add the mass matrix
  const FEMatrix& mass_block = *s.massMatrix_.get_blocks().at(0).get();
  s.stiffMatrix_.add_matrix_actives(mass_block, 1.0, {{0,0}}, {false});
  
  s.solution_.copy_nonactive(systems_[0].rhs_);
}

void CD3D_Temperature::reset_for_output()
{
  this->outputWriter = DataWriter3D(this->db);
  TFEFunction3D & feFunction = this->systems_.front().feFunction_;
  this->outputWriter.add_fe_function(&feFunction);
  auto& s = this->systems_.front();
  s.stiffMatrix_.reset();
  s.massMatrix_.reset();
  s.rhs_.reset();
}

