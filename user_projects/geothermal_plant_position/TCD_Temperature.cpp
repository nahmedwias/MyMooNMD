#include "TCD_Temperature.hpp"
#include "TimeConvectionDiffusion.h"
#include "Multigrid.h"
#include "LocalAssembling.h"

#ifdef __2D__
 #include "Assemble2D.h"
 #include "SquareMatrix2D.h"
 #include "AuxParam2D.h"
#include "FEVectFunct2D.h"
#include "FEFunction2D.h"
#else
 #include "Assemble3D.h"
 #include "SquareMatrix3D.h"
 #include "AuxParam3D.h"
#include "FEVectFunct3D.h"
#include "FEFunction3D.h"
#endif

#include <TimeDiscRout.h>

#ifdef _MPI
 #include "ParFECommunicator3D.h"
#endif

// see https://stackoverflow.com/a/8016853

#ifdef __2D__

//constexpr char TCD_Temperature<2>::required_database_name[];

/** ************************************************************************ */
template <int d>
TCD_Temperature<d>::TCD_Temperature(const TDomain& domain,
        const ParameterDatabase& param_db,
        Example_TimeCD example)
: TimeConvectionDiffusion<d>(domain, param_db, example)
{
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("assembling with a finite element function as convection is not "
            "implemented for multigrid");
  }

  /*
  // velocity solver database
  ParameterDatabase tcd2d_db = param_db;
  //auto db_name = std::string(TCD_Temperature<2>::required_database_name);
  // use the given database or one of its nested databases, depending on which
  // one has the correct name. Otherwise the default solver database is used.
  cout << "HIER !!!!!"<< param_db.has_nested_database(db_name) << endl;
  cout << db_name << endl;
 if(param_db.has_nested_database(db_name))
  {
    tcd2d_db.merge(param_db.get_nested_database(db_name), false);
    cout << "HELLO" << endl;
    exit(0);
  }
  */
}

/** ************************************************************************ */
void mapping_local_parameters(const double *in, double *out)
{
  // coordinates:  x at in[0], y at in[1]
  out[0] = in[2];
  out[1] = in[3];
}

/** ************************************************************************ */
void temperature_coefficients(int n_points, double *x, double *y,
        double **parameters, double **coeffs,
        double distance, double nu)
{
  for(int i = 0; i < n_points; ++i)
  {
    //another approx. for domain [0, 10] x [0, 6]
    double r_well = 0.2; // 20cm
    double epsDelta = 50*r_well;
    double T_in = 303.15; //injection temperature  = 30 + 273.15; 

    
    //std::array<double, 2> center_source = {{5000.0 - distance/2., 3000.0}};
    std::array<double, 2> center_source = {{4500.0, 3000.0}}; 
    
    double x_distance_to_source = std::pow(std::abs(x[i] - center_source[0]), 2);
    double y_distance_to_source = std::pow(std::abs(y[i] - center_source[1]), 2);
    bool at_source =(x_distance_to_source < epsDelta*epsDelta) *
	(y_distance_to_source < epsDelta*epsDelta);


    coeffs[i][0] = nu; // diffusion
    coeffs[i][1] = parameters[i][0]; // convection, x-direction
    coeffs[i][2] = parameters[i][1]; // convection, y-direction
    coeffs[i][3] = 0.; // reaction
    coeffs[i][4] = 0.; // right-hand side

    if(at_source)
    {
      
      double magnitude = cos(Pi*(x[i] - center_source[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_source[1])/epsDelta) + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][3] = magnitude; // reaction
      coeffs[i][4] = magnitude * T_in; // right-hand side
      
      //double penalty_factor = 1000.;
      //coeffs[i][3] = 1 * penalty_factor; // reaction
      //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
    }
  }
}

/** ************************************************************************ */
template <int d>
void TCD_Temperature<d>::assemble(const FEVectFunct& convection,
        const double * x, double nu)
{
  if(this->db["space_discretization_type"].is("supg"))
  {
    ErrThrow("SUPG is not yet supported");
  }

  
  double distance = x[0];
  auto u1 = convection.GetComponent(0);
  auto u2 = convection.GetComponent(1);
;
  
  

  std::array<TFEFunction2D*, 2> fe_functions_pointers{{u1, u2}};

  
  auto& s = this->TCD_Temperature<d>::systems.front();

  LocalAssembling_type la_type = LocalAssembling_type::TCDStiffRhs;

// New LB 20.11.18 start
  if(this->db["space_discretization_type"].is("supg"))
 {
    // In the SUPG case:
    // M = (u,v) + \tau (u,b.grad v)
    la_type = LocalAssembling_type::TCDStiffMassRhs;
 }

    LocalAssembling<d> la(this->db, la_type,
            fe_functions_pointers.data(),
            this->example.get_coeffs(), 0);

    //call assembling, including mass matrix (SUPG!)
    // OLD: call_assembling_routine(s, la_a_rhs, la_m_supg , true);

    using namespace std::placeholders;
    la.ResetCoeffFct(std::bind(temperature_coefficients,_1, _2, _3, _4, _5,
            distance, nu));
    la.setBeginParameter({0});
    la.SetN_Parameters(2);
    la.setN_ParamFct(1);
    la.setParameterFct({mapping_local_parameters});
    la.setN_FeValues(2);
    la.setFeValueFctIndex({0, 1});
    la.setFeValueMultiIndex({D00, D00});

    auto fe_space = &this->get_space();
    auto * boundary_conditions = fe_space->get_boundary_condition();
    BoundaryValuesFunction * non_const_bound_value[1] {this->example.get_bd()[0]};

/*    auto mass_blocks = s.mass_matrix.get_blocks_uniquely();// s.stiffness_matrix.get_blocks_uniquely();
    auto matrix_mass = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());

    s.rhs.reset();
#ifdef __3D__
    Assemble3D(
#else
    Assemble2D(
#endif
            1, &fe_space, 1, &matrix_mass, 0, nullptr, 0, nullptr,
               nullptr, &boundary_conditions, non_const_bound_value, la);
*/
    auto blocks_stiff = s.stiffness_matrix.get_blocks_uniquely();
    auto matrix_stiff = reinterpret_cast<TSquareMatrix2D*>(blocks_stiff.at(0).get());
    double * rhs_entries = s.rhs.get_entries();


#ifdef __3D__
 Assemble3D(
#else
 Assemble2D(
#endif
         1, &fe_space, 1, &matrix_stiff, 0, nullptr, 1, &rhs_entries,
        &fe_space, &boundary_conditions, non_const_bound_value, la);




  delete u1;
  delete u2;


  // here the modifications due to time discretization begin
  if (!this->db["algebraic_flux_correction"].is("none") )
  {
    this->do_algebraic_flux_correction();
    this->rhs_from_time_disc = this->systems.front().rhs;
    return; // modifications due to time discretization are per-
            // formed inside the afc scheme, so step out here!
  }

  // preparing the right hand side discretized by the used time
  // stepping scheme
  //auto& s = this->TCD_Temperature<d>::systems.front(); //TODO: for some reason this does not work so i plugged in the definition everywhere or used auto
  this->rhs_from_time_disc.reset();
  this->rhs_from_time_disc = s.rhs;
  // all matrices are available
  unsigned int n_sols = this->time_stepping_scheme.n_old_solutions();
  std::vector<BlockVector> old_sols(n_sols);
  old_sols[0] = s.solution_m1;
  if(old_sols.size() == 2)
    old_sols[1] = s.solution_m2;
  std::vector<BlockVector> rhs(2);
  rhs[0] = this->rhs_from_time_disc;
  rhs[1] = this->old_rhs;
  // prepare the right hand side from the previous time step
  this->time_stepping_scheme.prepare_rhs_from_time_disc(s.stiffness_matrix,
                                                  s.mass_matrix, rhs, old_sols);
  this->rhs_from_time_disc = rhs[0];
  this->old_rhs = s.rhs;
  this->rhs_from_time_disc.copy_nonactive(s.rhs);

  for(auto &s : this->systems)
    this->time_stepping_scheme.prepare_system_matrix(s.stiffness_matrix,
                                               s.mass_matrix);
  s.solution.copy_nonactive(s.rhs);
}

/** ************************************************************************ */
template <int d>
void TCD_Temperature<d>::reset_for_output()
{
  this->outputWriter = DataWriter2D(this->db);
  ///OLD this->timeDependentOutput = DataWriter2D(this->db);
  FEFunction & fe_function = this->systems.front().fe_function;
  this->outputWriter.add_fe_function(&fe_function);
  auto& s = this->systems.front();
  s.stiffness_matrix.reset();
  s.mass_matrix.reset();
  s.rhs.reset();
}

/* ****************************** 3D ***************************************** */
#else


template <int d>
TCD_Temperature<d>::TCD_Temperature(TDomain& domain,
                                   const ParameterDatabase& param_db,
                                    Example_TimeCD example
//#ifdef _MPI
//                                    ,int maxSubDomainPerDof
//#endif
                                   )
 :TimeConvectionDiffusion<d>( domain, param_db, example //(domain.refine_and_get_hierarchy_of_collections(param_db
//#ifdef _MPI
//                                    ,int maxSubDomainPerDof
//#endif
//                                    )
//                                    , param_db, example
//#ifdef _MPI
//                                    ,int maxSubDomainPerDof
//#endif
)
{
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("assembling with a finite element function as convection is not "
             "implemented for multigrid");
  }
}

/*
#ifdef _MPI
CD3D_Temperature::CD3D_Temperature(TDomain& domain,
                                   const ParameterDatabase& param_db,
                                    Example_TimeCD3D example,
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
                                    Example_TimeCD3D example)
:Time_CD3D(domain.refine_and_get_hierarchy_of_collections(param_db), param_db, Example_TimeCD3D(param_db))
{
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("assembling with a finite element function as convection is not "
             "implemented for multigrid");
  }
}
#endif
*/

/*
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
*/


/*
CD3D_Temperature::CD3D_Temperature(std::list< TCollection* > collections, //const TDomain& domain,
                                   const ParameterDatabase& param_db,
                                   const Example_TimeCD3D& example)
 : Time_CD3D(collections, //domain,
param_db, example)
{
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("assembling with a finite element function as convection is not "
             "implemented for multigrid");
  }
}
*/

//*************************************************************************//
void mapping_local_parameters(const double *in, double *out)
{
  // coordinates:  x at in[0], y at in[1], z at in[3]
  out[0] = in[2];
  out[1] = in[3];
  out[2] = in[4];
}

//*************************************************************************//
void temperature_coefficients(int n_points, double *x, double *y, double *z,
                              double **parameters, double **coeffs,
                              double distance, double nu)
{
  for(int i = 0; i < n_points; ++i)
  {
    //another approx. for domain [0, 10] x [0, 6] x [,. 6]
    double a = 0.05;
    double T_in = 50;
    std::array<double, 3> center_source = {{5.0 - distance/2., 3. , 3.}};
    double x_distance_to_source = std::pow(std::abs(x[i]-center_source[0]), 2);
    double y_distance_to_source = std::pow(std::abs(y[i]-center_source[1]), 2);
    double z_distance_to_source = std::pow(std::abs(z[i]-center_source[2]), 2);
    bool at_source = x_distance_to_source + y_distance_to_source + z_distance_to_source < a*a;
    coeffs[i][0] = nu; // diffusion
    coeffs[i][1] = parameters[i][0]; // convection, x-direction
    coeffs[i][2] = parameters[i][1]; // convection, y-direction
    coeffs[i][3] = parameters[i][1]; // convection, y-direction
    coeffs[i][4] = 0.; // reaction
    coeffs[i][5] = 0.; // right-hand side
    if(at_source)
    {
      //double magnitude = cos(Pi*x_distance_to_source/a) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/a) + 1;
      //magnitude /= 4.*a*a;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      double penalty_factor = 1000.;
      coeffs[i][4] = 1 * penalty_factor; // reaction
      coeffs[i][5] = T_in * penalty_factor; // right-hand side
    }
  }
}

//*************************************************************************//
template <int d>
void TCD_Temperature<d>::assemble(FEVectFunct& convection, const double * x, double nu)
{
  if(this->db["space_discretization_type"].is("supg"))
  {
    ErrThrow("SUPG is not yet supported");
  }
  double distance = x[0];
  auto u1 = convection.GetComponent(0);
  auto u2 = convection.GetComponent(1);
  auto u3 = convection.GetComponent(2);

  std::array<FEFunction*, 3> fe_functions_pointers{{u1, u2, u3}};
  auto& s = this->TCD_Temperature<d>::systems.front();

  LocalAssembling<d> la(this->db, LocalAssembling_type::TCDStiffRhs,
                       fe_functions_pointers.data(),
                       this->example.get_coeffs(), 1); // db["discretization_type"]);

  using namespace std::placeholders;
  la.ResetCoeffFct(std::bind(temperature_coefficients, _1, _2, _3, _4, _5, _6, distance, nu));

  la.setBeginParameter({0});
  la.SetN_Parameters(3);
  la.setN_ParamFct(1);
  la.setParameterFct({mapping_local_parameters});
  la.setN_FeValues(3);
  la.setFeValueFctIndex({0, 1, 2});
  la.setFeValueMultiIndex({D000, D000, D000});

  auto fe_space = &this->get_space();
  auto blocks = s.stiffness_matrix.get_blocks_uniquely();
  auto matrix = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
  double * rhs_entries = s.rhs.get_entries();
  auto * boundary_conditions = fe_space->get_boundary_condition();
  BoundaryValuesFunction * non_const_bound_value[1] {this->example.get_bd()[0]};
  s.rhs.reset();
  s.stiffness_matrix.reset();

#ifdef __3D__
    Assemble3D(
#else
    Assemble2D(
#endif
            1, &fe_space, 1, &matrix, 0, nullptr, 1, &rhs_entries,
             &fe_space, &boundary_conditions, non_const_bound_value, la);

  delete u1;
  delete u2;
  delete u3;


  // here the modifications due to time discretization begin
  if (!this->db["algebraic_flux_correction"].is("none") )
  {
    this->do_algebraic_flux_correction();
    this->rhs_from_time_disc = this->systems.front().rhs;
    return; // modifications due to time discretization are per-
            // formed inside the afc scheme, so step out here!
  }

  // preparing the right hand side discretized by the used time
  // stepping scheme
  //SystemPerGrid& s = this->systems.front();
  this->rhs_from_time_disc.reset();
  this->rhs_from_time_disc = s.rhs;
  // all matrices are available
  unsigned int n_sols = this->time_stepping_scheme.n_old_solutions();
  std::vector<BlockVector> old_sols(n_sols);
  old_sols[0] = s.solution_m1;
  if(old_sols.size() == 2)
    old_sols[1] = s.solution_m2;
  std::vector<BlockVector> rhs(2);
  rhs[0] = this->rhs_from_time_disc;
  rhs[1] = this->old_rhs;
  // prepare the right hand side from the previous time step
  this->time_stepping_scheme.prepare_rhs_from_time_disc(s.stiffness_matrix,
                                                  s.mass_matrix, rhs, old_sols);
  this->rhs_from_time_disc = rhs[0];
  this->old_rhs = s.rhs;
  this->rhs_from_time_disc.copy_nonactive(s.rhs);

  for(auto &s : this->systems)
    this->time_stepping_scheme.prepare_system_matrix(s.stiffness_matrix,
                                               s.mass_matrix);
  s.solution.copy_nonactive(s.rhs);
}

//*************************************************************************//
template <int d>
void TCD_Temperature<d>::reset_for_output()
{
  this->outputWriter = DataWriter3D(this->db);
  TFEFunction3D & feFunction = this->systems.front().fe_function;
  this->outputWriter.add_fe_function(&feFunction);
  auto& s = this->systems.front();
  s.stiffness_matrix.reset();
  s.mass_matrix.reset();
  s.rhs.reset();
}


#endif




#ifdef __3D__
template class TCD_Temperature<3>;
#else
template class TCD_Temperature<2>;
#endif
