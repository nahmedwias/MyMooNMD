#include "TCD2D_Temperature.hpp"
#include "Multigrid.h"
#include "LocalAssembling.h"
#include "Assemble2D.h"


// see https://stackoverflow.com/a/8016853
//constexpr char TCD2D_Temperature::required_database_name[];

/** ************************************************************************ */
TCD2D_Temperature::TCD2D_Temperature(const TDomain& domain,
        const ParameterDatabase& param_db,
        Example_TimeCD2D example)
: Time_CD2D(domain, param_db, example)
{
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("assembling with a finite element function as convection is not "
            "implemented for multigrid");
  }

  /*// velocity solver database
  ParameterDatabase tcd2d_db = param_db;
  auto db_name = std::string(TCD2D_Temperature::required_database_name);
  // use the given database or one of its nested databases, depending on which
  // one has the correct name. Otherwise the default solver database is used.
  cout << "HIER !!!!!"<< param_db.has_nested_database(db_name) << endl;
  cout << db_name << endl;
 if(param_db.has_nested_database(db_name))
  {
    tcd2d_db.merge(param_db.get_nested_database(db_name), false);
    cout << "HELLO" << endl;
    exit(0);
  }*/
}

/** ************************************************************************ */
void mapping_local_parameters(double *in, double *out)
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
    double a = 5.; //0.05;
    double T_in = 303.15; //= 30 + 273.15; //50;
    //cout <<"distance  in temperature!!!!"<< distance << endl;
    std::array<double, 2> center_source = {{5000.0 - distance/2., 3000.0}}; //{{5.0 - distance/2., 3}};

    //cout <<" center_source[0] in temperature!!!!"<<  center_source[0] << endl;
    //cout <<" center_source[1] in temperature!!!!"<<  center_source[1] << endl;

    double x_distance_to_source = std::pow(std::abs(x[i] - center_source[0]), 2);
    double y_distance_to_source = std::pow(std::abs(y[i] - center_source[1]), 2);
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

/** ************************************************************************ */
void TCD2D_Temperature::assemble(const TFEVectFunct2D& convection,
        const double * x, double nu)
{
  if(db["space_discretization_type"].is("supg"))
  {
    ErrThrow("SUPG is not yet supported");
  }

  double distance = x[0];
  auto u1 = convection.GetComponent(0);
  auto u2 = convection.GetComponent(1);
  std::array<TFEFunction2D*, 2> fe_functions_pointers{{u1, u2}};

  System_per_grid& s = this->systems.front();

  LocalAssembling2D la(this->db, LocalAssembling_type::TCD2D,
             fe_functions_pointers.data(),
             this->example.get_coeffs(), 0);

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
   auto blocks_stiff = s.stiff_matrix.get_blocks_uniquely();
   auto matrix_stiff = reinterpret_cast<TSquareMatrix2D*>(blocks_stiff.at(0).get());
   double * rhs_entries = s.rhs.get_entries();
   auto * boundary_conditions = fe_space->get_boundary_condition();
   BoundValueFunct2D * non_const_bound_value[1] {example.get_bd()[0]};

   s.rhs.reset();
   s.stiff_matrix.reset();
   Assemble2D(1, &fe_space, 1, &matrix_stiff, 0, nullptr, 1, &rhs_entries,
           &fe_space, &boundary_conditions, non_const_bound_value, la);


// New LB 20.11.18 start
  if(db["space_discretization_type"].is("supg"))
 {
    // In the SUPG case:
    // M = (u,v) + \tau (u,b.grad v)

    LocalAssembling2D la_m_supg(this->db, LocalAssembling_type::TCD2D_Mass,
            fe_functions_pointers.data(),
            this->example.get_coeffs(), 0);

    //call assembling, including mass matrix (SUPG!)
    // OLD: call_assembling_routine(s, la_a_rhs, la_m_supg , true);

    using namespace std::placeholders;
    la_m_supg.ResetCoeffFct(std::bind(temperature_coefficients,_1, _2, _3, _4, _5,
            distance, nu));
    la_m_supg.setBeginParameter({0});
    la_m_supg.SetN_Parameters(2);
    la_m_supg.setN_ParamFct(1);
    la_m_supg.setParameterFct({mapping_local_parameters});
    la_m_supg.setN_FeValues(2);
    la_m_supg.setFeValueFctIndex({0, 1});
    la_m_supg.setFeValueMultiIndex({D00, D00});

    auto mass_blocks = s.mass_matrix.get_blocks_uniquely();// s.stiff_matrix.get_blocks_uniquely();
    auto matrix_mass = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());

    s.rhs.reset();
    s.mass_matrix.reset();
    //matrix_mass[0]->reset();
    Assemble2D(1, &fe_space, 1, &matrix_mass, 0, nullptr, 0, nullptr,
               nullptr, &boundary_conditions, non_const_bound_value, la_m_supg);

 }

  // New LB 20.11.18 end



  delete u1;
  delete u2;

  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  if(TDatabase::TimeDB->THETA4)
  {
    // scale by time step length and theta4 (only active dofs)
    s.rhs.scaleActive(tau*TDatabase::TimeDB->THETA4);    
    // add old right hand side scaled by time step length and theta3 (only 
    // active dofs)
    if(TDatabase::TimeDB->THETA3 != 0.)
      s.rhs.addScaledActive((this->old_rhs), tau*TDatabase::TimeDB->THETA3);	

    // save old right hand side (only if THETA3 != 0)
    if(TDatabase::TimeDB->THETA3)
    {
      this->old_rhs.addScaledActive(s.rhs, -1./(tau*TDatabase::TimeDB->THETA3));
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
  s.mass_matrix.apply_scaled_add_actives(s.solution, s.rhs, 1.0);
  // rhs -= tau*theta2*A_old*uold
  s.rhs.addScaledActive(s.old_Au, -tau*TDatabase::TimeDB->THETA2);

  // preparing the left hand side, i.e., the system matrix
  // stiffness matrix is scaled by tau*THETA1, after solving 
  // the matrix needs to be descaled if the coeffs does not depends
  // on time

  //scale  stiffness matrix...
  const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
  s.stiff_matrix.scale_blocks_actives(tau*TDatabase::TimeDB->THETA1, cell_positions);
  // ...and add the mass matrix
  const FEMatrix& mass_block = *s.mass_matrix.get_blocks().at(0).get();
  s.stiff_matrix.add_matrix_actives(mass_block, 1.0, {{0,0}}, {false});

  s.solution.copy_nonactive(systems[0].rhs);
}

/** ************************************************************************ */
void TCD2D_Temperature::reset_for_output()
{
  this->timeDependentOutput = DataWriter2D(this->db);
  TFEFunction2D & fe_function = this->systems.front().fe_function;
  this->timeDependentOutput.add_fe_function(&fe_function);
  auto& s = this->systems.front();
  s.stiff_matrix.reset();
  s.mass_matrix.reset();
  s.rhs.reset();
}

