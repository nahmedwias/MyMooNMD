#include "GeothermalPlantsPositionOptimization.hpp"
#include "gppo_example.hpp"
#include "LocalAssembling.h"
#include "NSE_local_assembling_routines.h"
#include "Assemble2D.h"
#include "LoopInfo.h"
#include "Database.h"
#include "TimeDiscRout.h"
#include <algorithm>


ParameterDatabase GeothermalPlantsPositionOptimization::default_GPPO_database()
{
  // is only needed for the adjoint problem, but there is only one global 
  // parameter for this.
  TDatabase::ParamDB->NSTYPE = 4;
  auto db = ParameterDatabase::parmoon_default_database();
  db.merge(TDomain::default_domain_parameters());
  db.merge(Example2D::default_example_database());
  db.merge(LocalAssembling<2>::default_local_assembling_database());
  
  db.add("alpha_cost", 1., 
         "The scalar alpha in the functional J_hat, which is multiplied with "
         "the control term.", 0., 1.0e10);
  db.add("cost_functional", {1., 0.},
         "Switch between a few possible cost functionals via weights. The "
         "order of the functionals is: L2_norm_of_curl, backward_facing_step, ",
         0., 1.);        
  return db;
}

ParameterDatabase get_primal_flow_database(ParameterDatabase param_db)
{
  std::string basename = param_db["output_basename"];
  basename += std::string("flow");
  param_db["output_basename"].set(basename, false);
  return param_db;
}
ParameterDatabase get_primal_temperature_database(ParameterDatabase param_db)
{
  std::string basename = param_db["output_basename"];
  basename += std::string("temperature");
  param_db["output_basename"].set(basename, false);
  return param_db;
}

ParameterDatabase get_adjoint_database(const ParameterDatabase& param_db)
{
  ParameterDatabase adjoint_db(ParameterDatabase::parmoon_default_database());
  adjoint_db.merge(ParameterDatabase::default_output_database(), true);
  adjoint_db.merge(param_db, true, true);
  std::string adjoint_basename = param_db["output_basename"];
  adjoint_basename += std::string("gppo_adjoint");
  adjoint_db["output_basename"].set(adjoint_basename, false);
  adjoint_db["nonlinloop_damping_factor"] = 1.;
  return adjoint_db;
}

GeothermalPlantsPositionOptimization::GeothermalPlantsPositionOptimization(
  const TDomain& domain, const ParameterDatabase& param_db)
 : db(default_GPPO_database()), n_control(0), 
   brinkman2d_primal(domain, get_primal_flow_database(param_db),
                     get_gppo_flow_example(param_db)),
   cd2d_primal(domain, get_primal_temperature_database(param_db),
               get_gppo_temperature_example(param_db)),
   optimization_info("optimization", true, true, 1), // 1 -> full verbosity
   current_J_hat(std::numeric_limits<double>::infinity()), 
   control_old(), n_calls(0),
   n_computation_derivative(0)
{
  Output::print<5>("Creating the GeothermalPlantsPositionOptimization object");
  db.merge(param_db, false);
  
  // assemble all terms in matrix and right hand side. This matrix is not
  // changed afterwards
  brinkman2d_primal.assemble();
  
  n_control = 1;
  control_old = std::vector<double>(n_control, 0.0);
  Output::print<3>("Created the GeothermalPlantsPositionOptimization object, ",
                   "n_control = ", n_control);
}


double GeothermalPlantsPositionOptimization::compute_functional_and_derivative(
  unsigned int n, const double* x, double* grad)
{
  if(n != n_control)
  {
    ErrThrow("the given dimension ", n, " must match the dimension ", n_control,
             " of the control space");
  }
  ++n_calls;
  Output::print<2>("Call to 'compute_functional_and_derivative' ", n_calls);
  if(n_calls > 1)
  {
    double norm_diff_control = std::sqrt(
      std::inner_product(x, x+n, control_old.begin(), 0.0, std::plus<double>(),
                         [](double a, double b){ return (a-b) * (a-b); }));
    double norm_mean_control = std::sqrt(
      std::inner_product(x, x+n, control_old.begin(), 0.0, std::plus<double>(),
                         [](double a, double b){ return (a+b) * (a+b) / 4.; }));
    Output::print("diff of control (in l2): ", norm_diff_control,
                  "  mean: ", norm_mean_control);
  }
  std::copy(x, x+n, control_old.begin());
  
  apply_control_and_solve(x);
  current_J_hat = compute_functional(x);
  
  if(grad != nullptr)
  {
    ++n_computation_derivative;
    solve_adjoint_equation();
    compute_derivative(x, grad);
    double norm_of_grad = std::sqrt(std::accumulate(
      grad, grad+n, 0., [](double a, double b){ return a + b*b; }));
    Output::print("norm of gradient: ", norm_of_grad);
//     for(auto i = 0u; i < n_control; ++i)
//     {
//       Output::print("computed_gradient[", i, "] = ", grad[i]);
//     }
  }
  
  if(n_calls == 1)
  {
    optimization_info.restart("optimization", current_J_hat);
  }
  optimization_info.print(n_calls, current_J_hat);
  return current_J_hat;
}

void approximate_delta_functions(int n_points, double *x, double *y,
                                 double **parameters, double **coeffs,
                                 double distance)
{
  for(int i = 0; i < n_points; ++i)
  {
    //another approx. for domain [0, 10] x [0, 6]
    double a = 0.05;
    double r_Wells = 5*0.001;
    double Q_in = 0.03*Pi*r_Wells*r_Wells * 1000;
    std::array<double, 2> center_source = {{5.0 - distance/2., 3}};
    std::array<double, 2> center_sink = {{5.0 + distance/2., 3}};
    double x_distance_to_source = std::pow(std::abs(x[i]-center_source[0]), 2);
    double y_distance_to_source = std::pow(std::abs(y[i]-center_source[1]), 2);
    double x_distance_to_sink = std::pow(std::abs(x[i]-center_sink[0]), 2);
    double y_distance_to_sink = std::pow(std::abs(y[i]-center_sink[1]), 2);
    bool at_source = x_distance_to_source + y_distance_to_source < a*a;
    bool at_sink = x_distance_to_sink + y_distance_to_sink < a*a;
    coeffs[i][0] = 0.;
    coeffs[i][1] = 0.;
    coeffs[i][2] = 0.;
    coeffs[i][3] = 0.;
    if(at_source)
    {
      double magnitude = cos(Pi*x_distance_to_source/a) + 1;
      magnitude *= cos(Pi*y_distance_to_source/a) + 1;
      magnitude /= 4.*a*a;
      coeffs[i][3] += magnitude * Q_in; // source
    }
    if(at_sink)
    {
      double magnitude = cos(Pi*x_distance_to_sink/a) + 1;
      magnitude *= cos(Pi*y_distance_to_sink/a) + 1;
      magnitude /= 4.*a*a;
      coeffs[i][3] -=  magnitude * Q_in; // sink
    }
  }
}

void GeothermalPlantsPositionOptimization::apply_control_and_solve(const double* x)
{
  Chrono time;
  // apply control x
  double distance = x[0];
  Output::print("current control: ", std::setprecision(14), distance);
  using namespace std::placeholders;
  CoeffFct2D coeff = std::bind(approximate_delta_functions, _1, _2, _3, _4, _5,
                               distance);
  std::string disc_type = this->db["space_discretization_type"];
  bool nonsymm_gls = (disc_type == std::string("nonsymm_gls"));
  int rhs_div_sign = 1;
  if (nonsymm_gls)
  {
    rhs_div_sign = -1;
  }
  LocalAssembling<2> la(2, {D00, D00}, {0, 1}, {}, {}, {0, 0, 1}, coeff, 
			std::bind(NSRightHandSide<2>, _1, _2, _3, _4, _5, _6,
				  _7, _8, rhs_div_sign),
			nullptr, 0, 3, 0, {}, {}, 0, nullptr,
			0, {}, {});

  auto& v_space = brinkman2d_primal.get_velocity_space();
  auto& p_space = brinkman2d_primal.get_pressure_space();
  std::array<const TFESpace2D*, 3> fespaces = {{&v_space, &v_space, &p_space}};
  auto& rhs = brinkman2d_primal.get_rhs();
  rhs.reset();
  double *rhs_pointers[3] = {rhs.block(0), rhs.block(1), rhs.block(2)};
  BoundCondFunct2D * boundary_conditions[3] = {
      v_space.get_boundary_condition(), v_space.get_boundary_condition(),
      p_space.get_boundary_condition() };
  auto& example = brinkman2d_primal.get_example();
  std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];
  Assemble2D(2, fespaces.data()+1, 0, nullptr, 0, nullptr, 3, rhs_pointers,
             fespaces.data(), boundary_conditions,
             non_const_bound_values.data(), la);
  
  double mat_norm 
    = brinkman2d_primal.get_matrix().get_combined_matrix()->GetNorm();
  double rhs_norm = rhs.norm();
  Output::print("---- norms: ", mat_norm, "  ", rhs_norm);
  
  Output::print<2>("primal solve: flow");
  // solve the primal system
  brinkman2d_primal.solve();
  brinkman2d_primal.output(n_calls);
  time.restart_and_print("solving the stationary flow problem");
  Output::print<2>("primal solve: temperature");
  
  TDatabase::TimeDB->TIME_DISC = 2; // Crank-Nicolson
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  double end_time = TDatabase::TimeDB->ENDTIME;
  cd2d_primal.reset_for_output();
  cd2d_primal.assemble_initial_time();
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    SetTimeDiscParameters(1);
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    TDatabase::TimeDB->CURRENTTIME += tau;
      
    Output::print<1>("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    
    cd2d_primal.assemble(brinkman2d_primal.get_velocity(), x,
                         db["diffusion_coefficient"]);
    cd2d_primal.solve();
    cd2d_primal.descale_stiffness(tau, TDatabase::TimeDB->THETA1);
    cd2d_primal.output();
  }
  time.stop_and_print("solving the instationary temperature problem");
  time.print_total_time("solving coupled forward problem");
}


double GeothermalPlantsPositionOptimization::compute_functional(const double* x)
 const
{
  Output::print("computing current value of functional");
  double distance = x[0];
  auto u = brinkman2d_primal.get_velocity();
  auto u1 = u.GetComponent(0);
  auto u2 = u.GetComponent(1);
  //auto p = brinkman2d_primal.get_pressure();
  auto temperature = cd2d_primal.get_function();
  double alpha = db["alpha_cost"];
  std::vector<double> cost_functional = db["cost_functional"];
  
  
  auto temperature_at_sink = 
  [&](){
    double mean = temperature.compute_mean();
    return mean;
  };
  
  
  double functional_value = 0.5 * alpha * distance;
  functional_value += 0.5 * cost_functional[0] / distance;
  if(n_calls > 1)
  {
    Output::print("difference to previous functional ", std::setprecision(14),
                  current_J_hat - functional_value);
  }
  
  delete u1;
  delete u2;
  return functional_value;
}


void GeothermalPlantsPositionOptimization::solve_adjoint_equation()
{
  auto u = brinkman2d_primal.get_velocity();
  auto p = brinkman2d_primal.get_pressure();
  std::vector<double> cost_functional_weights = db["cost_functional"];
  bool restricted_curl_functional = db["restricted_curl_functional"];
  
  Output::print<2>("adjoint solve ", n_computation_derivative);
  ErrThrow("solving the adjoint equation is not yet implemented");
//   brinkman2d_primal.assemble(u, p, *stokes_sol, cost_functional_weights,
//                        restricted_curl_functional);
//   brinkman2d_primal.solve();
//   brinkman2d_primal.output(n_calls);
}

void GeothermalPlantsPositionOptimization::compute_derivative(const double* x, 
                                                        double* grad) const
{
  double alpha = db["alpha_cost"];
  auto & adjoint_solution = brinkman2d_primal.get_solution();
  auto adjoint_residual = brinkman2d_primal.get_rhs(); // copy
  auto & adjoint_mat = brinkman2d_primal.get_matrix();
  adjoint_mat.BlockMatrix::apply_scaled_add(adjoint_solution, adjoint_residual, -1.);
  auto length = adjoint_solution.length(0);
  ErrThrow("computing the derivative is not yet implemented");
//   for(auto i = 0u; i < adjoint_solution.length(); ++i)
//   {
//     Output::print("i ", i, " \t sol ", adjoint_solution[i], " \t rhs ", 
//                   nse_adjoint.get_rhs()[i], " \t residual ", 
//                   adjoint_residual[i]);
//   }
}
