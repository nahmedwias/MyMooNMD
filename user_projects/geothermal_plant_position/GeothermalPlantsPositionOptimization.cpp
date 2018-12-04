#include "GeothermalPlantsPositionOptimization.hpp"

#include "LocalAssembling.h"
#include "NSE_local_assembling_routines.h"

#ifdef __2D__
#include "Assemble2D.h"
#include "gppo_example.hpp"
#else
#include "Assemble3D.h"
#include "gppo_example3D.hpp"
#endif
#include "LoopInfo.h"
#include "Database.h"
#include "TimeDiscRout.h"
#include <algorithm>
#include <BaseCell.h>
#include <ParMooN_repository_info.h>

#include <TimeDiscRout.h>
#include <Domain.h>
#include <Database.h>
#include "TimeConvectionDiffusion.h"
#include "NSE_GPPO.hpp"


const std::string path = parmoon::source_directory;
const std::string path_to_repo = path + "/user_projects/geothermal_plant_position/";

template <int d>
constexpr char GeothermalPlantsPositionOptimization<d>::required_database_name_TCD_GPPO[];





/** ************************************************************************ */
template <int d>
ParameterDatabase GeothermalPlantsPositionOptimization<d>::default_GPPO_database()
{
  // is only needed for the adjoint problem, but there is only one global 
  // parameter for this.
  TDatabase::ParamDB->NSTYPE = 4;
  auto db = ParameterDatabase::parmoon_default_database();
  db.merge(TDomain::default_domain_parameters());

  #ifdef __2D__
  db.merge(Example2D::default_example_database());
  #else
  db.merge(Example3D::default_example_database());
  #endif
  
  db.merge(LocalAssembling<d>::default_local_assembling_database());

  db.add("n_control",1u,
	 " Dimension of the control space",0u,10000u);
	
  
  db.add("alpha_cost", 1., 
          "The scalar alpha in the functional J_hat, which is multiplied with "
          "the control term.", 0., 1.0e10);
  db.add("cost_functional", {1., 0.},
          "Switch between a few possible cost functionals via weights. The "
          "order of the functionals is: L2_norm_of_curl, backward_facing_step, ",
          0., 1.);
  return db;
}



/** ************************************************************************ */
template <int d>
ParameterDatabase GeothermalPlantsPositionOptimization<d>::
get_primal_flow_database(ParameterDatabase param_db)
{
  std::string basename = param_db["output_basename"];
  basename += std::string("_flow");
  param_db["output_basename"].set(basename, false);

  ParameterDatabase NSE_GPPO_db =  NavierStokes<d>::default_nse_database();

  NSE_GPPO_db.add("variable_sigma_fct_type", false,
		  "Assemble the Brinkman/Darcy problem with space dependent"
		  " permeability.",
		  {true, false});

  NSE_GPPO_db.add("variable_sigma_fct_file", path_to_repo +
		  "default_coeff_fct.txt" ,
		  "Provide inverse of permeability from an input file. "
		  "Note: The File has to fit with the mesh (refinement level). "
		  "A default file  is contained in input_files/ .");

  NSE_GPPO_db.merge(param_db, false);

  return  NSE_GPPO_db;
}


#ifdef __2D__
/** ************************************************************************ */
template <int d>
ParameterDatabase GeothermalPlantsPositionOptimization<d>::
get_primal_temperature_database(ParameterDatabase param_db)
{
  std::string basename = param_db["output_basename"];
  basename += std::string("temperature");
  param_db["output_basename"].set(basename, false);

  // velocity solver database
  ParameterDatabase tcd_GPPO_db = param_db;
  auto db_name = std::string(GeothermalPlantsPositionOptimization<d>::required_database_name_TCD_GPPO);
  // use the given database or one of its nested databases, depending on which
  // one has the correct name. Otherwise the default solver database is used.
  cout << "HIER !!!!!"<< param_db.has_nested_database(db_name) << endl;
  cout << db_name << endl;
  
  if(param_db.has_nested_database(db_name))
    {
      tcd_GPPO_db.merge(param_db.get_nested_database(db_name), false);
    }
  return tcd_GPPO_db;
}
#else
/** ************************************************************************ */
template <int d>
ParameterDatabase GeothermalPlantsPositionOptimization<d>::get_primal_temperature_database(ParameterDatabase param_db)
{
  std::string basename = param_db["output_basename"];
  basename += std::string("temperature");
  param_db["output_basename"].set(basename, false);

  return param_db;
}
#endif

/** ************************************************************************ */
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

/** ************************************************************************ */
template <int d>
GeothermalPlantsPositionOptimization<d>::
GeothermalPlantsPositionOptimization(
				     const TDomain& domain, const ParameterDatabase& param_db)
  : db(default_GPPO_database()), n_control(0),
#ifdef __2D__
    brinkman_mixed(domain, get_primal_flow_database(param_db),
		      get_gppo_flow_example(param_db)),
    tcd_primal(domain, get_primal_temperature_database(param_db),
		 get_gppo_temperature_example(param_db)),
#else
    brinkman_mixed(domain, get_primal_flow_database(param_db),
		      get_3D_gppo_flow_example(param_db)),
    tcd3d_primal(domain, get_primal_temperature_database(param_db),
		 get_3D_gppo_temperature_example(param_db)),
#endif
    optimization_info("optimization", true, true, 1), // 1 -> full verbosity
    current_J_hat(std::numeric_limits<double>::infinity()),
    control_old(), n_calls(0),
    n_computation_derivative(0)
{
  Output::print<5>("Creating the GeothermalPlantsPositionOptimization object");
  db.merge(param_db, false);

#ifdef __2D__
  bool variable_sigma = get_primal_flow_database(param_db)["variable_sigma_fct_type"];
  brinkman_mixed.assemble_with_coefficient_fct(variable_sigma);
#else
  brinkman_mixed.assemble_linear_terms();
#endif
  
  /*
    // remove this and redefine the function brinkman_mixed.assemble_with_coefficient_fct
  if (get_primal_flow_database(param_db)["variable_sigma_fct_type"]) {
    brinkman_mixed.
      assemble_with_coefficient_fct(& brinkman_mixed.get_coefficient_function());
  } else {
    brinkman_mixed.assemble_with_coefficient_fct(nullptr);
  }
  */

  if(!(param_db["problem_type"].is(3)) && !(param_db["problem_type"].is(7)))
  {
    Output::print<1>("ERROR: NSE3D_GPPO::assemble() is not considering nonlinear terms ",
		     " but your chosen problem_type ",
		     param_db["problem_type"] , " does! ");
    return;
  }
  
  n_control = param_db["n_control"];
  control_old = std::vector<double>(n_control, 0.0);
  Output::print<3>("Created the GeothermalPlantsPositionOptimization object, ",
          "n_control = ", n_control);
}

/** ************************************************************************ */
template <int d>
double GeothermalPlantsPositionOptimization<d>::compute_functional_and_derivative(
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

/** ************************************************************************ */
void approximate_delta_functions(int n_points, double *x, double *y,
        double **parameters, double **coeffs,
        double distance)
{
  for(int i = 0; i < n_points; ++i)
  {

    coeffs[i][0] = 0.;
    coeffs[i][1] = 0.;
    coeffs[i][2] = 0.;
    coeffs[i][3] = 0.;
    
    // Note: this function should be consistent with the used domain/mesh
    double r_well = 0.2; // 20cm
    double epsDelta = 25*r_well;
    double Qin = 150./3600.;


    std::vector<double> singular_x,singular_y,singular_sign;
    // source in (xi,yi)
    singular_x.push_back(4500.);
    singular_y.push_back(3000.);
    singular_sign.push_back(1.);
    // sink in (xe,ye)    
    singular_x.push_back(5500.);
    singular_y.push_back(3000.);
    singular_sign.push_back(-1.);
    
    /*
    //another approx. for domain [0, 10] x [0, 6]
    double a = 5.; //0.05;
    double r_Wells = 5*0.001;
    double Q_in = 0.03*Pi*r_Wells*r_Wells * 0.1;
    //cout <<"distance  in brinkman!!!!"<< distance << endl;
    */
    
    //std::array<double, 2> center_source = {{5000.0 - distance/2., 3000.0}}; //{{5.0 - distance/2., 3}};
    //std::array<double, 2> center_sink = {{5000.0 + distance/2., 3000.0}}; //{{5.0 + distance/2., 3}};
    //cout <<" center_source[0]) in brinkman!!!!"<<  center_source[0] << endl;
    //cout <<" center_source[1]) in brinkman!!!!"<<  center_source[1]<< endl;

    //cout <<" center_sink[0]) in brinkman!!!!"<<  center_sink[0] << endl;
    //cout <<" center_sink[1]) in brinkman!!!!"<<  center_sink[1]<< endl;

    for (unsigned int m=0; m<singular_x.size(); m++) {
      double x_center_source = singular_x[m];
      double y_center_source = singular_y[m];
      double x_distance_to_source = std::pow(std::abs(x[i] - x_center_source), 2);
      double y_distance_to_source = std::pow(std::abs(y[i] - y_center_source), 2);
      bool at_source =(x_distance_to_source < epsDelta*epsDelta) *
	(y_distance_to_source < epsDelta*epsDelta);
    

	if(at_source)
	  {
	    Output::print<4>(" adding a singular source/sink - point ", m);
	    double magnitude = cos(Pi*(x[i] - x_center_source)/epsDelta) + 1;
	    magnitude *= cos(Pi*(y[i] - y_center_source)/epsDelta) + 1;
	    magnitude /= 4.*epsDelta*epsDelta;
	    coeffs[i][3] += singular_sign[m] * magnitude * Qin;
	  }
    }
  }
}

/** ************************************************************************ */
/*// For 3 wells:
void approximate_delta_functions_3(int n_points, double *x, double *y,
        double **parameters, double **coeffs,
        double x_i, double x_p1, double x_p2, double y_i, double y_p1, double y_p2 )
{
  for(int i = 0; i < n_points; ++i)
  {
    //another approx. for domain [0, 10] x [0, 6]
    double a = 0.25;
    double r_Wells = 5*0.001;
    double Q_in = 0.03*Pi*r_Wells*r_Wells * 0.1;
    std::array<double, 2> center_source = {{x_i, y_i}};
    std::array<double, 2> center_sink_1 = {{x_p1, y_p1}};
    std::array<double, 2> center_sink_2 = {{x_p2, y_p2}};
    double x_distance_to_source = std::pow(std::abs(x[i]-center_source[0]), 2);
    double y_distance_to_source = std::pow(std::abs(y[i]-center_source[1]), 2);
    double x_distance_to_sink_1 = std::pow(std::abs(x[i]-center_sink_1[0]), 2);
    double y_distance_to_sink_1 = std::pow(std::abs(y[i]-center_sink_1[1]), 2);
    double x_distance_to_sink_2 = std::pow(std::abs(x[i]-center_sink_2[0]), 2);
    double y_distance_to_sink_2 = std::pow(std::abs(y[i]-center_sink_2[1]), 2);
    bool at_source = x_distance_to_source + y_distance_to_source < a*a;
    bool at_sink_1 = x_distance_to_sink_1 + y_distance_to_sink_1 < a*a;
    bool at_sink_2 = x_distance_to_sink_2 + y_distance_to_sink_2 < a*a;
    coeffs[i][0] = 0.;
    coeffs[i][1] = 0.;
    coeffs[i][2] = 0.;
    coeffs[i][3] = 0.;

    if(at_source)
    {
      double magnitude = cos(Pi*x_distance_to_source/a) + 1;
      magnitude *= cos(Pi*y_distance_to_source/a) + 1;
      magnitude /= 4.*a*a;
      //coeffs[i][3] += magnitude * Q_in; // source
      coeffs[i][3] += 0.001; // source
    }
    if(at_sink_1)
    {
      double magnitude_1 = cos(Pi*x_distance_to_sink_1/a) + 1;
      magnitude_1 *= cos(Pi*y_distance_to_sink_1/a) + 1;
      magnitude_1 /= 4.*a*a;
      //coeffs[i][3] -=  magnitude_1 * Q_in; // sink
      coeffs[i][3] -=  0.001; // sink
    }
    if(at_sink_2)
    {
      double magnitude_2 = cos(Pi*x_distance_to_sink_2/a) + 1;
      magnitude_2 *= cos(Pi*y_distance_to_sink_2/a) + 1;
      magnitude_2 /= 4.*a*a;
      //coeffs[i][3] -=  magnitude_2 * Q_in; // sink
      coeffs[i][3] -=  0.001; // sink
    }
  }
}
*/

/** ************************************************************************ */
template <int d>
void GeothermalPlantsPositionOptimization<d>::apply_control_and_solve(const double* x)
{
  Chrono time;
  double t_start = GetTime();
  // apply control x
  double distance = x[0];
  Output::print("current control: ", std::setprecision(14), distance);
  using namespace std::placeholders;  
#ifdef __2D__
  CoeffFct2D coeff = std::bind(approximate_delta_functions,
			       _1, _2, _3, _4, _5,distance);
  // try to generalize to more wells
  //   CoeffFct2D coeff = std::bind(approximate_delta_functions_3,
  //_1, _2, _3, _4, _5, x_i, x_p1, x_p2, y_i, y_p1, y_p2);
#else
  CoeffFct3D coeff = std::bind(approximate_delta_functions,
			       _1, _2, _3, _4, _5, _6, distance);
#endif
  
  std::string disc_type = this->db["space_discretization_type"];
  bool nonsymm_gls = (disc_type == std::string("nonsymm_gls"));
  int rhs_div_sign = 1;
  if (nonsymm_gls)
  {
    rhs_div_sign = -1;
  }

#ifdef __2D__
  LocalAssembling<d> la(2, {D00, D00}, {0, 1}, {}, {},
			{0, 0, 1}, coeff,
			std::bind(NSRightHandSide<2>, _1, _2, _3, _4, _5, _6, _7, _8,
				  rhs_div_sign),nullptr,
			0, 3, 0, {}, {}, 0, nullptr,
			0, {}, {});
#else
  LocalAssembling<d> la(2, {D000, D000},  {0, 1},  {},  {},
			{0, 0, 0, 1}, coeff,
			std::bind(NSRightHandSide<3>, _1, _2, _3, _4, _5, _6, _7, _8,
				  rhs_div_sign),nullptr,
			0, 4, 0,  {},  {},   0,   nullptr,  0,  {}, {});
#endif

  auto& v_space = brinkman_mixed.get_velocity_space();
  auto& p_space = brinkman_mixed.get_pressure_space();
  auto& rhs = brinkman_mixed.get_rhs();
  rhs.reset();

#ifdef __2D__
  std::array<const TFESpace2D*, 3> fespaces = {{&v_space, &v_space, &p_space}};
  double *rhs_pointers[3] = {rhs.block(0), rhs.block(1), rhs.block(2)};
  BoundCondFunct2D * boundary_conditions[3] = {
          v_space.get_boundary_condition(), v_space.get_boundary_condition(),
          p_space.get_boundary_condition() };
  auto& example = brinkman_mixed.get_example();
  std::array<BoundValueFunct2D*, 3> non_const_bound_values;
  non_const_bound_values[0] = example.get_bd()[0];
  non_const_bound_values[1] = example.get_bd()[1];
  non_const_bound_values[2] = example.get_bd()[2];
  
#else
  std::array<const TFESpace3D*, 4> fespaces = {{&v_space, &v_space, &v_space, &p_space}};
  double *rhs_pointers[4] = {rhs.block(0), rhs.block(1), rhs.block(2), rhs.block(3)};
  BoundCondFunct3D * boundary_conditions[4] = {
      v_space.get_boundary_condition(), v_space.get_boundary_condition(),
      v_space.get_boundary_condition(), p_space.get_boundary_condition() };
  auto& example = brinkman_mixed.get_example();
  std::array<BoundValueFunct3D*, 4> non_const_bound_values;
  non_const_bound_values[0] = example.get_bd()[0];
  non_const_bound_values[1] = example.get_bd()[1];
  non_const_bound_values[2] = example.get_bd()[2];
  non_const_bound_values[3] = example.get_bd()[3];

#endif
  
#ifdef __2D__
  Assemble2D(
#else
  Assemble3D(
#endif
	     2, fespaces.data()+1, 0, nullptr, 0, nullptr, 3, rhs_pointers,
	     fespaces.data(), boundary_conditions,
	     non_const_bound_values.data(), la);
  
#ifdef __2D__
  double mat_norm 
  = brinkman_mixed.get_matrix().get_combined_matrix()->GetNorm();
  double rhs_norm = rhs.norm();
  Output::print("---- norms: ", mat_norm, "  ", rhs_norm);
  
  Output::print<2>("primal solve: flow");
  // solve the primal system
  brinkman_mixed.solve();
  // brinkman_mixed.get_pressure().WriteSol("/Users/blank/ParMooN/Tests/Geothermal_Plants_Position_Optimization/", "written_fe_function_pressure");
  
  brinkman_mixed.output(n_calls);
  time.restart_and_print("solving the stationary flow problem");
  Output::print<2>("primal solve: temperature");
  
  TDatabase::TimeDB->TIME_DISC = 2; // Crank-Nicolson
  //double end_time = TDatabase::TimeDB->ENDTIME;
  tcd_primal.reset_for_output();
  
  TimeDiscretization& tss = tcd_primal.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  tcd_primal.assemble_initial_time();
  
  double start_time = this->db["time_start"];
  double end_time   = this->db["time_end"];
  TDatabase::TimeDB->CURRENTTIME = start_time; //TDatabase::TimeDB->STARTTIME;
  
  tcd_primal.output();

#else
  
  double mat_norm
    = brinkman_mixed.get_matrix().get_combined_matrix()->GetNorm();
  double rhs_norm = rhs.norm();
  Output::print("---- norms: ", mat_norm, "  ", rhs_norm);

  Output::print<2>("primal solve: flow");
  // solve the primal system
  brinkman_mixed.solve();
  brinkman_mixed.output(n_calls);
  time.restart_and_print("solving the stationary flow problem");
  Output::print<2>("primal solve: temperature");

  cout << "Before reset_for_output()"<<endl;
  tcd_primal.reset_for_output();
  TimeDiscretization& tss = tcd_primal.get_time_stepping_scheme();
     tss.current_step_ = 0;
     tss.set_time_disc_parameters();
   tcd_primal.assemble_initial_time();

   double start_time = this->db["time_start"];
   double end_time   = this->db["time_end"];
   TDatabase::TimeDB->CURRENTTIME = start_time; //TDatabase::TimeDB->STARTTIME;

   tcd_primal.output();
#endif

  
   // ======================================================================
   // time iteration
   // ======================================================================
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    tss.current_step_++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    tss.set_time_disc_parameters();
    double tau = this->db["time_step_length"];
    TDatabase::TimeDB->CURRENTTIME += tau;

    Output::print<1>("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    SetTimeDiscParameters(1);

#ifdef __2D__
    tcd_primal.assemble(brinkman_mixed.get_velocity(), x,
			  db["diffusion_coefficient"]);
    tcd_primal.solve();
    if((tss.current_step_-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      tcd_primal.output();
#else
    tcd_primal.assemble(brinkman_mixed.get_velocity(), x,
			  db["diffusion_coefficient"]);
    tcd_primal.solve();
    if((tss.current_step_-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      tcd_primal.output();
#endif
  }
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();

}

/** ************************************************************************ */
template <int d>
double GeothermalPlantsPositionOptimization<d>::compute_functional(const double* x)
const
{
  Output::print("computing current value of functional");
  double distance = x[0];
  // For three wells:
  // double x_i =   ;
  // double x_p1 =  ;
  // double x_p2 =  ;
  // double y_i =  ;
  // double y_p1 =   ;
  // double y_p2 =  ;
  auto u = brinkman_mixed.get_velocity();
  auto u1 = u.GetComponent(0);
  auto u2 = u.GetComponent(1);
#ifdef __3D__
  auto u3 = u.GetComponent(2);
#endif
  auto pressure = brinkman_mixed.get_pressure();

  auto temperature = tcd_primal.get_function();

  double alpha = db["alpha_cost"];

  /*
  //New LB 19.11.18 START
  double pressure_prod[1], pressure_inj[1], temperature_prod[1];
  double x_prod = 5000.0 + distance/2.; //5.5;
  double y_prod = 3000.; //3.;
  double x_inj = 5000.0 - distance/2.; //4.5;
  double y_inj = 3000.; //3.;


  TBaseCell *cell;
  TBaseCell *cell_prod, *cell_inj;
  int cell_prod_num, cell_inj_num;


  TCollection *coll_pressure, *coll_temperature;
  coll_pressure = brinkman_mixed.get_pressure_space().GetCollection();
  coll_temperature = tcd_primal.get_space().GetCollection();

  for (int j = 0; j < coll_pressure->GetN_Cells(); j++)
    {
      cell = coll_pressure->GetCell(j);
      if (cell->PointInCell(x_prod, y_prod))
      {
        cell_prod = cell;
        cell_prod_num = j;
      }
      else if (cell->PointInCell(x_inj, y_inj))
      {
       cell_inj = cell;
       cell_inj_num = j;
      }
    }


  pressure.FindValueLocal(cell_prod, cell_prod_num, x_prod, y_prod, pressure_prod);
  pressure.FindValueLocal(cell_inj, cell_inj_num, x_inj, y_inj, pressure_inj);

  temperature.FindValueLocal(cell_prod, cell_prod_num, x_prod, y_prod, temperature_prod);
  //temperature->FindValueLocal (const TBaseCell *cell, int cell_no, x_inj, inj, values);


  double temperature_inj = 303.15; // = 30 + 273.15;
  double Q = 150/360;//24 * 50; // 50 - 300
  double dens = 1050;
  double Cap = 4200; // 4200/(24*3600*24*3600);
  double functional_value_new = Q/(0.6) * alpha * (pressure_prod[1] - pressure_inj[1])  -  Q * dens * Cap * (temperature_prod[1] - temperature_inj); // Net energy AFTER 50 years
 //New LB 19.11.18 END

 */

  std::vector<double> cost_functional = db["cost_functional"];


  auto temperature_at_sink = 
          [&](){
    double mean = temperature.compute_mean();
    return mean;
  };

  /*auto temperature_at_sink1 =
          [&](){
    double mean = temperature.compute_mean();
    return mean;
  };
  auto temperature_at_sink2 =
          [&](){
    double mean = temperature.compute_mean();
    return mean;
  };
  */

  // 3Wells
  //  double functional_value = 0.5 * alpha * distance;
  // functional_value += 0.5 * cost_functional[0] / distance;
  double functional_value = (1000/1000001)  * alpha * distance; //0.5 * alpha * distance;
  functional_value +=  (1000/1000001)  * cost_functional[0] / distance; //0.5 * cost_functional[0] / distance;
  if(n_calls > 1)
  {
    Output::print("difference to previous functional ", std::setprecision(14),
            current_J_hat - functional_value);
  }

  delete u1;
  delete u2;
#ifdef __3D__
  delete u3;
#endif  
  return functional_value; //functional_value_new;
}

/** ************************************************************************ */
template <int d>
void GeothermalPlantsPositionOptimization<d>::solve_adjoint_equation()
{
  auto u = brinkman_mixed.get_velocity();
  auto p = brinkman_mixed.get_pressure();
  std::vector<double> cost_functional_weights = db["cost_functional"];
  bool restricted_curl_functional = db["restricted_curl_functional"];

  Output::print<2>("adjoint solve ", n_computation_derivative);
  ErrThrow("solving the adjoint equation is not yet implemented");
  //   brinkman_mixed.assemble(u, p, *stokes_sol, cost_functional_weights,
  //                        restricted_curl_functional);
  //   brinkman_mixed.solve();
  //   brinkman_mixed.output(n_calls);
}

/** ************************************************************************ */
template <int d>
void GeothermalPlantsPositionOptimization<d>::compute_derivative(const double* x,
        double* grad) const
{
  double alpha = db["alpha_cost"];
  auto & adjoint_solution = brinkman_mixed.get_solution();
  auto adjoint_residual = brinkman_mixed.get_rhs(); // copy
  auto & adjoint_mat = brinkman_mixed.get_matrix();
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






#ifdef __3D__
template class GeothermalPlantsPositionOptimization<3>;
#else
template class GeothermalPlantsPositionOptimization<2>;
#endif
