#include "GeothermalPlantsPositionOptimization.hpp"

#include "LocalAssembling.h"
#include "NSE_local_assembling_routines.h"

#ifdef __2D__
#include "Assemble2D.h"
#else
#include "Assemble3D.h"
#endif

#include "gppo_example.hpp"
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

#include <fstream>
#include <iterator>

#include <BaseCell.h>
#include <Collection.h>

const std::string path = parmoon::source_directory;
const std::string path_to_repo = path + "/user_projects/geothermal_plant_position/";

template <int d>
constexpr char GeothermalPlantsPositionOptimization<d>::required_database_name_TCD_GPPO[];


/** ************************************************************************ */
template <int d>
ParameterDatabase GeothermalPlantsPositionOptimization<d>::default_GPPO_database()
{
  auto db = ParameterDatabase::parmoon_default_database();
  db.merge(TDomain::default_domain_parameters());

#ifdef __2D__
  db.merge(Example2D::default_example_database());
#else
  db.merge(Example3D::default_example_database());
#endif
  db.merge(TimeConvectionDiffusion<d>::default_tcd_database());
  db.merge(LocalAssembling<d>::default_local_assembling_database());

  db.add("n_control",1u,
	 " Dimension of the control space.",
	 0u, 10000u);
  
  db.add("alpha_cost", 1., 
          "The scalar alpha in the functional J_hat, which is multiplied with "
          "the control term.",
          0., 1.0e10);
  
  db.add("cost_functional", {1., 0.},
          "Switch between a few possible cost functionals via weights. The "
          "order of the functionals is: L2_norm_of_curl, backward_facing_step, ",
          0., 1.);
  
  db.add("longitudinal_dispersion_factor", 0.,
       "Use thermal dispersion tensor in the thermal conductivity",
       -10000000., 10000000.);

  db.add("transversal_dispersion_factor", 0.,
       "Use thermal dispersion tensor in the thermal conductivity",
       -10000000., 10000000.);
  
  db.add("temperature_injection_well", 303.15, // = 30 °C + 273.15;
          "The temperature of the injected fluid in Kelvin",
          -400., 400.);
  
  db.add("fluid_density", 1050.,
          "The density of the fluid (e.g., brine)",
          0., 10000.);
  
  db.add("fluid_heat_capacity", 4200.,
          "The heat capacity of the injected fluid (e.g., brine)",
         0., 10000.);
  
  db.add("pump_efficiency", 1.,
          "The efficiency of the pumps",
          0., 1.);
  
  db.add("minimum_temperature_production_well",  343.15, // = 70 °C + 273.15;
          "The temperature of the fluid at the production well necessery to "
          "continue with production",
          0., 500.);
  
  //todo: this should be rather the minimum temperature on a circle with
  //      certain distance to the center, since the cold water fronts might
  //      be non-symmetric
  db.add("x_distance_from_well_center",  20., // meters
          "The distance from the center of the injection (+) and "
          "production (-) well, used to compute temperatures and pressures.",
          0., 500.);
  
  db.add("well_radius",  0.1, // meters
          "The radius of the bore holes.",
          0., 100.);
  
  db.add("delta_fct_eps_factor", 1., // epsDelta = delta_fct_eps_factor * well_radius
          "The scaling of the delta-distribution, given by "
          "epsDelta = delta_fct_eps_factor * db[well_radius].",
          0., 10000000.);

  //2 Doublets & 3 Double Doublets
  db.add("well_distance", 1000.,
          "Determines the distance between the injection and production well "
          "for both doublets",
          0., 10000000.);

  db.add("u_in", 1e-5,
          "Sets the normal velocity (outflow) w.r.t. a circular source (2D)"
          " representing the injection well",
          0., 10000000.);

  db.add("scenario", "1doublet_optimize_distance",
           "Determine the scenario to be optimized.",
           {"1doublet_optimize_distance", "2doublets_fixed_well_distance",
                   "3_rows_of_double_doublets_fixed_well_distance",
                   "3_rows_of_double_doublets_varying_row_distance",
                   "lattice", "hexagon"});
  db.add("lattice_n_wells_per_row", 4u , "Parameter for the scenario 'lattice'."
        " Determines the number of wells per row and column. "
        "Has to be an even number since (lattice_n_wells_per_row*lattice_n_wells_per_row)/2 production "
        "resp. injection wells are assumed to be present.", 0u, 1000000u);

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

  /*
  NSE_GPPO_db.add("variable_sigma_fct_type", false,
		  "Assemble the Brinkman/Darcy problem with space dependent"
		  " permeability.",
		  {true, false});

  NSE_GPPO_db.add("variable_sigma_fct_file", path_to_repo +
		  "default_coeff_fct.txt" ,
		  "Provide inverse of permeability from an input file. "
		  "Note: The File has to fit with the mesh (refinement level). "
		  "A default file  is contained in input_files/ .");
  */
  
  NSE_GPPO_db.merge(param_db, false);

  return  NSE_GPPO_db;
}

/** ************************************************************************ */
template <int d>
ParameterDatabase GeothermalPlantsPositionOptimization<d>::
get_primal_temperature_database(ParameterDatabase param_db)
{
  std::string basename = param_db["output_basename"];
  basename += std::string("_temperature");
  param_db["output_basename"].set(basename, false);

  ParameterDatabase tcd_GPPO_db = TimeConvectionDiffusion<d>::default_tcd_database();

  // temperature solver database
  // ParameterDatabase tcd_GPPO_db = param_db;
  auto db_name = std::string(GeothermalPlantsPositionOptimization<d>::
			     required_database_name_TCD_GPPO);
  // use the given database or one of its nested databases, depending on which
  // one has the correct name. Otherwise the default solver database is used.
  cout << "Nested database, name: " << db_name << endl;

  tcd_GPPO_db.merge(param_db, true); //OLD 15.01.19:tcd_GPPO_db.merge(param_db,false);

  // replace parameter values if some parameter is specific for tcd problem
  if(param_db.has_nested_database(db_name))
    {
      tcd_GPPO_db.merge(param_db.get_nested_database(db_name), false);      
    }

  return tcd_GPPO_db;
}

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
GeothermalPlantsPositionOptimization(const TDomain& domain,
				     const ParameterDatabase& param_db)
  : db(default_GPPO_database()), n_control(0),
    brinkman_mixed(domain, get_primal_flow_database(param_db),
		   get_gppo_flow_example(param_db)),
    tcd_primal(domain, get_primal_temperature_database(param_db),
	       get_gppo_temperature_example(param_db)),
    optimization_info("optimization", true, true, 1), // 1 -> full verbosity
    current_J_hat(std::numeric_limits<double>::infinity()),
    control_old(), n_calls(0),
    n_computation_derivative(0)
{
  Output::print<5>("Creating the GeothermalPlantsPositionOptimization object");
  db.merge(param_db, false);
  
#ifdef __2D__
  bool variable_sigma = false;
  //get_primal_flow_database(param_db)["variable_sigma_fct_type"];
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
  }

  if(n_calls == 1)
  {
    optimization_info.restart("optimization", current_J_hat);
  }
  optimization_info.print(n_calls, current_J_hat);
  return current_J_hat;
}

/** ************************************************************************ */
void approximate_delta_functions_hexagon(int n_points,
					 double *x, double *y,
#ifdef __3D__
					 double *z,
#endif
					 double **parameters, double **coeffs,
					 // problem dependent
					 double u_in,double distance,
					 double well_radius,
					 double delta_fct_eps_factor)
{
  int dim = 2;
#ifdef __3D__
  dim = 3;
#endif
    
  double domain_Lx = 10000.;
  double domain_Ly = 6000.;
  size_t n_wells = 6;
  
  // set well position
  std::vector<double> singular_x, singular_y, singular_sign, flow_rate;

  singular_sign.push_back(-1.);
  singular_sign.push_back(1.);
  singular_sign.push_back(1.);
  singular_sign.push_back(-1.);
  singular_sign.push_back(1);
  singular_sign.push_back(1);

  double Qin = u_in * (2./well_radius); //(dim/well_radius);
  flow_rate.push_back(Qin);
  flow_rate.push_back(Qin/2.);
  flow_rate.push_back(Qin/2.);
  flow_rate.push_back(Qin);
  flow_rate.push_back(Qin/2.);
  flow_rate.push_back(Qin/2.);
  
  double epsDelta = delta_fct_eps_factor * well_radius;
  
  double x0 = domain_Lx/2.;
  double y0 = domain_Ly/2.;
  double pi = acos(-1.);
  int count = 0;
  for (unsigned int k1=0; k1 < n_wells; k1++) {
    double xk = x0 + distance*cos(2.*pi*k1/n_wells);
    double yk = y0 + distance*sin(2.*pi*k1/n_wells);
    singular_x.push_back(xk);
    singular_y.push_back(yk);
  }
  
    
  for (int i = 0; i < n_points; ++i)
  {
    for (size_t k = 0; k <= dim+1; k++)
    {
      coeffs[i][k] = 0.;
    }

    

    
    for (unsigned int m = 0; m < singular_x.size(); m++)
    {
      double x_center_source = singular_x[m];
      double y_center_source = singular_y[m];
      
      double x_distance_to_source = std::pow(std::abs(x[i] - x_center_source), 2);
      double y_distance_to_source = std::pow(std::abs(y[i] - y_center_source), 2);
      
      bool at_source = (x_distance_to_source < epsDelta*epsDelta) *
	(y_distance_to_source < epsDelta*epsDelta);

      if (at_source)
      {
        double magnitude = cos(Pi*(x[i] - x_center_source)/epsDelta) + 1;
        magnitude *= cos(Pi*(y[i] - y_center_source)/epsDelta) + 1;
	magnitude /= 4.*epsDelta*epsDelta;
	
	coeffs[i][dim+1] += singular_sign[m] * magnitude * flow_rate[m];
        Output::print<4>(" adding a singular source/sink - point ", m,
                " coeff[", dim+1, "] = ", coeffs[i][dim+1]);
      }
    }
  }
  
}
/** ************************************************************************ */
void approximate_delta_functions_lattice(int n_points,
					 double *x, double *y,
#ifdef __3D__
					 double *z,
#endif
					 double **parameters, double **coeffs,
					 // problem dependent
					 double u_in,double distance,
					 double well_radius,
					 double delta_fct_eps_factor,
					 size_t lattice_n_wells_per_row)
{
  int dim = 2;
#ifdef __3D__
  dim = 3;
#endif
    
  double domain_Lx = 10000.;
  double domain_Ly = 6000.;

  // has to be an even number
  size_t n_wells_per_row = lattice_n_wells_per_row; //equals number of wells per column

  // set well position
  std::vector<double> singular_x, singular_y, singular_sign;
  double x0 = domain_Lx/2.- (n_wells_per_row-1)/2.*distance;
  double y0 = domain_Ly/2.- (n_wells_per_row-1)/2.*distance;
  int count = 0;
  for (unsigned int k1=0; k1 < n_wells_per_row; k1++) {
    for (unsigned int k2=0; k2 < n_wells_per_row; k2++) {
      double xk = x0 + distance*k2;
      double yk = y0 + distance*k1;
      singular_sign.push_back(pow(-1,k1)*pow(-1,k2));
      singular_x.push_back(xk);
      singular_y.push_back(yk);
      Output::print<4>("well= ", singular_x.size(), " ", xk,",",yk,
		       " sign = ",pow(-1,k1)*pow(-1,k2));
    }
  }
  
    
  for (int i = 0; i < n_points; ++i)
  {
    for (size_t k = 0; k <= dim+1; k++)
    {
      coeffs[i][k] = 0.;
    }

    double epsDelta = delta_fct_eps_factor * well_radius;
    double Qin = u_in * (2./well_radius); //(dim/well_radius);

    
    for (unsigned int m = 0; m < singular_x.size(); m++)
    {
      double x_center_source = singular_x[m];
      double y_center_source = singular_y[m];
      
      double x_distance_to_source = std::pow(std::abs(x[i] - x_center_source), 2);
      double y_distance_to_source = std::pow(std::abs(y[i] - y_center_source), 2);
      
      bool at_source = (x_distance_to_source < epsDelta*epsDelta) *
	(y_distance_to_source < epsDelta*epsDelta);

      if (at_source)
      {
        double magnitude = cos(Pi*(x[i] - x_center_source)/epsDelta) + 1;
        magnitude *= cos(Pi*(y[i] - y_center_source)/epsDelta) + 1;
	magnitude /= 4.*epsDelta*epsDelta;
	
	coeffs[i][dim+1] += singular_sign[m] * magnitude * Qin;
        Output::print<4>(" adding a singular source/sink - point ", m,
                " coeff[", dim+1, "] = ", coeffs[i][dim+1]);
      }
    }
  }
  
}

void approximate_delta_functions(int n_points, double *x, double *y,
#ifdef __3D__
        double *z,
#endif
        double **parameters, double **coeffs, double u_in,
        double distance, double well_radius, double delta_fct_eps_factor)
{
int dim = 2;
#ifdef __3D__
  dim = 3;
#endif

  for (int i = 0; i < n_points; ++i)
  {
    for (size_t k = 0; k <= dim+1; k++)
    {
      coeffs[i][k] = 0.;
    }

    double epsDelta = delta_fct_eps_factor * well_radius;  // ~h
    /*double H = 75.;
    double Volume = epsDelta*epsDelta*Pi*H;  //well_radius*well_radius*Pi*H;
    double Qin =  150./3600./Volume; // m^3/h thickness = 1000m
     */
    double Qin = u_in * (2/well_radius); //(dim/well_radius);

#ifdef __3D__
    Qin = u_in * (2/well_radius/well_radius);
#endif
    
    std::vector<double> singular_x, singular_y, singular_z,singular_sign;

    double z0 = 50;
    double dz = 50;
    // source in (xi,yi)
#ifdef __3D__
    //note: this should match the depth of the actual mesh
    // [0,500]
    for (unsigned int k=0; k<9; k++) {
      singular_x.push_back(5000. - (distance)/2.);
      singular_y.push_back(3000.);
      singular_sign.push_back(1.);singular_z.push_back(z0);
      z0 = z0 + dz;
    }


#else
    // in 2D, we only need a (dummy) point along z for each well
    singular_x.push_back(5000. - (distance)/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(1.);
    singular_z.push_back(0.);
#endif


    // sink in (xe,ye)    
#ifdef __3D__
    // note: this should match the depth of the actual mesh
    // [0,500]
    z0 = 50.;
    for (unsigned int k=0; k<9; k++) {
      singular_x.push_back(5000.+ (distance)/2.);
      singular_y.push_back(3000.);
      singular_sign.push_back(-1.);
      singular_z.push_back(z0);
      z0 = z0 + dz;
    }
    
#else
    singular_x.push_back(5000.+ (distance)/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(-1.);
    // in 2D, we only need a (dummy) point along z for each well
    singular_z.push_back(0.);
#endif

    for (unsigned int m = 0; m < singular_x.size(); m++)
    {
      double x_center_source = singular_x[m];
      double y_center_source = singular_y[m];
      double z_center_source = singular_z[m];//0-500

      double x_distance_to_source = std::pow(std::abs(x[i] - x_center_source), 2);
      double y_distance_to_source = std::pow(std::abs(y[i] - y_center_source), 2);
      
      double z_distance_to_source = 1.;
#ifdef __3D__
      z_distance_to_source = std::pow(std::abs(z[i] - z_center_source), 2);
#endif

      bool at_source = (x_distance_to_source < epsDelta*epsDelta) *
	(y_distance_to_source < epsDelta*epsDelta)*
	(z_distance_to_source < epsDelta*epsDelta);

      if (at_source)
      {
        double magnitude = cos(Pi*(x[i] - x_center_source)/epsDelta) + 1;
        magnitude *= cos(Pi*(y[i] - y_center_source)/epsDelta) + 1;
	magnitude /= 4.*epsDelta*epsDelta;
	
#ifdef __3D__
	magnitude *= cos(Pi*(z[i] - z_center_source)/epsDelta) + 1;
	magnitude /= epsDelta;
#endif

        coeffs[i][dim+1] += singular_sign[m] * magnitude * Qin;
        Output::print<4>(" adding a singular source/sink - point ", m,
                " coeff[", dim+1, "] = ", coeffs[i][dim+1]);
      }
    }
  }
}


/** ************************************************************************ */
void approximate_delta_functions_2doublets(int n_points, double *x, double *y,
#ifdef __3D__
        double *z,
#endif
        double **parameters, double **coeffs, double u_in,
        double center_x_moving_doublet, double well_radius,
        double delta_fct_eps_factor, double well_distance)
{
int dim = 2;
#ifdef __3D__
  dim = 3;
#endif

  for (int i = 0; i < n_points; ++i)
  {
    for (size_t k = 0; k <= dim+1; k++)
    {
      coeffs[i][k] = 0.;
    }

    double epsDelta = delta_fct_eps_factor * well_radius;  // ~h
    /*double H = 75.;
    double Volume = epsDelta*epsDelta*Pi*H;  //well_radius*well_radius*Pi*H;
    double Qin =  150./3600./Volume; // m^3/h thickness = 1000m
     */
    double Qin = u_in * (2/well_radius); //(dim/well_radius);

    std::vector<double> singular_x, singular_y,
#ifdef __3D__
    singular_z,
#endif
    singular_sign;

    // source in (xi,yi) for fixed doublet
    singular_x.push_back(5000. - (well_distance)/2.);
    singular_y.push_back(2500.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    // sink in (xe,ye) for fixed doublet
    singular_x.push_back(5000.+ (well_distance)/2.);
    singular_y.push_back(2500.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    // source in (xi,yi) for moving doublet
    singular_x.push_back(center_x_moving_doublet - (well_distance)/2.);
    singular_y.push_back(3500.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    // sink in (xe,ye) for moving doublet
    singular_x.push_back(center_x_moving_doublet + (well_distance)/2.);
    singular_y.push_back(3500.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    for (unsigned int m = 0; m < singular_x.size(); m++)
    {
      double x_center_source = singular_x[m];
      double y_center_source = singular_y[m];
#ifdef __3D__
      double z_center_source = singular_z[m];//0-500
#endif
      double x_distance_to_source = std::pow(std::abs(x[i] - x_center_source), 2);
      double y_distance_to_source = std::pow(std::abs(y[i] - y_center_source), 2);
      /*
#ifdef __3D__
      double z_distance_to_source = std::pow(std::abs(z[i] - z_center_source), 2);
#endif
       */

      bool at_source = (x_distance_to_source < epsDelta*epsDelta) *
              (y_distance_to_source < epsDelta*epsDelta)
              /*
#ifdef __3D__
               * (z_distance_to_source < epsDelta*epsDelta)
#endif
               */
              ;

      if (at_source
#ifdef __3D__
              && z[i] < 490 && z[i] > 10
#endif
      )
      {
        double magnitude = cos(Pi*(x[i] - x_center_source)/epsDelta) + 1;
        magnitude *= cos(Pi*(y[i] - y_center_source)/epsDelta) + 1;

#ifdef __3D__
        //// magnitude *= 1 + 1;
        /*
        magnitude *= cos(Pi*(z[i] - z_center_source)/epsDelta) + 1;
         */
#endif

        magnitude /= 4.*epsDelta*epsDelta;
#ifdef __3D__
        magnitude /= 480.;
#endif

        coeffs[i][dim+1] += singular_sign[m] * magnitude * Qin;
        Output::print<4>(" adding a singular source/sink - point ", m,
                " coeff[", dim+1, "] = ", coeffs[i][dim+1]);
      }
    }
  }
}

/** ************************************************************************ */
void approximate_delta_functions_3doubledoublets(int n_points, double *x, double *y,
#ifdef __3D__
        double *z,
#endif
        double **parameters, double **coeffs, double u_in,
        double center_x_moving_doublet_bottom_row, double center_x_moving_doublet_top_row, double well_radius,
        double delta_fct_eps_factor, double well_distance)
{
int dim = 2;
#ifdef __3D__
  dim = 3;
#endif

  for (int i = 0; i < n_points; ++i)
  {
    for (size_t k = 0; k <= dim+1; k++)
    {
      coeffs[i][k] = 0.;
    }

    double epsDelta = delta_fct_eps_factor * well_radius;  // ~h
    /*double H = 75.;
    double Volume = epsDelta*epsDelta*Pi*H;  //well_radius*well_radius*Pi*H;
    double Qin =  150./3600./Volume; // m^3/h thickness = 1000m
     */
    double Qin = u_in * (2/well_radius);//(dim/well_radius);

    std::vector<double> singular_x, singular_y,
#ifdef __3D__
    singular_z,
#endif
    singular_sign;

    //%%%%%%%%%%%%%%%%%% FIXED DOUBLE DOUBLET (center) %%%%%%%%%%%%%%%%%%%%%%%
    // source1 in (xi,yi) for fixed doublet
    singular_x.push_back(5000. - (well_distance)*3/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink1 in (xe,ye) for fixed doublet
    singular_x.push_back(5000.- (well_distance)/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // source2 in (xi,yi) for fixed doublet
    singular_x.push_back(5000. + (well_distance)/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink2 in (xe,ye) for fixed doublet
    singular_x.push_back(5000.+ (well_distance)*3/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    //%%%%%%%%%%%%%%%%%% MOVING DOUBLE DOUBLET 1 (lower) %%%%%%%%%%%%%%%%%%%%%%%
    // source1 in (xi,yi) for moving doublet
    singular_x.push_back(center_x_moving_doublet_bottom_row - (well_distance)*3/2.);
    singular_y.push_back(2000.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink1 in (xe,ye) for moving doublet
    singular_x.push_back(center_x_moving_doublet_bottom_row - (well_distance)/2.);
    singular_y.push_back(2000.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // source2 in (xi,yi) for moving doublet
    singular_x.push_back(center_x_moving_doublet_bottom_row + (well_distance)/2.);
    singular_y.push_back(2000.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink2 in (xe,ye) for moving doublet
    singular_x.push_back(center_x_moving_doublet_bottom_row + (well_distance)*3/2.);
    singular_y.push_back(2000.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    //%%%%%%%%%%%%%% MOVING DOUBLE DOUBLET 2 (upper) %%%%%%%%%%%%%%%%%%%%%%%%%%%
    // source1 in (xi,yi) for moving doublet
    singular_x.push_back(center_x_moving_doublet_top_row - (well_distance)*3/2.);
    singular_y.push_back(4000.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink1 in (xe,ye) for moving doublet
    singular_x.push_back(center_x_moving_doublet_top_row - (well_distance)/2.);
    singular_y.push_back(4000.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // source2 in (xi,yi) for moving doublet
    singular_x.push_back(center_x_moving_doublet_top_row + (well_distance)/2.);
    singular_y.push_back(4000.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink2 in (xe,ye) for moving doublet
    singular_x.push_back(center_x_moving_doublet_top_row + (well_distance)*3/2.);
    singular_y.push_back(4000.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    for (unsigned int m = 0; m < singular_x.size(); m++)
    {
      double x_center_source = singular_x[m];
      double y_center_source = singular_y[m];
#ifdef __3D__
      double z_center_source = singular_z[m];//0-500
#endif
      double x_distance_to_source = std::pow(std::abs(x[i] - x_center_source), 2);
      double y_distance_to_source = std::pow(std::abs(y[i] - y_center_source), 2);
      /*
#ifdef __3D__
      double z_distance_to_source = std::pow(std::abs(z[i] - z_center_source), 2);
#endif
       */

      bool at_source = (x_distance_to_source < epsDelta*epsDelta) *
              (y_distance_to_source < epsDelta*epsDelta)
              /*
#ifdef __3D__
               * (z_distance_to_source < epsDelta*epsDelta)
#endif
               */
              ;

      if (at_source)
      {
     double magnitude = cos(Pi*(x[i] - x_center_source)/epsDelta) + 1;
        magnitude *= cos(Pi*(y[i] - y_center_source)/epsDelta) + 1;

#ifdef __3D__
        //// magnitude *= 1 + 1;
        /*
        magnitude *= cos(Pi*(z[i] - z_center_source)/epsDelta) + 1;
         */
#endif

        magnitude /= 4.*epsDelta*epsDelta;
#ifdef __3D__
        magnitude /= 480.;
#endif

        coeffs[i][dim+1] += singular_sign[m] * magnitude * Qin;
        Output::print<4>(" adding a singular source/sink - point ", m,
                " coeff[", dim+1, "] = ", coeffs[i][dim+1]);
      }
    }
  }
}

/** ************************************************************************ */
void approximate_delta_functions_3doubledoublets_4controls(int n_points, double *x, double *y,
#ifdef __3D__
        double *z,
#endif
        double **parameters, double **coeffs, double u_in,
        double center_x_moving_doublet_bottom_row, double center_x_moving_doublet_top_row, double well_radius,
        double delta_fct_eps_factor, double well_distance,
        double height_of_bottom_row, double height_of_top_row)
{
int dim = 2;
#ifdef __3D__
  dim = 3;
#endif

  for (int i = 0; i < n_points; ++i)
  {
    for (size_t k = 0; k <= dim+1; k++)
    {
      coeffs[i][k] = 0.;
    }

    double epsDelta = delta_fct_eps_factor * well_radius;  // ~h
    /*double H = 75.;
    double Volume = epsDelta*epsDelta*Pi*H;  //well_radius*well_radius*Pi*H;
    double Qin =  150./3600./Volume; // m^3/h thickness = 1000m
     */
    double Qin = u_in * (2/well_radius);//(dim/well_radius);

    std::vector<double> singular_x, singular_y,
#ifdef __3D__
    singular_z,
#endif
    singular_sign;

    //%%%%%%%%%%%%%%%%%% FIXED DOUBLE DOUBLET (center) %%%%%%%%%%%%%%%%%%%%%%%
    // source1 in (xi,yi) for fixed doublet
    singular_x.push_back(5000. - (well_distance)*3/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink1 in (xe,ye) for fixed doublet
    singular_x.push_back(5000.- (well_distance)/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // source2 in (xi,yi) for fixed doublet
    singular_x.push_back(5000. + (well_distance)/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink2 in (xe,ye) for fixed doublet
    singular_x.push_back(5000.+ (well_distance)*3/2.);
    singular_y.push_back(3000.);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    //%%%%%%%%%%%%%%%%%% MOVING DOUBLE DOUBLET 1 (lower) %%%%%%%%%%%%%%%%%%%%%%%
    // source1 in (xi,yi) for moving doublet
    singular_x.push_back(center_x_moving_doublet_bottom_row - (well_distance)*3/2.);
    singular_y.push_back(height_of_bottom_row);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink1 in (xe,ye) for moving doublet
    singular_x.push_back(center_x_moving_doublet_bottom_row - (well_distance)/2.);
    singular_y.push_back(height_of_bottom_row);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // source2 in (xi,yi) for moving doublet
    singular_x.push_back(center_x_moving_doublet_bottom_row + (well_distance)/2.);
    singular_y.push_back(height_of_bottom_row);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink2 in (xe,ye) for moving doublet
    singular_x.push_back(center_x_moving_doublet_bottom_row + (well_distance)*3/2.);
    singular_y.push_back(height_of_bottom_row);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    //%%%%%%%%%%%%%% MOVING DOUBLE DOUBLET 2 (upper) %%%%%%%%%%%%%%%%%%%%%%%%%%%
    // source1 in (xi,yi) for moving doublet
    singular_x.push_back(center_x_moving_doublet_top_row - (well_distance)*3/2.);
    singular_y.push_back(height_of_top_row);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink1 in (xe,ye) for moving doublet
    singular_x.push_back(center_x_moving_doublet_top_row - (well_distance)/2.);
    singular_y.push_back(height_of_top_row);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // source2 in (xi,yi) for moving doublet
    singular_x.push_back(center_x_moving_doublet_top_row + (well_distance)/2.);
    singular_y.push_back(height_of_top_row);
    singular_sign.push_back(1.);
#ifdef __3D__
    //singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif
    // sink2 in (xe,ye) for moving doublet
    singular_x.push_back(center_x_moving_doublet_top_row + (well_distance)*3/2.);
    singular_y.push_back(height_of_top_row);
    singular_sign.push_back(-1.);
#ifdef __3D__
    // singular_z.push_back(3000.); //0-500
    singular_z.push_back(250.);
#endif

    for (unsigned int m = 0; m < singular_x.size(); m++)
    {
      double x_center_source = singular_x[m];
      double y_center_source = singular_y[m];
#ifdef __3D__
      double z_center_source = singular_z[m];//0-500
#endif
      double x_distance_to_source = std::pow(std::abs(x[i] - x_center_source), 2);
      double y_distance_to_source = std::pow(std::abs(y[i] - y_center_source), 2);
      /*
#ifdef __3D__
      double z_distance_to_source = std::pow(std::abs(z[i] - z_center_source), 2);
#endif
       */

      bool at_source = (x_distance_to_source < epsDelta*epsDelta) *
              (y_distance_to_source < epsDelta*epsDelta)
              /*
#ifdef __3D__
               * (z_distance_to_source < epsDelta*epsDelta)
#endif
               */
              ;

      if (at_source)
      {
     double magnitude = cos(Pi*(x[i] - x_center_source)/epsDelta) + 1;
        magnitude *= cos(Pi*(y[i] - y_center_source)/epsDelta) + 1;

#ifdef __3D__
        //// magnitude *= 1 + 1;
        /*
        magnitude *= cos(Pi*(z[i] - z_center_source)/epsDelta) + 1;
         */
#endif

        magnitude /= 4.*epsDelta*epsDelta;
#ifdef __3D__
        magnitude /= 480.;
#endif

        coeffs[i][dim+1] += singular_sign[m] * magnitude * Qin;
        Output::print<4>(" adding a singular source/sink - point ", m,
                " coeff[", dim+1, "] = ", coeffs[i][dim+1]);
      }
    }
  }
}

/** ************************************************************************ */
template <int d>
void GeothermalPlantsPositionOptimization<d>::apply_control_and_solve(const double* x)
{
  Chrono time;
  double t_start = GetTime();
  // apply control x

  // for 1 doublet (varying distance between source and sink)
  double distance = x[0];

  // for 2 doublets (1 is moving to the left and to the right)
  double center_x_moving_doublet = x[0];

  // for 3 double doublets (2 rows are moving to the left and to the right)
  double center_x_moving_doublet_bottom_row = x[0];
  double center_x_moving_doublet_top_row = x[1];

  // // for 3 double doublets (2 rows are moving to the left/right and up/down)
  //double center_x_moving_doublet_bottom_row = x[0];
  //double center_x_moving_doublet_top_row = x[1];
  double height_of_bottom_row = x[2];
  double height_of_top_row = x[3];

  using namespace std::placeholders;
  
/*
 #ifdef __2D__
  CoeffFct2D coeff = std::bind(approximate_delta_functions,
			       _1, _2, _3, _4, _5, distance, (double) this->db["well_radius"], (double) this->db["delta_fct_eps_factor"]);
  // try to generalize to more wells
  //   CoeffFct2D coeff = std::bind(approximate_delta_functions_3,
  //_1, _2, _3, _4, _5, x_i, x_p1, x_p2, y_i, y_p1, y_p2);
#else
  CoeffFct3D coeff = std::bind(approximate_delta_functions,
			       _1, _2, _3, _4, _5, _6, distance, (double) this->db["well_radius"], (double) this->db["delta_fct_eps_factor"]);
#endif
  */

  // Set coeffs according to scenario
  CoeffFct coeff;
  if (this->db["scenario"].is("3_rows_of_double_doublets_varying_row_distance") )
    {
      if(n_control != 4)
        ErrThrow("The chosen 'scenario' and the chosen 'n_control' do not fit. "
                "This scenario necessitates n_control=4.");

      Output::print<1>("CURRENT CONTROL center_x_moving_doublet_bottom_row: ",
              std::setprecision(14), center_x_moving_doublet_bottom_row);
      Output::print<1>("CURRENT CONTROL center_x_moving_doublet_top_row: ",
              std::setprecision(14), center_x_moving_doublet_top_row);
      Output::print<1>("CURRENT CONTROL height_of_bottom_row: ",
              std::setprecision(14), height_of_bottom_row);
      Output::print<1>("CURRENT CONTROL height_of_top_row: ",
              std::setprecision(14), height_of_top_row);

      coeff = std::bind(approximate_delta_functions_3doubledoublets_4controls, _1, _2, _3, _4, _5,
  #ifdef __3D__
              _6,
  #endif
              this->db["u_in"], center_x_moving_doublet_bottom_row,
              center_x_moving_doublet_top_row, (double) this->db["well_radius"],
              (double) this->db["delta_fct_eps_factor"],
              (double) this->db["well_distance"],
              height_of_bottom_row, height_of_top_row);
    }
  else if (this->db["scenario"].is("3_rows_of_double_doublets_fixed_well_distance") )
  {
    if(n_control != 2)
      ErrThrow("The chosen 'scenario' and the chosen 'n_control' do not fit. "
              "This scenario necessitates n_control=2.");

    Output::print<1>("CURRENT CONTROL center_x_moving_doublet_bottom_row: ",
            std::setprecision(14), center_x_moving_doublet_bottom_row);
    Output::print<1>("CURRENT CONTROL center_x_moving_doublet_top_row: ",
            std::setprecision(14), center_x_moving_doublet_top_row);
    coeff = std::bind(approximate_delta_functions_3doubledoublets, _1, _2, _3, _4, _5,
#ifdef __3D__
            _6,
#endif
            this->db["u_in"], center_x_moving_doublet_bottom_row,
            center_x_moving_doublet_top_row, (double) this->db["well_radius"],
            (double) this->db["delta_fct_eps_factor"],
            (double) this->db["well_distance"]);
  }
  else if (this->db["scenario"].is("2doublets_fixed_well_distance") )
  {
    if(n_control != 1)
      ErrThrow("The chosen 'scenario' and the chosen 'n_control' do not fit. "
              "This scenario necessitates n_control=1.");

    Output::print<1>("CURRENT CONTROL center_x_moving_doublet: ",
            std::setprecision(14), center_x_moving_doublet);
    coeff = std::bind(approximate_delta_functions_2doublets, _1, _2, _3, _4, _5,
#ifdef __3D__
            _6,
#endif
            this->db["u_in"], center_x_moving_doublet,
            (double) this->db["well_radius"],
            (double) this->db["delta_fct_eps_factor"],
            (double) this->db["well_distance"]);
  }
  else if (this->db["scenario"].is("1doublet_optimize_distance"))
  {
    if(n_control != 1)
      ErrThrow("The chosen 'scenario' and the chosen 'n_control' do not fit. "
              "This scenario necessitates n_control=1.");

    Output::print<1>("CURRENT CONTROL well distance: ",
            std::setprecision(14), distance);
    coeff = std::bind(approximate_delta_functions, _1, _2, _3, _4, _5,
#ifdef __3D__
            _6,
#endif
            this->db["u_in"], distance, (double) this->db["well_radius"],
            (double) this->db["delta_fct_eps_factor"]);
  }
  else if (this->db["scenario"].is("lattice"))
  {
    if(n_control != 1)
      ErrThrow("The chosen 'scenario' and the chosen 'n_control' do not fit. "
              "This scenario necessitates n_control=1.");

    Output::print<1>("CURRENT CONTROL well distance: ",
            std::setprecision(14), distance);
    coeff = std::bind(approximate_delta_functions_lattice,
		      _1, _2, _3, _4, _5,
#ifdef __3D__
            _6,
#endif
		      this->db["u_in"], distance, (double) this->db["well_radius"],
		      (double) this->db["delta_fct_eps_factor"],
		      (size_t) this->db["lattice_n_wells_per_row"]);
  }
  else if (this->db["scenario"].is("hexagon"))
  {
    if(n_control != 1)
      ErrThrow("The chosen 'scenario' and the chosen 'n_control' do not fit. "
              "This scenario necessitates n_control=1.");

    Output::print<1>("CURRENT CONTROL well distance: ",
            std::setprecision(14), distance);
    coeff = std::bind(approximate_delta_functions_hexagon,
		      _1, _2, _3, _4, _5,
#ifdef __3D__
            _6,
#endif
		      this->db["u_in"], distance, (double) this->db["well_radius"],
		      (double) this->db["delta_fct_eps_factor"]);
  }
      

  std::string disc_type = this->db["space_discretization_type"];
  bool nonsymm_gls = (disc_type == std::string("nonsymm_gls"));
  int rhs_div_sign = 1;
  if (nonsymm_gls)
  {
    rhs_div_sign = -1;
  }

#ifdef __2D__
  LocalAssembling<d> la(2, {D00, D00}, {0, 1}, {}, {},	{0, 0, 1},
			coeff, {std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8,
					rhs_div_sign)} ,nullptr,
			0, d+1, 0, {}, {}, 0, nullptr, 0, {}, {});
#else
  LocalAssembling<d> la(2, {D000, D000},{0, 1}, {}, {},	{0, 0, 0, 1},
			coeff, {std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8,
					rhs_div_sign)} ,nullptr,
			0, d+1, 0,  {},  {},   0,   nullptr,  0,  {}, {});
#endif

  auto v_space = brinkman_mixed.get_velocity_space();
  auto p_space = brinkman_mixed.get_pressure_space();
  auto& rhs = brinkman_mixed.get_rhs();
  rhs.reset(); //????

#ifdef __2D__
  std::array<const FESpace*, d+1> fespaces = {{v_space.get(), v_space.get(), p_space.get()}};
  double *rhs_pointers[d+1] = {rhs.block(0), rhs.block(1), rhs.block(2)};
  BoundaryConditionFunction* boundary_conditions[d+1] = {
          v_space->get_boundary_condition(), v_space->get_boundary_condition(),
          p_space->get_boundary_condition() };
  auto& example = brinkman_mixed.get_example();
  std::array<BoundaryValuesFunction*, d+1> non_const_bound_values;
  non_const_bound_values[0] = example.get_bd()[0];
  non_const_bound_values[1] = example.get_bd()[1];
  non_const_bound_values[2] = example.get_bd()[2];
#else
  std::array<const FESpace*, d+1> fespaces = {{v_space.get(), v_space.get(), v_space.get(), p_space.get()}};
  double *rhs_pointers[d+1] = {rhs.block(0), rhs.block(1), rhs.block(2), rhs.block(3)};
  BoundaryConditionFunction* boundary_conditions[d+1] = {
      v_space->get_boundary_condition(), v_space->get_boundary_condition(),
      v_space->get_boundary_condition(), p_space->get_boundary_condition() };
  auto& example = brinkman_mixed.get_example();
  std::array<BoundaryValuesFunction*, d+1> non_const_bound_values;
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
	     2, fespaces.data()+1, 0, nullptr, 0, nullptr, d+1, rhs_pointers,
	     fespaces.data(), boundary_conditions,
	     non_const_bound_values.data(), la);

  double mat_norm 
  = brinkman_mixed.get_matrix().get_combined_matrix()->GetNorm();
  double rhs_norm = rhs.norm();
  Output::print("---- norms: ", mat_norm, "  ", rhs_norm);

  // solve the flow system
  Output::print<2>(" ******   flow   *******");

  brinkman_mixed.solve();
  // brinkman_mixed.get_pressure().WriteSol("/Users/blank/ParMooN/Tests/Geothermal_Plants_Position_Optimization/", "written_fe_function_pressure");
  brinkman_mixed.output(n_calls);
  time.restart_and_print("solving the stationary flow problem");

  Output::print<2>(" ******   temperature   *******");
  
  TDatabase::TimeDB->TIME_DISC = 2; // Crank-Nicolson
  tcd_primal.reset_for_output();
  
  TimeDiscretization& tss = tcd_primal.get_time_stepping_scheme();
  tss.current_step_ = 0;
  
  this->temperature_production_well_at_time_steps.clear();
  
  tss.set_time_disc_parameters();
  Output::print<2>("  --> assemble ");
  tcd_primal.assemble_initial_time();

  double start_time = db["time_start"];
  double end_time   = db["time_end"];
  Output::print<2>("  --> solve temperature, t0 = ", start_time, ", tEnd = ", end_time);
  TDatabase::TimeDB->CURRENTTIME = start_time;

  // output initial condition
  tcd_primal.output();

  // ======================================================================
  // time iteration
  // ======================================================================
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    tss.current_step_++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    cout << "db[scaling_time_derivative]: " <<  db["scaling_time_derivative"] <<endl;
    tss.set_time_disc_parameters();

    double tau = db["time_step_length"];

    TDatabase::TimeDB->CURRENTTIME += tau;

    Output::print<1>("\nCURRENT time: ", TDatabase::TimeDB->CURRENTTIME);

    /// OLD: SetTimeDiscParameters(1);

    tcd_primal.assemble(brinkman_mixed.get_velocity(), x,  db["diffusion_coefficient"]);

    tcd_primal.solve();
    if((tss.current_step_-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      tcd_primal.output();

    // Temperature at production well
    auto temperature = tcd_primal.get_function();

#ifdef __2D__
      const TCollection* Coll = temperature.GetFESpace2D()->GetCollection();
#else
      const TCollection* Coll = temperature.GetFESpace3D()->GetCollection();
#endif

      if (this->db["scenario"].is("lattice"))
      {
          double domain_Lx = 10000.;
          double domain_Ly = 6000.;

        // has to be an even number
        size_t n_wells_per_row = (size_t) db["lattice_n_wells_per_row"]; //equals number of wells per column

        // lower left well
        double x0 = domain_Lx/2.- (n_wells_per_row-1)/2.*distance;
        double y0 = domain_Ly/2.- (n_wells_per_row-1)/2.*distance;

        size_t n_wells = n_wells_per_row*n_wells_per_row;

        if (n_wells % 2 != 0)
          ErrThrow("n_wells_per_row should be an even number with equally many production "
                  "and injection wells, but we have ", n_wells, " many wells.");

        std::vector<double> x_production_well, y_production_well;

        for (unsigned int k1 = 0; k1 < n_wells_per_row; k1++) {
            for (unsigned int k2 = 0; k2 < n_wells_per_row; k2++) {
              double xk = x0 + distance*k2;
              double yk = y0 + distance*k1;
              if ( pow(-1,k1)*pow(-1,k2) == -1 )
              {
              x_production_well.push_back(xk);
              y_production_well.push_back(yk);
              }
            }
          }

        x_production_well.resize(n_wells/2.);
        y_production_well.resize(n_wells/2.);

  #ifdef __3D__
        std::vector<double> z_production_well(n_wells/2.,0.);//3000.;
        for (int i = 0; i < n_wells/2.; i++)
          z_production_well[i] = 250.;
  #endif

        std::vector<double> average(n_wells/2., 0.), min(n_wells/2., 0.);
        bool minimum_temperature_reached = false;
        size_t Num_circle_points = 10;
        double total_average = 0.;

        for (int i = 0; i < x_production_well.size(); i++)
        {
          std::array<double, d> centers = {x_production_well[i], y_production_well[i]
  #ifdef __3D__
             , z_production_well[i]
  #endif
          };

          sinks sink((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
          sink.find_average_and_min_along_circle(&temperature, average[i], min[i]);

          total_average += average[i];

          if (min[i] < (double) db["minimum_temperature_production_well"])
          {
            minimum_temperature_reached = true;
            Output::print("Minimum production temperature obtained in the production well at position ",
                          x_production_well[i], ", ", y_production_well[i],
#ifdef __3D__
                          ", ", z_production_well[i],
#endif
                  " at time step ",
                  (double) tss.current_step_, " .");
            break;
          }
        }

        if ( !minimum_temperature_reached )
             {
               total_average = total_average/(x_production_well.size());
               Output::print(" *** T(average) = ", total_average);
               this->temperature_production_well_at_time_steps.push_back(total_average);
             }
             else
               break;

      }
      else if (this->db["scenario"].is("hexagon"))
      {
        // todo
        Output::print(" Temperature along well not yet implemented");
      }
      else if (this->db["scenario"].is("3_rows_of_double_doublets_varying_row_distance"))
    {
      std::vector<double> x_production_well(6,0.), y_production_well(6,0.);
      x_production_well[0] = center_x_moving_doublet_top_row - ((double)  db["well_distance"])/2.;
      x_production_well[1] = center_x_moving_doublet_top_row + ((double)  db["well_distance"])*3/2.;
      x_production_well[2] = 5000. - ((double)  db["well_distance"])/2.;
      x_production_well[3] = 5000. + ((double)  db["well_distance"])*3/2.;
      x_production_well[4] = center_x_moving_doublet_bottom_row - ((double)  db["well_distance"])/2.;
      x_production_well[5] = center_x_moving_doublet_bottom_row + ((double)  db["well_distance"])*3/2.;
      y_production_well[0] = height_of_top_row;
      y_production_well[1] = height_of_top_row;
      y_production_well[2] = 3000.;
      y_production_well[3] = 3000.;
      y_production_well[4] = height_of_bottom_row;
      y_production_well[5] = height_of_bottom_row;
#ifdef __3D__
      std::vector<double> z_production_well(6,0.);//3000.;
      for (int i = 0; i < 6; i++)
        z_production_well[i] = 250.;
#endif

      std::vector<double> average(6, 0.), min(6, 0.);
      bool minimum_temperature_reached = false;
      size_t Num_circle_points = 10;

      for (int i = 0; i < x_production_well.size(); i++)
      {
        std::array<double, d> centers = {x_production_well[i], y_production_well[i]
#ifdef __3D__
           , z_production_well[i]
#endif
        };

        sinks sink((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
        sink.find_average_and_min_along_circle(&temperature, average[i], min[i]);

        if (min[i] < (double) db["minimum_temperature_production_well"])
        {
          minimum_temperature_reached = true;
          break;
        }
      }

      Output::print(" *** T(moving_well_upper_row, average) = ", (average[0] + average[1])/2.);
      Output::print(" *** T(fixed_row, average) = ", (average[2] + average[3])/2.);
      Output::print(" *** T(moving_well_lower_row, average) = ", (average[4] + average[5])/2.);

      if ( !minimum_temperature_reached )
      {
        this->temperature_production_well_at_time_steps.push_back(
                average[0] + average[1] + average[2] + average[3] + average[4] + average[5]);
      }
      else
      {
        Output::print("Minimum production temperature obtained in one of the production wells at time step ", (double) tss.current_step_, " .");
        break;
      }

    }
    else if (this->db["scenario"].is("3_rows_of_double_doublets_fixed_well_distance"))
    {
      std::vector<double> x_production_well(6,0.), y_production_well(6,0.);
      x_production_well[0] = center_x_moving_doublet_top_row - ((double)  db["well_distance"])/2.;
      x_production_well[1] = center_x_moving_doublet_top_row + ((double)  db["well_distance"])*3/2.;
      x_production_well[2] = 5000. - ((double)  db["well_distance"])/2.;
      x_production_well[3] = 5000. + ((double)  db["well_distance"])*3/2.;
      x_production_well[4] = center_x_moving_doublet_bottom_row - ((double)  db["well_distance"])/2.;
      x_production_well[5] = center_x_moving_doublet_bottom_row + ((double)  db["well_distance"])*3/2.;
      y_production_well[0] = 4000.;
      y_production_well[1] = 4000.;
      y_production_well[2] = 3000.;
      y_production_well[3] = 3000.;
      y_production_well[4] = 2000.;
      y_production_well[5] = 2000.;
#ifdef __3D__
      std::vector<double> z_production_well(6,0.);//3000.;
      for (int i = 0; i < 6; i++)
        z_production_well[i] = 250.;
#endif

      std::vector<double> average(6, 0.), min(6, 0.);
      bool minimum_temperature_reached = false;
      size_t Num_circle_points = 10;

      for (int i = 0; i < x_production_well.size(); i++)
      {
        std::array<double, d> centers = {x_production_well[i], y_production_well[i]
#ifdef __3D__
              , z_production_well[i]
#endif
        };

        sinks sink((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
        sink.find_average_and_min_along_circle(&temperature, average[i], min[i]);

        if (min[i] < (double) db["minimum_temperature_production_well"])
        {
          minimum_temperature_reached = true;
          break;
        }
      }
      Output::print(" *** T(moving_well_upper_row, average) = ", (average[0] + average[1])/2.);
      Output::print(" *** T(fixed_row, average) = ", (average[2] + average[3])/2.);
      Output::print(" *** T(moving_well_lower_row, average) = ", (average[4] + average[5])/2.);

      if ( !minimum_temperature_reached )
      {
        this->temperature_production_well_at_time_steps.push_back(
                average[0] + average[1] + average[2] + average[3] + average[4] + average[5]);
      }
      else
      {
        Output::print("Minimum production temperature obtained in one of the production wells at time step ", (double) tss.current_step_, " .");
        break;
      }
    }
    else if (this->db["scenario"].is("2doublets_fixed_well_distance") )
    {
      std::vector<double> x_production_well(2,0.), y_production_well(2,0.);
      x_production_well[0] = 5000. + ((double)  db["well_distance"])/2.;
      x_production_well[1] = center_x_moving_doublet + ((double)  db["well_distance"])/2.;
      y_production_well[0] = 2500.;
      y_production_well[1] = 3500.;
#ifdef __3D__
      std::vector<double> z_production_well(2,0.);//3000.;
      for (int i = 0; i < 2; i++)
        z_production_well[i] = 250.;
#endif

#ifdef __2D__
      double temperature_values_fixed_production_well[d+1], temperature_values_moving_production_well[d+1];
#else
      std::vector<double> temperature_values_fixed_production_well(d+1), temperature_values_moving_production_well(d+1);
#endif

      std::vector<double> average(2, 0.), min(2, 0.);
      bool minimum_temperature_reached = false;
      size_t Num_circle_points = 10;

      for (int i = 0; i < x_production_well.size(); i++)
      {
        std::array<double, d> centers = {x_production_well[i], y_production_well[i]
#ifdef __3D__
                                                                                 , z_production_well[i]
#endif
        };

        sinks sink((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
        sink.find_average_and_min_along_circle(&temperature, average[i], min[i]);

        if (min[i] < (double) db["minimum_temperature_production_well"])
        {
          minimum_temperature_reached = true;
          break;
        }
      }
      Output::print(" *** T(fixed_well) = ", average[0]);
      Output::print(" *** T(moving_wells) = ", average[1]);

      if ( !minimum_temperature_reached )
      {
        this->temperature_production_well_at_time_steps.push_back(
                average[0] + average[1]);
      }
      else
      {
        Output::print("Minimum production temperature obtained in one of the production wells at time step ", (double) tss.current_step_, " .");
        break;
      }
    }
    else if (this->db["scenario"].is("1doublet_optimize_distance") )
    {
      std::vector<double> x_production_well(1,0.), y_production_well(1,0.);
      x_production_well[0] = 5000. + (distance)/2.;
      y_production_well[0] = 3000.;
#ifdef __3D__
      std::vector<double> z_production_well(1,0.);
      z_production_well[0] = 250.;//3000.;
#endif

#ifdef __2D__
      double temperature_values[d+1];
#else
      std::vector<double> temperature_values(d+1);
#endif

      std::vector<double> average(1, 0.), min(1, 0.);
      bool minimum_temperature_reached = false;

      for (int i = 0; i < x_production_well.size(); i++)
      {
        std::array<double, d> centers = {x_production_well[i], y_production_well[i]
#ifdef __3D__
                                                                                 , z_production_well[i]
#endif
        };

        size_t Num_circle_points = 10;
        sinks sink((double) db["delta_fct_eps_factor"], (double) db["well_radius"],
                centers, Num_circle_points, Coll);
        sink.find_average_and_min_along_circle(&temperature, average[i], min[i]);

        if (min[i] < (double) db["minimum_temperature_production_well"])
        {
          minimum_temperature_reached = true;
          break;
        }
      }
      Output::print(" *** T(moving_well) = ", average[0]);

      if ( !minimum_temperature_reached )
      {
        this->temperature_production_well_at_time_steps.push_back(
                average[0] + average[1]);
      }
      else
      {
        Output::print("Minimum production temperature obtained in one of the "
                "production wells at time step ", (double) tss.current_step_, " .");
        break;
      }
    }
  }
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================


}

/** ************************************************************************ */
template <int d>
double GeothermalPlantsPositionOptimization<d>::compute_functional(const double* x)
const
{
  Output::print("computing current value of functional");
  double distance = x[0];
  double center_x_moving_doublet = x[0]; // for 2 doublets (1 is moving)

  double center_x_moving_doublet_bottom_row = x[0]; // for 3 double doublets (2 are moving left/right)
  double center_x_moving_doublet_top_row = x[1];

  double height_of_bottom_row = 2000.; // for 3 double doublets (2 are moving left/right and top/bottom)
  double height_of_top_row = 4000.;

  auto u = brinkman_mixed.get_velocity();
  auto u1 = u.GetComponent(0);
  auto u2 = u.GetComponent(1);
#ifdef __3D__
  auto u3 = u.GetComponent(2);
#endif
  auto pressure = brinkman_mixed.get_pressure();

  double functional_value;

#ifdef __2D__
      const TCollection* Coll = pressure.GetFESpace2D()->GetCollection();
#else
      const TCollection* Coll = pressure.GetFESpace3D()->GetCollection();
#endif

      if (this->db["scenario"].is("lattice"))
      {
        double domain_Lx = 10000.;
        double domain_Ly = 6000.;

      // has to be an even number
      size_t n_wells_per_row = (size_t) db["lattice_n_wells_per_row"]; //equals number of wells per column

      // lower left well
      double x0 = domain_Lx/2.- (n_wells_per_row-1)/2.*distance;
      double y0 = domain_Ly/2.- (n_wells_per_row-1)/2.*distance;

      size_t n_wells = n_wells_per_row*n_wells_per_row;

      if (n_wells % 2 != 0)
        ErrThrow("n_wells_per_row should be an even number with equally many production "
                "and injection wells, but we have ", n_wells, " many wells.");

      std::vector<double> x_production_wells, y_production_wells, x_injection_wells, y_injection_wells;

      for ( unsigned int k1 = 0; k1 < n_wells_per_row; k1++)
      {
          for (unsigned int k2 = 0; k2 < n_wells_per_row; k2++)
          {
            double xk = x0 + distance*k2;
            double yk = y0 + distance*k1;
            if ( pow(-1,k1)*pow(-1,k2) == -1 )
            {
            x_production_wells.push_back(xk);
            y_production_wells.push_back(yk);
            }
            else
            {
              x_injection_wells.push_back(xk);
              y_injection_wells.push_back(yk);
            }
          }
        }

      x_production_wells.resize(n_wells/2.);
      y_production_wells.resize(n_wells/2.);
      x_injection_wells.resize(n_wells/2.);
      y_injection_wells.resize(n_wells/2.);

#ifdef __3D__
      std::vector<double> z_production_wells(n_wells/2.,0.), z_injection_wells(n_wells/2.,0.);//3000.;
      for (int i = 0; i < n_wells/2.; i++)
      {
        z_production_wells[i] = 250.;
        z_injection_wells[i] = 250.;
      }
#endif


        std::vector<double> average_injection(n_wells/2., 0.), min_injection(n_wells/2., 0.),
                average_production(n_wells/2., 0.), min_production(n_wells/2., 0.);
        double pressure_values_production_wells = 0.;
        double pressure_values_injection_wells = 0.;
            size_t Num_circle_points = 10;

            for (int i = 0; i < x_production_wells.size(); i++)
            {
              std::array<double, d> centers = {x_production_wells[i], y_production_wells[i]
      #ifdef __3D__
                , z_production_wells[i]
      #endif
              };

              sinks sink((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
              sink.find_average_and_min_along_circle(&pressure, average_production[i], min_production[i]);

              pressure_values_production_wells += average_production[i];
            }
            for (int i = 0; i < x_injection_wells.size(); i++)
            {
              std::array<double, d> centers = {x_injection_wells[i], y_injection_wells[i]
      #ifdef __3D__
                 , z_injection_wells[i]
      #endif
              };

              sinks source((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
              source.find_average_and_min_along_circle(&pressure, average_injection[i], min_injection[i]);

              pressure_values_injection_wells += average_injection[i];
            }


        //todo: compute Q from u_in
        /// double Q = 150/360;//24 * 50; // 50 - 300

        int number_of_time_steps_for_production = 0;
        double Delta_Temp = 0;

        for (int i = 0; i < this->temperature_production_well_at_time_steps.size(); i++)
        {
          if (this->temperature_production_well_at_time_steps.at(i) >= ((double) db["minimum_temperature_production_well"])/6.)
          {
            Delta_Temp += this->temperature_production_well_at_time_steps.at(i) -  (n_wells/2.)*(double) db["temperature_injection_well"];
            cout <<" temperature_production_well_at_time_steps: "<<  this->temperature_production_well_at_time_steps.at(i) << ", step: "<< i <<endl;
          }
          else
            break;

          number_of_time_steps_for_production = i+1;
        }

        double alpha = db["alpha_cost"];
        /*
         * double functional_value_new // = Q/(0.6)* Delta t * (pressure_prod[1] - pressure_inj[1])  -  Q * Delta t * fluid_density * fluid_heat_capacity * (temperature_prod[1] - temperature_inj); // Net energy AFTER 50 years
                                       //= Q * Delta t * (   1/(0.6) * (pressure_prod[1] - pressure_inj[1])  - fluid_density * fluid_heat_capacity * (temperature_prod[1] - temperature_inj)   );
       //since Q_i=const, Delta t = const we can minimize
         */

        functional_value =  (number_of_time_steps_for_production * 1/(double)db["pump_efficiency"] *
                   ( pressure_values_injection_wells - pressure_values_production_wells)
                  - (double) db["fluid_density"] * (double) db["fluid_heat_capacity"] * Delta_Temp)
          +  alpha * std::abs(distance - 1000.);

        Output::print("functional_value: ", functional_value);

        //write to stream
        std::ofstream outputFile("temperature_values_at_production_well_and_net_energy_for_disctance" + std::to_string(distance) +  ".txt");
        std::copy(this->temperature_production_well_at_time_steps.begin(), this->temperature_production_well_at_time_steps.end(), std::ostream_iterator<int>(outputFile, "\n"));
        //std::copy(functional_value, std::ostream_iterator<int>(outputFile, "\n"));
        outputFile << functional_value;


        std::vector<double> cost_functional = db["cost_functional"];

        auto temperature = tcd_primal.get_function();
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
      }
 else if (this->db["scenario"].is("3_rows_of_double_doublets_fixed_well_distance")
        || this->db["scenario"].is("3_rows_of_double_doublets_varying_row_distance"))
{
  if ( this->db["scenario"].is("3_rows_of_double_doublets_varying_row_distance"))
  {
    height_of_bottom_row = x[2];
    height_of_top_row = x[3];
  }
  std::vector<double> x_production_wells(6, 0.), x_injection_wells(6, 0.), y_production_wells(6, 0.), y_injection_wells(6, 0.)
#ifdef __3D__
 , z_production_wells(6, 0.), z_injection_wells(6, 0.)
#endif
  ;
//%%%%%%%%%%%%% UPPER ROW (moving) %%%%%%%%%%%%%%%%%%
  x_production_wells[0] = center_x_moving_doublet_bottom_row - (double) this->db["well_distance"]/2.; //5.5;
  y_production_wells[0] = height_of_top_row; //3.;
  x_production_wells[1] = center_x_moving_doublet_bottom_row + (double) this->db["well_distance"]*3/2.; //5.5;
  y_production_wells[1] = height_of_top_row; //3.;
  x_injection_wells[0] = center_x_moving_doublet_bottom_row - (double) this->db["well_distance"]*3/2.; //4.5;
  y_injection_wells[0] = height_of_top_row; //3.;
  x_injection_wells[1] = center_x_moving_doublet_bottom_row + (double) this->db["well_distance"]/2.; //4.5;
  y_injection_wells[1] = height_of_top_row; //3.;
  //%%%%%%%%%%%%% CENTER ROW (fixed) %%%%%%%%%%%%%%%%%%
  x_production_wells[2] = 5000. - (double) this->db["well_distance"]/2.; //5.5;
  y_production_wells[2] = 3000.; //3.;
  x_production_wells[3] =  5000. + (double) this->db["well_distance"]*3/2.; //5.5;
  y_production_wells[3] = 3000.; //3.;
  x_injection_wells[2] = 5000. - (double) this->db["well_distance"]*3/2.; //4.5;
  y_injection_wells[2] = 3000.; //3.;
  x_injection_wells[3] =  5000. + (double) this->db["well_distance"]/2.; //4.5;
  y_injection_wells[3] = 3000.; //3.;
  //%%%%%%%%%%%%% LOWER ROW (moving) %%%%%%%%%%%%%%%%%%
  x_production_wells[4] = center_x_moving_doublet_top_row - (double) this->db["well_distance"]/2.; //5.5;
  y_production_wells[4] = height_of_bottom_row; //3.;
  x_production_wells[5] = center_x_moving_doublet_top_row + (double) this->db["well_distance"]*3/2.; //5.5;
  y_production_wells[5] = height_of_bottom_row; //3.;
  x_injection_wells[4] = center_x_moving_doublet_top_row - (double) this->db["well_distance"]*3/2.; //4.5;
  y_injection_wells[4] = height_of_bottom_row; //3.;
  x_injection_wells[5] = center_x_moving_doublet_top_row + (double) this->db["well_distance"]/2.; //4.5;
  y_injection_wells[5] = height_of_bottom_row; //3.;
#ifdef __3D__
  for (int i = 0; i < 6; i++)
  {
    z_production_wells[i] = 250.;
    z_injection_wells[i] = 250.;
  }
#endif


  std::vector<double> average_injection(6, 0.), min_injection(6, 0.),
          average_production(6, 0.), min_production(6, 0.);
  double pressure_values_production_wells = 0.;
  double pressure_values_injection_wells = 0.;
      size_t Num_circle_points = 10;

      for (int i = 0; i < x_production_wells.size(); i++)
      {
        std::array<double, d> centers = {x_production_wells[i], y_production_wells[i]
#ifdef __3D__
          , z_production_wells[i]
#endif
        };

        sinks sink((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
        sink.find_average_and_min_along_circle(&pressure, average_production[i], min_production[i]);

        pressure_values_production_wells += average_production[i];
      }
      for (int i = 0; i < x_injection_wells.size(); i++)
      {
        std::array<double, d> centers = {x_injection_wells[i], y_injection_wells[i]
#ifdef __3D__
           , z_injection_wells[i]
#endif
        };

        sinks source((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
        source.find_average_and_min_along_circle(&pressure, average_injection[i], min_injection[i]);

        pressure_values_injection_wells += average_injection[i];
      }


  //todo: compute Q from u_in
  /// double Q = 150/360;//24 * 50; // 50 - 300

  int number_of_time_steps_for_production = 0;
  double Delta_Temp = 0;

  for (int i = 0; i < this->temperature_production_well_at_time_steps.size(); i++)
  {
    if (this->temperature_production_well_at_time_steps.at(i) >= ((double) db["minimum_temperature_production_well"])/6.)
    {
      Delta_Temp += this->temperature_production_well_at_time_steps.at(i) -  6.*(double) db["temperature_injection_well"];
      cout <<" temperature_production_well_at_time_steps: "<<  this->temperature_production_well_at_time_steps.at(i) << ", step: "<< i <<endl;
    }
    else
      break;

    number_of_time_steps_for_production = i+1;
  }

  double alpha = db["alpha_cost"];
  /*
   * double functional_value_new // = Q/(0.6)* Delta t * (pressure_prod[1] - pressure_inj[1])  -  Q * Delta t * fluid_density * fluid_heat_capacity * (temperature_prod[1] - temperature_inj); // Net energy AFTER 50 years
                                 //= Q * Delta t * (   1/(0.6) * (pressure_prod[1] - pressure_inj[1])  - fluid_density * fluid_heat_capacity * (temperature_prod[1] - temperature_inj)   );
 //since Q_i=const, Delta t = const we can minimize
   */

  functional_value =  (number_of_time_steps_for_production * 1/(double)db["pump_efficiency"] *
             ( pressure_values_injection_wells - pressure_values_production_wells)
            - (double) db["fluid_density"] * (double) db["fluid_heat_capacity"] * Delta_Temp)
    +  alpha * std::abs(center_x_moving_doublet_bottom_row - 2000.);

  Output::print("functional_value: ", functional_value);

  //write to stream
  std::ofstream outputFile("temperature_values_at_production_well_and_net_energy_for_disctance" + std::to_string(distance) +  ".txt");
  std::copy(this->temperature_production_well_at_time_steps.begin(), this->temperature_production_well_at_time_steps.end(), std::ostream_iterator<int>(outputFile, "\n"));
  //std::copy(functional_value, std::ostream_iterator<int>(outputFile, "\n"));
  outputFile << functional_value;


  std::vector<double> cost_functional = db["cost_functional"];

  auto temperature = tcd_primal.get_function();
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
}
else if (this->db["scenario"].is("2doublets_fixed_well_distance") )
  {
    std::vector<double> x_production_wells(2, 0.), x_injection_wells(2, 0.), y_production_wells(2, 0.), y_injection_wells(2, 0.)
#ifdef __3D__
   , z_production_wells(2, 0.), z_injection_wells(2, 0.)
#endif
    ;

    x_production_wells[0] = 5000. + (double) this->db["well_distance"]/2.; //5.5;
    y_production_wells[0] = 2500.; //3.;
    x_production_wells[1] = center_x_moving_doublet + (double) this->db["well_distance"]/2.; //5.5;
    y_production_wells[1] = 3500.; //3.;
    x_injection_wells[0] = 5000. - (double) this->db["well_distance"]/2.; //4.5;
    y_injection_wells[0] = 2500.; //3.;
    x_injection_wells[1] = center_x_moving_doublet - (double) this->db["well_distance"]/2.; //4.5;
    y_injection_wells[1] = 3500.; //3.;
#ifdef __3D__
    z_production_wells[0] = 250.; //3000.;
    z_injection_wells[0] = 250.; //3000.;
    z_production_wells[1] = 250.; //3000.;
    z_injection_wells[1] = 250.; //3000.;
#endif

    std::vector<double> average_injection(2, 0.), min_injection(2, 0.),
            average_production(2, 0.), min_production(2, 0.);
    double pressure_values_production_wells = 0.;
    double pressure_values_injection_wells = 0.;
         size_t Num_circle_points = 10;

         for (int i = 0; i < x_production_wells.size(); i++)
         {
           std::array<double, d> centers = {x_production_wells[i], y_production_wells[i]
   #ifdef __3D__
              , z_production_wells[i]
   #endif
           };

           sinks sink((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
           sink.find_average_and_min_along_circle(&pressure, average_production[i], min_production[i]);

           pressure_values_production_wells += average_production[i];
         }
         for (int i = 0; i < x_injection_wells.size(); i++)
         {
           std::array<double, d> centers = {x_injection_wells[i], y_injection_wells[i]
   #ifdef __3D__
              , z_injection_wells[i]
   #endif
           };

           sinks source((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
           source.find_average_and_min_along_circle(&pressure, average_injection[i], min_injection[i]);

           pressure_values_injection_wells += average_injection[i];
         }

    //todo: compute Q from u_in
    /// double Q = 150/360;//24 * 50; // 50 - 300

    int number_of_time_steps_for_production = 0;
    double Delta_Temp = 0;

    for (int i = 0; i < this->temperature_production_well_at_time_steps.size(); i++)
    {
      if (this->temperature_production_well_at_time_steps.at(i) >= ((double) db["minimum_temperature_production_well"])/2.)
      {
        Delta_Temp += this->temperature_production_well_at_time_steps.at(i) -  2.*(double) db["temperature_injection_well"];
        cout <<" temperature_production_well_at_time_steps: "<<  this->temperature_production_well_at_time_steps.at(i) << ", step: "<< i <<endl;
      }
      else
        break;

      number_of_time_steps_for_production = i+1;
    }

    double alpha = db["alpha_cost"];
    /*
     * double functional_value_new // = Q/(0.6)* Delta t * (pressure_prod[1] - pressure_inj[1])  -  Q * Delta t * fluid_density * fluid_heat_capacity * (temperature_prod[1] - temperature_inj); // Net energy AFTER 50 years
                                   //= Q * Delta t * (   1/(0.6) * (pressure_prod[1] - pressure_inj[1])  - fluid_density * fluid_heat_capacity * (temperature_prod[1] - temperature_inj)   );
   //since Q_i=const, Delta t = const we can minimize
     */
    functional_value =  (number_of_time_steps_for_production * 1/(double)db["pump_efficiency"] *
            (pressure_values_injection_wells - pressure_values_production_wells)
            - (double) db["fluid_density"] * (double) db["fluid_heat_capacity"] * Delta_Temp)
      +  alpha * std::abs(center_x_moving_doublet - 2000.);

    Output::print("functional_value: ", functional_value);

    //write to stream
    std::ofstream outputFile("temperature_values_at_production_well_and_net_energy_for_disctance" + std::to_string(distance) +  ".txt");
    std::copy(this->temperature_production_well_at_time_steps.begin(), this->temperature_production_well_at_time_steps.end(), std::ostream_iterator<int>(outputFile, "\n"));
    //std::copy(functional_value, std::ostream_iterator<int>(outputFile, "\n"));
    outputFile << functional_value;

    std::vector<double> cost_functional = db["cost_functional"];

    auto temperature = tcd_primal.get_function();
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


    //  double functional_value = 1.;
  }
else if (this->db["scenario"].is("1doublet_optimize_distance") )
{
  std::vector<double> x_production_wells, y_production_wells, x_injection_wells, y_injection_wells;
  x_production_wells[0] = 5000. + (distance)/2.; //5.5;
  y_production_wells[0] = 3000.; //3.;
  x_injection_wells[0] = 5000. - (distance)/2.; //4.5;
  y_injection_wells[0] = 3000.; //3.;
#ifdef __3D__
  std::vector<double> z_production_wells, z_injection_wells;
  z_production_wells[0] = 250.; //3000.;
  z_injection_wells[0] = 250.; //3000.;
#endif

       std::vector<double> average_injection(1, 0.), min_injection(1, 0.),
               average_production(1, 0.), min_production(1, 0.);
       double pressure_values_production_wells = 0.;
       double pressure_values_injection_wells = 0.;
       size_t Num_circle_points = 10;

       for (int i = 0; i < x_production_wells.size(); i++)
       {
         std::array<double, d> centers = {x_production_wells[i], y_production_wells[i]
 #ifdef __3D__
            , z_production_wells[i]
 #endif
         };

         sinks sink((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
         sink.find_average_and_min_along_circle(&pressure, average_production[i], min_production[i]);

         pressure_values_production_wells = average_production[i];
       }
       for (int i = 0; i < x_injection_wells.size(); i++)
       {
         std::array<double, d> centers = {x_injection_wells[i], y_injection_wells[i]
 #ifdef __3D__
            , z_injection_wells[i]
 #endif
         };

         sinks source((double) db["delta_fct_eps_factor"], (double) db["well_radius"], centers, Num_circle_points, Coll);
         source.find_average_and_min_along_circle(&pressure, average_injection[i], min_injection[i]);

         pressure_values_injection_wells = average_injection[i];
       }

  //todo: compute Q from u_in
  /// double Q = 150/360;//24 * 50; // 50 - 300

  int number_of_time_steps_for_production = 0;
  double Delta_Temp = 0;

  for (int i = 0; i < this->temperature_production_well_at_time_steps.size(); i++)
  {
    if (this->temperature_production_well_at_time_steps.at(i) >= (double) db["minimum_temperature_production_well"])
    { 
      Delta_Temp += this->temperature_production_well_at_time_steps.at(i) -  (double) db["temperature_injection_well"];
      cout <<" temperature_production_well_at_time_steps: "<<  this->temperature_production_well_at_time_steps.at(i) << ", step: "<< i <<endl;
    }
    else 
      break;

    number_of_time_steps_for_production = i+1;
  }

  double alpha = db["alpha_cost"];
  /* 
   * double functional_value_new // = Q/(0.6)* Delta t * (pressure_prod[1] - pressure_inj[1])  -  Q * Delta t * fluid_density * fluid_heat_capacity * (temperature_prod[1] - temperature_inj); // Net energy AFTER 50 years
                                 //= Q * Delta t * (   1/(0.6) * (pressure_prod[1] - pressure_inj[1])  - fluid_density * fluid_heat_capacity * (temperature_prod[1] - temperature_inj)   );
 //since Q_i=const, Delta t = const we can minimize
   */
  functional_value =  (number_of_time_steps_for_production * 1/(double)db["pump_efficiency"] *
          (pressure_values_injection_wells - pressure_values_production_wells)
          - (double) db["fluid_density"] * (double) db["fluid_heat_capacity"] * Delta_Temp) 
  +  alpha * std::abs(distance - 2000.);

  Output::print("functional_value: ", functional_value);

  //write to stream
  std::ofstream outputFile("temperature_values_at_production_well_and_net_energy_for_disctance" + std::to_string(distance) +  ".txt");
  std::copy(this->temperature_production_well_at_time_steps.begin(), this->temperature_production_well_at_time_steps.end(), std::ostream_iterator<int>(outputFile, "\n"));
  //std::copy(functional_value, std::ostream_iterator<int>(outputFile, "\n"));
  outputFile << functional_value;


  std::vector<double> cost_functional = db["cost_functional"];

  auto temperature = tcd_primal.get_function();
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


  //  double functional_value = 1.;
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


/** ************************************************************************ */
template <int d>
GeothermalPlantsPositionOptimization<d>::sinks::sinks(double eps_delta_fct, double well_radius,
                       std::array<double, d> center_point, size_t Num_circle_points, const TCollection* Coll)
:center(center_point), meine_punkte(Num_circle_points)
{
  double radius = eps_delta_fct*well_radius;
  for (size_t k = 0; k < Num_circle_points; k++)
  {
    this->meine_punkte[k].coordinates[0] = radius*cos(k*2.*Pi/(Num_circle_points + 1)) + center_point[0];
    this->meine_punkte[k].coordinates[1] = radius*sin(k*2.*Pi/(Num_circle_points + 1)) + center_point[1];
    if (d==3)
    {
      this->meine_punkte[k].coordinates[2] = center_point[2];
    }

    int N_Cells = Coll->GetN_Cells();
    for (int i = 0; i < N_Cells; i++)
    {
      auto cell = Coll->GetCell(i);
      {
        if (cell->PointInCell(this->meine_punkte[k].coordinates[0], this->meine_punkte[k].coordinates[1]
#ifdef __3D__
        ,this->meine_punkte[k].coordinates[2]
#endif
                                           ))
        {
          this->meine_punkte[k].cell = cell;
          this->meine_punkte[k].cell_index = i;
        }
      }
    }
  }
}

/** ************************************************************************ */
template <int d>
void GeothermalPlantsPositionOptimization<d>::sinks::find_average_and_min_along_circle(
        const FEFunction* function, double & average, double & min)
{
  min = 1.e9;
  average = 0.;
  int Num_circle_points = this->meine_punkte.size();
  double val;

  for (size_t k = 0; k < Num_circle_points; k++)
  {
    auto current_point = this->meine_punkte[k];
 function->FindValueLocal(current_point.cell, current_point.cell_index, current_point.coordinates[0],
         current_point.coordinates[1],
#ifdef __3D__
                                  current_point.coordinates[2],
#endif
                                  &val);
  if (val < min)
     min = val;

  average += val;
  }

  if (min == 1.e9)
    ErrThrow("no minimum found");

average /= Num_circle_points;
}





#ifdef __3D__
template class GeothermalPlantsPositionOptimization<3>;
#else
template class GeothermalPlantsPositionOptimization<2>;
#endif
