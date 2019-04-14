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

/** ************************************************************************ */
template <int d>
TCD_Temperature<d>::TCD_Temperature(const TDomain& domain,
        const ParameterDatabase& param_db,
        Example_TimeCD example
        //#ifdef _MPI
        //                                    ,int maxSubDomainPerDof
        //#endif
        )
: TimeConvectionDiffusion<d>(domain, param_db, example
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

  this->TimeConvectionDiffusion<d>::db.merge(param_db, true);

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
template <int d>
void mapping_local_parameters(const double *in, double *out)
{
  // coordinates:  x at in[0], y at in[1]
  out[0] = in[d];
  out[1] = in[d+1];
  if (d == 3)
  {
    //z at in[2]
    out[2] = in[d+2];
  }
}


/** ************************************************************************ */
void temperature_coefficients(int n_points, double *x, double *y,
#ifdef __3D__
           double *z,
#endif
        double **parameters, double **coeffs,
        double distance, double nu, double transversal_dispersion_factor, double longitudinal_dispersion_factor,
        double well_radius, double temperature_injection_well, double delta_fct_eps_factor,
        double fluid_density, double fluid_heat_capacity)
{
  for(int i = 0; i < n_points; ++i)
  {
    //another approx. for domain [0, 10] x [0, 6]
    // bubble diameter = epsDelta
    // bubble height = 1/(epsDelta*epsDelta)
    //double r_well = 0.2; // 20cm
    double epsDelta = delta_fct_eps_factor * well_radius; //10*r_well; // 25*r_well; //50*r_well;
    //double T_in = 303.15; //injection temperature  = 30 + 273.15;


#ifdef __2D__
      std::array<double, 2> center_source = {{5000. - (distance)/2., 3000.}};
#else
      std::array<double, 3> center_source = {{5000. - (distance/2.), 3000., 3000.}};
      double z_distance_to_source = std::pow(std::abs(z[i]-center_source[2]), 2);
#endif
      double x_distance_to_source = std::pow(std::abs(x[i]-center_source[0]), 2);
      double y_distance_to_source = std::pow(std::abs(y[i]-center_source[1]), 2);


      /*     bool at_source = x_distance_to_source + y_distance_to_source
     #ifdef __3D__
                   + z_distance_to_source
     #endif
                   < epsDelta*epsDelta;
                   */

      bool at_source = (x_distance_to_source < epsDelta*epsDelta) * (y_distance_to_source < epsDelta*epsDelta)
    #ifdef __3D__
     // OLD 05.02.19             * (z_distance_to_source < epsDelta*epsDelta)
    #endif
     ;

double norm_u = sqrt(parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]
#ifdef __3D__
              + parameters[i][2]*parameters[i][2]
#endif
);
    coeffs[i][0] = nu + transversal_dispersion_factor * norm_u;  // diffusion
    coeffs[i][1] = parameters[i][0]; // convection, x-direction
    coeffs[i][2] = parameters[i][1]; // convection, y-direction

#ifdef __2D__
      coeffs[i][3] = 0.; // reaction
      coeffs[i][4] = 0.; // right-hand side
      if(norm_u)
      coeffs[i][5] = //fluid_density * fluid_heat_capacity *
              (longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
      else
      coeffs[i][5] = 0.;

     /* cout << "!!!!      !!!!     !!!! parameters[i][0]: " << parameters[i][0] <<"parameters[i][1]: "<< parameters[i][1] <<endl;
      cout << "!!!!      !!!!     !!!!parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]: " << parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]<<endl;
      cout << "!!!!      !!!!     !!!! fluid_density * fluid_heat_capacity * (longitudinal_dispersion_factor - transversal_dispersion_factor): " <<  fluid_density * fluid_heat_capacity * (longitudinal_dispersion_factor - transversal_dispersion_factor)<<endl;
      cout << "!!!!      !!!!     !!!! coeffs[i][5]: " <<  coeffs[i][5] <<endl;
      */
#else
      coeffs[i][3] = parameters[i][2]; // convection, z-direction
      coeffs[i][4] = 0.; // reaction
      coeffs[i][5] = 0.; // right-hand side

      if(norm_u)
      coeffs[i][6] = //fluid_density * fluid_heat_capacity *
              (longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
      else
      coeffs[i][6] = 0.;
#endif

      if(at_source)
    {
#ifdef __2D__
        double magnitude = cos(Pi*(x[i] - center_source[0])/epsDelta) + 1;
        magnitude *= cos(Pi*(y[i] - center_source[1])/epsDelta) + 1;
        magnitude /= 4.*epsDelta*epsDelta;

        coeffs[i][3] = magnitude; // reaction
        coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

        //double penalty_factor = 1000.;
        //coeffs[i][3] = 1 * penalty_factor; // reaction
        //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
#else
        //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
        //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
        //magnitude /= 4.*epsDelta*epsDelta;
        //coeffs[i][3] += magnitude; // reaction
        //coeffs[i][4] -= magnitude * T_in; // source
      /*  double penalty_factor = 1000.;
        coeffs[i][4] = 1 * penalty_factor; // reaction
        coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
        */
        double magnitude = cos(Pi*(x[i] - center_source[0])/epsDelta) + 1;
        magnitude *= cos(Pi*(y[i] - center_source[1])/epsDelta) + 1;
       //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
        ////magnitude *= 1 + 1;
        magnitude /= 4.*epsDelta*epsDelta;

        coeffs[i][4] = magnitude; // reaction
        coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
#endif

    }
  }
}

/** ************************************************************************ */
void temperature_coefficients_hexagon(int n_points, double *x, double *y,
#ifdef __3D__
				      double *z,
#endif
				      double **parameters, double **coeffs,
				      double distance, double nu, double transversal_dispersion_factor,
				      double longitudinal_dispersion_factor,
				      double well_radius, double temperature_injection_well, double delta_fct_eps_factor,
				      double fluid_density, double fluid_heat_capacity)
{

   double domain_Lx = 10000.;
  double domain_Ly = 6000.;
  size_t n_wells = 6;
  
  // set well position
  std::vector<double> singular_x, singular_y, singular_sign, flow_rate;

  singular_sign.push_back(-1.); // injection
  singular_sign.push_back(1.);  // production
  singular_sign.push_back(1.);
  singular_sign.push_back(-1.);
  singular_sign.push_back(1);
  singular_sign.push_back(1);
  
  double epsDelta = delta_fct_eps_factor * well_radius;
  
  double x0 = domain_Lx/2.;
  double y0 = domain_Ly/2.;
  double pi = acos(-1.);
  int count = 0;
  for (unsigned int k1 = 0; k1 < n_wells; k1++) {
    double xk = x0 + distance*cos(2.*pi*k1/n_wells);
    double yk = y0 + distance*sin(2.*pi*k1/n_wells);
    if ( singular_sign[k1] == -1) {
      singular_x.push_back(xk);
      singular_y.push_back(yk);
    }
  }
  
  
  for(int i = 0; i < n_points; ++i)
  {
    double norm_u = sqrt(parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]);
    coeffs[i][0] = nu + transversal_dispersion_factor * norm_u;  // diffusion
    coeffs[i][1] = parameters[i][0]; // convection, x-direction
    coeffs[i][2] = parameters[i][1]; // convection, y-direction

    coeffs[i][3] = 0.; // reaction
    coeffs[i][4] = 0.; // right-hand side
    if(norm_u)
      coeffs[i][5] = //fluid_density * fluid_heat_capacity *
	(longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
    else
      coeffs[i][5] = 0.;
    
    // Add contribution of singular sources
    double epsDelta = delta_fct_eps_factor * well_radius; //10*r_well; // 25*r_well; //50*r_well;
    
    for (unsigned int m = 0; m < singular_x.size(); m++)
    {
      double x_distance_to_source = std::pow(std::abs(x[i]-singular_x[m]), 2);
      double y_distance_to_source = std::pow(std::abs(y[i]-singular_y[m]), 2);
      bool at_source = (x_distance_to_source < epsDelta*epsDelta) * (y_distance_to_source < epsDelta*epsDelta);
      if(at_source)
      {
	double magnitude = cos(Pi*(x[i] - singular_x[m])/epsDelta) + 1;
	magnitude *= cos(Pi*(y[i] - singular_y[m])/epsDelta) + 1;
	magnitude /= 4.*epsDelta*epsDelta;
	
        coeffs[i][3] = magnitude; // reaction
        coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side
      }
    }
    
  }
}


/** ************************************************************************ */
void temperature_coefficients_lattice(int n_points, double *x, double *y,
#ifdef __3D__
				      double *z,
#endif
				      double **parameters, double **coeffs,
				      double distance, double nu, double transversal_dispersion_factor,
				      double longitudinal_dispersion_factor,
				      double well_radius, double temperature_injection_well, double delta_fct_eps_factor,
				      double fluid_density, double fluid_heat_capacity)
{

  double domain_Lx = 10000.;
  double domain_Ly = 6000.;
  size_t n_wells_per_row = 4;
  // set well position
  std::vector<double> singular_x, singular_y;

  double x0 = domain_Lx/2.- (n_wells_per_row-1)/2.*distance;
  double y0 = domain_Ly/2.- (n_wells_per_row-1)/2.*distance;
  int count = 0;
  for (unsigned int k1=0; k1 < n_wells_per_row; k1++) {
    for (unsigned int k2=0; k2 < n_wells_per_row; k2++) {
      double xk = x0 + distance*k2;
      double yk = y0 + distance*k1;
      int singular_sign = pow(-1,k1)*pow(-1,k2);
      // add singular source only at injection wells
      if (singular_sign > 0)
      {
	     singular_x.push_back(xk);
	     singular_y.push_back(yk);
      }

    }
  }
  
  for(int i = 0; i < n_points; ++i)
  {
    double norm_u = sqrt(parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]);
    coeffs[i][0] = nu + transversal_dispersion_factor * norm_u;  // diffusion
    coeffs[i][1] = parameters[i][0]; // convection, x-direction
    coeffs[i][2] = parameters[i][1]; // convection, y-direction

    coeffs[i][3] = 0.; // reaction
    coeffs[i][4] = 0.; // right-hand side
    if(norm_u)
      coeffs[i][5] = //fluid_density * fluid_heat_capacity *
	(longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
    else
      coeffs[i][5] = 0.;
    
    // Add contribution of singular sources
    double epsDelta = delta_fct_eps_factor * well_radius; //10*r_well; // 25*r_well; //50*r_well;
    
    for (unsigned int m = 0; m < singular_x.size(); m++)
    {
      double x_distance_to_source = std::pow(std::abs(x[i]-singular_x[m]), 2);
      double y_distance_to_source = std::pow(std::abs(y[i]-singular_y[m]), 2);
      bool at_source = (x_distance_to_source < epsDelta*epsDelta) * (y_distance_to_source < epsDelta*epsDelta);
      if(at_source)
      {
	double magnitude = cos(Pi*(x[i] - singular_x[m])/epsDelta) + 1;
	magnitude *= cos(Pi*(y[i] - singular_y[m])/epsDelta) + 1;
	magnitude /= 4.*epsDelta*epsDelta;
	
        coeffs[i][3] = magnitude; // reaction
        coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side
      }
    }
    
  }
}


/** ************************************************************************ */
void temperature_coefficients_2doublets(int n_points, double *x, double *y,
#ifdef __3D__
           double *z,
#endif
        double **parameters, double **coeffs,
        double center_x_moving_doublet, double nu,
        double transversal_dispersion_factor, double longitudinal_dispersion_factor,
        double well_radius, double temperature_injection_well, double delta_fct_eps_factor,
        double fluid_density, double fluid_heat_capacity, double well_distance)
{
  for(int i = 0; i < n_points; ++i)
  {
    //another approx. for domain [0, 10] x [0, 6]
    // bubble diameter = epsDelta
    // bubble height = 1/(epsDelta*epsDelta)
    //double r_well = 0.2; // 20cm
    double epsDelta = delta_fct_eps_factor * well_radius; //10*r_well; // 25*r_well; //50*r_well;
    //double T_in = 303.15; //injection temperature  = 30 + 273.15;


#ifdef __2D__
    std::array<double, 2> center_fixed_source = {{5000. - (well_distance)/2., 2500.}};
    std::array<double, 2> center_moving_source = {{center_x_moving_doublet - (well_distance)/2., 3500.}};
#else
    std::array<double, 3> center_fixed_source = {{5000. - (well_distance/2.), 2500., 250.}};
    std::array<double, 3> center_moving_source = {{center_x_moving_doublet - (well_distance/2.), 3500., 250.}};
    double z_distance_to_fixed_source = std::pow(std::abs(z[i]-center_fixed_source[2]), 2);
    double z_distance_to_moving_source = std::pow(std::abs(z[i]-center_moving_source[2]), 2);
#endif
    double x_distance_to_fixed_source = std::pow(std::abs(x[i]-center_fixed_source[0]), 2);
    double y_distance_to_fixed_source = std::pow(std::abs(y[i]-center_fixed_source[1]), 2);
    double x_distance_to_moving_source = std::pow(std::abs(x[i]-center_moving_source[0]), 2);
    double y_distance_to_moving_source = std::pow(std::abs(y[i]-center_moving_source[1]), 2);


    /*     bool at_source = x_distance_to_source + y_distance_to_source
     #ifdef __3D__
                   + z_distance_to_source
     #endif
                   < epsDelta*epsDelta;
     */

    bool at_fixed_source = (x_distance_to_fixed_source < epsDelta*epsDelta) * (y_distance_to_fixed_source < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_fixed_source < epsDelta*epsDelta)
#endif
             ;

    bool at_moving_source = (x_distance_to_moving_source < epsDelta*epsDelta) * (y_distance_to_moving_source < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_moving_source < epsDelta*epsDelta)
#endif
             ;

    double norm_u = sqrt(parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]
#ifdef __3D__
                                                                                            + parameters[i][2]*parameters[i][2]
#endif
    );
    coeffs[i][0] = nu + transversal_dispersion_factor * norm_u;  // diffusion
    coeffs[i][1] = parameters[i][0]; // convection, x-direction
    coeffs[i][2] = parameters[i][1]; // convection, y-direction

#ifdef __2D__
    coeffs[i][3] = 0.; // reaction
    coeffs[i][4] = 0.; // right-hand side
    if(norm_u)
      coeffs[i][5] = //fluid_density * fluid_heat_capacity *
              (longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
    else
      coeffs[i][5] = 0.;

    /* cout << "!!!!      !!!!     !!!! parameters[i][0]: " << parameters[i][0] <<"parameters[i][1]: "<< parameters[i][1] <<endl;
      cout << "!!!!      !!!!     !!!!parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]: " << parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]<<endl;
      cout << "!!!!      !!!!     !!!! fluid_density * fluid_heat_capacity * (longitudinal_dispersion_factor - transversal_dispersion_factor): " <<  fluid_density * fluid_heat_capacity * (longitudinal_dispersion_factor - transversal_dispersion_factor)<<endl;
      cout << "!!!!      !!!!     !!!! coeffs[i][5]: " <<  coeffs[i][5] <<endl;
     */
#else
    coeffs[i][3] = parameters[i][2]; // convection, z-direction
    coeffs[i][4] = 0.; // reaction
    coeffs[i][5] = 0.; // right-hand side

    if(norm_u)
      coeffs[i][6] = //fluid_density * fluid_heat_capacity *
              (longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
    else
      coeffs[i][6] = 0.;
#endif

    if(at_fixed_source)
    {
#ifdef __2D__
      double magnitude = cos(Pi*(x[i] - center_fixed_source[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_fixed_source[1])/epsDelta) + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][3] = magnitude; // reaction
      coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

      //double penalty_factor = 1000.;
      //coeffs[i][3] = 1 * penalty_factor; // reaction
      //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
#else
      //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
      //magnitude /= 4.*epsDelta*epsDelta;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      /*  double penalty_factor = 1000.;
        coeffs[i][4] = 1 * penalty_factor; // reaction
        coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
       */
      double magnitude = cos(Pi*(x[i] - center_fixed_source[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_fixed_source[1])/epsDelta) + 1;
      //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
      ////magnitude *= 1 + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][4] = magnitude; // reaction
      coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
#endif
    }
    if(at_moving_source)
    {
#ifdef __2D__
      double magnitude = cos(Pi*(x[i] - center_moving_source[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_source[1])/epsDelta) + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][3] = magnitude; // reaction
      coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

      //double penalty_factor = 1000.;
      //coeffs[i][3] = 1 * penalty_factor; // reaction
      //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
#else
      //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
      //magnitude /= 4.*epsDelta*epsDelta;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      /*  double penalty_factor = 1000.;
        coeffs[i][4] = 1 * penalty_factor; // reaction
        coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
       */
      double magnitude = cos(Pi*(x[i] - center_moving_source[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_source[1])/epsDelta) + 1;
      //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
      ////magnitude *= 1 + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][4] = magnitude; // reaction
      coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
#endif
    }
  }
}

/** ************************************************************************ */
void temperature_coefficients_3doubledoublets(int n_points, double *x, double *y,
#ifdef __3D__
           double *z,
#endif
        double **parameters, double **coeffs,
        double center_x_moving_doublet_bottom_row, double center_x_moving_doublet_top_row, double nu,
        double transversal_dispersion_factor, double longitudinal_dispersion_factor,
        double well_radius, double temperature_injection_well, double delta_fct_eps_factor,
        double fluid_density, double fluid_heat_capacity, double well_distance)
{
  for(int i = 0; i < n_points; ++i)
  {
    //another approx. for domain [0, 10] x [0, 6]
    // bubble diameter = epsDelta
    // bubble height = 1/(epsDelta*epsDelta)
    //double r_well = 0.2; // 20cm
    double epsDelta = delta_fct_eps_factor * well_radius; //10*r_well; // 25*r_well; //50*r_well;
    //double T_in = 303.15; //injection temperature  = 30 + 273.15;


#ifdef __2D__
    std::array<double, 2> center_fixed_source_1 = {{5000. - (well_distance)*3/2., 3000.}};
    std::array<double, 2> center_fixed_source_2 = {{5000. + (well_distance)/2., 3000.}};
    std::array<double, 2> center_moving_upper_source_1 = {{center_x_moving_doublet_top_row - (well_distance)*3/2., 4000.}};
    std::array<double, 2> center_moving_upper_source_2 = {{center_x_moving_doublet_top_row + (well_distance)/2., 4000.}};
    std::array<double, 2> center_moving_lower_source_1 = {{center_x_moving_doublet_bottom_row - (well_distance)*3/2., 2000.}};
    std::array<double, 2> center_moving_lower_source_2 = {{center_x_moving_doublet_bottom_row + (well_distance)/2., 2000.}};
#else
    std::array<double, 3> center_fixed_source_1 = {{5000. - (well_distance)*3/2., 3000., 250.}};
    std::array<double, 3> center_fixed_source_2 = {{5000. + (well_distance)/2., 3000., 250.}};
    std::array<double, 3> center_moving_upper_source_1 = {{center_x_moving_doublet_top_row - (well_distance)*3/2., 4000., 250.}};
    std::array<double, 3> center_moving_upper_source_2 = {{center_x_moving_doublet_top_row + (well_distance)/2., 4000., 250.}};
    std::array<double, 3> center_moving_lower_source_1 = {{center_x_moving_doublet_bottom_row - (well_distance)*3/2., 2000., 250.}};
    std::array<double, 3> center_moving_lower_source_2 = {{center_x_moving_doublet_bottom_row + (well_distance)/2., 2000., 250.}};
    double z_distance_to_fixed_source_1 = std::pow(std::abs(z[i]-center_fixed_source_1[2]), 2);
    double z_distance_to_fixed_source_2 = std::pow(std::abs(z[i]-center_fixed_source_2[2]), 2);
    double z_distance_to_moving_upper_source_1 = std::pow(std::abs(z[i]-center_moving_upper_source_1[2]), 2);
    double z_distance_to_moving_upper_source_2 = std::pow(std::abs(z[i]-center_moving_upper_source_2[2]), 2);
    double z_distance_to_moving_lower_source_1 = std::pow(std::abs(z[i]-center_moving_lower_source_1[2]), 2);
    double z_distance_to_moving_lower_source_2 = std::pow(std::abs(z[i]-center_moving_lower_source_2[2]), 2);
#endif
    double x_distance_to_fixed_source_1 = std::pow(std::abs(x[i]-center_fixed_source_1[0]), 2);
    double y_distance_to_fixed_source_1 = std::pow(std::abs(y[i]-center_fixed_source_1[1]), 2);
    double x_distance_to_fixed_source_2 = std::pow(std::abs(x[i]-center_fixed_source_2[0]), 2);
    double y_distance_to_fixed_source_2 = std::pow(std::abs(y[i]-center_fixed_source_2[1]), 2);
    double x_distance_to_moving_upper_source_1 = std::pow(std::abs(x[i]-center_moving_upper_source_1[0]), 2);
    double y_distance_to_moving_upper_source_1 = std::pow(std::abs(y[i]-center_moving_upper_source_1[1]), 2);
    double x_distance_to_moving_upper_source_2 = std::pow(std::abs(x[i]-center_moving_upper_source_2[0]), 2);
    double y_distance_to_moving_upper_source_2 = std::pow(std::abs(y[i]-center_moving_upper_source_2[1]), 2);
    double x_distance_to_moving_lower_source_1 = std::pow(std::abs(x[i]-center_moving_lower_source_1[0]), 2);
    double y_distance_to_moving_lower_source_1 = std::pow(std::abs(y[i]-center_moving_lower_source_1[1]), 2);
    double x_distance_to_moving_lower_source_2 = std::pow(std::abs(x[i]-center_moving_lower_source_2[0]), 2);
    double y_distance_to_moving_lower_source_2 = std::pow(std::abs(y[i]-center_moving_lower_source_2[1]), 2);


    /*     bool at_source = x_distance_to_source + y_distance_to_source
     #ifdef __3D__
                   + z_distance_to_source
     #endif
                   < epsDelta*epsDelta;
     */

    bool at_fixed_source_1 = (x_distance_to_fixed_source_1 < epsDelta*epsDelta) * (y_distance_to_fixed_source_1 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_fixed_source < epsDelta*epsDelta)
#endif
             ;
    bool at_fixed_source_2 = (x_distance_to_fixed_source_2 < epsDelta*epsDelta) * (y_distance_to_fixed_source_2 < epsDelta*epsDelta)
  #ifdef __3D__
               // OLD 05.02.19             * (z_distance_to_fixed_source < epsDelta*epsDelta)
  #endif
               ;

    bool at_moving_source_upper_1 = (x_distance_to_moving_upper_source_1 < epsDelta*epsDelta) * (y_distance_to_moving_upper_source_1 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_moving_source < epsDelta*epsDelta)
#endif
             ;

    bool at_moving_source_upper_2 = (x_distance_to_moving_upper_source_2 < epsDelta*epsDelta) * (y_distance_to_moving_upper_source_2 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_moving_source < epsDelta*epsDelta)
#endif
             ;

    bool at_moving_source_lower_1 = (x_distance_to_moving_lower_source_1 < epsDelta*epsDelta) * (y_distance_to_moving_lower_source_1 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_moving_source < epsDelta*epsDelta)
#endif
             ;

    bool at_moving_source_lower_2 = (x_distance_to_moving_lower_source_2 < epsDelta*epsDelta) * (y_distance_to_moving_lower_source_2 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_moving_source < epsDelta*epsDelta)
#endif
             ;

    double norm_u = sqrt(parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]
#ifdef __3D__
                                                                                            + parameters[i][2]*parameters[i][2]
#endif
    );
    coeffs[i][0] = nu + transversal_dispersion_factor * norm_u;  // diffusion
    coeffs[i][1] = parameters[i][0]; // convection, x-direction
    coeffs[i][2] = parameters[i][1]; // convection, y-direction

#ifdef __2D__
    coeffs[i][3] = 0.; // reaction
    coeffs[i][4] = 0.; // right-hand side
    if(norm_u)
      coeffs[i][5] = //fluid_density * fluid_heat_capacity *
              (longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
    else
      coeffs[i][5] = 0.;

    /* cout << "!!!!      !!!!     !!!! parameters[i][0]: " << parameters[i][0] <<", parameters[i][1]: "<< parameters[i][1] <<endl;
     */
#else
    coeffs[i][3] = parameters[i][2]; // convection, z-direction
    coeffs[i][4] = 0.; // reaction
    coeffs[i][5] = 0.; // right-hand side

    if(norm_u)
      coeffs[i][6] = //fluid_density * fluid_heat_capacity *
              (longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
    else
      coeffs[i][6] = 0.;
#endif

    if(at_fixed_source_1)
    {
#ifdef __2D__
      double magnitude = cos(Pi*(x[i] - center_fixed_source_1[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_fixed_source_1[1])/epsDelta) + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][3] = magnitude; // reaction
      coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

      //double penalty_factor = 1000.;
      //coeffs[i][3] = 1 * penalty_factor; // reaction
      //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
#else
      //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
      //magnitude /= 4.*epsDelta*epsDelta;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      /*  double penalty_factor = 1000.;
        coeffs[i][4] = 1 * penalty_factor; // reaction
        coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
       */
      double magnitude = cos(Pi*(x[i] - center_fixed_source_1[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_fixed_source_1[1])/epsDelta) + 1;
      //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
      ////magnitude *= 1 + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][4] = magnitude; // reaction
      coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
#endif
    }
    if(at_fixed_source_2)
       {
   #ifdef __2D__
         double magnitude = cos(Pi*(x[i] - center_fixed_source_2[0])/epsDelta) + 1;
         magnitude *= cos(Pi*(y[i] - center_fixed_source_2[1])/epsDelta) + 1;
         magnitude /= 4.*epsDelta*epsDelta;

         coeffs[i][3] = magnitude; // reaction
         coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

         //double penalty_factor = 1000.;
         //coeffs[i][3] = 1 * penalty_factor; // reaction
         //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
   #else
         //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
         //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
         //magnitude /= 4.*epsDelta*epsDelta;
         //coeffs[i][3] += magnitude; // reaction
         //coeffs[i][4] -= magnitude * T_in; // source
         /*  double penalty_factor = 1000.;
           coeffs[i][4] = 1 * penalty_factor; // reaction
           coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
          */
         double magnitude = cos(Pi*(x[i] - center_fixed_source_2[0])/epsDelta) + 1;
         magnitude *= cos(Pi*(y[i] - center_fixed_source_2[1])/epsDelta) + 1;
         //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
         ////magnitude *= 1 + 1;
         magnitude /= 4.*epsDelta*epsDelta;

         coeffs[i][4] = magnitude; // reaction
         coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
   #endif
       }
    if(at_moving_source_upper_1)
    {
#ifdef __2D__
      double magnitude = cos(Pi*(x[i] - center_moving_upper_source_1[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_upper_source_1[1])/epsDelta) + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][3] = magnitude; // reaction
      coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

      //double penalty_factor = 1000.;
      //coeffs[i][3] = 1 * penalty_factor; // reaction
      //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
#else
      //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
      //magnitude /= 4.*epsDelta*epsDelta;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      /*  double penalty_factor = 1000.;
        coeffs[i][4] = 1 * penalty_factor; // reaction
        coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
       */
      double magnitude = cos(Pi*(x[i] - center_moving_upper_source_1[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_upper_source_1[1])/epsDelta) + 1;
      //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
      ////magnitude *= 1 + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][4] = magnitude; // reaction
      coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
#endif
    } if(at_moving_source_upper_2)
    {
#ifdef __2D__
      double magnitude = cos(Pi*(x[i] - center_moving_upper_source_2[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_upper_source_2[1])/epsDelta) + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][3] = magnitude; // reaction
      coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

      //double penalty_factor = 1000.;
      //coeffs[i][3] = 1 * penalty_factor; // reaction
      //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
#else
      //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
      //magnitude /= 4.*epsDelta*epsDelta;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      /*  double penalty_factor = 1000.;
        coeffs[i][4] = 1 * penalty_factor; // reaction
        coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
       */
      double magnitude = cos(Pi*(x[i] - center_moving_upper_source_2[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_upper_source_2[1])/epsDelta) + 1;
      //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
      ////magnitude *= 1 + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][4] = magnitude; // reaction
      coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
#endif
    }
    if(at_moving_source_lower_1)
       {
   #ifdef __2D__
         double magnitude = cos(Pi*(x[i] - center_moving_lower_source_1[0])/epsDelta) + 1;
         magnitude *= cos(Pi*(y[i] - center_moving_lower_source_1[1])/epsDelta) + 1;
         magnitude /= 4.*epsDelta*epsDelta;

         coeffs[i][3] = magnitude; // reaction
         coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

         //double penalty_factor = 1000.;
         //coeffs[i][3] = 1 * penalty_factor; // reaction
         //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
   #else
         //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
         //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
         //magnitude /= 4.*epsDelta*epsDelta;
         //coeffs[i][3] += magnitude; // reaction
         //coeffs[i][4] -= magnitude * T_in; // source
         /*  double penalty_factor = 1000.;
           coeffs[i][4] = 1 * penalty_factor; // reaction
           coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
          */
         double magnitude = cos(Pi*(x[i] - center_moving_lower_source_1[0])/epsDelta) + 1;
         magnitude *= cos(Pi*(y[i] - center_moving_lower_source_1[1])/epsDelta) + 1;
         //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
         ////magnitude *= 1 + 1;
         magnitude /= 4.*epsDelta*epsDelta;

         coeffs[i][4] = magnitude; // reaction
         coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
   #endif
       }
    if(at_moving_source_lower_2)
          {
      #ifdef __2D__
            double magnitude = cos(Pi*(x[i] - center_moving_lower_source_2[0])/epsDelta) + 1;
            magnitude *= cos(Pi*(y[i] - center_moving_lower_source_2[1])/epsDelta) + 1;
            magnitude /= 4.*epsDelta*epsDelta;

            coeffs[i][3] = magnitude; // reaction
            coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

            //double penalty_factor = 1000.;
            //coeffs[i][3] = 1 * penalty_factor; // reaction
            //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
      #else
            //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
            //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
            //magnitude /= 4.*epsDelta*epsDelta;
            //coeffs[i][3] += magnitude; // reaction
            //coeffs[i][4] -= magnitude * T_in; // source
            /*  double penalty_factor = 1000.;
              coeffs[i][4] = 1 * penalty_factor; // reaction
              coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
             */
            double magnitude = cos(Pi*(x[i] - center_moving_lower_source_2[0])/epsDelta) + 1;
            magnitude *= cos(Pi*(y[i] - center_moving_lower_source_2[1])/epsDelta) + 1;
            //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
            ////magnitude *= 1 + 1;
            magnitude /= 4.*epsDelta*epsDelta;

            coeffs[i][4] = magnitude; // reaction
            coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
      #endif
          }
  }
}


/** ************************************************************************ */
void temperature_coefficients_3doubledoublets_4controls(int n_points, double *x, double *y,
#ifdef __3D__
           double *z,
#endif
        double **parameters, double **coeffs,
        double center_x_moving_doublet_bottom_row, double center_x_moving_doublet_top_row, double nu,
        double transversal_dispersion_factor, double longitudinal_dispersion_factor,
        double well_radius, double temperature_injection_well, double delta_fct_eps_factor,
        double fluid_density, double fluid_heat_capacity, double well_distance,
        double height_of_bottom_row, double height_of_top_row)
{
for(int i = 0; i < n_points; ++i)
  {
    //another approx. for domain [0, 10] x [0, 6]
    // bubble diameter = epsDelta
    // bubble height = 1/(epsDelta*epsDelta)
    //double r_well = 0.2; // 20cm
    double epsDelta = delta_fct_eps_factor * well_radius; //10*r_well; // 25*r_well; //50*r_well;
    //double T_in = 303.15; //injection temperature  = 30 + 273.15;


#ifdef __2D__
    std::array<double, 2> center_fixed_source_1 = {{5000. - (well_distance)*3/2., 3000.}};
    std::array<double, 2> center_fixed_source_2 = {{5000. + (well_distance)/2., 3000.}};
    std::array<double, 2> center_moving_upper_source_1 = {{center_x_moving_doublet_top_row - (well_distance)*3/2., height_of_top_row}};
    std::array<double, 2> center_moving_upper_source_2 = {{center_x_moving_doublet_top_row + (well_distance)/2., height_of_top_row}};
    std::array<double, 2> center_moving_lower_source_1 = {{center_x_moving_doublet_bottom_row - (well_distance)*3/2., height_of_bottom_row}};
    std::array<double, 2> center_moving_lower_source_2 = {{center_x_moving_doublet_bottom_row + (well_distance)/2., height_of_bottom_row}};
#else
    std::array<double, 3> center_fixed_source_1 = {{5000. - (well_distance)*3/2., 3000., 250.}};
    std::array<double, 3> center_fixed_source_2 = {{5000. + (well_distance)/2., 3000., 250.}};
    std::array<double, 3> center_moving_upper_source_1 = {{center_x_moving_doublet_top_row - (well_distance)*3/2., height_of_top_row, 250.}};
    std::array<double, 3> center_moving_upper_source_2 = {{center_x_moving_doublet_top_row + (well_distance)/2., height_of_top_row, 250.}};
    std::array<double, 3> center_moving_lower_source_1 = {{center_x_moving_doublet_bottom_row - (well_distance)*3/2., height_of_bottom_row, 250.}};
    std::array<double, 3> center_moving_lower_source_2 = {{center_x_moving_doublet_bottom_row + (well_distance)/2., height_of_bottom_row, 250.}};
    double z_distance_to_fixed_source_1 = std::pow(std::abs(z[i]-center_fixed_source_1[2]), 2);
    double z_distance_to_fixed_source_2 = std::pow(std::abs(z[i]-center_fixed_source_2[2]), 2);
    double z_distance_to_moving_upper_source_1 = std::pow(std::abs(z[i]-center_moving_upper_source_1[2]), 2);
    double z_distance_to_moving_upper_source_2 = std::pow(std::abs(z[i]-center_moving_upper_source_2[2]), 2);
    double z_distance_to_moving_lower_source_1 = std::pow(std::abs(z[i]-center_moving_lower_source_1[2]), 2);
    double z_distance_to_moving_lower_source_2 = std::pow(std::abs(z[i]-center_moving_lower_source_2[2]), 2);
#endif
    double x_distance_to_fixed_source_1 = std::pow(std::abs(x[i]-center_fixed_source_1[0]), 2);
    double y_distance_to_fixed_source_1 = std::pow(std::abs(y[i]-center_fixed_source_1[1]), 2);
    double x_distance_to_fixed_source_2 = std::pow(std::abs(x[i]-center_fixed_source_2[0]), 2);
    double y_distance_to_fixed_source_2 = std::pow(std::abs(y[i]-center_fixed_source_2[1]), 2);
    double x_distance_to_moving_upper_source_1 = std::pow(std::abs(x[i]-center_moving_upper_source_1[0]), 2);
    double y_distance_to_moving_upper_source_1 = std::pow(std::abs(y[i]-center_moving_upper_source_1[1]), 2);
    double x_distance_to_moving_upper_source_2 = std::pow(std::abs(x[i]-center_moving_upper_source_2[0]), 2);
    double y_distance_to_moving_upper_source_2 = std::pow(std::abs(y[i]-center_moving_upper_source_2[1]), 2);
    double x_distance_to_moving_lower_source_1 = std::pow(std::abs(x[i]-center_moving_lower_source_1[0]), 2);
    double y_distance_to_moving_lower_source_1 = std::pow(std::abs(y[i]-center_moving_lower_source_1[1]), 2);
    double x_distance_to_moving_lower_source_2 = std::pow(std::abs(x[i]-center_moving_lower_source_2[0]), 2);
    double y_distance_to_moving_lower_source_2 = std::pow(std::abs(y[i]-center_moving_lower_source_2[1]), 2);


    /*     bool at_source = x_distance_to_source + y_distance_to_source
     #ifdef __3D__
                   + z_distance_to_source
     #endif
                   < epsDelta*epsDelta;
     */

    bool at_fixed_source_1 = (x_distance_to_fixed_source_1 < epsDelta*epsDelta) * (y_distance_to_fixed_source_1 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_fixed_source < epsDelta*epsDelta)
#endif
             ;
    bool at_fixed_source_2 = (x_distance_to_fixed_source_2 < epsDelta*epsDelta) * (y_distance_to_fixed_source_2 < epsDelta*epsDelta)
  #ifdef __3D__
               // OLD 05.02.19             * (z_distance_to_fixed_source < epsDelta*epsDelta)
  #endif
               ;

    bool at_moving_source_upper_1 = (x_distance_to_moving_upper_source_1 < epsDelta*epsDelta) * (y_distance_to_moving_upper_source_1 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_moving_source < epsDelta*epsDelta)
#endif
             ;

    bool at_moving_source_upper_2 = (x_distance_to_moving_upper_source_2 < epsDelta*epsDelta) * (y_distance_to_moving_upper_source_2 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_moving_source < epsDelta*epsDelta)
#endif
             ;

    bool at_moving_source_lower_1 = (x_distance_to_moving_lower_source_1 < epsDelta*epsDelta) * (y_distance_to_moving_lower_source_1 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_moving_source < epsDelta*epsDelta)
#endif
             ;

    bool at_moving_source_lower_2 = (x_distance_to_moving_lower_source_2 < epsDelta*epsDelta) * (y_distance_to_moving_lower_source_2 < epsDelta*epsDelta)
#ifdef __3D__
             // OLD 05.02.19             * (z_distance_to_moving_source < epsDelta*epsDelta)
#endif
             ;

    double norm_u = sqrt(parameters[i][0]*parameters[i][0] + parameters[i][1]*parameters[i][1]
#ifdef __3D__
                                                                                            + parameters[i][2]*parameters[i][2]
#endif
    );
    coeffs[i][0] = nu + transversal_dispersion_factor * norm_u;  // diffusion
    coeffs[i][1] = parameters[i][0]; // convection, x-direction
    coeffs[i][2] = parameters[i][1]; // convection, y-direction

#ifdef __2D__
    coeffs[i][3] = 0.; // reaction
    coeffs[i][4] = 0.; // right-hand side
    if(norm_u)
      coeffs[i][5] = //fluid_density * fluid_heat_capacity *
              (longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
    else
      coeffs[i][5] = 0.;

    /* cout << "!!!!      !!!!     !!!! parameters[i][0]: " << parameters[i][0] <<", parameters[i][1]: "<< parameters[i][1] <<endl;
     */
#else
    coeffs[i][3] = parameters[i][2]; // convection, z-direction
    coeffs[i][4] = 0.; // reaction
    coeffs[i][5] = 0.; // right-hand side

    if(norm_u)
      coeffs[i][6] = //fluid_density * fluid_heat_capacity *
              (longitudinal_dispersion_factor - transversal_dispersion_factor) * 1/norm_u;
    else
      coeffs[i][6] = 0.;
#endif

    if(at_fixed_source_1)
    {
#ifdef __2D__
      double magnitude = cos(Pi*(x[i] - center_fixed_source_1[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_fixed_source_1[1])/epsDelta) + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][3] = magnitude; // reaction
      coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

      //double penalty_factor = 1000.;
      //coeffs[i][3] = 1 * penalty_factor; // reaction
      //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
#else
      //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
      //magnitude /= 4.*epsDelta*epsDelta;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      /*  double penalty_factor = 1000.;
        coeffs[i][4] = 1 * penalty_factor; // reaction
        coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
       */
      double magnitude = cos(Pi*(x[i] - center_fixed_source_1[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_fixed_source_1[1])/epsDelta) + 1;
      //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
      ////magnitude *= 1 + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][4] = magnitude; // reaction
      coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
#endif
    }
    if(at_fixed_source_2)
       {
   #ifdef __2D__
         double magnitude = cos(Pi*(x[i] - center_fixed_source_2[0])/epsDelta) + 1;
         magnitude *= cos(Pi*(y[i] - center_fixed_source_2[1])/epsDelta) + 1;
         magnitude /= 4.*epsDelta*epsDelta;

         coeffs[i][3] = magnitude; // reaction
         coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

         //double penalty_factor = 1000.;
         //coeffs[i][3] = 1 * penalty_factor; // reaction
         //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
   #else
         //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
         //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
         //magnitude /= 4.*epsDelta*epsDelta;
         //coeffs[i][3] += magnitude; // reaction
         //coeffs[i][4] -= magnitude * T_in; // source
         /*  double penalty_factor = 1000.;
           coeffs[i][4] = 1 * penalty_factor; // reaction
           coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
          */
         double magnitude = cos(Pi*(x[i] - center_fixed_source_2[0])/epsDelta) + 1;
         magnitude *= cos(Pi*(y[i] - center_fixed_source_2[1])/epsDelta) + 1;
         //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
         ////magnitude *= 1 + 1;
         magnitude /= 4.*epsDelta*epsDelta;

         coeffs[i][4] = magnitude; // reaction
         coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
   #endif
       }
    if(at_moving_source_upper_1)
    {
#ifdef __2D__
      double magnitude = cos(Pi*(x[i] - center_moving_upper_source_1[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_upper_source_1[1])/epsDelta) + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][3] = magnitude; // reaction
      coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

      //double penalty_factor = 1000.;
      //coeffs[i][3] = 1 * penalty_factor; // reaction
      //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
#else
      //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
      //magnitude /= 4.*epsDelta*epsDelta;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      /*  double penalty_factor = 1000.;
        coeffs[i][4] = 1 * penalty_factor; // reaction
        coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
       */
      double magnitude = cos(Pi*(x[i] - center_moving_upper_source_1[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_upper_source_1[1])/epsDelta) + 1;
      //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
      ////magnitude *= 1 + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][4] = magnitude; // reaction
      coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
#endif
    } if(at_moving_source_upper_2)
    {
#ifdef __2D__
      double magnitude = cos(Pi*(x[i] - center_moving_upper_source_2[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_upper_source_2[1])/epsDelta) + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][3] = magnitude; // reaction
      coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

      //double penalty_factor = 1000.;
      //coeffs[i][3] = 1 * penalty_factor; // reaction
      //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
#else
      //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
      //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
      //magnitude /= 4.*epsDelta*epsDelta;
      //coeffs[i][3] += magnitude; // reaction
      //coeffs[i][4] -= magnitude * T_in; // source
      /*  double penalty_factor = 1000.;
        coeffs[i][4] = 1 * penalty_factor; // reaction
        coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
       */
      double magnitude = cos(Pi*(x[i] - center_moving_upper_source_2[0])/epsDelta) + 1;
      magnitude *= cos(Pi*(y[i] - center_moving_upper_source_2[1])/epsDelta) + 1;
      //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
      ////magnitude *= 1 + 1;
      magnitude /= 4.*epsDelta*epsDelta;

      coeffs[i][4] = magnitude; // reaction
      coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
#endif
    }
    if(at_moving_source_lower_1)
       {
   #ifdef __2D__
         double magnitude = cos(Pi*(x[i] - center_moving_lower_source_1[0])/epsDelta) + 1;
         magnitude *= cos(Pi*(y[i] - center_moving_lower_source_1[1])/epsDelta) + 1;
         magnitude /= 4.*epsDelta*epsDelta;

         coeffs[i][3] = magnitude; // reaction
         coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

         //double penalty_factor = 1000.;
         //coeffs[i][3] = 1 * penalty_factor; // reaction
         //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
   #else
         //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
         //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
         //magnitude /= 4.*epsDelta*epsDelta;
         //coeffs[i][3] += magnitude; // reaction
         //coeffs[i][4] -= magnitude * T_in; // source
         /*  double penalty_factor = 1000.;
           coeffs[i][4] = 1 * penalty_factor; // reaction
           coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
          */
         double magnitude = cos(Pi*(x[i] - center_moving_lower_source_1[0])/epsDelta) + 1;
         magnitude *= cos(Pi*(y[i] - center_moving_lower_source_1[1])/epsDelta) + 1;
         //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
         ////magnitude *= 1 + 1;
         magnitude /= 4.*epsDelta*epsDelta;

         coeffs[i][4] = magnitude; // reaction
         coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
   #endif
       }
    if(at_moving_source_lower_2)
          {
      #ifdef __2D__
            double magnitude = cos(Pi*(x[i] - center_moving_lower_source_2[0])/epsDelta) + 1;
            magnitude *= cos(Pi*(y[i] - center_moving_lower_source_2[1])/epsDelta) + 1;
            magnitude /= 4.*epsDelta*epsDelta;

            coeffs[i][3] = magnitude; // reaction
            coeffs[i][4] = magnitude * temperature_injection_well; // right-hand side

            //double penalty_factor = 1000.;
            //coeffs[i][3] = 1 * penalty_factor; // reaction
            //coeffs[i][4] = T_in * penalty_factor; // right-hand side*/
      #else
            //double magnitude = cos(Pi*x_distance_to_source/epsDelta) + 1;
            //magnitude *= cos(Pi*y_distance_to_source/epsDelta) + 1;
            //magnitude /= 4.*epsDelta*epsDelta;
            //coeffs[i][3] += magnitude; // reaction
            //coeffs[i][4] -= magnitude * T_in; // source
            /*  double penalty_factor = 1000.;
              coeffs[i][4] = 1 * penalty_factor; // reaction
              coeffs[i][5] = temperature_injection_well * penalty_factor; // right-hand side
             */
            double magnitude = cos(Pi*(x[i] - center_moving_lower_source_2[0])/epsDelta) + 1;
            magnitude *= cos(Pi*(y[i] - center_moving_lower_source_2[1])/epsDelta) + 1;
            //Old 05.02.19  magnitude *= cos(Pi*(z[i] - center_source[2])/epsDelta) + 1;
            ////magnitude *= 1 + 1;
            magnitude /= 4.*epsDelta*epsDelta;

            coeffs[i][4] = magnitude; // reaction
            coeffs[i][5] = magnitude * temperature_injection_well; // right-hand side
      #endif
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
  double center_x_moving_doublet = x[0];

  double center_x_moving_doublet_bottom_row = x[0];
  double center_x_moving_doublet_top_row = x[1];

  double height_of_bottom_row = 2000.;
  double height_of_top_row = 4000.;

  auto u1 = convection.GetComponent(0);
  auto u2 = convection.GetComponent(1);

#ifdef __3D__
  auto u3 = convection.GetComponent(2);
  std::array<FEFunction*, d> fe_functions_pointers{{u1, u2, u3}};
#else
  std::array<FEFunction*, d> fe_functions_pointers{{u1, u2}};
  #endif
  
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

    /*3D ???
        LocalAssembling<d> la(this->db, LocalAssembling_type::TCDStiffRhs,
                             fe_functions_pointers.data(),
                             this->example.get_coeffs(), 1); // db["discretization_type"]);
    */


    //call assembling, including mass matrix (SUPG!)
    // OLD: call_assembling_routine(s, la_a_rhs, la_m_supg , true);

    using namespace std::placeholders;

    if (this->db["scenario"].is("3_rows_of_double_doublets_varying_row_distance") )
    {
      height_of_bottom_row = x[2];
      height_of_top_row = x[3];

      la.ResetCoeffFct(std::bind(temperature_coefficients_3doubledoublets_4controls,_1, _2, _3, _4, _5,
  #ifdef __3D__
              _6,
  #endif
              center_x_moving_doublet_bottom_row, center_x_moving_doublet_top_row, nu,
              (double) this->db["transversal_dispersion_factor"], (double) this->db["longitudinal_dispersion_factor"],
              (double) this->db["well_radius"], (double) this->db["temperature_injection_well"],
              (double) this->db["delta_fct_eps_factor"], (double) this->db["fluid_density"],
              (double) this->db["fluid_heat_capacity"], (double) this->db["well_distance"],
              height_of_bottom_row, height_of_top_row));
    }
    else if (this->db["scenario"].is("3_rows_of_double_doublets_fixed_well_distance") )
    {
      la.ResetCoeffFct(std::bind(temperature_coefficients_3doubledoublets,_1, _2, _3, _4, _5,
  #ifdef __3D__
              _6,
  #endif
              center_x_moving_doublet_bottom_row, center_x_moving_doublet_top_row, nu,
              (double) this->db["transversal_dispersion_factor"], (double) this->db["longitudinal_dispersion_factor"],
              (double) this->db["well_radius"], (double) this->db["temperature_injection_well"],
              (double) this->db["delta_fct_eps_factor"], (double) this->db["fluid_density"],
              (double) this->db["fluid_heat_capacity"], (double) this->db["well_distance"] ));
    }
    else if (this->db["scenario"].is("2doublets_fixed_well_distance") )
    {
      la.ResetCoeffFct(std::bind(temperature_coefficients_2doublets,_1, _2, _3, _4, _5,
  #ifdef __3D__
              _6,
  #endif
              center_x_moving_doublet, nu,
              (double) this->db["transversal_dispersion_factor"], (double) this->db["longitudinal_dispersion_factor"],
              (double) this->db["well_radius"], (double) this->db["temperature_injection_well"],
              (double) this->db["delta_fct_eps_factor"], (double) this->db["fluid_density"],
              (double) this->db["fluid_heat_capacity"], (double) this->db["well_distance"] ));
    }
    else  if (this->db["scenario"].is("1doublet_optimize_distance") )
    {
    la.ResetCoeffFct(std::bind(temperature_coefficients,_1, _2, _3, _4, _5,
#ifdef __3D__
            _6,
#endif
            distance, nu,
            (double) this->db["transversal_dispersion_factor"], (double) this->db["longitudinal_dispersion_factor"],
            (double) this->db["well_radius"], (double) this->db["temperature_injection_well"],
            (double) this->db["delta_fct_eps_factor"], (double) this->db["fluid_density"],
            (double) this->db["fluid_heat_capacity"]));
    }
    else  if (this->db["scenario"].is("lattice") )
    {
    la.ResetCoeffFct(std::bind(temperature_coefficients_lattice,_1, _2, _3, _4, _5,
#ifdef __3D__
            _6,
#endif
            distance, nu,
            (double) this->db["transversal_dispersion_factor"], (double) this->db["longitudinal_dispersion_factor"],
            (double) this->db["well_radius"], (double) this->db["temperature_injection_well"],
            (double) this->db["delta_fct_eps_factor"], (double) this->db["fluid_density"],
            (double) this->db["fluid_heat_capacity"]));
    }
    else  if (this->db["scenario"].is("hexagon") )
    {
    la.ResetCoeffFct(std::bind(temperature_coefficients_hexagon,_1, _2, _3, _4, _5,
#ifdef __3D__
            _6,
#endif
            distance, nu,
            (double) this->db["transversal_dispersion_factor"], (double) this->db["longitudinal_dispersion_factor"],
            (double) this->db["well_radius"], (double) this->db["temperature_injection_well"],
            (double) this->db["delta_fct_eps_factor"], (double) this->db["fluid_density"],
            (double) this->db["fluid_heat_capacity"]));
    }
    
    la.setBeginParameter({0});
    la.SetN_Parameters(d);
    la.setN_ParamFct(1);
    la.setParameterFct({mapping_local_parameters<d>});
    la.setN_FeValues(d);

#ifdef __2D__
    la.setFeValueFctIndex({0, 1});
    la.setFeValueMultiIndex({D00, D00});
#else
    la.setFeValueFctIndex({0, 1, 2});
    la.setFeValueMultiIndex({D000, D000, D000});
#endif

    auto fe_space =  s.fe_space.get(); //this->get_space();
    double * rhs_entries = s.rhs.get_entries();
    auto * boundary_conditions = fe_space->get_boundary_condition();
    BoundaryValuesFunction * non_const_bound_value[1] {this->example.get_bd()[0]};
    auto blocks_stiff = s.stiffness_matrix.get_blocks_uniquely();

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

auto matrix_stiff = reinterpret_cast<SquareMatrixD*>(blocks_stiff.at(0).get());


#ifdef __3D__
 Assemble3D(
#else
 Assemble2D(
#endif
         1, &fe_space, 1, &matrix_stiff, 0, nullptr, 1, &rhs_entries,
        &fe_space, &boundary_conditions, non_const_bound_value, la);



/* 3D ??
     s.rhs.reset();
     s.stiffness_matrix.reset();
 */

  delete u1;
  delete u2;

#ifdef __3D__
     delete u3;
#endif

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
#ifdef __2D__
    this->outputWriter = DataWriter2D(this->db);
  ///OLD this->timeDependentOutput = DataWriter2D(this->db);
  FEFunction & fe_function = this->systems.front().fe_function;
#else
    this->outputWriter = DataWriter3D(this->db);
  TFEFunction3D & fe_function = this->systems.front().fe_function;
#endif

  this->outputWriter.add_fe_function(&fe_function);
  auto& s = this->systems.front();

  //reload the initial condition as fe_function for the next optimization loop
  s.fe_function.Interpolate(this->example.get_initial_cond(0));
  s.stiffness_matrix.reset();
  s.mass_matrix.reset();
  s.rhs.reset();
}


/* ************************************************************************* */
/*template <int d>
void TCD_Temperature<d>::output(int i)
{
  SystemPerGrid &s = this->systems.front();
  s.fe_function.PrintMinMax();
  bool i_am_root = true;
#ifdef _MPI
  int my_rank;
  int root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == root);
  // computing errors as well as writing vtk files requires a minimum
  // consistency level of 1
  s.fe_space->get_communicator().consistency_update(s.solution.get_entries(), 1);
#endif // _MPI

  //write solution for visualization
  outputWriter.write(time_stepping_scheme.current_time_);

  // compute errors
  if(db["output_compute_errors"])
  {
    const int n_errors = 4;
    std::array<double, n_errors+1> locError;
#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D all_derivatives[4] = { D000, D100, D010, D001 };
#else
    TAuxParam2D aux;
    MultiIndex2D all_derivatives[3] = { D00, D10, D01 };
#endif
    const FESpace* space = s.fe_space.get();

    s.fe_function.GetErrors(example.get_exact(0), d+1, all_derivatives,
                            n_errors, conv_diff_l2_h1_linf_error<d>,
                            example.get_coeffs(), &aux, 1, &space,
                            locError.data());
#ifdef _MPI
    /// @todo the GetErrors method in TFEFunction3D should already do the
    /// communication, it's surprising that in the mpi case the errors are
    /// squared, while the square root has been taken already in the sequential
    /// case.
    // global (across all processes) error, L2, H1, Linf, SD-error (for SUPG)
    std::vector<double> errorsReduced(n_errors);
    MPI_Reduce(locError.data(), errorsReduced.data(), n_errors-1, MPI_DOUBLE,
               MPI_SUM, root, MPI_COMM_WORLD);
    MPI_Reduce(&locError[n_errors-1], &errorsReduced[n_errors-1], 1,
               MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
    // correct values only on the root process!
    for(int i = 0; i < n_errors-1; ++i)
      locError[i] = std::sqrt(errorsReduced[i]);
    locError[n_errors-1] = errorsReduced[n_errors-1];
#endif // _MPI

    if(i_am_root)
    {
      Output::print<1>("time   : ", time_stepping_scheme.current_time_);
      Output::print<1>("  L2   : ", std::setprecision(14), locError[0]);
      Output::print<1>("  H1   : ", std::setprecision(14), locError[1]);
      Output::print<1>("  L_inf: ", std::setprecision(14), locError.at(3));
    }
    double tau = time_stepping_scheme.get_step_length();

    // errors[1] is the previous L2 error squared
    errors[0] += (locError[0] * locError[0] + errors[1] * errors[1]) * tau*0.5;
    errors[1] = locError[0];
    // errors[3] is the previous H1 error squared
    errors[2] += (locError[1] * locError[1] + errors[3] * errors[3]) * tau*0.5;
    errors[3] = locError[1];

    if(errors.at(4) < locError[3])
      errors[4] = locError[3];

    if(i_am_root)
    {
      Output::print<1>("  L2(0,T;L2)      : ", sqrt(errors[0]));
      Output::print<1>("  L2(0,T;H1)      : ", sqrt(errors[2]));
      Output::print<1>("  L_inf(0,T,L_inf): " , errors[4]);
    }
  }
}
*/





#ifdef __3D__
template class TCD_Temperature<3>;
#else
template class TCD_Temperature<2>;
#endif
