#include "ParameterDatabase.h"
#include "Example_NSE2D.h"
#ifdef __2D__
#include "Example_TimeCD2D.h"
#else
#include "Example_TimeCD3D.h"
#endif

constexpr double surrounding_temperature = 348.15; //=75 + 273.15; //150.;

#ifdef __2D__
// we consider an approximated exact solution for the 2D case
void doublet_ux_solution(double x, double y, double *values)
{

  double r_0 = 50.; // epsilon used for approximate delta function
  double r_1 = 3000.; // minimal distance well-boundary
  double sigma = 10000.;
  double Qin = 150./3600.;
  double u0 = Qin/(2.*Pi*r_0);
  double xi = 4500.;
  double yi = 3000.;
  double xe = 5500.;
  double ye = 3000.;
  
  double x_1 = x - xi;
  double y_1 = y - yi;
  double r2_1 = x_1*x_1 + y_1*y_1;
  values[0] = u0 * r_0 * x_1/r2_1; 
  values[1] = u0 * r_0 * ( (-x_1*x_1+y_1*y_1)/(r2_1 * r2_1) );
  values[2] = u0 * r_0 * ( -2*x_1*y_1/(r2_1 * r2_1) );

  double x_2 = x - xe;
  double y_2 = y - ye;
  double r2_2 = x_2*x_2 + y_2*y_2;
  values[0] -= u0 * r_0 * x_2/r2_2; 
  values[1] -= u0 * r_0 * ( (-x_2*x_2+y_2*y_2)/(r2_2*r2_2) );
  values[2] -= u0 * r_0 * ( -2*x_2*y_2/(r2_2*r2_2) );

  values[3] = 0.;
}

void doublet_uy_solution(double x, double y, double *values)
{
  double r_0 = 50.; // epsilon used for approximate delta function
  double r_1 = 3000.; // minimal distance well-boundary
  double sigma = 10000.;
  double Qin = 150./3600.;
  double u0 = Qin/(2.*Pi*r_0);
  double xi = 4500.;
  double yi = 3000.;
  double xe = 5500.;
  double ye = 3000.;
  
  double x_1 = x - xi;
  double y_1 = y - yi;
  double r2_1 = x_1*x_1 + y_1*y_1;
  values[0] = u0 * r_0 * y_1/r2_1; 
  values[1] = u0 * r_0 * (-2*y_1*x_1/(r2_1*r2_1));
  values[2] = u0 * r_0 * (-y_1*y_1 + x_1*x_1)/( r2_1*r2_1);
    
  double x_2 = x - xe;
  double y_2 = y - ye;
  double r2_2 = x_2*x_2 + y_2*y_2;
  values[0] = u0 * r_0 * y_2/r2_2; 
  values[1] = u0 * r_0 * (-2*y_2*x_2/(r2_2*r2_2));
  values[2] = u0 * r_0 * (-y_2*y_2 + x_2*x_2)/( r2_2*r2_2);

  values[3] = 0.;
}

void doublet_p_solution(double x, double y, double *values)
{
  double r_0 = 50.; // epsilon used for approximate delta function
  double r_1 = 3000.; // minimal distance well-boundary
  double sigma = 10000.;
  double Qin = 150./3600.;
  double u0 = Qin/(2.*Pi*r_0);
  double xi = 4500.;
  double yi = 3000.;
  double xe = 5500.;
  double ye = 3000.;
  
  double x_1 = x - xi;
  double y_1 = y - yi;
  double r2_1 = x_1*x_1 + y_1*y_1;
  values[0] = -sigma * u0 * r_0 * 0.5 * log( r2_1/(r_1*r_1) );     
  values[1] = -sigma * u0 * r_0 * x_1/r2_1;
  values[2] = -sigma * u0 * r_0 * y_1/r2_1;

    
  double x_2 = x - xe;
  double y_2 = y - ye;
  double r2_2 = x_2*x_2 + y_2*y_2;
  values[0] -= -sigma * u0 * r_0 * 0.5 * log( r2_2/(r_1*r_1) );     
  values[1] -= -sigma * u0 * r_0 * x_2/r2_2;
  values[2] -= -sigma * u0 * r_0 * y_2/r2_2;

  values[3] = 0.;
}
#endif

void unknown_solution(double, double, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
#ifdef __3D__
  values[4] = 0.;
#endif
}

void all_Dirichlet_boundary_condition(
#ifdef __2D__
				      int , double , BoundCond &cond
#else
				      double, double, double , BoundCond &cond
#endif
				      )
{
  cond = DIRICHLET;
}

void all_Neumann_boundary_condition(
#ifdef __2D__
				    int , double , BoundCond &cond
#else
				    double, double, double , BoundCond &cond
#endif
				    )
{
  cond = NEUMANN;
}

void zero_boundary_value(
#ifdef __2D__
			 int , double , double& value
#else
			 double, double, double , double& value
#endif
			 )
{
  value = 0;
}

void temperature_boundary_value(
#ifdef __2D__
				int , double , double& value
#else
				double, double, double , double& value
#endif
				)
{
  value = surrounding_temperature;
}

void initial_condition_temperature(double x, double y,
#ifdef __3D__
				   double z,
#endif
				   double *values)
{
  values[0] = surrounding_temperature;
}


void pde_coefficients_flow(int n_points, double *x, double *y,
#ifdef __3D__
			   double *z,
#endif
			   double ** parameters,
			   double **coeffs, double nu, double sigma,
			   bool use_parameters = false)
{
  for(int i = 0; i < n_points; i++)
  {
    int dimension;
#ifdef __2D__
    dimension = 2;
#else
    dimension = 3;
#endif
    // physical parameters
    coeffs[i][0] = nu;
    coeffs[i][dimension+2] = sigma;
    // momentum source
    for (int k=1; k<=dimension; k++)
      coeffs[i][k] = 0.; // f_k

    //divergence soucr
    coeffs[i][dimension+1] = 0.; // divergence

    /*
    ///@todo implement space dependent parameters
    if (!use_parameters) // initial assembling // if (parameters[i] == nullptr)
    {
      if (abs(x[i] - 1.5 - y[i] + 0.01) < 0.05) // crack through the sink
      {
	coeffs[i][4] = sigma; //1e10;
	//cout<< "erste Gerade" <<endl;
      }
      else if (abs(x[i] - 2.5 - y[i] - 0.01) < 0.05) // crack through the source
      {
	coeffs[i][4] = sigma; //1e10;
      }
      else
	coeffs[i][4] = sigma;
    }
    else
    {
      coeffs[i][4] = parameters[i][0];
      
    }
    */
    
  }

  
}


void pde_coefficients_temperature(int n_points, double *, double *,
#ifdef __3D__
				  double *,
#endif
                                  double **parameters, double **coeffs,
                                  double nu)
{
  int dim=2;
#ifdef __3D__
  dim=3;
#endif
  for (int i = 0; i < n_points; i++)
  {
    // physical parameters
    coeffs[i][0] = nu;
    if(parameters[i] == nullptr) // initial assembling
    {
      coeffs[i][1] = 0.; // convection, x-direction
      coeffs[i][2] = 0.; // convection, y-direction
#ifdef __3D__
      coeffs[i][3] = 0.; // convection, z-direction
#endif
    }
    else
    {
      coeffs[i][1] = parameters[i][0]; // convection, x-direction
      coeffs[i][2] = parameters[i][1]; // convection, y-direction
#ifdef __3D__
      coeffs[i][3] = parameters[i][2]; // convection, z-direction
#endif
    }
    coeffs[i][dim+1] = 0.; // reaction
    coeffs[i][dim+2] = 0.; // f
  }
}

#ifdef __2D__
Example_NSE2D get_gppo_flow_example(const ParameterDatabase & db)
{
  int example_code = db["example"];
  Output::print<1>(" Example code: ", example_code);
  switch (example_code) {
  case 1:
  {
    //std::vector<DoubleFunct2D *> exact(3, unknown_solution);
    std::vector<DoubleFunct2D *> exact;
    exact.push_back(doublet_ux_solution);
    exact.push_back(doublet_uy_solution);
    exact.push_back(doublet_p_solution);
    
    std::vector<BoundCondFunct2D *> bc{{all_Neumann_boundary_condition,
	  all_Neumann_boundary_condition, all_Neumann_boundary_condition}};
    std::vector<BoundValueFunct2D *> bd(3, zero_boundary_value); 
    double reynolds_number = db["reynolds_number"];
    double sigma = db["inverse_permeability"];
    using namespace std::placeholders;
    bool use_coeff_fct = false; // db["variable_sigma_fct_type"];
    CoeffFct2D coeffs = std::bind(pde_coefficients_flow, _1, _2, _3, _4, _5,
				  1./reynolds_number, sigma, use_coeff_fct);
    return Example_NSE2D(exact, bc, bd, coeffs, 1./reynolds_number);
    break;
  }
  default:
  {
    Output::print(" ** ERROR, example ", example_code , " not implemented");
    exit(1);
  }
  }
}

Example_TimeCD2D get_gppo_temperature_example(const ParameterDatabase & db)
{
  //int example = db["example"];
  std::vector<DoubleFunct2D *> exact(1, unknown_solution);
  std::vector<BoundCondFunct2D *> bc(1, all_Dirichlet_boundary_condition);
  std::vector<BoundValueFunct2D *> bd(1, temperature_boundary_value);
  std::vector <DoubleFunct2D*> ic(1, initial_condition_temperature);
  double nu = db["diffusion_coefficient"];
  using namespace std::placeholders;
  CoeffFct2D coeffs = std::bind(pde_coefficients_temperature, _1, _2, _3, _4, _5, nu);
  return Example_TimeCD2D(exact, bc, bd, coeffs, false, false, ic);
}
 
#else

 Example_NSE3D get_gppo_flow_example(const ParameterDatabase & db)
{
  //int example = db["example"];
  std::vector<DoubleFunct3D *> exact(4, unknown_solution_3D);
  std::vector<BoundCondFunct3D *> bc{{all_Dirichlet_boundary_condition, all_Dirichlet_boundary_condition, 
                                      all_Dirichlet_boundary_condition, all_Neumann_boundary_condition}};
  std::vector<BoundValueFunct3D *> bd_3D(4, zero_boundary_value_3D); 
  
  double reynolds_number = db["reynolds_number"];
  double sigma = db["inverse_permeability"];
  using namespace std::placeholders;
  CoeffFct3D coeffs = std::bind(pde_coefficients_flow, _1, _2, _3, _4, _5, _6, 1./reynolds_number, sigma);
  
  //cout << " **********INSIDE 3D_gppo_flow_example************"<< endl;
  return Example_NSE3D(exact, bc, bd_3D, coeffs, 1./reynolds_number);
}


Example_TimeCD3D get_gppo_temperature_example(const ParameterDatabase & db)
{
  //int example = db["example"];
  std::vector<DoubleFunct3D *> exact(1, unknown_solution_3D);
  std::vector<BoundCondFunct3D *> bc(1, all_Dirichlet_boundary_condition);
  std::vector<BoundValueFunct3D *> bd(1, temperature_boundary_value_3D);
  std::vector <DoubleFunct3D*> ic(1, initial_condition_temperature_3D);
  
  double nu = db["diffusion_coefficient"];
  using namespace std::placeholders;
  CoeffFct3D coeffs = std::bind(pde_coefficients_temperature, _1, _2, _3, _4, _5, _6, nu);
  //cout << " **********INSIDE 3D_gppo_temperature_example************"<< endl;
  return Example_TimeCD3D(exact, bc, bd, coeffs, false, false, ic);
}
#endif
