#include "ParameterDatabase.h"
#include "Example_NSE2D.h"
#ifdef __2D__
#include "Example_TimeCD2D.h"
#else
#include "Example_TimeCD3D.h"
#endif

constexpr double surrounding_temperature = 348.15; //=75 + 273.15; //150.;

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
    // physical parameters
    coeffs[i][0] = nu;
    coeffs[i][1] = 0.; // f1
    coeffs[i][2] = 0.; // f2
#ifdef __2D__
    coeffs[i][3] = 0.; // divergence
#else
    coeffs[i][3] = 0.; // f3
    coeffs[i][4] = 0.; // div
#endif

#ifdef __2D__
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
#else
    coeffs[i][5] = sigma;
#endif
    
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
  //int example = db["example"];
  std::vector<DoubleFunct2D *> exact(3, unknown_solution);
  std::vector<BoundCondFunct2D *> bc{{all_Dirichlet_boundary_condition,
    all_Dirichlet_boundary_condition, all_Neumann_boundary_condition}};
  std::vector<BoundValueFunct2D *> bd(3, zero_boundary_value); 
  double reynolds_number = db["reynolds_number"];
  double sigma = db["inverse_permeability"];
  using namespace std::placeholders;
  bool use_coeff_fct = false; // db["variable_sigma_fct_type"];
  CoeffFct2D coeffs = std::bind(pde_coefficients_flow, _1, _2, _3, _4, _5,
                                1./reynolds_number, sigma, use_coeff_fct);
  return Example_NSE2D(exact, bc, bd, coeffs, 1./reynolds_number);
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
Example_NSE3D get_3D_gppo_flow_example(const ParameterDatabase & db)
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


Example_TimeCD3D get_3D_gppo_temperature_example(const ParameterDatabase & db)
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
