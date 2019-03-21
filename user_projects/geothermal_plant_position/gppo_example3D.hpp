#include "ParameterDatabase.h"
#include "Example_NSE3D.h"
#include "Example_TimeCD3D.h"
//#include "Example3D.h"

constexpr double surrounding_temperature_3D = 150.;

void unknown_solution_3D(double, double, double, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
  values[4] = 0.;
}

void all_Dirichlet_boundary_condition(double, double, double , BoundCond &cond)
{
  cond = DIRICHLET;
}

void all_Neumann_boundary_condition(double, double, double , BoundCond &cond)
{
  cond = NEUMANN;
}

void zero_boundary_value_3D(double, double, double, double &value)
{
  value = 0.;
}

void temperature_boundary_value_3D(double, double, double, double& value)
{
  value = surrounding_temperature_3D; 
}

void initial_condition_temperature_3D(double x, double y, double z, double *values)
{
  values[0] = surrounding_temperature_3D;
}

void pde_coefficients_flow(int n_points, double *, double *, double *, double **,
                           double **coeffs, double nu, double sigma)
{
  for(int i = 0; i < n_points; i++)
  {
    // physical parameters
    coeffs[i][0] = nu;
    coeffs[i][1] = 0.; // f1
    coeffs[i][2] = 0.; // f2
    coeffs[i][3] = 0.; // f3
    coeffs[i][4] = 0.; // divergence
    coeffs[i][5] = sigma;
  }
}
void pde_coefficients_temperature(int n_points, double *, double *, double *, 
                                  double **parameters, double **coeffs,
                                  double nu)
{
  for(int i = 0; i < n_points; i++)
  {
    // physical parameters
    coeffs[i][0] = nu;
   if(parameters[i] == nullptr) // initial assembling
    {
      coeffs[i][1] = 0.; // convection, x-direction
      coeffs[i][2] = 0.; // convection, y-direction
      coeffs[i][3] = 0.; // convection, z-direction
    }
     else
    {
      coeffs[i][1] = parameters[i][0]; // convection, x-direction
      coeffs[i][2] = parameters[i][1]; // convection, y-direction
      coeffs[i][3] = parameters[i][2]; // convection, z-direction
    }
    coeffs[i][4] = 0.; // reaction
    coeffs[i][5] = 0.; // f
  }
}

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
