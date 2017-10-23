#ifndef USER_PROJECTS_EXAMPLES_AXISYMMETRIC_ASA_CRYSTALLIZER_H_
#define USER_PROJECTS_EXAMPLES_AXISYMMETRIC_ASA_CRYSTALLIZER_H_

/**
 * Axisymmetric_ASA_Crystallizer.h
 *
 * A ParMooN example to be used with the Example_TimeCoupledCDR2D class.
 * Simulation of a coninuously operated ASA flow crystallizer, as reported
 * in:
 * Eder et al. (2010): Continuously seeded (...) tubular crystallizer for
 * the production of active pharmaceutical ingredients. Crystal Growth & Design 10,
 * pp. 2247--2257.
 *
 * Contains a coupling of temperature, ASA concentration and crystal population
 * balance, all subject to a flow in a long, thin tube.
 *
 * This is the axisymmetric version of the very example.
 *
 * @date 2017/06/21
 * @author Clemens Bartsch
 */

#include "vector"
#include "functional"
#include "string"

// some integer constants used to identify boundary parts
int bdry_wall = 0;
int bdry_outflow = 1;
int bdry_rotational = 2;
int bdry_inflow = 3;


// hard coded parameter sets for four different inflow velocities
//  - can be controlled via input database
int VELOCITY_CODE = 0;

double tube_length = 15; //length of the tube in m - must be in accordance to the geometry!

double inflow_T[4] = {307.6, 312.9, 313.1, 313.7}; //unit is K
double inflow_c = 1511.11; // unit is mol/m^3
double surrounding_T = 297.5;

// physical parameters
namespace Physics
{
  double lambda_E = 0.1676; 	// W/(m * K), thermal conductivity (of Ethanol)
  double C_E = 2441.3; 		// J/(kg * K), specific heat capacity (of Ethanol)
  double rho_E = 790; 		// kg/m^3, density (of Ethanol)
  double rho_ASA = 1350; 	// kg/m^3, density (of ASA)
  double rho_susp = 916.76; // kg/m^3, the assumed density of the suspension (see my modelling tables)
  double delta_h_cryst = 1.6541e5; //specific heat of fusion [J/kg]
  // should be diffusion coefficient of ASA in ethanol, is actually diffusion
  // coefficient of urea in EtOH as reported by Anker et al. 2015...
  double D = 1.35e-9; // m^2/s

  double M_Ethanol = 0.04607; //molar mass of Ethanol (kg/mol)
  double M_ASA = 0.18016; //molar mass of ASA (kg/mol)
}

namespace Temperature
{
double T_amb = 297.5;   // ambient temperature (room temperature)
double T_feed[4] = {307.6 , 312.9 , 313.1 , 313.7};  // feed stream temperature in K
double r_outer = 0.002; // 2mm
double r_inner = 0.001;  // 1mm

// calculate suspension heat capacity
double m_dot_asa [4] ={0.000058, 0.000088, 0.000116, 0.000128};  // kg / s
double c_p_ASA = 1260;      // J/(kg K)
double m_dot_EtOH [4] = { 0.000116, 0.000175, 0.000232, 0.000257}; // kg / s
double c_p_EtOH = 2400; // J/(kg K)
//FIXME I have the suspicion, that the following line takes ALWAYS  VELOCITY_CODE=0!
double m_tc_t = m_dot_asa[VELOCITY_CODE] * c_p_ASA + m_dot_EtOH[VELOCITY_CODE] * c_p_EtOH;

// calculate total heat transfer coefficient
double alpha_inner = 306; // W/(m^2 K)
double alpha_outer =  70; // W/(m^2 K)
double lambda_tube = 0.3; // W/(m K)
double k = 1 / (r_outer * ( 1/(r_inner*alpha_inner) + log(r_outer/r_inner)/lambda_tube + 1/(r_outer*alpha_outer)) ) ;
}

double temperature_bound_cond_exp( double x )
{
    using namespace Temperature;
    return T_amb + (T_feed[VELOCITY_CODE]-T_amb) * std::exp( -k*2*r_outer*M_PI/(m_tc_t) * x );
}

double couplingTerm_T( const double* const params )
{
  return 0;   // currently there is no direct reaction between fluid quantities -
  	  	  	  // they are only couled via growth, i.e. Brush takes care of it
}

double couplingTerm_C_ASA( const double* const params )
{
  return 0; // currently there is no direct reaction between fluid quantities -
  	  	  	// they are only couled via growth, i.e. Brush takes care of it
}


//Just print some informations on the example.
void ExampleFile()
{
  Output::info("Example", "ASA crystallizer coupling example in axisymmetric formulation.");
}

// ///////////// Temperature uncoupled part ///////////// //
void BoundCond_T(int BdComp, double t, BoundCond &cond)
{
	  if(BdComp == bdry_inflow)
	    cond = DIRICHLET;
	  else if(BdComp == bdry_outflow)
	    cond = NEUMANN;
	  else if(BdComp == bdry_wall) //wall boundary
	    cond = DIRICHLET;
	  else if(BdComp == bdry_rotational) //the spurious boundary, Neumann bdry (see theoretical work Ganesan & Tobiska 2008)
	    cond = NEUMANN;
}

void BoundValue_T(int BdComp, double Param, double &value)
{
  if(BdComp == bdry_inflow){
    value = inflow_T[VELOCITY_CODE]; //constant, depends on parameter set
  }
  else if(BdComp == bdry_outflow)
  {
    value = 0;
  }
  else if(BdComp == bdry_wall)//wall boundary
  {
    double x = Param*tube_length;
    value = temperature_bound_cond_exp(x);
  }
  else if(BdComp == bdry_rotational)
	  value = 0;

}

void Coefficients_T(int n_points, double *x, double *y,
                      double **parameters, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = Physics::lambda_E / (Physics::rho_susp * Physics::C_E); //heat diffusion coefficient
    coeffs[i][1] = parameters[i][1]; //convection in z direction
    coeffs[i][2] = parameters[i][2]; //convection in r direction
    coeffs[i][3] = 0; //no reaction.

    coeffs[i][4] = parameters[i][3]; //rhs, interpolated sources and sinks from Brush

    coeffs[i][5] = parameters[i][0]; //y=r coordinate
  }
}

void InitialCondition_T(double x,  double y, double *values)
{
  values[0] = surrounding_T;
}

// /////////// ASA Concentration uncoupled part /////////// //
void BoundCond_C_ASA(int BdComp, double t, BoundCond &cond)
{
  if ( BdComp == bdry_inflow )
    cond = DIRICHLET;
  else if ( BdComp == bdry_outflow )
    cond = NEUMANN;
  else //wall boundary and rotational (spurious) boundary
    cond = NEUMANN; //impermeability/natural condition
}

void BoundValue_C_ASA(int BdComp, double Param, double &value)
{
  if(BdComp == bdry_inflow){
    double t = TDatabase::TimeDB->CURRENTTIME;
    if (t < 0.1) //standard trick: linearly raise inflow
      value = 10 * t * inflow_c; //constant, depends on parameter set
    else
      value = inflow_c; //constant, depends on parameter set
  }
  else if ( BdComp == bdry_outflow )
    value = 0;
  else //wall boundary and rotational (spurious) boundary
    value = 0; //impermeability/natural condition
}

void Coefficients_C_ASA(int n_points, double *x, double *y,
                      double **parameters, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = Physics::D; //diffusion coefficient
    coeffs[i][1] =  parameters[i][1];//convection in z direction
    coeffs[i][2] =  parameters[i][2];//convection in r direction
    coeffs[i][3] = 0; //no reaction.

    coeffs[i][4] = parameters[i][3]; //rhs, interpolated sources and sinks from Brush

    coeffs[i][5] = parameters[i][0]; //y=r coordinate
  }
}

void InitialCondition_C_ASA(double x,  double y, double *values)
{
  values[0] = 0; //initially, there is a total absence of ASA concentration
}

/****************************************************************************************
 * ParameterFunction and AssemblingFunctions used in the "Linearized Decoupled" strategy.
 ****************************************************************************************/
// FIXME All the following are not actually needed here and just used as dummies,
// since there is no direct reactive coupling in this example.
void ParameterFunction(double* in, double* out){
  // Pass on as many values as 1 + nCoupled_ + number of further_functions.
  // - for this example that would be 1 + 2 + 1 = 4.
  // The first entry is the y coordinate (=r), needed in axisymmetric case due to integral trafo.
  for (size_t i = 0; i<4 ; ++i){
    out[i] = in [i+1];
  }
}

void RhsAssemblingFunction_C_ASA(
    double Mult, double *coeff, double *param,
    double hK, double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs)
{
  // For convenience: rename the place to write to (Rhs) and the place to read from (Orig,
  // which suposedly contains values of test functions)
  double* Rhs = LocRhs[0];
  double* Orig = OrigValues[0];

  // Calculate the value of the coupling term at the quad point.
  // This relies on the correct interaction with the ParameterFunction and the
  // input order of fe functions to the aux object.

  double coupledTerm = couplingTerm_C_ASA(param);

  //Loop over all local base functions.
  for(int i=0;i<N_BaseFuncts[0];i++)
  {
    Rhs[i] -= Mult*coupledTerm*Orig[i]; //(Mult contains quad weigth, and maybe even the determinant due to trafo)
  }
}

void RhsAssemblingFunction_T(
    double Mult, double *coeff, double *param,
    double hK, double **OrigValues, int *N_BaseFuncts,
    double ***LocMatrices, double **LocRhs)
{
  // For convenience: rename the place to write to (Rhs) and the place to read from (Orig,
  // which suposedly contains values of test functions)
  double* Rhs = LocRhs[0];
  double* Orig = OrigValues[0];

  // Calculate the value of the coupling term at the quad point.
  // This relies on the correct interaction with the ParameterFunction and the
  // input order of fe functions to the aux object.

  // coupled term evaluated at the quad point this
  // AssembleFctParam2D is called upon
  double coupledTerm = couplingTerm_T(param);

  //Loop over all local base functions.
  for(int i=0;i<N_BaseFuncts[0];i++)
  {
    Rhs[i] -= Mult*coupledTerm*Orig[i]; //(Mult contains quad weigth, and maybe even the determinant due to trafo)
  }
}

/* ****************************************************************************
 * This is example information which is of relevance to the Brush Wrapper that
 * must be used for this example.
 * ****************************************************************************/

/* SOURCE AND SINK information: information needed to control Brush to ParMooN*/

std::vector<Exmpl::SourceAndSinkTerms> source_and_sink_fct_requests =
  { Exmpl::SourceAndSinkTerms::ASACrystEnergyRelease,
    Exmpl::SourceAndSinkTerms::ASACrystConcConsumption };
std::vector<std::string> source_and_sink_term_names = {"T_sources", "c_ASA_sinks" };

/* PARAMETER information: information needed to control from ParMooN to Brush*/
size_t parameter_spatial_dimension = 2;
size_t parameter_n_specs_primary = 2; //temperature and dissolved ASA concentration
size_t parameter_n_specs_derived = 2; //EtOH and ASASUP

std::vector<std::string> parameter_term_names = {"uz","ur","p","T","ASA","CH3CH2OH","ASASUP"};

// With given vaues uz, ur, p, T and ASA conc in a certain point,
// will evaluate the molar concentration of EtOH, which is a derived quantity.
double derived_concentration_EtOH(const std::vector<double>& data)
{
  if (data.size() != 7)
    throw std::runtime_error("derived_concentration_EtOH: expected 7 data points."
        " ur, uz, p, T, ASA, 0(CH3CH2OH) , 0(ASASUP)");

  return 1; // TODO The concentration of EtOH is not needed in Brush,
            // yet it would be nice to have a correct value there.
}

// With given values uz, ur, p, T and ASA conc in a certain point,
// will evaluate the supersaturation concentration of ASA, which is a
// derived quantity.
// Temperature must be in K and ASA concentration in mol/m^3. Output
// will also be in mol/m^3.
double derived_concentration_ASASUP(const std::vector<double>& data)
{
  if (data.size() != 7)
    throw std::runtime_error("derived_concentration_EtOH: expected 7 data points."
        " uz, ur, p, T, ASA, CH3CH2OH , 0(ASASUP)");

    double T = data[3];     // grab temperature
    double c_asa = data[4]; // grab ASA concentration

    //theoretical supersaturation in terms of ASA mole fraction, as given by Eder et al.
    double chi_sat = pow(10, 27.769 - (2500.906/T) - 8.323 * log10(T));

    //now this has to be transformed to supersaturation in
    // terms of molar concentration, which is a bit cumbersome

    //molar mass of saturated solution
    double M_solsat = chi_sat*Physics::M_ASA + (1-chi_sat)*Physics::M_Ethanol;

    // Var I: constant density assumption, here applied to the saturated solution
    double rho_solsat = Physics::rho_susp;
    // Var II: would be ideal mixture assumption - but means that the calculation
    // of the supersaturation is the only place (besides to the computation of initial values)
    // where we do not assume constant density.
    // ... TODO

    //now here comes supersaturation in terms of molar concentration
    double c_sat = chi_sat * rho_solsat / M_solsat;
    double c_supsat = std::max(0.0, c_asa - c_sat);

//    if (c_supsat > 0)
//      Output::print("c_asa: ", c_asa, ", c_sat: ", c_sat , ", supsat: ", c_supsat, " mol/m^3");

    return c_supsat;
}

std::vector<std::function<double(const std::vector<double>&)>> parameter_specs_derived_fcts =
{ derived_concentration_EtOH, derived_concentration_ASASUP};

#endif /* USER_PROJECTS_EXAMPLES_AXISYMMETRIC_ASA_CRYSTALLIZER_H_ */
