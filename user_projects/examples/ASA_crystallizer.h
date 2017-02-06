/**
 * ASA_crystallizer.h
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
 * @date 2016/06/08
 * @author Clemens Bartsch
 */

#ifndef USER_PROJECTS_EXAMPLES_ASA_CRYSTALLIZER_H_
#define USER_PROJECTS_EXAMPLES_ASA_CRYSTALLIZER_H_

// some integer constants used to identify boundary parts
int bdry_inflow = 3;
int bdry_upper = 0;
int bdry_outflow = 1;
int bdry_lower = 2;

// hard coded parameter sets for four different inflow velocities
//  - can be controlled via input database
int VELOCITY_CODE = 0;

double inflow_T[4] = {307.6, 312.9, 313.1, 313.7}; //unit is K
double inflow_c = 2005; // unit is mol/m^3
double surrounding_T = 297.5;

// physical parameters
namespace Physics
{
  double lambda = 0.1676; // W/(m * K), thermal conductivity (of Ethanol)
  double C_E = 2441.3; // J/(kg * K), specific heat capacity (of Ethanol)
  double rho_E = 790; // kg/m^3, density (of Ethanol)
  double delta_h_cryst = 1.6541e5; //fusion enthalpy [J/kg]
  // should be diffusion coefficient of ASA in ethanol, is actually some
  // diffusion of some kind of ASA in NaOH reported in
  // http://www.iiste.org/Journals/index.php/JHMN/article/viewFile/10040/10256
  double D = 1.69e-8; // m^2/s

  double M_Ethanol = 0.04607; //molar mass of Ethanol (kg/mol)
  double M_ASA = 0.18016; //molar mass of ASA (kg/mol)
}

//term responsible for the coupling on the right hand side in temperature equation
double couplingTerm_T( const double* const params ) //must contain all that is necessary: T,c,F (zeroth moment)
{
//  //term responsible for the coupling on the right hand side in temperature equation
//  double T = params[0];
//  double c = params[1];
//  double F = params[2]; //parameter functions takes care of that
//
//  double constant = Physics::delta_h_cryst / (Physics::rho_E * Physics::C_E);
//
//  return constant * F_growth(T,c,F); //watch out! sign differs from the way the equation was scetched
  return 0; // no more coupling - the coupling terms will be reactivated
            // if there is an actual reaction taking place in the fluid
}

double couplingTerm_C_ASA( const double* const params ) //must contain all that is necessary: T,c,F (zeroth moment)
{
//  double T = params[0];
//  double c = params[1];
//  double F = params[2]; //parameter functions takes care of that
//
//  double constant = 1 / Physics::M_ASA;
//
//  return - constant * F_growth(T,c,F); //watch out! sign differs from the way the equation was scetched
  return 0; // no more coupling - the coupling terms will be reactivated
            // if there is an actual reaction taking place in the fluid
}


//Just print some informations on the example.
void ExampleFile()
{
  Output::info("Example", "This is the ASA crystallizer coupling example, "
      "it is work in progress.");
}

// ///////////// Temperature uncoupled part ///////////// //
void BoundCond_T(int BdComp, double t, BoundCond &cond)
{
  if(BdComp == bdry_inflow)
    cond = DIRICHLET;
  else if(BdComp == bdry_outflow)
    cond = NEUMANN;
  else //lower and upper boundary
    cond = DIRICHLET;
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
  else //wall boundary
  {//TODO this modelling assumption is supposedly really bad.
    value = surrounding_T; // surrounding temperature
  }
}

void Coefficients_T(int n_points, double *x, double *y,
                      double **parameters, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = Physics::lambda / (Physics::rho_E * Physics::C_E); //diffusion coefficient
    coeffs[i][1] = parameters[i][0];//convection in x direction
    coeffs[i][2] = parameters[i][1];//convection in y direction
    coeffs[i][3] = 0; //no reaction.

    coeffs[i][4] = parameters[i][2]; //rhs, interpolated sources and sinks from Brush
  }
}

void InitialCondition_T(double x,  double y, double *values)
{
  values[0] = surrounding_T; //TODO must be replaced by precomputed steady-state temp field
}

// /////////// ASA Concentration uncoupled part /////////// //
void BoundCond_C_ASA(int BdComp, double t, BoundCond &cond)
{
  if ( BdComp == bdry_inflow )
    cond = DIRICHLET;
  else if ( BdComp == bdry_outflow )
    cond = NEUMANN;
  else //wall boundary
    cond = NEUMANN; //impermeability
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
  else //wall boundary
    value = 0; //impermeability
}

void Coefficients_C_ASA(int n_points, double *x, double *y,
                      double **parameters, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = Physics::D; //diffusion coefficient
    coeffs[i][1] =  parameters[i][0];//convection in x direction
    coeffs[i][2] =  parameters[i][1];//convection in y direction
    coeffs[i][3] = 0; //no reaction.

    coeffs[i][4] = parameters[i][2]; //rhs, interpolated sources and sinks from Brush
  }
}

void InitialCondition_C_ASA(double x,  double y, double *values)
{
  values[0] = 0;
}

/****************************************************************************************
 * ParameterFunction and AssemblingFunctions used in the "Linearized Decoupled" strategy.
 ****************************************************************************************/

void ParameterFunction(double* in, double* out){
  // Just skip the first two entries - this is where TAuxParam2D->GetParameters() places the x and y value,
  // and pass on as many values as nCoupled_ + number of further_functions.
  // - for this example that would be 2 + 1 = 3.
  for (size_t i = 0; i<3 ; ++i){
    out[i] = in [i+2];
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

////A helper function relevant in the coupling of concentration, temperature and psd
//// TODO Turn this into a function that translates ASA and T into ASASUP (supersaturation)
//double F_growth(double T, double c, double F)
//{
//  //theoretical supersaturation in terms of ASA mole fraction, as given by Eder et al.
//  double chi_sat = pow(10, 27.769 - (2500.906/T) - 8.323 * log10(T));
//
//  //now this has to be transformed to supersaturation in
//  // terms of molar concentration, which is a bit cumbersome
//  double M_S = chi_sat*Physics::M_ASA + (1-chi_sat)*Physics::M_Ethanol; //molar mass of solution
//  double w_Ethanol = 1 - chi_sat*(Physics::M_ASA/M_S); //Ethanol mass fraction
//  double rho_S = Physics::rho_E / w_Ethanol; //density of solution calculated under "ideal solution" asumption
//
//  //now here comes supersaturation in terms of molar concentration
//  double c_sat = chi_sat * rho_S / M_S;
//
//  // calculate the actual growth coefficient
//  double k_g = 1.2e-5; //some coefficient, chosen as in Carina's diss
//  double g = 1; //some exponent, chosen as in Carina's diss
//  double G = c > c_sat ? k_g * ( (c - c_sat) / c_sat) : 0; //^g - do this only if g is not 1!
//
//  double F_growth = - G * F;
//
//  return F_growth;
//}

#endif /* USER_PROJECTS_EXAMPLES_ASA_CRYSTALLIZER_H_ */
