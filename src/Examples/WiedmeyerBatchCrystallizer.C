/**
 * @file Batch crystallizer example - complete.
 * Experiment was reported by Viktoria Wiedmeyer (VW) in 2017.
 */

#include <Database.h>
#include <MooNMD_Io.h>
#include <WiedmeyerBatchCrystallizer.h>

#include <cmath>

namespace wiedmeyer_batch_crystallizer
{
enum class OutCondition{DO_NOTHING, CONSTANT, PARABOLIC};

OutCondition out_condition;

namespace FluidProperties
{

double u_infty = 1;    // (m/s) the characteristic velocity of the fluid
double l_infty = 1;    // (m) the characteristic length scale of the tube

double r_in = 0.01;   //m the radius of the inlet
double r_out = 0.075; //m the radius of the outlet
double z_out = 0.5;   //m the distance from inflow to outflow

double mass_flow_rate = 0; // (kg/s) the mass flow rate at in- and outflow
double u_avg_in  = 0;      //m/s Those quantitites are computed according
double u_max_in  = 0;      //m/s to the entered mass flow, assuming HP-inflow.
double u_avg_out  = 0;     //m/s

double M_ALUM = 0.4743884;       // kg/mol, the molar mass of potash alum dodecahydrate
double M_A = 0.258205;           // kg/mol, the molar mass of potash alum anhydrate
double beta_A = 0.54429029;      // 0.258205 kg/mol / 0.4743884 kg/mol - the mass fraction
// of KAL(SO_4)_2 (anhydrate) in the potash alum dodecahydrate

double eta = 0.0014;     // ( kg /(m*s) ) the dynamic viscosity, (of a Kalialaun solution, VW)
double rho = 1050; // ( kg / m^3  ) the density (of a Kalialaun solution, VW) (assumed to be constant)

double get_r_in(){return r_in;};
double get_eta(){return eta;};
double get_rho(){return rho;};
double get_u_max_in(){return u_max_in;};

void set_r_out(double new_r_out)
{
  r_out = new_r_out;
}

void set_z_out(double new_z_out)
{
  z_out = new_z_out;
}

void set_mass_flow_rate(double mfr)
{
  mass_flow_rate = mfr/3600; // (kg/s) the mass flow rate at in- and outflow
  u_avg_in = mass_flow_rate / (rho * M_PI * r_in * r_in ); //0.047157;   //m/s
  u_max_in = 2*u_avg_in; //m/s, assuming HP-inflow
  u_avg_out = mass_flow_rate / (rho * M_PI * r_out*r_out ); //m/s
}

void set_out_condition(const std::string& cond)
{
  if(cond == "do_nothing")
    out_condition = OutCondition::DO_NOTHING;
  else if(cond == "constant")
    out_condition = OutCondition::CONSTANT;
  else if(cond == "parabolic")
    out_condition = OutCondition::PARABOLIC;
  else
    ErrThrow("Unrecognized outflow option. Choose between 'do_nothing',"
        "'constant' and 'parabolic'.");
}

std::string out_condition_string()
{
  if(out_condition == OutCondition::DO_NOTHING)
    return "do_nothing";
  if(out_condition == OutCondition::CONSTANT)
    return "constant";
  if(out_condition == OutCondition::PARABOLIC)
    return "parabolic";
  else
    return "terribly wrong";
}

// note: in the coefficients function the de-dimensionalized diffusion
// coefficient will be calculated as:
//      eps = (eta/rho) / (u_infty*l_infty);
}

namespace TemperatureConditions
{
double T_start;
double T_end;
double t_start;
double t_end;
void set_T_start(double new_T_start){T_start = new_T_start;}
void set_T_end(double new_T_end){T_end = new_T_end;}
void set_t_start(double new_t_start){t_start = new_t_start;}
void set_t_end(double new_t_end){t_end = new_t_end;}

double lambda = 0.6;  // W/(m  * K), thermal  conductivity  (of the fluid), from Volker's kali alaun paper stub
double c_p = 3841;    // J/(kg * K), specific heat capacity (of the fluid), from Volker's kali alaun paper stub
}

namespace ConcentrationProperties
{
double diffusion_ALUM = 5.4e-10; // m^2/s, from Volker's kali alaun paper stub
double M_ALUManhydrate = 0.258192; // [kg/mol]
double cALUM_start = 0;  //mol anhydrate/m^3 This will be set from the outside.

std::vector<double> inlet_values;
size_t current_index; // The place in the vector 'inlet_values' were
                      //   a) the current inflow is read from and
                      //   b) the current use of concentration is substracted from

double get_M_ALUManhydrate(){return M_ALUManhydrate;}

void subtract_material(double delta_c)
{
  inlet_values.at(current_index) -= delta_c;
//  Output::print("new value at ",current_index," is ", std::setprecision(16),
//                inlet_values.at(current_index),". (Subtracted ",delta_c,").");
}

void update_index()
{
  ++current_index;
  current_index = current_index % inlet_values.size();
}

void initialise_inlet_values(double c_start, double t_cycle, double delta_t)
{
  cALUM_start = c_start;
  int n_values = t_cycle/delta_t;
  inlet_values = std::vector<double>(n_values, cALUM_start);
  current_index = 0;

  Output::info("INLET CONCENTRATIONS LIST", "cALUM_start was set to ",
               cALUM_start, " mol anhydrate/m^3.");
  Output::info("INLET CONCENTRATIONS LIST",
               "Initialized ", n_values, " values , cycle time is ",
               t_cycle, " s and delta_t is ", delta_t," s.");
}
}


// This computes the "concentration" of alum supersaturation, which is needed
// by Brush in order to compute the growth rate. The growth model we use here
// was communicated by V.Wiedmeyer, and it is formulated in terms of the
// relative mass-fraction supersaturation. The growth model comes from
//    Temmel et.al.(2016): A Short-Cut Method for the Quantification of
//    Crystallization Kinetics. 1. Method Development. Crystal Growth &
//    Design, 16 (12), 6743-6755,
// and it is formulated in a different supersaturation measure, i.e.,
// in kg hydrate / kg added solvent.
// To be precise, the model requires the term:
//    (S-1)^g
// with exponent g = 1.4 and the supersaturation measure
//    S = ... TODO describe the reformulation due to different concentration measures here
// w_A [kg/kg]is the mass fraction of anhydrate
// w_eq,A [kg/kg] is the equilibrium anhydrate mass fraction, which is given by
// a fitted quadratic polynomial (in temperature T) here.
// This dimensionless number is computed by this function, and then handed over
// to Brush by the name "PAL_SUPSAT_POWG".

double derived_concentration_PAL_SUPSAT_POWG(const std::vector<double>& data)
{
  if(data.size() != 7)
    throw std::runtime_error("derived_concentration_PAL_SUPSAT_POWG: "
        "expected 7 data points. ux, uy, uz, p, T, POTASHALUM, 0(PAL_SUPSAT_POWG)");

  double T = data[4] - 273.15;     // grab temperature, convert to degree C
  double c_A = data[5];   // grab alum anhydrate concentration [mol/m^3]

  double m_A = c_A * FluidProperties::M_A; // anhydrate, [kg/m^3]
  double m_W = FluidProperties::rho - m_A; // water    , [kg/m^3]

  //Just hack in Viktorias formula.
  double w_A = m_A / (FluidProperties::beta_A * m_W - m_A * (1 - FluidProperties::beta_A ));
  // Temmel model for solubility.
  double A = 5.06;
  double B = 0.23;
  double C = 7.76e-3;
  double D = -2.43e-4;
  double E = 4.86e-6;
  double w_A_eq = (A + B * T + C * T * T
      + D * T *T * T + E * T * T * T * T) / 100 ;

  double rel_sup_sat = (w_A / w_A_eq) - 1;

  // return the supersaturation to the power of 1.4, as the model requires, or 0 if there is no supersaturation
  return rel_sup_sat > 0 ? pow(rel_sup_sat, 1.4) : 0;

}

namespace BrushInfo
{
//! The spatial dimension is obviously 3.
size_t parameter_spatial_dimension = 3;
//! There are two "species" to be communicated: temperature and dissolved potash alum conc
size_t parameter_n_specs_primary = 2;
// Potash alum saturation concentration "to the power of 2.1" is a derived quantity.
size_t parameter_n_specs_derived = 1;
// Names of the parameters to be sent to Brush, just for double-checking.
std::vector<std::string> parameter_term_names = {"ux","uy","uz","p","T","POTASHALUM","PAL_SUPSAT_POWG"};
// Function to compute potash alum saturation concentration "to the power of 2.1" is a derived quantity.
std::vector<std::function<double(const std::vector<double>&)>> parameter_specs_derived_fcts =
{derived_concentration_PAL_SUPSAT_POWG};
//! Names of the functions that are requested from Brush.
std::vector<Exmpl::SourceAndSinkTerms> source_and_sink_fct_requests =
{ Exmpl::SourceAndSinkTerms::PotashAlumCrystEnergyRelease,
  Exmpl::SourceAndSinkTerms::PotashAlumCrystConcConsumption };
//! Names of the functions that are expected from Brush in return.
std::vector<std::string> source_and_sink_term_names = {"T_sources", "c_POTASHALUM_sinks" };

/// The getters.
size_t get_parameter_spatial_dimension(){return parameter_spatial_dimension;};
size_t get_parameter_n_specs_primary(){return parameter_n_specs_primary;};
size_t get_parameter_n_specs_derived(){return parameter_n_specs_derived;};
std::vector<std::string> get_parameter_term_names(){return parameter_term_names;};
std::vector<std::function<double(const std::vector<double>&)>> get_parameter_specs_derived_fcts(){return parameter_specs_derived_fcts;};
std::vector<Exmpl::SourceAndSinkTerms> get_source_and_sink_fct_requests(){return source_and_sink_fct_requests;};
std::vector<std::string> get_source_and_sink_term_names(){return source_and_sink_term_names;};
}


void ExampleFile()
{
  Output::info<1>("EXAMPLE","WiedmeyerBatchCrystallizer.h");
  Output::info<1>("EXAMPLE","With mass flow rate ", FluidProperties::mass_flow_rate * 3600, " kg/h.");
  Output::info<1>("EXAMPLE","and outflow condition '", FluidProperties::out_condition_string() ,"'.");
}

//Boundary parts of the Geometry, numbered "bottom up"
enum class BoundaryPart {BOTTOM, WALL, TOP};

//Find out which boundary part a boundary point lies on. Measures are in m.
BoundaryPart determine_boundary_part(double x, double y, double z)
{
  double tol = 1e-8;
  //catch first those values that lie on the wall parts of top and bottom
  if(    (fabs(sqrt(x*x + y*y) - FluidProperties::r_in) < tol && fabs(z-0) < tol)      //bottom wall
      || (fabs(sqrt(x*x + y*y) - FluidProperties::r_out) < tol && fabs(z-FluidProperties::z_out) < tol) ) //top wall
  {
    return BoundaryPart::WALL;
  }
  if (fabs(z-0) < tol )
  {
    return BoundaryPart::BOTTOM;
  }
  else if (fabs(z-FluidProperties::z_out) < tol )
  {
    return BoundaryPart::TOP;
  }
  else //everything "complicated" is on the wall boundary
    return BoundaryPart::WALL;
}

// kind of boundary condition
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  if (determine_boundary_part(x,y,z) == BoundaryPart::BOTTOM)
    cond = DIRICHLET;
  else if (determine_boundary_part(x,y,z) == BoundaryPart::WALL)
    cond = DIRICHLET;
  else //TOP
    cond = out_condition == OutCondition::DO_NOTHING ? NEUMANN : DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  using namespace FluidProperties;
  if (determine_boundary_part(x,y,z) == BoundaryPart::BOTTOM) //inflow
  {//HP inflow profile
    double R = 0.01;
    double r = sqrt(x*x + y*y);

    value = u_max_in * (1 - (r*r)/(R*R));
    double t = TDatabase::TimeDB->CURRENTTIME;
    if(t < 0.1) //within first second of the simulated time
      //multiply inflow with t ("anstroemen")
      value  *= (10*t);
  }
  else if (determine_boundary_part(x,y,z) == BoundaryPart::WALL)	//no-slip
  {
    value = 0;
  }
  else //TOP - outflow
  {
    if(out_condition == OutCondition::DO_NOTHING)
    {//do-nothing outflow condition
      value = 0;
    }
    else if(out_condition == OutCondition::CONSTANT)
    {//constant outflow condition
      value = FluidProperties::u_avg_out;
      double t = TDatabase::TimeDB->CURRENTTIME;
      if(t < 0.1) //within first second of the simulated time
        //multiply outflow with t ("anstroemen")
        value  *= (10*t);
    }
    else if(out_condition == OutCondition::PARABOLIC)
    {//parabolic outflow profile
      double R = r_out;
      double r = sqrt(x*x + y*y);
      value = (2*u_avg_out) * (1 - (r*r)/(R*R));
      double t = TDatabase::TimeDB->CURRENTTIME;
      if(t < 0.1) //within first second of the simulated time
        //multiply outflow with t ("anstroemen")
        value  *= (10*t);
    }
  }
}

// ========================================================================
// coefficients for Stokes form
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  using namespace FluidProperties;

  static double eps = (eta/rho) / (u_infty*l_infty);
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = -9.81; // f3 - gravity forcing
    coeff[4] = 0; // g
  }
}

// ========================================================================
// exact solution is unknown, everything set to 0.
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void InitialU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// ALUM mass balance equation
// ========================================================================
// exact solution is unknown
void Exact_cALUM(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_cALUM(double x, double y, double z, BoundCond &cond)
{
  if (determine_boundary_part(x,y,z) == BoundaryPart::BOTTOM)
  {
    cond=DIRICHLET;
  }
  else if (determine_boundary_part(x,y,z) == BoundaryPart::WALL)
  {
    cond=NEUMANN;
  }
  else if (determine_boundary_part(x,y,z) == BoundaryPart::TOP)
  {
    cond=NEUMANN;
  }
}

// value of boundary condition
void BoundValue_cALUM(double x, double y, double z, double &value)
{
  using namespace FluidProperties;
  if (determine_boundary_part(x,y,z) == BoundaryPart::BOTTOM)
  {
    value=ConcentrationProperties::inlet_values.at(ConcentrationProperties::current_index);
  }
  else if (determine_boundary_part(x,y,z) == BoundaryPart::WALL)
  {
    value=0;
  }
  else if (determine_boundary_part(x,y,z) == BoundaryPart::TOP)
  {
    value=0;
  }
}

// initial conditon
void InitialCondition_cALUM(double x, double y, double z, double *values)
{
  values[0] = ConcentrationProperties::cALUM_start;
}

void BilinearCoeffs_cALUM(int n_points, double *X, double *Y, double *Z,
                          double **parameters, double **coeffs)
{
  double *coeff;

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = ConcentrationProperties::diffusion_ALUM; //diffusion coefficient
    coeff[1] = parameters[i][0];   // ux
    coeff[2] = parameters[i][1];   // uy
    coeff[3] = parameters[i][2];   // uz
    coeff[4] = 0;                  // no reaction in this system
    coeff[5] = parameters[i][3];   //rhs, interpolated sources and sinks from Brush
  }
}

// ========================================================================
// (Thermal) Energy balance equation
// ========================================================================
// exact solution is unknown
void Exact_T(double x, double y, double z, double *values)
{
  //  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_T(double x, double y, double z, BoundCond &cond)
{
  if (   determine_boundary_part(x,y,z) == BoundaryPart::BOTTOM
      || determine_boundary_part(x,y,z) == BoundaryPart::WALL)
  {
    cond=DIRICHLET;
  }
  else if (determine_boundary_part(x,y,z) == BoundaryPart::TOP)
  {
    cond=NEUMANN;
  }
}

// value of boundary condition
void BoundValue_T(double x, double y, double z, double &value)
{
  using namespace TemperatureConditions;
  if (   determine_boundary_part(x,y,z) == BoundaryPart::BOTTOM
      || determine_boundary_part(x,y,z) == BoundaryPart::WALL)
  {
    double t = TDatabase::TimeDB->CURRENTTIME;
    //current 'value' is a linear interpolation of start and end temperature
    value = (T_start * (t_end - t) + T_end * (t- t_start)) / (t_end - t_start);
  }
  else if (determine_boundary_part(x,y,z) == BoundaryPart::TOP)
  {
    value=0;
  }
}

// initial conditon
void InitialCondition_T(double x, double y, double z, double *values)
{
  using namespace TemperatureConditions;
  values[0] = T_start;
}

void BilinearCoeffs_T(int n_points, double *X, double *Y, double *Z,
                      double **parameters, double **coeffs)
{
  using namespace TemperatureConditions;
  double *coeff;
  //  double t=TDatabase::TimeDB->CURRENTTIME;

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = TemperatureConditions::lambda/(FluidProperties::rho * TemperatureConditions::c_p);
    coeff[1] = parameters[i][0];   // ux
    coeff[2] = parameters[i][1];   // uy
    coeff[3] = parameters[i][2];   // uz
    coeff[4] = 0;                  // no reaction in this system
    coeff[5] = parameters[i][3];   //rhs, interpolated sources and sinks from Brush
  }
}
}
