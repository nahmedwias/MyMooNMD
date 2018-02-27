/*
 * WiedmeyerBatchCrystallizer.h
 *
 *  Created on: 27 Feb 2018
 *      Author: bartsch
 */

#ifndef INCLUDE_EXAMPLES_WIEDMEYERBATCHCRYSTALLIZER_H_
#define INCLUDE_EXAMPLES_WIEDMEYERBATCHCRYSTALLIZER_H_

#include <functional>
#include <string>
#include <vector>

#include <parmoon_source_and_sink_terms.h>

namespace wiedmeyer_batch_crystallizer
{

namespace FluidProperties
{
double get_r_in();
double get_eta();
double get_rho();
double get_u_max_in();
void set_r_out(double new_r_out);
void set_z_out(double new_z_out);
void set_mass_flow_rate(double mfr);
void set_out_condition(const std::string& cond);
std::string out_condition_string();
}

namespace TemperatureConditions
{
void set_T_start(double new_T_start);
void set_T_end(double new_T_end);
void set_t_start(double new_t_start);
void set_t_end(double new_t_end);
}

namespace ConcentrationProperties
{
void subtract_material(double delta_c);
void update_index();
void initialise_inlet_values(double c_start, double t_cycle, double delta_t);
}

double derived_concentration_PAL_SUPSAT_POWG(const std::vector<double>& data);

namespace BrushInfo
{
size_t get_parameter_spatial_dimension();
size_t get_parameter_n_specs_primary();
size_t get_parameter_n_specs_derived();
std::vector<std::string> get_parameter_term_names();
std::vector<std::function<double(const std::vector<double>&)>> get_parameter_specs_derived_fcts();
std::vector<Exmpl::SourceAndSinkTerms> get_source_and_sink_fct_requests();
std::vector<std::string> get_source_and_sink_term_names();
}


void ExampleFile();
enum class BoundaryPart;
BoundaryPart determine_boundary_part(double x, double y, double z);
void BoundCondition(double x, double y, double z, BoundCond &cond);
void U1BoundValue(double x, double y, double z, double &value);
void U2BoundValue(double x, double y, double z, double &value);
void U3BoundValue(double x, double y, double z, double &value);
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs);
void ExactU1(double x, double y,  double z, double *values);
void ExactU2(double x, double y,  double z, double *values);
void ExactU3(double x, double y,  double z, double *values);
void ExactP(double x, double y,  double z, double *values);
void InitialU1(double x, double y,  double z, double *values);
void InitialU2(double x, double y,  double z, double *values);
void InitialU3(double x, double y,  double z, double *values);
void InitialP(double x, double y,  double z, double *values);
void Exact_cALUM(double x, double y, double z, double *values);
void BoundCondition_cALUM(double x, double y, double z, BoundCond &cond);
void BoundValue_cALUM(double x, double y, double z, double &value);
void InitialCondition_cALUM(double x, double y, double z, double *values);
void BilinearCoeffs_cALUM(int n_points, double *X, double *Y, double *Z,
                          double **parameters, double **coeffs);
void Exact_T(double x, double y, double z, double *values);
void BoundCondition_T(double x, double y, double z, BoundCond &cond);
void BoundValue_T(double x, double y, double z, double &value);
void InitialCondition_T(double x, double y, double z, double *values);
void BilinearCoeffs_T(int n_points, double *X, double *Y, double *Z,
                      double **parameters, double **coeffs);
}

#endif /* INCLUDE_EXAMPLES_WIEDMEYERBATCHCRYSTALLIZER_H_ */
