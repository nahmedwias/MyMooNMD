#ifndef GEOTHERMALPLANTSPOSITIONOPTIMIZATION_H
#define GEOTHERMALPLANTSPOSITIONOPTIMIZATION_H

#include "Domain.h"
#include "ParameterDatabase.h"
#include "TCD_Temperature.hpp"
#include "LoopInfo.h"

#include "NSE_GPPO.hpp"


class ParameterDatabase;

/** ************************************************************************ */
template <int d>
class GeothermalPlantsPositionOptimization
{

public:
  constexpr static char required_database_name_TCD_GPPO[] = "TCD parameter database";

#ifdef __2D__
  GeothermalPlantsPositionOptimization(const TDomain& domain,
          const ParameterDatabase& param_db);
#else
  GeothermalPlantsPositionOptimization( TDomain& domain,
                                       const ParameterDatabase& param_db);
#endif

  /// @brief compute the functional \f$ \hat J \f$ and, if necessary, its
  /// gradient.
  /// This includes a primal solve with the given control 'x'. If the
  /// gradient is required, an adjoint equation is solved.
  double compute_functional_and_derivative(unsigned n, const double *x,
          double *grad);
  /// @brief return the dimension of the control space
  unsigned get_n_control() const { return n_control; }

  static ParameterDatabase default_GPPO_database();

  ParameterDatabase get_primal_flow_database(ParameterDatabase param_db);
  ParameterDatabase get_primal_temperature_database(ParameterDatabase param_db);

protected:

  /// @brief keep all parameters for this optimization in one database
  ParameterDatabase db;
  /// @brief the size (dimension) of the control space
  unsigned n_control;

#ifdef __2D__
  /// @brief the Brinkman and TCD2D  objects representing the primal solve
  NSE_GPPO<d> brinkman2d_primal;
  TCD_Temperature<d> tcd2d_primal;
#else
  /// @brief the Brinkman and CD2D  objects representing the primal solve
  NSE_GPPO<d> brinkman3d_primal;
  TCD_Temperature<d> tcd3d_primal;

#endif

  /// @brief the Brinkman2D and TCD2D objects representing the adjoint solve
  //Brinkman2D_Adjoint brinkman2d_adjoint;
  //TCD2D_Adjoint tcd2d_adjoint;

  /// variables during the optimization loop
  /// @brief keeping track of the optimization loop:
  LoopInfo optimization_info;
  /// @brief value of current functional \f$ \hat J \f$
  double current_J_hat;
  /// @brief the control from the previous iteration
  std::vector<double> control_old;
  /// @brief counting the calls to 'compute_functional_and_derivative'
  unsigned n_calls;
  /// @brief counting how often derivatives are computed
  unsigned n_computation_derivative;

  /// @brief apply the control and solve the (primal) pde
  ///
  /// This essentially sets copies the values to the right-hand side
  void apply_control_and_solve(const double * x);

  /// @brief compute \f$ \hat J \f$ using the (primal) solution and control
  double compute_functional(const double * x) const;

  /// @brief using the (primal) solution, compute the adjoint
  void solve_adjoint_equation();

  /// @brief compute \f$ \hat J' \f$ using the adjoint solution and control
  void compute_derivative(const double * x, double* grad) const;

};

#endif // GEOTHERMALPLANTSPOSITIONOPTIMIZATION_H
