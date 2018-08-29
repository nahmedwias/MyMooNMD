#ifndef TIMEBOUNDARYCONTROLLEDOPTIMIZATION_H
#define TIMEBOUNDARYCONTROLLEDOPTIMIZATION_H

#include "ParameterDatabase.h"
#include "Time_NSE2D.h"
#include "Time_NSE2D_Adjoint.hpp"
#include "LoopInfo.h"


class Time_BoundaryControlledOptimization
{
  public:
    Time_BoundaryControlledOptimization(const TDomain& domain,
                                   const ParameterDatabase& param_db);
    
    /// @brief compute the functional \f$ \hat J \f$ and, if necessary, its 
    ///        gradient
    /// 
    /// This includes a primal solve with the given control 'x'. If the 
    /// gradient is required, an adjoint equation is solved.
    double compute_functional_and_derivative(unsigned n, const double *x, 
                                             double *grad);
    /// @brief return the dimension of the control space
    unsigned get_n_control() const { return n_control; }
    
    static ParameterDatabase default_BCO_database();
    
  protected:
    /// @brief keep all parameters for this optimization in one database
    ParameterDatabase db;
    /// @brief the size (dimension) of the control space
    unsigned n_control;
    /// @brief the dofs which are to be controlled.
    /// These are the dofs of the space for each component, it is therefore
    /// n_control = 2*control_dofs.size();
    std::vector<int> control_dofs;
    /// @brief the Navier--Stokes object representing the primal solve
    Time_NSE2D tnse_primal;
    /// @brief the Navier--Stokes object representing the adjoint solve
    Time_NSE2D_Adjoint tnse_adjoint;
    /// @brief Stokes solution as a 'desired state'
    std::shared_ptr<BlockVector> stokes_fe_vector;
    std::shared_ptr<TFEVectFunct2D> stokes_sol;
    
    
    /// variables during the optimization loop
    /// @brief keeping track of the optimization loop:
    LoopInfo optimization_info;
    /// @brief keeping track of the nonlinear loops:
    LoopInfo nonlinear_info;
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
    
    void apply_control_in_rhs_and_sol(const double * x);

    /// @brief compute \f$ \hat J \f$ using the (primal) solution and control
    double compute_functional() const;
    
    /// @brief using the (primal) solution, compute the adjoint
    void solve_adjoint_equation();
    
    /// @brief compute \f$ \hat J' \f$ using the adjoint solution and control
    void compute_derivative(const double * x, double* grad) const;
};

#endif // TIMEBOUNDARYCONTROLLEDOPTIMIZATION_H
