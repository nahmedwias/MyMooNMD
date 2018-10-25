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
    
    static bool check_input_parameter_consistency(const ParameterDatabase& db);

  protected:
    /// @brief keep all parameters for this optimization in one database
    ParameterDatabase db;
    /// @brief the size (dimension) of the control space
    unsigned n_control;
    /// @brief the dofs which are to be controlled.
    /// These are the dofs of the space for each component, it is therefore
    /// n_control = 2*control_dofs.size() multiplied by the number of time steps
    std::vector<int> control_dofs;
    /// @brief the Navier--Stokes object representing the primal solve
    Time_NSE2D tnse_primal;
    /// @brief the Navier--Stokes object representing the adjoint solve
    Time_NSE2D_Adjoint tnse_adjoint;
    /// @brief Stokes solution as a 'desired state'
    std::shared_ptr<BlockVector> stokes_fe_vector;
    std::shared_ptr<TFEVectFunct2D> stokes_sol;

    /// @brief Number of time steps of the time-dependent primal (and adjoint) NSE
    int n_time_steps_;
    /// @brief a struct which stores one solution of tnse_primal at a given step t
    struct tnse_primal_solution{
        BlockVector vector_at_timestep_t_;
        int timestep_t_;
        TFEVectFunct2D u_at_timestep_t_;
        TFEFunction2D p_at_timestep_t_;

        /** @brief default constructor*/
        tnse_primal_solution();

        /** @brief constructor*/
        tnse_primal_solution(const BlockVector& solution_vector, int step,
                             const TFEVectFunct2D& u_sol,const TFEFunction2D& p_sol);
    };
    
    /// @brief a collection which stores all the solutions of tnse_primal in time
     std::deque<tnse_primal_solution> tnse_primal_solutions_;

    
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
    
    void impose_control_in_rhs_and_sol(const double * x, int current_time_step);

    /// @brief compute \f$ \hat J \f$ using the (primal) solution and control
    double compute_functional_at_t(int time_step) const;
    
    /// @brief compute \int_t0^t1 \f$ \hat J \f$ using the (primal) solution and control
    double compute_functional_in_time(int t0, int t1) const;

    /// @brief using the (primal) solution, compute the adjoint
    void solve_adjoint_equation();
    
    /// @brief compute \f$ \hat J' \f$ using the adjoint solution and control
    void compute_derivative(const double * x, double* grad) const;

    /// @brief write solution vectors of all time steps into files (only for
    /// checking and debugging purposes)
    void write_all_solutions();

    ///
    ParameterDatabase get_primal_database(const ParameterDatabase& param_db);

    ///
    ParameterDatabase get_adjoint_database(const ParameterDatabase& param_db);
};

#endif // TIMEBOUNDARYCONTROLLEDOPTIMIZATION_H
