#ifndef __SYSTEM_TIMENAVIERSTOKES_ADJOINT_H__
#define __SYSTEM_TIMENAVIERSTOKES_ADJOINT_H__

#include "TimeNavierStokes.h"

template<int d>
class TimeNavierStokes_Adjoint : public TimeNavierStokes<d>
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FEVectFunct = typename Template_names<d>::FEVectFunct;
    using DoubleFunction = typename Template_names<d>::DoubleFunction;
    using BoundaryValuesFunction = typename Template_names<d>::BoundaryValuesFunction;
    using Example_TimeNSE = typename Template_names<d>::Example_TimeNSE;
    using FESpace = typename Template_names<d>::FESpace;
    
    TimeNavierStokes_Adjoint(const TimeNavierStokes<d>& nse,
                             const ParameterDatabase& param_db);
  
    /// @brief assemble all terms needed for the first time step (from t=T to 
    /// t=T-dt).
    void assemble_initial_time(const FEVectFunct& u, const FEFunction& p);
    
    /// @brief assemble all terms in the matrix and right-hand side which
    /// contain the primal solutions u or p.
    void assemble_matrices_rhs(const FEVectFunct& u,
                               const FEFunction& p);
    
    /// @brief solve the system.
    /// This is essentially a call to Time_NSE2D::solve, however afterwards scaling
    /// on the diagonal in the Dirichlet rows is reversed
    void solve();
};

#endif // __SYSTEM_TIMENAVIERSTOKES_ADJOINT_H__
