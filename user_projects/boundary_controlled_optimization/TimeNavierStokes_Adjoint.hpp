#ifndef TIMENSE2D_ADJOINT_H
#define TIMENSE2D_ADJOINT_H

#include "TimeNavierStokes.h"

template<int d>
class TimeNavierStokes_Adjoint : public TimeNavierStokes<d>
{
  public:
    TimeNavierStokes_Adjoint(const TimeNavierStokes<d>& nse, const ParameterDatabase& param_db);
  
    /// @brief assemble all terms in the matrix and right-hand side which
    /// contain the primal solutions u or p.
    void assemble(const TFEVectFunct2D& u, const TFEFunction2D& p,
                  const TFEVectFunct2D& stokes_u, 
                  std::vector<double> cost_functional_weights, 
                  bool restricted_curl_functional);
    
    /// @brief solve the system.
    /// This is essentially a call to Time_NSE2D::solve, however afterwards scaling
    /// on the diagonal in the Dirichlet rows is reversed
    void solve();
};

#endif // TIMENSE2D_ADJOINT_H
