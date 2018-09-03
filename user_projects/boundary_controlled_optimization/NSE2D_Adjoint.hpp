#ifndef NSE2D_ADJOINT_H
#define NSE2D_ADJOINT_H

#include "NSE2D.h"

class NSE2D_Adjoint : public NSE2D
{
  public:
    NSE2D_Adjoint(const NSE2D& nse2d, const ParameterDatabase& param_db);
  
    /// @brief assemble all terms in the matrix and right-hand side which
    /// contain the primal solutions u or p.
    void assemble(const TFEVectFunct2D& u, const TFEFunction2D& p,
                  const TFEVectFunct2D& stokes_u, 
                  std::vector<double> cost_functional_weights, 
                  bool restricted_curl_functional);
    
    /// @brief solve the system.
    /// This is essentially a call to NSE2D::solve, however afterwards scaling
    /// on the diagonal in the Dirichlet rows is reversed
    void solve();
};

#endif // NSE2D_ADJOINT_H
