#ifndef __SYSTEM_NAVIERSTOKES_ADJOINT_H__
#define __SYSTEM_NAVIERSTOKES_ADJOINT_H__

#include "NavierStokes.h"

template <int d>
class NavierStokes_Adjoint : public NavierStokes<d>
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FEVectFunct = typename Template_names<d>::FEVectFunct;
    using DoubleFunction = typename Template_names<d>::DoubleFunction;
    using BoundaryValuesFunction = typename Template_names<d>::BoundaryValuesFunction;
    using FESpace = typename Template_names<d>::FESpace;
    
    NavierStokes_Adjoint(const NavierStokes<d>& nse, const ParameterDatabase& param_db);
  
    /// @brief assemble all terms in the matrix and right-hand side which
    /// contain the primal solutions u or p.
    void assemble(const FEVectFunct& u, const FEFunction& p,
                  const FEVectFunct& stokes_u, 
                  std::vector<double> cost_functional_weights, 
                  bool restricted_curl_functional);
    
    /// @brief solve the system.
    /// This is essentially a call to NSE2D::solve, however afterwards scaling
    /// on the diagonal in the Dirichlet rows is reversed
    void solve();
};

#endif // __SYSTEM_NAVIERSTOKES_ADJOINT_H__
