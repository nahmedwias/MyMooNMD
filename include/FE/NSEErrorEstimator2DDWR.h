// =============================================================================
//
// Purpose:
//
// Author: Moritz Hoffmann
//
// =============================================================================

#ifndef NSEERRORESTIMATOR2DDWR_H
#define NSEERRORESTIMATOR2DDWR_H

#include <NSEErrorEstimator2D.h>
#include <functional>
#include <Domain.h>
#include <NSE2D.h>

class NSEErrorEstimator2DDWR : public NSEErrorEstimator2D {
protected:
    // the functional of interest
    typedef std::function<std::vector<double>(const TFEVectFunct2D*, const TFEFunction2D*, double, double, TBaseCell&)> NSE2DDwrFunctional;
    //J(solution, x, y, cell) = value
    NSE2DDwrFunctional dwrFunctional;
    // the domain
    TDomain domain;

public:

    NSEErrorEstimator2DDWR(Example2D &ex, NSE2DDwrFunctional &functional, TDomain &domain);

    void estimate(const TFEVectFunct2D &u, const TFEFunction2D &p, TAuxParam2D &Aux);

    void manipulate(const TFEVectFunct2D *sol_u, const TFEFunction2D *sol_p, int N_Points, double* X, double* Y, double **Coeff, double **parameters, TBaseCell *cell);
};

class NSE2DDual : public NSE2D {
public:
    NSE2DDual(const TDomain& domain, Example_NSE2D example, unsigned int reference_id = -4711);
    void assemble(const ManipulateAssemblingFct2D &manipulate);
    void output_dual(size_t level);
    void solve_transposed();
};

#endif //NSEERRORESTIMATOR2DDWR_H
