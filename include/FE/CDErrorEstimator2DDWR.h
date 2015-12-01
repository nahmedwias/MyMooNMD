// =============================================================================
// 
// Purpose: 
// 
// Author: mho, 27/11/15.
//
// =============================================================================

#ifndef CDERRORESTIMATOR2DDWR_H
#define CDERRORESTIMATOR2DDWR_H

#include <CDErrorEstimator2D.h>
#include <functional>
#include <Domain.h>
#include <CD2D.h>

class CDErrorEstimator2DDWR : public CDErrorEstimator2D {
protected:
    // the functional of interest
    typedef std::function<double(const TFEFunction2D*, double, double, TBaseCell&)> CD2DDwrFunctional;
    //J(solution, x, y, cell) = value
    CD2DDwrFunctional dwrFunctional;
    // the domain
    TDomain domain;

    void estimate(const std::vector<MultiIndex2D> &derivatives, const TFEFunction2D &fe_function2D);
public:

    CDErrorEstimator2DDWR(Example2D &ex, CD2DDwrFunctional &functional, TDomain &domain);

    void estimate(const TFEFunction2D &tfeFunction2D) {
        CDErrorEstimator2D::estimate(tfeFunction2D);
    };

    void manipulate(const TFEFunction2D *sol, int N_Points, double* X, double* Y, double **Coeff, double **parameters, TBaseCell *cell);
};

class CD2DDual : public CD2D {
public:
    CD2DDual(const TDomain& domain, Example_CD2D example, int reference_id = -4711);
    void assemble(const ManipulateAssemblingFct2D &manipulate);
    void output_dual(size_t level);
    void solve_transposed();
};

#endif //CDERRORESTIMATOR2DDWR_H
