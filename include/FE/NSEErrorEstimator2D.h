#ifndef NSEERRORESTIMATOR2D_H
#define NSEERRORESTIMATOR2D_H

#include <ErrorEstimator2D.h>
#include <Example_NSE2D.h>

enum class NSE2DErrorEstimatorType {
    // 0 - gradient indicator
    gradient_indicator = 0,
    // 1 - residual estimator h1
    residual_estimator_h1,
    // 2 - residual estimator l2
    residual_estimator_l2,
    // 3 - residual estimator energy norm quasi robust
    residual_estimator_energy_quasi_robust,
    // 4 - gradient recovery
    gradient_recovery,
    // 5 - implicit estimator neumann
    implicit_estimator_neumann
};

#define N_NSE2D_ESTIMATOR_TYPES 6

// enable output of NSE error estimator type
std::ostream &operator<<(std::ostream &os, NSE2DErrorEstimatorType &type);

class NSEErrorEstimator2D : public ErrorEstimator2D {
protected:
    // the selected estimator type
    NSE2DErrorEstimatorType estimatorType;
    struct EdgeData;
    struct EdgeRefData;
    // parameter indicating whether the grid is conforming or not
    int conform_grid;
    // indicating if navier stokes or not
    bool is_nse;
    // needed derivatives
    const std::vector<MultiIndex2D> derivatives_u {D00, D01, D10, D02, D20};
    const std::vector<MultiIndex2D> derivatives_p {D00, D01, D10};
    // estimation method
    void estimate(const std::vector<MultiIndex2D> &derivatives, const TFEFunction2D &fe_function2D) {};

    void calculateEtaK(TFEVectFunct2D &fe_function2D_u, TFEFunction2D &fe_function2D_p, TBaseCell *cell, unsigned int N_Points, unsigned int N_Points1D, double *AbsDetjk, double *weights, double **Derivatives,
                           double **coeffs, Example2D &example, EdgeData &edgeData, EdgeRefData &edgeRefData, int *global_numbers_u, int *begin_index_u, double *values_u, int *global_numbers_p, int *begin_index_p,
                           double *values_p, double *estimated_local_error);
    std::vector<double> getWeights(const double hK, const double delta, const double *coeff);
    unsigned int get_max_n_base_functions(const TFESpace2D &fe_space);
public:
    NSEErrorEstimator2D(Example2D &ex);
    NSEErrorEstimator2D(Example2D &ex, int type, bool is_nse);

    void estimate(TFEVectFunct2D &fe_function2D_u, TFEFunction2D &fe_function2D_p, TAuxParam2D &Aux);

    NSE2DErrorEstimatorType GetEstimatorType() const {
        return estimatorType;
    }
};

#endif //NSEERRORESTIMATOR2D_H