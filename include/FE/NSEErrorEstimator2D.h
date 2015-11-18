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
    // parameter indicating whether the grid is conforming or not
    int conform_grid;
    // needed derivatives
    const std::vector<MultiIndex2D> derivatives_u {D00, D01, D10, D02, D20};
    const std::vector<MultiIndex2D> derivatives_p {D00, D01, D10};
    // estimation method
    void estimate(const std::vector<MultiIndex2D> &derivatives, const TFEFunction2D &fe_function2D) {};

    void calculateEtaK(const TFESpace2D &fe_space_u, const TFESpace2D &fe_space_p, TBaseCell *cell,
                            int N_Points, double *X, double *Y, double *AbsDetjk, double *weights, double **Derivatives, double **coeffs,
                            Example2D &example,
                            EdgeData &edgeData,
                            int *global_numbers_u, int *begin_index_u, int *DOF_u, double *values_u,
                            int *global_numbers_p, int *begin_index_p, int *DOF_p, double *values_p,
                            double *estimated_local_error);
public:
    NSEErrorEstimator2D(Example2D &ex, TDomain &domain, int type);

    void estimate(TFEVectFunct2D &fe_function2D_u, TFEFunction2D &fe_function2D_p, TAuxParam2D &Aux);

};

#endif //NSEERRORESTIMATOR2D_H