//
// Created by Moritz Hoffmann on 05/07/15.
//

#ifndef CDERRORESTIMATOR2D_H
#define CDERRORESTIMATOR2D_H


#include <ErrorEstimator2D.h>
#include <Example_CD2D.h>

enum class CDErrorEstimatorType {
    // 0 - gradient indicator
            GradientIndicator = 0,
    // 1 - H^1 estimator
            H1_ResidualEstimator,
    // 2 - L^2 estimator
            L2_ResidualEstimator,
    // 3 - energy norm + dual norm, Verf"urth 2005
            Energy_ResidualEstimatorQuasiRobust,
    // 4 - energy norm estimator without jumps
            Energy_ResidualEstimatorWithoutJumps,
    // 5 - supg estimator John/Novo, upper estimate
            SUPG_Upper,
    // 6 - supg estimator John/Novo, lower estimate
            SUPG_Lower
};

#define N_CD2D_ESTIMATOR_TYPES 7

// enable output of CD error estimator type
std::ostream &operator<<(std::ostream &os, CDErrorEstimatorType &type);

class CDErrorEstimator2D : public ErrorEstimator2D {
protected:
    // needed derivatives
    const std::vector<MultiIndex2D> derivatives{D10, D01, D00, D20, D02};
    // the selected estimator
    CDErrorEstimatorType estimatorType;

    // opaque struct holding relevant data of the edges
    struct EdgeData;
    // opaque struct holding relevant data of the ref edges in jump calculation
    struct EdgeRefData;

    // member variable indicating if the grid is conform or not
    int conform_grid;

    // internal function calculating jumps across edges at the boundary
    bool handleJump_BoundaryEdge(double *result, Example2D &example2D, const int estimatorType, const int N_QuadraturePoints1D, double *const &weights1D, const CDErrorEstimator2D::EdgeData &edgeData, BoundCond &Cond0, const double meas,
                                 const double *coeff, double linfb, const std::vector<double> &alpha, int edgeIdx, const TJoint *joint) const;

    // estimate!
    void estimate(const std::vector<MultiIndex2D> &derivatives, const TFEFunction2D &fe_function2D);

    // calculates eta_K for a single cell K
    double calculateEtaK(TBaseCell *cell, const TFEFunction2D &fe_function, TCollection &coll, std::vector<double *> &derivativesPerQuadPoint,
                         double (&AbsDetjk)[MaxN_QuadPoints_2D], double *(&weights), std::vector<double *> &coefficientsPerQuadPoint,
                         int n_quadrature_points, const unsigned int N_QuadraturePoints1D, double *(&weights1D), const EdgeData &edgeData, const double* zeta, EdgeRefData &edgeRefData) const;

public:
    // constructor
    CDErrorEstimator2D(Example_CD2D &ex, TCollection &collection, int type);

    void estimate(const TFEFunction2D &fe_function2D) {
        estimate(derivatives, fe_function2D);
    }


    int isConformGrid() const;

    void setConformGrid(int conform_grid);

    ~CDErrorEstimator2D() {};

};


#endif //CDERRORESTIMATOR2D_H
