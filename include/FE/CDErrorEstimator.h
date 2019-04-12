//
// Created by Moritz Hoffmann on 05/07/15.
//

#ifndef CDERRORESTIMATOR2D_H
#define CDERRORESTIMATOR2D_H

#include "templateNames.h"
#include "ErrorEstimator.h"
#ifdef __2D__
#include "Example_CD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_CD3D.h"
#include "FEFunction3D.h"
#endif
#include "Joint.h"

enum class CDErrorEstimatorType
{
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

constexpr int N_CD2D_ESTIMATOR_TYPES = 7;

// enable output of CD error estimator type
std::ostream &operator<<(std::ostream &os, CDErrorEstimatorType type);

/**
 * @brief Error estimator specialization for convection-diffusion problems.
 * 
 * This (template) class essentially provides only one method 'estimate' whose
 * main purpose it is to fill the vector 'eta_K' along with the members 
 * 'maximal_local_error', 'estimated_global_error', and 'currentCollection' in
 * its base (template) class.
 * 
 * @note This is only a start, it is not thoroughly tested and currently only 
 *       works in 2D.
 */
template <int d>
class CDErrorEstimator : public ErrorEstimator<d>
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_CD = typename Template_names<d>::Example_CD;
    using MultiIndex_vector = typename Template_names<d>::MultiIndex_vector;

    // opaque struct holding relevant data of the edges
    struct JointData;
    // opaque struct holding relevant data of the ref edges in jump calculation
    struct JointRefData;
  
  protected:
    // needed derivatives
#ifdef __3D__
    const MultiIndex_vector derivatives{D000, D100, D010, D001, D200, D110, D101, D020, D011, D002};
#else
    const MultiIndex_vector derivatives{D00, D10, D01, D20, D02};
#endif
    // the selected estimator
    CDErrorEstimatorType estimatorType;

    // member variable indicating if the grid is conform or not
    bool conform_grid;
    
    Parameter space_discretization;

    // internal function calculating jumps across edges at the boundary
    bool handleJump_BoundaryEdge(double *result, const Example_CD &example,
                                 const int estimatorType,
                                 const int N_QuadraturePoints1D,
                                 const double *const &weights1D,
                                 const CDErrorEstimator::JointData &edgeData,
                                 const double meas, const double *coeff,
                                 double linfb, const std::vector<double> &alpha,
                                 int edgeIdx, const TJoint *joint) const;


    // calculates eta_K for a single cell K
    double calculateEtaK(const TBaseCell *cell,
                         const FEFunction &fe_function,
                         const TCollection &coll, const Example_CD & example,
                         std::vector<double *> &derivativesPerQuadPoint,
                         double (&AbsDetjk)[MaxN_QuadPoints_2D],
                         const double *(&weights),
                         std::vector<double *> &coefficientsPerQuadPoint,
                         int n_quadrature_points,
                         const unsigned int N_QuadraturePoints1D,
                         const double *(&weights1D), const JointData &edgeData,
                         const double* zeta, JointRefData &edgeRefData) const;

  public:
    // constructor
    CDErrorEstimator(int type, bool conform_grid, Parameter space_disc);

    void estimate(const Example_CD &ex, const FEFunction &fe_function);
    
    virtual void info() override;
};


#endif //CDERRORESTIMATOR2D_H
