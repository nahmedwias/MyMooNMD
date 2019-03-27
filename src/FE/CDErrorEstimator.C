//
// Created by Moritz Hoffmann on 05/07/15.
//

#include <Database.h>
#include "CDErrorEstimator.h"
#include <ConvDiff.h>
#ifdef __3D__
#include "FEDatabase3D.h"
#else
#endif
#include "FEDatabase2D.h"
#include "BaseCell.h"
#include <array>
#include <cmath>
#include <BoundEdge.h>
#include <memory>

/**
 * Data structure holding relevant data of the edges.
 *
 * XEdge1D, YEdge1D: coordinates of the edges
 * xi1D, eta1D: arrays holding the line quadrature points on reference cell of the form arr[baseFunction][edgeIdx][quadraturePoint] = coordinate
 * xietaval, xideriv, etaderiv: arrays holding values of the derivatives at the quadrature points, same structure as xi1D and eta1D
 * AbsDetjk1D: determinant of the affine mapping for each edge of the form arr[edgeIdx][quadPoint] = value
 * xyval_ref1D: values in reference cell
 * xderiv_ref1D, yderiv_ref1D: derivative values in reference cell
 * xyval_1D: values in original cell
 * xderiv_1D, yderiv_1D: values of derivatives in original cell
 */
template <int d>
struct CDErrorEstimator<d>::JointData
{
  protected:
    std::vector<double> xieta_ref1D_data;
  public:
    // constants
    constexpr static int MaxN_QuadPoints_1D_loc = MaxN_QuadPoints_1D;
    constexpr static int N_BaseFuncts2D_loc = N_BaseFuncts2D;

    // coordinates of the edges
    std::array<std::vector<double>, 4> XEdge1D{}, YEdge1D{};
    // arrays holding the line quadrature points on reference cell
    double xi1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc]{};
    double eta1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc]{};
    // mapping (base_function2D, edge, quadPoint, BaseFunction) -> value
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, 4>> xietaval_ref1D{N_BaseFuncts2D_loc};
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, 4>> xideriv_ref1D{N_BaseFuncts2D_loc};
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, 4>> etaderiv_ref1D{N_BaseFuncts2D_loc};
    // determinant of the affine mapping for each edge
    std::array<std::array<double, MaxN_QuadPoints_2D>, 4> AbsDetjk1D{};
    // values and derivative values in reference cell
    std::vector<std::array<std::array<double, MaxN_QuadPoints_1D_loc>, 4>> xyval_ref1D{};
    std::vector<std::array<std::array<double, MaxN_QuadPoints_1D_loc>, 4>> xderiv_ref1D{};
    std::vector<std::array<std::array<double, MaxN_QuadPoints_1D_loc>, 4>> yderiv_ref1D{};
    // values and derivative values in original cell
    std::array<std::vector<double>, 4> xyval_1D{};
    std::array<std::vector<double>, 4> xderiv_1D{};
    std::array<std::vector<double>, 4> yderiv_1D{};

    JointData()
    {
      // initialize structures holding values and derivatives on reference edge
      xieta_ref1D_data.resize(N_BaseFuncts2D_loc * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc * 3);
      {
        // back xietaval_ref1D, xideriv_ref1D, etaderiv_ref1D by xieta_ref1D_data
        auto *ptr = xieta_ref1D_data.data();
        xietaval_ref1D[0][0][0] = &ptr[0];
        xideriv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D_loc * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc];
        etaderiv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D_loc * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc * 2];
        for(size_t ii = 0; ii < N_BaseFuncts2D_loc; ii++)
        {
          for(size_t jj = 0; jj < 4; jj++)
          {
            for(size_t kk = 0; kk < MaxN_QuadPoints_1D_loc; kk++)
            {
              xietaval_ref1D[ii][jj][kk] = xietaval_ref1D[0][0][0] + (kk * N_BaseFuncts2D_loc + jj * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc + ii * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc);
              xideriv_ref1D[ii][jj][kk] = xideriv_ref1D[0][0][0] + (kk * N_BaseFuncts2D_loc + jj * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc + ii * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc);
              etaderiv_ref1D[ii][jj][kk] = etaderiv_ref1D[0][0][0] + (kk * N_BaseFuncts2D_loc + jj * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc + ii * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc);
            }
          }
        }
      }
    }
};

/**
 * Edge data holding relevant data to the referring edges in the process of jump calculation.
 */
template <int d>
struct CDErrorEstimator<d>::JointRefData
{
  private:
    std::vector<double> refNeigh1D_data;
    unsigned int refDataQuadPoints1D;
  public:
    constexpr static int MaxN_QuadPoints_1D_loc = MaxN_QuadPoints_1D;

    // these data structures depend on the maximal number of 2d base functions
    unsigned int max_n_base_funct_2d;
    std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, MaxN_BaseFunctions2D> xietaval_refNeigh1D;
    std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, MaxN_BaseFunctions2D> xideriv_refNeigh1D;
    std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, MaxN_BaseFunctions2D> etaderiv_refNeigh1D;

    // these data structures depend on the number of quadrature points (and therefore may be updated in each cell)
    std::vector<double> xi1DNeigh;
    std::vector<double> eta1DNeigh;
    std::vector<double> X1DNeigh;
    std::vector<double> Y1DNeigh;
    std::vector<double> X1DCell;
    std::vector<double> Y1DCell;
    std::vector<double> xderiv_Neigh1D;
    std::vector<double> yderiv_Neigh1D;
    std::vector<double> xyval_Neigh1D;
    std::vector<double> xderiv_Cell1D;
    std::vector<double> yderiv_Cell1D;
    std::vector<double> xyval_Cell1D;
    std::array<double *, MaxN_QuadPoints_1D_loc> xyval_refNeigh1D{};
    std::array<double *, MaxN_QuadPoints_1D_loc> xderiv_refNeigh1D{};
    std::array<double *, MaxN_QuadPoints_1D_loc> yderiv_refNeigh1D{};
    std::vector<double> FEFunctValuesNeigh;

    JointRefData (unsigned int max_n_base_funct_2d)
     : refDataQuadPoints1D{0}, max_n_base_funct_2d{max_n_base_funct_2d}
    {
      unsigned int k = MaxN_BaseFunctions2D * MaxN_QuadPoints_1D_loc * max_n_base_funct_2d;
      refNeigh1D_data = std::vector<double>(3 * k);
      xietaval_refNeigh1D[0][0] = &refNeigh1D_data[0];
      xideriv_refNeigh1D[0][0] = &refNeigh1D_data[0] + k;
      etaderiv_refNeigh1D[0][0] = &refNeigh1D_data[0] + 2 * k;
      for(unsigned int i = 0; i < MaxN_BaseFunctions2D; i++)
      {
        unsigned int n = i * MaxN_BaseFunctions2D * max_n_base_funct_2d;
        for(unsigned int j = 0; j < MaxN_QuadPoints_1D_loc; j++)
        {
          unsigned int l = n + j * max_n_base_funct_2d;
          xietaval_refNeigh1D[i][j] = xietaval_refNeigh1D[0][0] + l;
          xideriv_refNeigh1D[i][j] = xideriv_refNeigh1D[0][0] + l;
          etaderiv_refNeigh1D[i][j] = etaderiv_refNeigh1D[0][0] + l;
        }
      }
    }

    void updateQuadPointData(unsigned int refDataQuadPoints1D)
    {
      if (this->refDataQuadPoints1D != refDataQuadPoints1D)
      {
        this->refDataQuadPoints1D = refDataQuadPoints1D;
        xi1DNeigh.resize(refDataQuadPoints1D);
        eta1DNeigh.resize(refDataQuadPoints1D);
        X1DNeigh.resize(refDataQuadPoints1D);
        Y1DNeigh.resize(refDataQuadPoints1D);
        X1DCell.resize(refDataQuadPoints1D);
        Y1DCell.resize(refDataQuadPoints1D);
        xderiv_Neigh1D.resize(refDataQuadPoints1D);
        yderiv_Neigh1D.resize(refDataQuadPoints1D);
        xyval_Neigh1D.resize(refDataQuadPoints1D);
        xderiv_Cell1D.resize(refDataQuadPoints1D);
        yderiv_Cell1D.resize(refDataQuadPoints1D);
        xyval_Cell1D.resize(refDataQuadPoints1D);

        FEFunctValuesNeigh.resize((3 * refDataQuadPoints1D + 1) * MaxN_BaseFunctions2D);
        unsigned int k = MaxN_BaseFunctions2D;
        for(size_t i = 0; i < refDataQuadPoints1D; i++)
        {
          xyval_refNeigh1D[i] = &FEFunctValuesNeigh.data()[0] + k;
          k += MaxN_BaseFunctions2D;
          xderiv_refNeigh1D[i] = &FEFunctValuesNeigh.data()[0] + k;
          k += MaxN_BaseFunctions2D;
          yderiv_refNeigh1D[i] = &FEFunctValuesNeigh.data()[0] + k;
          k += MaxN_BaseFunctions2D;
        }
      }
    }
};

//
// Helper functions
//
namespace
{
  // returns the maximal number of used base functions 2D per element for a fe space
  template<int d>
  unsigned int getMaxN_BaseFunctions(const typename CDErrorEstimator<d>::FESpace &fe_space)
  {
    unsigned int MaxN_BaseFunctions_loc = 0;
    // # used finite elements
    unsigned int n{(unsigned int) fe_space.GetN_UsedElements()};
    // used finite elements
    auto *UsedElements = fe_space.GetUsedElements();
    // for all finite elements
    for(unsigned int j = 0; j < n; j++)
    {
      // if number of base functions is higher, update value
#ifdef __3D__
      unsigned int k = (unsigned int) TFEDatabase3D::GetN_BaseFunctFromFE3D(UsedElements[j]);
#else
      unsigned int k = (unsigned int) TFEDatabase2D::GetN_BaseFunctFromFE2D(UsedElements[j]);
#endif
      MaxN_BaseFunctions_loc = std::max(MaxN_BaseFunctions_loc, k);
    }
    return MaxN_BaseFunctions_loc;
  }

  // this returns weights for the cell residuals and calculates linfb (pass by reference)
  // alpha[0] for H^1 estimator
  // alpha[1] for L^2 estimator
  // alpha[2] for energy norm estimator
  // alpha[3] for energy norm estimator without jumps
  // alpha[4| for the SUPG norm, version 1, John/Novo
  // alpha[5] for the SUPG norm, version 2, John/Novo
  template<int d>
  std::vector<double> getCellWeights(const double hK, bool supg,
                                     const double *coeff, const double hK_tilde,
                                     double &linfb, double &delta_K)
  {
    double nu = coeff[0];
    std::vector<double> alpha;
    alpha.resize(6);
    alpha[0] = hK * hK;
    alpha[1] = hK * hK * hK * hK;
    alpha[2] = hK * hK / nu;
    if(TDatabase::ParamDB->INTERNAL_COERCIVITY > 0)
    {
      if (1.0 / TDatabase::ParamDB->INTERNAL_COERCIVITY < alpha[2])
      {
        // update weight for energy norm estimator
        alpha[2] = 1.0 / TDatabase::ParamDB->INTERNAL_COERCIVITY;
      }
    }
    alpha[3] = alpha[2];
    alpha[4] = alpha[3];
    linfb = std::max(std::abs(coeff[1]), std::abs(coeff[2]));
    if(d == 3)
      linfb = std::max(linfb, std::abs(coeff[3]));
    if(supg)
    {
      // compute stabilization parameter
      std::array<double, d> convection;
      std::copy(coeff+1, coeff+1+d, convection.begin());
      delta_K = Compute_SDFEM_delta<d>(hK, nu, convection, coeff[d+1], linfb);
      if(alpha[4] > 24 * delta_K)
      {
        alpha[4] = 24 * delta_K;
      }
      // second contribution
      alpha[4] += 24 * delta_K;
    }
    alpha[5] = hK;
    if(supg)
    {
      alpha[5] = hK * hK / (3 * std::sqrt(10.0) * nu);
      if(coeff[d+1] > 0)
      {
        if(1.0 / coeff[d+1] < alpha[5])
        {
          alpha[5] = 1.0 / coeff[d+1];
        }
      }
      // compute stabilization parameter
      if(hK_tilde / (std::sqrt(2.0) * linfb) < alpha[5])
      {
        alpha[5] = hK / (std::sqrt(2.0) * linfb);
      }
      alpha[5] *= 2 * alpha[5];
    }
    return alpha;
  }
}


template<int d>
CDErrorEstimator<d>::CDErrorEstimator(int type, bool conform_grid,
                                      Parameter space_disc)
: ErrorEstimator<d>(), estimatorType{CDErrorEstimatorType(type)},
  conform_grid{conform_grid}, space_discretization{space_disc}
{
  if(type < 0 || type >= N_CD2D_ESTIMATOR_TYPES)
    ErrThrow("unknown CDErrorEstimator type ", type, 
             ". Please choose a value between 0 and ", N_CD2D_ESTIMATOR_TYPES);
  this->estimated_global_error.resize(N_CD2D_ESTIMATOR_TYPES);
}

template <int d>
void eta_to_edge_ref(BF2DRefElements bf2DrefelementsNeigh, int neigh_edge,
                     unsigned int N_QuadraturePoints1D, const double * zeta,
                     typename CDErrorEstimator<d>::JointRefData& edgeData)
{
  switch (bf2DrefelementsNeigh)
  {
    case BFUnitSquare:
    {
      // edge 0
      if(neigh_edge == 0)
      {
        // for all quadrature points
        for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
        {
          edgeData.xi1DNeigh[i] = -zeta[i];
          edgeData.eta1DNeigh[i] = -1;
        }
      }
      // edge 1
      if(neigh_edge == 1)
      {
        for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
        {
          // for all quadrature points
          edgeData.xi1DNeigh[i] = 1;
          edgeData.eta1DNeigh[i] = -zeta[i];
        }
      }
      // edge 2
      if(neigh_edge == 2)
      {
        // for all quadrature points
        for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
        {
          edgeData.xi1DNeigh[i] = zeta[i];
          edgeData.eta1DNeigh[i] = 1;
        }
      }

      // edge 3
      if(neigh_edge == 3)
      {
        // for all quadrature points
        for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
        {
          edgeData.xi1DNeigh[i] = -1;
          edgeData.eta1DNeigh[i] = zeta[i];
        }
      }
      break;
    }
    case BFUnitTriangle:
    {
      if(neigh_edge == 0)
      {
        for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
        {
          // for all quadrature points
          edgeData.xi1DNeigh[i] = (-zeta[i] + 1) / 2;
          edgeData.eta1DNeigh[i] = 0;
        }
      }
      if(neigh_edge == 1)
      {
        for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
        {
          // for all quadrature points
          edgeData.xi1DNeigh[i] = (zeta[i] + 1) / 2;
          edgeData.eta1DNeigh[i] = (-zeta[i] + 1) / 2;
        }
      }
      if(neigh_edge == 2)
      {
        for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
        {
          // for all quadrature points
          edgeData.xi1DNeigh[i] = 0;
          edgeData.eta1DNeigh[i] = (zeta[i] + 1) / 2;
        }
      }
      break;
    }
  }
}

template <int d>
void eta_to_edge_ref(bool conform_grid, int part,
                     BF2DRefElements bf2DrefelementsNeigh, int neigh_edge,
                     unsigned int N_QuadraturePoints1D, const double * zeta,
                     typename CDErrorEstimator<d>::JointRefData& edgeRefData)
{
  if(conform_grid)
  {
    // compute coordinates of line quadrature points in reference cell
    eta_to_edge_ref<d>(bf2DrefelementsNeigh, neigh_edge, N_QuadraturePoints1D,
                       zeta, edgeRefData);
  }
  else
  {
    // compute coordinates of line quadrature
    // this is only for 1-regular triangulations
    switch(bf2DrefelementsNeigh)
    {
      // points in reference cell
      case BFUnitSquare:
      {
        // edge 0
        if(neigh_edge == 0)
        {
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            // for all quadrature points
            edgeRefData.xi1DNeigh[i] = (-zeta[i] + part) / 2;
            edgeRefData.eta1DNeigh[i] = -1;
          }
        }
        // edge 1
        if(neigh_edge == 1)
        {
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            // for all quadrature points
            edgeRefData.xi1DNeigh[i] = 1;
            edgeRefData.eta1DNeigh[i] = (-zeta[i] + part) / 2;
          }
        }
        // edge 2
        if(neigh_edge == 2)
        {
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            // for all quadrature points
            edgeRefData.xi1DNeigh[i] = (zeta[i] - part) / 2;
            edgeRefData.eta1DNeigh[i] = 1;
          }
        }
        // edge 3
        if(neigh_edge == 3)
        {
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            // for all quadrature points
            edgeRefData.xi1DNeigh[i] = -1;
            edgeRefData.eta1DNeigh[i] = (zeta[i] - part) / 2;
          }
        }
        break;
      }
      case BFUnitTriangle:
      {
        if(neigh_edge == 0)
        {
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            // for all quadrature points
            if (part == -1) part = 0;
            edgeRefData.xi1DNeigh[i] = ((-zeta[i] + 1) / 2 + part) / 2;
            edgeRefData.eta1DNeigh[i] = 0;
          }
        }
        if(neigh_edge == 1)
        {
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            // for all quadrature points
            if (part == 1) part = 0;
            edgeRefData.xi1DNeigh[i] = ((zeta[i] + 1) / 2 - part) / 2;
            if (part == 0) part = 1;
            if (part == -1) part = 0;
            edgeRefData.eta1DNeigh[i] = ((-zeta[i] + 1) / 2 + part) / 2;
            if (part == 0) part = -1;
          }
        }
        if (neigh_edge == 2)
        {
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            // for all quadrature points
            if (part == 1) part = 0;
            edgeRefData.xi1DNeigh[i] = 0;
            edgeRefData.eta1DNeigh[i] = ((zeta[i] + 1) / 2 - part) / 2;
          }
        }
        break;
      }
    }
  }
}

template <int d>
void get_values_ref(BF2DRefElements bf2Drefelements, int N_Points1D,
                    const double* zeta, const TBaseFunct2D* bf,
                    typename CDErrorEstimator<d>::JointData& edgeData,
                    int baseFunct_loc)
{
  switch (bf2Drefelements)
  {
    // quadrilateral cell
    case BFUnitSquare:
    {
      // edge 0
      for(auto j = 0; j < N_Points1D; j++)
      {
        edgeData.xi1D[baseFunct_loc][0][j] = zeta[j];
        edgeData.eta1D[baseFunct_loc][0][j] = -1;
        bf->GetDerivatives(D00, zeta[j], -1,
                           edgeData.xietaval_ref1D[baseFunct_loc][0][j]);
        bf->GetDerivatives(D10, zeta[j], -1,
                           edgeData.xideriv_ref1D[baseFunct_loc][0][j]);
        bf->GetDerivatives(D01, zeta[j], -1,
                           edgeData.etaderiv_ref1D[baseFunct_loc][0][j]);
      }
      // edge 1
      for(auto j = 0; j < N_Points1D; j++)
      {
        edgeData.xi1D[baseFunct_loc][1][j] = 1;
        edgeData.eta1D[baseFunct_loc][1][j] = zeta[j];
        bf->GetDerivatives(D00, 1, zeta[j],
                           edgeData.xietaval_ref1D[baseFunct_loc][1][j]);
        bf->GetDerivatives(D10, 1, zeta[j],
                           edgeData.xideriv_ref1D[baseFunct_loc][1][j]);
        bf->GetDerivatives(D01, 1, zeta[j],
                           edgeData.etaderiv_ref1D[baseFunct_loc][1][j]);
      }
      // edge 2
      for(auto j = 0; j < N_Points1D; j++)
      {
        edgeData.xi1D[baseFunct_loc][2][j] = -zeta[j];
        edgeData.eta1D[baseFunct_loc][2][j] = 1;
        bf->GetDerivatives(D00, -zeta[j], 1,
                           edgeData.xietaval_ref1D[baseFunct_loc][2][j]);
        bf->GetDerivatives(D10, -zeta[j], 1,
                           edgeData.xideriv_ref1D[baseFunct_loc][2][j]);
        bf->GetDerivatives(D01, -zeta[j], 1,
                           edgeData.etaderiv_ref1D[baseFunct_loc][2][j]);
      }
      // edge 3
      for(auto j = 0; j < N_Points1D; j++)
      {
        edgeData.xi1D[baseFunct_loc][3][j] = -1;
        edgeData.eta1D[baseFunct_loc][3][j] = -zeta[j];
        bf->GetDerivatives(D00, -1, -zeta[j],
                           edgeData.xietaval_ref1D[baseFunct_loc][3][j]);
        bf->GetDerivatives(D10, -1, -zeta[j],
                           edgeData.xideriv_ref1D[baseFunct_loc][3][j]);
        bf->GetDerivatives(D01, -1, -zeta[j],
                           edgeData.etaderiv_ref1D[baseFunct_loc][3][j]);
      }
      break;
    }
    // triangular cell
    case BFUnitTriangle:
    {
      // edge 0
      for(auto j = 0; j < N_Points1D; j++)
      {
        edgeData.xi1D[baseFunct_loc][0][j] = (zeta[j] + 1) / 2;
        edgeData.eta1D[baseFunct_loc][0][j] = 0;
        bf->GetDerivatives(D00, (zeta[j] + 1) / 2, 0,
                           edgeData.xietaval_ref1D[baseFunct_loc][0][j]);
        bf->GetDerivatives(D10, (zeta[j] + 1) / 2, 0,
                           edgeData.xideriv_ref1D[baseFunct_loc][0][j]);
        bf->GetDerivatives(D01, (zeta[j] + 1) / 2, 0,
                           edgeData.etaderiv_ref1D[baseFunct_loc][0][j]);
      }
      // edge 1
      for(auto j = 0; j < N_Points1D; j++)
      {
        edgeData.xi1D[baseFunct_loc][1][j] = (-zeta[j] + 1) / 2;
        edgeData.eta1D[baseFunct_loc][1][j] = (zeta[j] + 1) / 2;
        bf->GetDerivatives(D00, (-zeta[j] + 1) / 2, (zeta[j] + 1) / 2,
                           edgeData.xietaval_ref1D[baseFunct_loc][1][j]);
        bf->GetDerivatives(D10, (-zeta[j] + 1) / 2, (zeta[j] + 1) / 2,
                           edgeData.xideriv_ref1D[baseFunct_loc][1][j]);
        bf->GetDerivatives(D01, (-zeta[j] + 1) / 2, (zeta[j] + 1) / 2,
                           edgeData.etaderiv_ref1D[baseFunct_loc][1][j]);
      }
      // edge 2
      for(auto j = 0; j < N_Points1D; j++)
      {
        edgeData.xi1D[baseFunct_loc][2][j] = 0;
        edgeData.eta1D[baseFunct_loc][2][j] = (-zeta[j] + 1) / 2;
        bf->GetDerivatives(D00, 0, (-zeta[j] + 1) / 2,
                           edgeData.xietaval_ref1D[baseFunct_loc][2][j]);
        bf->GetDerivatives(D10, 0, (-zeta[j] + 1) / 2,
                           edgeData.xideriv_ref1D[baseFunct_loc][2][j]);
        bf->GetDerivatives(D01, 0, (-zeta[j] + 1) / 2,
                           edgeData.etaderiv_ref1D[baseFunct_loc][2][j]);
      }
      break;
    }
  }
}


template<int d>
void CDErrorEstimator<d>::estimate(const Example_CD & example,
                                   const FEFunction &fe_function)
{
  this->eta_K.clear();
  std::fill(this->estimated_global_error.begin(),
            this->estimated_global_error.end(), 0.);

  // initialization
  auto coll = fe_function.GetFESpace2D()->GetCollection();
  // this call is obligatory,
  // otherwise the collection is unknown to the refinement strategy
  this->currentCollection = coll;

  this->eta_K.resize(coll->GetN_Cells());
  this->maximal_local_error = 0.;

  // array of pointers holding quad point -> derivatives
  std::vector<double *> derivativesPerQuadPoint;

  // actual data for quad point -> derivatives in a 1D structure
  std::vector<double> derivativesPerQuadPointData;
  derivativesPerQuadPointData.resize(MaxN_QuadPoints_2D * derivatives.size());

  // map this structure into the 2D array
  for(unsigned int j = 0; j < MaxN_QuadPoints_2D; j++)
  {
    derivativesPerQuadPoint.push_back(&derivativesPerQuadPointData[j * derivatives.size()]);
  }

  // vector containing quadpoint index -> coefficients
  std::vector<double *> coefficientsPerQuadPoint;

  const unsigned int stride = 20;
  std::vector<double> coefficientVectorData;
  coefficientVectorData.resize(MaxN_QuadPoints_2D * stride);
  for(unsigned int j = 0; j < MaxN_QuadPoints_2D; j++)
  {
    coefficientsPerQuadPoint.push_back(&coefficientVectorData[j * stride]);
  }


  // the FE basis functions
  BaseFunct2D *baseFunctions = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();

  // the fe space
  auto fe_space = fe_function.GetFESpace2D();
  // attributes of the reference cell (coordinates and determinant)
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];

  JointData edgeData{};
  // initialize and resize vectors / arrays accordingly
  const unsigned int max_n_base_funct_2d = getMaxN_BaseFunctions<d>(*fe_space);
  JointRefData edgeRefData{max_n_base_funct_2d};
  {
    edgeData.xietaval_ref1D.resize(max_n_base_funct_2d);
    edgeData.xideriv_ref1D.resize(max_n_base_funct_2d);
    edgeData.etaderiv_ref1D.resize(max_n_base_funct_2d);
    edgeData.xyval_ref1D.resize(max_n_base_funct_2d);
    edgeData.xderiv_ref1D.resize(max_n_base_funct_2d);
    edgeData.yderiv_ref1D.resize(max_n_base_funct_2d);
    for(unsigned int x = 0; x < 4; x++)
    {
      edgeData.XEdge1D[x] = std::vector<double>(max_n_base_funct_2d);
      edgeData.YEdge1D[x] = std::vector<double>(max_n_base_funct_2d);
      edgeData.xyval_1D[x] = std::vector<double> {};
      edgeData.xderiv_1D[x] = std::vector<double> {};
      edgeData.yderiv_1D[x] = std::vector<double> {};
    }
  }

  {
    /**
     * prepare error estimates
     */
    // do for all mesh cells
    for(int i = 0; i < coll->GetN_Cells(); i++)
    {
      // on the finest level
      auto cell = coll->GetCell(i);
      // for all edges
      for(int j = 0; j < cell->GetN_Edges(); j++)
      {
        // neighbour cell
        auto neigh = cell->GetJoint(j)->GetNeighbour(cell);
        // set clipboard to -1
        if (neigh) neigh->SetClipBoard(-1);
      }
      cell->SetClipBoard(-1);
    }
    // now non finest neighbours of finest cells have clipboard -1
    for(int i = 0; i < coll->GetN_Cells(); i++)
    {
      // set clipboard of cells on finest
      auto cell = coll->GetCell(i);
      cell->SetClipBoard(i);
    }
  }

  int N_Points;
  const double *xi, *eta, *weights;
  // we do need second derivatives on the reference cell
  bool secondDerivatives[] = {true};

  // ########################################################################
  // calculate values of base functions and derivatives on ref element
  // ########################################################################

  auto n_used_elements = fe_space->GetN_UsedElements();

  // store used finite elements
  std::vector<FE2D> UsedElements(n_used_elements);
  {
    std::vector<int> Used = std::vector<int>(N_FEs2D);
    auto n = fe_space->GetN_UsedElements();
    auto use = fe_space->GetUsedElements();
    for(auto j = 0; j < n; j++)
    {
      FE2D CurrentElement = use[j];
      Used[CurrentElement] = 1;
    }
    int N_UsedElements = 0;
    for(auto i = 0; i < N_FEs2D; i++)
      if (Used[i]) N_UsedElements++;

    size_t j = 0;
    for(auto i = 0; i < N_FEs2D; i++)
    {
      if(Used[i])
      {
        UsedElements[j] = (FE2D) i;
        j++;
      }
    }
    if(n_used_elements != N_UsedElements)
    {
      std::cerr << "neineineineineineineinein" << std::endl;
    }
  }

  int N_Points1D;
  for(auto i = 0; i < n_used_elements; i++)
  {
    FE2D usedElement = UsedElements.at(i);//fe_space->GetUsedElements()[i];
    QuadFormula1D qf = TFEDatabase2D::GetQFLineFromDegree(2 * TFEDatabase2D::GetPolynomialDegreeFromFE2D(usedElement));
    TQuadFormula1D *quadFormula = TFEDatabase2D::GetQuadFormula1D(qf);
    const double *weights1D, *zeta;
    quadFormula->GetFormulaData(N_Points1D, weights1D, zeta);
    if(N_Points1D > edgeData.MaxN_QuadPoints_1D_loc)
    {
      ErrThrow("CD2DErrorEstimator: too many 1D quadrature points ",
               N_Points1D, "Increase  MaxN_QuadPoints_1D_loc !!!");
    }

    {
      for(auto edge = 0; edge < 4; edge++)
      {
        edgeData.xyval_1D[edge].resize((unsigned long) N_Points1D);
        edgeData.xderiv_1D[edge].resize((unsigned long) N_Points1D);
        edgeData.yderiv_1D[edge].resize((unsigned long) N_Points1D);
      }
    }

    BaseFunct2D baseFunct2D = baseFunctions[usedElement];
    TBaseFunct2D *bf = TFEDatabase2D::GetBaseFunct2D(baseFunct2D);// get base functions
    BF2DRefElements bf2Drefelements = bf->GetRefElement();
    int baseFunct_loc = 0;
    {
      int j;
      for(j = 0; j < n_used_elements; j++)
      {
        if((int) baseFunct2D == (int) UsedElements[j])
        {
          break;
        }
      }
      baseFunct_loc = j;
    }

    // compute coordinates of and values of basis functions at line quadrature 
    // points in reference cell
    get_values_ref<d>(bf2Drefelements, N_Points1D, zeta, bf, edgeData,
                      baseFunct_loc);
  }


  // for each cell
  for(int cellIdx = 0; cellIdx < coll->GetN_Cells(); cellIdx++)
  {
    // the cell
    const TBaseCell *cell = coll->GetCell(cellIdx);

    auto fe = fe_space->get_fe(cellIdx);
    // fe2d element on cell
    auto element_id = fe_space->GetFE2D(cellIdx, cell);
    TBaseFunct2D * basis_functions = fe.GetBaseFunct2D();
    int n_base_functs = fe.GetN_DOF();
    
    // get reference transformation RefTrans2D
    RefTrans2D refTrans = TFEDatabase2D::GetOrig(1, &element_id, coll, cell,
                                                 secondDerivatives, N_Points,
                                                 xi, eta, weights, X, Y,
                                                 AbsDetjk);

    // get fe function values for every basis function of the element
    std::vector<double> FEFunctValues(n_base_functs);
    for(int l = 0; l < n_base_functs; l++)
    {
      // fe values of dofs
      int *DOF = fe_space->GetGlobalDOF(cellIdx);
      FEFunctValues[l] = fe_function.GetValues()[DOF[l]];
    }

    // compute values for all derivatives in all quadrature points in original
    // mesh cell for all derivatives
    for(size_t k = 0; k < derivatives.size(); k++)
    {
      // get values in original cell by dof of current mesh cell
      double **OrigFEValues = TFEDatabase2D::GetOrigElementValues(
        basis_functions->GetID(), derivatives[k]);

      // for all quadrature points
      for(long j = 0; j < N_Points; j++)
      {
        // value in original cell
        double *Orig = OrigFEValues[j];
        double value = 0;
        // for all basis functions
        for(int l = 0; l < n_base_functs; l++)
        {
          // accumulate value of derivative in point j
          value += FEFunctValues[l] * Orig[l];
        }
        // for k-th derivative
        derivativesPerQuadPoint[j][k] = value;
      }
    }

    // problem's coefficients
    if(example.get_coeffs())
    {
      auto coeff_fct = example.get_coeffs();
      coeff_fct(N_Points, X, Y, NULL, &coefficientsPerQuadPoint[0]);
    }

    // 1D quadrature formula
    QuadFormula1D qf = TFEDatabase2D::GetQFLineFromDegree(
      2 * TFEDatabase2D::GetPolynomialDegreeFromFE2D(element_id));
    TQuadFormula1D *quadFormula = TFEDatabase2D::GetQuadFormula1D(qf);
    int N_QuadraturePoints1D;
    const double *weights1D, *zeta;
    quadFormula->GetFormulaData(N_QuadraturePoints1D, weights1D, zeta);


    // update data base
    basis_functions->MakeRefElementData(qf);


    int BaseFunct_loc = 0;
    {
      int j = 0;
      for(j = 0; j < n_used_elements; j++)
      {
        if((int) basis_functions->GetID() == (int) fe_space->GetUsedElements()[j])
        {
          break;
        }
      }
      BaseFunct_loc = j;
    }

    // loop over all edges of cell
    for(auto edgeIdx = 0; edgeIdx < cell->GetN_Edges(); edgeIdx++)
    {
      // get original coordinates of edge quad. points
      TFEDatabase2D::GetOrigFromRef(refTrans, N_QuadraturePoints1D,
                                    edgeData.xi1D[BaseFunct_loc][edgeIdx],
                                    edgeData.eta1D[BaseFunct_loc][edgeIdx],
                                    edgeData.XEdge1D[edgeIdx].data(),
                                    edgeData.YEdge1D[edgeIdx].data(),
                                    edgeData.AbsDetjk1D[edgeIdx].data());
      // get values and derivatives in original cell
      for(int k = 0; k < N_QuadraturePoints1D; k++)
      {
        TFEDatabase2D::GetOrigValues(refTrans, edgeData.xi1D[BaseFunct_loc][edgeIdx][k],
                                     edgeData.eta1D[BaseFunct_loc][edgeIdx][k],
                                     basis_functions, coll, (TGridCell *) cell,
                                     edgeData.xietaval_ref1D[BaseFunct_loc][edgeIdx][k],
                                     edgeData.xideriv_ref1D[BaseFunct_loc][edgeIdx][k],
                                     edgeData.etaderiv_ref1D[BaseFunct_loc][edgeIdx][k],
                                     &edgeData.xyval_ref1D[edgeIdx][k][0],
                                     &edgeData.xderiv_ref1D[edgeIdx][k][0],
                                     &edgeData.yderiv_ref1D[edgeIdx][k][0]);
      }
      double val[3];
      // for all quadrature points
      for(auto k = 0; k < N_QuadraturePoints1D; k++)
      {
        val[0] = val[1] = val[2] = 0;
        // for all basis functions
        for(auto l = 0; l < n_base_functs; l++)
        {
          // accumulate value of derivative
          val[0] += FEFunctValues[l] * edgeData.xyval_ref1D[edgeIdx][k][l];
          // accumulate value of derivative
          val[1] += FEFunctValues[l] * edgeData.xderiv_ref1D[edgeIdx][k][l];
          // accumulate value of derivative
          val[2] += FEFunctValues[l] * edgeData.yderiv_ref1D[edgeIdx][k][l];
        }
        /*std::cout << "cell=" << cellIdx << ", " << "xyval_1d[" << edgeIdx << "][" << k << "]=" << val[0] << std::endl;
        for(auto l = 0; l < n_base_functs; l++) {
            std::cout << "\tl="<<l<<", Fe_val="<<FEFunctValues[l] <<", " << "xyval="<<edgeData.xyval_ref1D[edgeIdx][k][l]<<std::endl;
        }*/
        edgeData.xyval_1D[edgeIdx][k] = val[0];
        edgeData.xderiv_1D[edgeIdx][k] = val[1];
        edgeData.yderiv_1D[edgeIdx][k] = val[2];
      }
    }

    if(space_discretization.is("supg"))
    {
      TDatabase::ParamDB->INTERNAL_LOCAL_DOF = cellIdx;
      int cellEdges = cell->GetN_Edges();
      for(int ij = 0; ij < cellEdges; ij++)
      {
        TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
        TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
      }
      if(cellEdges == 3)
      {
        TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      }
    }
    TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1; // reset for every cell


    double result = calculateEtaK(cell, fe_function, *coll, example,
                                  derivativesPerQuadPoint, AbsDetjk,
                                  weights, coefficientsPerQuadPoint,
                                  N_Points,
                                  (const unsigned int) N_QuadraturePoints1D,
                                  weights1D, edgeData, zeta, edgeRefData);

    this->eta_K[cellIdx] = result;
    this->estimated_global_error[int(estimatorType)] += result;
    this->maximal_local_error = std::max(result, this->maximal_local_error);
  } // end loop over cells
  this->estimated_global_error[int(estimatorType)] 
    = std::sqrt(this->estimated_global_error[int(estimatorType)]);
  this->maximal_local_error = std::sqrt(this->maximal_local_error);
}


bool is_Dirichlet_cell(const TBaseCell *cell,
                       BoundCondFunct2D* boundary_condition)
{
  int N_Edges = cell->GetN_Edges();
  // loop over all edges of cell
  for(int j = 0; j < N_Edges; j++)
  {
    auto joint = cell->GetJoint(j);
    if(joint->InnerJoint())
      continue;
    const TBoundEdge *bdryEdge = (const TBoundEdge *) joint;
    int boundary_comp = bdryEdge->GetBoundComp()->GetID();
    double t0, t1;
    bdryEdge->GetParameters(t0, t1);
    // get boundary condition
    BoundCond Cond0;
    //Output::print<1>(comp << " " << w);
    boundary_condition(boundary_comp, (t0 + t1) / 2.0, Cond0);
    // Dirichlet
    if(Cond0 == DIRICHLET)
      return true;
  }
  return false;
}

double compute_estimator_weight(int estimatorType, double hE, double nu,
                                double linfb, const std::vector<double> &alpha,
                                double meas)
{
  std::vector<double> beta(N_CD2D_ESTIMATOR_TYPES - 1, 0);
  // weight for H^1 estimator
  beta[0] = hE;
  // weight for L^2 estimator
  beta[1] = hE * hE * hE;
  // weight for energy norm estimator
  beta[2] = hE / nu;
  if(TDatabase::ParamDB->INTERNAL_COERCIVITY > 0)
  {
    double w = 1 / std::sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * nu);
    if(w < beta[2])
        beta[2] = w;
  }
  // beta[3] remains zero, because Energy_ResidualEstimatorWithoutJumps
  if(24.0 / linfb < beta[2])
  {
    beta[4] = 24.0 / linfb;
  }
  else
  {
    beta[4] = beta[2];
  }
  beta[5] = alpha[5] * hE / (4.0 * meas);
  return beta[estimatorType - 1];
}

template<int d>
double CDErrorEstimator<d>::calculateEtaK(
  const TBaseCell *cell, const FEFunction &fe_function, const TCollection &coll,
  const Example_CD & example,
  std::vector<double *> &derivativesPerQuadPoint,
  double (&AbsDetjk)[MaxN_QuadPoints_2D], const double *(&weights),
  std::vector<double *> &coefficientsPerQuadPoint, int n_quadrature_points,
  const unsigned int N_QuadraturePoints1D, const double *(&weights1D),
  const JointData &edgeData, const double *zeta, JointRefData &edgeRefData) const
{
  // SUPG parameter
  double delta_K = 0.0;

  auto fe_space = fe_function.GetFESpace2D();

  if(TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS
     && is_Dirichlet_cell(cell, example.get_bc(0)))
  {
    return 0;
  }

  double result = 0.0;
  if(estimatorType == CDErrorEstimatorType::GradientIndicator)
  {
    // calculate the eta_K
    for(int j = 0; j < n_quadrature_points; j++)
    {
      // all derivatives in quadrature points
      const double *deriv = derivativesPerQuadPoint[j];
      // weight of the gradient obtained by transformation of the element
      const double w = weights[j] * AbsDetjk[j];

      // x derivative
      const double e1 = deriv[0];
      // y derivative
      const double e2 = deriv[1];
      result += w * (e1 * e1 + e2 * e2);
    }
    return result;
  }
  else
  {
    /**
     * allocate space
     */
    double absdet1D[MaxN_QuadPoints_2D];

    // values of the fe function
    const double *Values = fe_function.GetValues();

    // weights of the jumps
    double beta[N_CD2D_ESTIMATOR_TYPES - 1];
    /**
     * initialize space
     */
    {
      edgeRefData.updateQuadPointData(N_QuadraturePoints1D);
    }

    /**
     * other initializations
     */
    bool check_cont;
    if(TDatabase::ParamDB->ANSATZ_ORDER > 0)
      check_cont = 1;
    else
      check_cont = 0;

    const double hK = cell->GetDiameter();
    const int N_Edges = cell->GetN_Edges();
    // INTERNAL_VERTEX_X and ..._Y are used in Mesh_size_in_convection_direction
    for(int i = 0; i < N_Edges; i++)
    {
      TDatabase::ParamDB->INTERNAL_VERTEX_X[i] = cell->GetVertex(i)->GetX();
      TDatabase::ParamDB->INTERNAL_VERTEX_Y[i] = cell->GetVertex(i)->GetY();
    }
    if(N_Edges == 3)
    {
      TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
    }

    /**
     * Calculate strong residual
     */
    const double meas = cell->GetMeasure();
    double strong_residual = 0;
    // for all quadrature points
    for(int i = 0; i < n_quadrature_points; i++)
    {
      double *coeff = coefficientsPerQuadPoint[i];
      // all derivatives in quadrature points
      const double *deriv = derivativesPerQuadPoint[i];
      const double w = weights[i] * AbsDetjk[i];
      // strong residual
      double e1 = -coeff[0] * (deriv[3] + deriv[4]) + coeff[1] * deriv[0] + coeff[2] * deriv[1] + coeff[3] * deriv[2] - coeff[4];

      // time-dependent problem
//       if(TDatabase::ParamDB->P15 == -4711.0)
//       {
//         e1 += coeff[6];
//       }
      // L^2 norm
      strong_residual += w * e1 * e1;
    }
    double *coeff = coefficientsPerQuadPoint[n_quadrature_points - 1];

    double hK_tilde = Mesh_size_in_convection_direction<d>(
      hK, {coeff[1], coeff[2]});

    // getting the cell weight in a vector [weight_0,...,weight_5], see description of the method
    double linfb;
    std::vector<double> alpha = getCellWeights<d>(
      hK, space_discretization.is("supg"), coeff, hK_tilde, linfb, delta_K);
    double estimator_cell_weight = alpha[int(estimatorType) - 1];

    // this is the accordingly weighted strong residual without jumps
    result = estimator_cell_weight * strong_residual;

    /**
     * Calculate jumps
     */
    for(int edgeIdx = 0; edgeIdx < cell->GetN_Edges(); edgeIdx++)
    {
      const TJoint *joint = cell->GetJoint(edgeIdx);

      // boundary edge
      if(joint->GetType() == BoundaryEdge
        || joint->GetType() == IsoBoundEdge)
      {
        // This is a boundary edge
        double boundaryJump = 0.0;
        bool shouldNotAbortEstimation = handleJump_BoundaryEdge(
          &boundaryJump, example, int(estimatorType), N_QuadraturePoints1D,
          weights1D, edgeData, meas, coeff, linfb, alpha, edgeIdx, joint);
        if(shouldNotAbortEstimation)
        {
          result += boundaryJump;
        }
        else
        {
          return 0;
        }
      }
      else
      {
        // this is an interior edge
        const TRefDesc *refdesc = cell->GetRefDesc();
        // get refinement descriptor
        const int *TmpEdVer;
        refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
        // get vertices of edge edgeIdx
        const TVertex *ver0 = cell->GetVertex(TmpEdVer[2 * edgeIdx]);
        const TVertex *ver1 = cell->GetVertex(TmpEdVer[2 * edgeIdx + 1]);
        // coordinates of edge "edgeIdx"
        double cell_x0 = ver0->GetX();
        double cell_y0 = ver0->GetY();
        double cell_x1 = ver1->GetX();
        double cell_y1 = ver1->GetY();
        double jump = 0;
        const TBaseCell *neigh = joint->GetNeighbour(cell);
        // compute normal
        double nx = cell_y1 - cell_y0;
        double ny = cell_x0 - cell_x1;
        // length of edge
        double hE = std::sqrt(nx * nx + ny * ny);
        // normalized normal vector
        nx /= hE;
        ny /= hE;
        
        //Output::print<2>(" A ", cell_x0, " ", cell_y0);
        //Output::print<2>(" B ", cell_x1, " ", cell_y1);
        //Output::print<2>(" n ", nx, " ", ny);

        /*************************************************************************/
        /*  no neighbour, find neighbour of parent                               */
        /*************************************************************************/
        if (!neigh)
        {
          // there is no neighbour on the same level
          //  => finer cell in 1 regularity

          // parent edge vertex, child edge, newEdgeOldEdge
          const int *TmpEdVerParent, *TmpCE, *TmpoEnlE;
          // maxlen of child edge
          int MaxLen1;
          // parent cell
          const TBaseCell *parent = cell->GetParent();
          refdesc = parent->GetRefDesc();
          refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerParent);
          refdesc->GetChildEdge(TmpCE, MaxLen1);
          refdesc->GetNewEdgeOldEdge(TmpoEnlE);
          // the current cell's local child number of parent cell
          int neighIdx = 0;
          {
            // local child number
            while (parent->GetChild(neighIdx) != cell) neighIdx++;
          }
          // number of father edge
          int parent_edge = TmpCE[neighIdx * MaxLen1 + edgeIdx];
          // number of father edge
          parent_edge = TmpoEnlE[parent_edge];

          const TJoint *parent_joint = parent->GetJoint(parent_edge);
          // neighbour to parent
          neigh = parent_joint->GetNeighbour(parent);
          if (!neigh)
          {
            std::cerr << "cell was no boundary cell and had no neighbour, "
                         "thus its parent cell must have a neighbor\n";
          }

          // edge index of the parents neighbor
          int neigh_edge = 0;
          {
            // iterate through the neighbours edges until the neighbor's neighbor is the parent cell itself again
            while (neigh->GetJoint(neigh_edge)->GetNeighbour(neigh) != parent) neigh_edge++;
          }
          // vertices of edge
          const TVertex *ver2 = neigh->GetVertex(TmpEdVerParent[2 * neigh_edge]);
          const TVertex *ver3 = neigh->GetVertex(TmpEdVerParent[2 * neigh_edge + 1]);
          // variable indicating if it is the first or the second part of the long edge, i.e.,
          // comparison if ver1 (the second vertex of the cell) matches ver2 (first vertex of parent edge)
          // or ver0 (the first vertex of the cell) matches ver3 (the second vertex of the parent edge)
          int part = 0;
          if(ver1 == ver2)
          {
            // first part of long edge
            part = -1;
          }
          else if (ver0 == ver3)
          {
            // second part of long edge
            part = 1;
          }
          else
          {
            std::cerr << "\"2nd vertex of cell edge\" != \"1st vertex of parent edge\" and \"1st vertex of cell edge\" != \"2nd vertex of parent edge\"" << std::endl
            << "this should not happen since the vertices are ordered" << std::endl;
          }

          // number of neighbour in iterator
          int neigh_N_ = neigh->GetClipBoard();
          if (neigh_N_ == -1)
          {
            std::cerr << "neighbors clipboard (i.e., number in iterator) was -1 which should not happen" << std::endl;
          }

          // finite element on neighbour
          FE2D CurrEleNeigh = fe_space->GetFE2D(neigh_N_, neigh);
          TFE2D *eleNeigh = TFEDatabase2D::GetFE2D(CurrEleNeigh);

          // basis functions on neighbour
          BaseFunct2D BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();
          // number of basis functions
          int N_Neigh = eleNeigh->GetN_DOF();

          TBaseFunct2D *bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);
          // reference cell of neighbour
          BF2DRefElements bf2DrefelementsNeigh = bfNeigh->GetRefElement();

          eta_to_edge_ref<d>(conform_grid, part, bf2DrefelementsNeigh,
                             neigh_edge, N_QuadraturePoints1D, zeta,
                             edgeRefData);
          
          // compute gradients in reference cell of the neighbour
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            bfNeigh->GetDerivatives(D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
            bfNeigh->GetDerivatives(D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
            bfNeigh->GetDerivatives(D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
          }

          // reftrafo of neighbour
          RefTrans2D RefTransNeigh = eleNeigh->GetRefTransID();
          TFEDatabase2D::SetCellForRefTrans(neigh, RefTransNeigh);

          int *GlobalNumbers = fe_space->GetGlobalNumbers();
          int *BeginIndex = fe_space->GetBeginIndex();

          int *DOF = GlobalNumbers + BeginIndex[neigh_N_];
          for(int i = 0; i < N_Neigh; i++)
          {
            edgeRefData.FEFunctValuesNeigh[i] = Values[DOF[i]];
            //Output::print<2>(" value ", edgeRefData.FEFunctValuesNeigh[i]);
          }

          // get values and derivatives in original cell
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            TFEDatabase2D::GetOrigValues(RefTransNeigh, edgeRefData.xi1DNeigh[i],
                                         edgeRefData.eta1DNeigh[i],
                                         TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                                         &coll, (TGridCell *) neigh,
                                         edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                                         edgeRefData.xideriv_refNeigh1D.data()[BaseFunctNeigh][i],
                                         edgeRefData.etaderiv_refNeigh1D.data()[BaseFunctNeigh][i],
                                         edgeRefData.xyval_refNeigh1D.data()[i],
                                         edgeRefData.xderiv_refNeigh1D.data()[i],
                                         edgeRefData.yderiv_refNeigh1D.data()[i]);
          }

          double val[3];
          // for all quadrature points
          for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
          {
            val[0] = val[1] = val[2] = 0;
            // for all basis functions
            for(neighIdx = 0; neighIdx < N_Neigh; neighIdx++)
            {
              // accumulate value of derivative
              val[0] += edgeRefData.FEFunctValuesNeigh[neighIdx] * edgeRefData.xderiv_refNeigh1D[i][neighIdx];
              // accumulate value of derivative
              val[1] += edgeRefData.FEFunctValuesNeigh[neighIdx] * edgeRefData.yderiv_refNeigh1D[i][neighIdx];
              // accumulate value of derivative
              val[2] += edgeRefData.FEFunctValuesNeigh[neighIdx] * edgeRefData.xyval_refNeigh1D[i][neighIdx];
              //Output::print<2>(neighIdx, "  ",
              //                 edgeRefData.xderiv_refNeigh1D[i][neighIdx],
              //                 "  ",
              //                 edgeRefData.yderiv_refNeigh1D[i][neighIdx],
              //                 "  ",
              //                 edgeRefData.FEFunctValuesNeigh[neighIdx]);
            }
            edgeRefData.xderiv_Neigh1D[i] = val[0];
            edgeRefData.yderiv_Neigh1D[i] = val[1];
            edgeRefData.xyval_Neigh1D[i] = val[2];
          }

          TFEDatabase2D::GetOrigFromRef(RefTransNeigh, N_QuadraturePoints1D, edgeRefData.xi1DNeigh.data(), edgeRefData.eta1DNeigh.data(), edgeRefData.X1DNeigh.data(), edgeRefData.Y1DNeigh.data(), absdet1D);
          jump = 0.0;
          double absdetjk1D = hE / 2.0;
          // compute jump
          for(size_t i = 0; i < N_QuadraturePoints1D; i++)
          {
            if((std::abs(edgeData.XEdge1D[edgeIdx][i]
                       - edgeRefData.X1DNeigh[i])
              + std::abs(edgeData.YEdge1D[edgeIdx][i]
                       - edgeRefData.Y1DNeigh[i])) > 1e-8)
            {
              Output::print(" wrong quad points 1 ",
                            edgeData.XEdge1D[edgeIdx][i], " , ",
                            edgeData.YEdge1D[edgeIdx][i], "   ",
                            edgeRefData.X1DNeigh[i], " , ",
                            edgeRefData.Y1DNeigh[i]);
            }
            if(check_cont 
               && std::abs(edgeRefData.xyval_Neigh1D[i]
                           - edgeData.xyval_1D[edgeIdx][i]) > 1e-8)
            {
                Output::print("quad points a ", edgeData.XEdge1D[edgeIdx][i],
                              " , ", edgeData.YEdge1D[edgeIdx][i]);
                Output::print(" i ", i, " vala ",
                              edgeData.xyval_1D[edgeIdx][i], " neigha ",
                              edgeRefData.xyval_Neigh1D[i], " ",
                              std::abs(edgeData.xyval_1D[edgeIdx][i]
                                       - edgeRefData.xyval_Neigh1D[i]));
            }
            double e1 = coeff[0] * ((edgeData.xderiv_1D[edgeIdx][i]
                                     - edgeRefData.xderiv_Neigh1D[i]) * nx
                                  + (edgeData.yderiv_1D[edgeIdx][i]
                                     - edgeRefData.yderiv_Neigh1D[i]) * ny);
            //Output::print<2>(i, " jumpx ", edgeData.xderiv_1D[edgeIdx][i],
            //                 " ", edgeRefData.xderiv_Neigh1D[i]);
            double w = weights1D[i] * absdetjk1D;
            jump += w * e1 * e1;                       // integral on the edge
          }
          //Output::print<2>("jump ", jump);
          
          result += compute_estimator_weight((int)estimatorType, hE, coeff[0],
                                             linfb, alpha, meas) * jump / 2.0;
        }
        /*************************************************************************/
        /*  neighbour is not on the finest level, find children of neighbour     */
        /*************************************************************************/
        else // there is a neighbour on the same level
        {
          int n = neigh->GetClipBoard();
          if(n == -1)
          {
            // the neighbour is no member of the collection
            // check whether the children of neigh are in collection
            // find the local edge of neigh on which cell is -> l

            int edge2neigh = 0;
            while (neigh->GetJoint(edge2neigh)->GetNeighbour(neigh) != cell)
              edge2neigh++;                       // find connections between cells
            refdesc = neigh->GetRefDesc();          // ref desc of neighbour
            int MaxLen1, MaxLen2, MaxLen3;
            const int *TmpEdVerNeigh, *TmpoEnE, *TmpLen1, *TmpEC, *TmpLen2, *TmpoEnlE, *TmpECI, *TmpLen3;
            // get edges
            refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerNeigh);
            // get connection to child edges
            refdesc->GetOldEdgeNewEdge(TmpoEnE, TmpLen1, MaxLen1);
            // get cell belonging to child edge (TmpEC)
            refdesc->GetEdgeChild(TmpEC, TmpLen2, MaxLen2);
            // get local no.s of child edge
            refdesc->GetOldEdgeNewLocEdge(TmpoEnlE);
            // not general !!!
            const unsigned int N_child = conform_grid ? 1 : 2;

            // find children of neigh on face l -> child
            for(unsigned int r = 0; r < N_child; r++)
            {
              if(!TmpoEnE)
              {
                break;
              }
              // edge child, not general !!!
              int edge1 = TmpoEnE[edge2neigh * MaxLen1 + r];
              // local number of child cell
              int chnum1 = TmpEC[edge1 * MaxLen2];
              // child cell
              const TBaseCell *child = neigh->GetChild(chnum1);
              // get local indices of child edge
              refdesc->GetEdgeChildIndex(TmpECI, TmpLen3, MaxLen3);
              // local index of child edge
              int l_child = TmpECI[edge1 * MaxLen3];

              // ref desc of child
              const TRefDesc *refdesc_child = child->GetRefDesc();
              // conn. edge -> vertices
              refdesc_child->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
              // vertices of edge
              const TVertex *ver2 = child->GetVertex(TmpEdVer[2 * l_child]);
              const TVertex *ver3 = child->GetVertex(TmpEdVer[2 * l_child + 1]);

              //Output::print<2>("ver 0 ", ver0->GetX(), "  ", ver0->GetY());
              //Output::print<2>("ver 1 ", ver1->GetX(), "  ", ver1->GetY());
              //Output::print<2>("ver 2 ", ver2->GetX(), "  ", ver2->GetY());
              //Output::print<2>("ver 3 ", ver3->GetX(), "  ", ver3->GetY());

              int part = 0;
              if (ver1 == ver2)
              {
                part = 1;
              }
              else if (ver0 == ver3)
              {
                part = -1;
              }
              else
              {
                cout << " something wrong 5 " << endl;
                cout << "ver 0 " << ver0->GetX() << "  " << ver0->GetY() << endl;
                cout << "ver 1 " << ver1->GetX() << "  " << ver1->GetY() << endl;
                cout << "ver 2 " << ver2->GetX() << "  " << ver2->GetY() << endl;
                cout << "ver 3 " << ver3->GetX() << "  " << ver3->GetY() << endl;
              }
              // now from point of view of child cell -> cell becomes the neighbour
              // prepare intergration for the half part of edge j

              // number of original cell  in iterator
              int neigh_N_ = cell->GetClipBoard();
              if(neigh_N_ == -1)
              {
                cout << "Clipboard was overwritten, cell's clipboard value was -1." << endl;
              }
              // finite element on neighbour
              FE2D CurrEleNeigh = fe_space->GetFE2D(neigh_N_, cell);
              TFE2D *eleNeigh = TFEDatabase2D::GetFE2D(CurrEleNeigh);

              // basis functions on neighbout
              BaseFunct2D BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();
              int N_Neigh = eleNeigh->GetN_DOF();     // number of basis functions

              TBaseFunct2D *bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);
              // referenz cell of neighbour
              BF2DRefElements bf2DrefelementsNeigh = bfNeigh->GetRefElement();

              int neigh_edge = edgeIdx;
              eta_to_edge_ref<d>(conform_grid, part, bf2DrefelementsNeigh,
                                 neigh_edge, N_QuadraturePoints1D, zeta,
                                 edgeRefData);
              
              // compute gradients in reference cell of the neighbour
              for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)          // for all quadrature points
              {
                bfNeigh->GetDerivatives(D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
              }
              // reftrafo of neighbour
              RefTrans2D RefTransNeigh = eleNeigh->GetRefTransID();
              TFEDatabase2D::SetCellForRefTrans(cell, RefTransNeigh);

              int *DOF = fe_space->GetGlobalNumbers() + fe_space->GetBeginIndex()[neigh_N_];
              for(int i = 0; i < N_Neigh; i++)
              {
                edgeRefData.FEFunctValuesNeigh[i] = Values[DOF[i]];
                //Output::print<2>(" value ", edgeRefData.FEFunctValuesNeigh[i]);
              }

              for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)          // get values and derivatives in original cell
              {
                TFEDatabase2D::GetOrigValues(RefTransNeigh, edgeRefData.xi1DNeigh[i],
                                             edgeRefData.eta1DNeigh[i],
                                             TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                                             &coll, (TGridCell *) neigh,
                                             edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                                             edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
                                             edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
                                             edgeRefData.xyval_refNeigh1D[i],
                                             edgeRefData.xderiv_refNeigh1D[i],
                                             edgeRefData.yderiv_refNeigh1D[i]);
              }

              double val[3];
              for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)           // for all quadrature points
              {
                val[0] = val[1] = val[2] = 0;
                for(int l = 0; l < N_Neigh; l++)            // for all basis functions
                {
                  // accumulate value of derivative
                  val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l];
                  // accumulate value of derivative
                  val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l];
                  // accumulate value of derivative
                  val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l];
                  //Output::print<2>(l, "  ",
                  //                 edgeRefData.xderiv_refNeigh1D[i][l], "  ",
                  //                 edgeRefData.yderiv_refNeigh1D[i][l], "  ",
                  //                 edgeRefData.FEFunctValuesNeigh[l]);
                }
                edgeRefData.xderiv_Cell1D[i] = val[0];         // for k-th
                edgeRefData.yderiv_Cell1D[i] = val[1];         // for k-th
                edgeRefData.xyval_Cell1D[i] = val[2];          // for k-th
              }

              TFEDatabase2D::GetOrigFromRef(RefTransNeigh, N_QuadraturePoints1D, edgeRefData.xi1DNeigh.data(),
                                            edgeRefData.eta1DNeigh.data(),
                                            edgeRefData.X1DCell.data(), edgeRefData.Y1DCell.data(), absdet1D);
              // prepare integration for the child of the neighbour belong to the half part
              // of edge edgeIdx

              // number of neighbour in iterator
              neigh_N_ = child->GetClipBoard();
              // finite element on neighbour
              CurrEleNeigh = fe_space->GetFE2D(neigh_N_, child);
              eleNeigh = TFEDatabase2D::GetFE2D(CurrEleNeigh);

              // basis functions on neighbout
              BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();
              // number of basis functions
              N_Neigh = eleNeigh->GetN_DOF();

              bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);
              // referenz cell of neighbour
              bf2DrefelementsNeigh = bfNeigh->GetRefElement();

              neigh_edge = l_child;
              // compute coordinates of line quadrature
              // points in reference cell
              switch (bf2DrefelementsNeigh)
              {
                case BFUnitSquare :
                  if (neigh_edge == 0) {
                    for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                      edgeRefData.xi1DNeigh[i] = zeta[i];
                      edgeRefData.eta1DNeigh[i] = -1;
                    }
                  }
                  if (neigh_edge == 1) {
                    for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                      edgeRefData.xi1DNeigh[i] = 1;
                      edgeRefData.eta1DNeigh[i] = zeta[i];
                    }
                  }
                  if (neigh_edge == 2) {
                    for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                      edgeRefData.xi1DNeigh[i] = -zeta[i];
                      edgeRefData.eta1DNeigh[i] = 1;
                    }
                  }

                  if (neigh_edge == 3) {
                    for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                      edgeRefData.xi1DNeigh[i] = -1;
                      edgeRefData.eta1DNeigh[i] = -zeta[i];
                    }
                  }
                  break;

                case BFUnitTriangle :
                  if (neigh_edge == 0) {
                    for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                      edgeRefData.xi1DNeigh[i] = (zeta[i] + 1) / 2;
                      edgeRefData.eta1DNeigh[i] = 0;
                    }
                  }
                  if (neigh_edge == 1) {
                    for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                      edgeRefData.xi1DNeigh[i] = (-zeta[i] + 1) / 2;
                      edgeRefData.eta1DNeigh[i] = (zeta[i] + 1) / 2;
                    }
                  }
                  if (neigh_edge == 2) {
                    for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                      edgeRefData.xi1DNeigh[i] = 0;
                      edgeRefData.eta1DNeigh[i] = (-zeta[i] + 1) / 2;
                    }
                  }
                  break;
              }

              //for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
              //{
              //  Output::print<2>("xiN ", edgeRefData.xi1DNeigh[i], " etaN ",
              //                   edgeRefData.eta1DNeigh[i]);
              //}

              // compute gradients in reference cell of the neighbour
              for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)          // for all quadrature points
              {
                bfNeigh->GetDerivatives(
                  D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i],
                  edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(
                  D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i],
                  edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(
                  D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i],
                  edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
              }
              // reftrafo of neighbour
              RefTransNeigh = eleNeigh->GetRefTransID();
              TFEDatabase2D::SetCellForRefTrans(child, RefTransNeigh);

              DOF = fe_space->GetGlobalNumbers() + fe_space->GetBeginIndex()[neigh_N_];
              for(int i = 0; i < N_Neigh; i++)
              {
                edgeRefData.FEFunctValuesNeigh[i] = Values[DOF[i]];
                //Output::print<2>(" value ", edgeRefData.FEFunctValuesNeigh[i]);
              }

              // get values and derivatives in original cell
              for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
              {
                TFEDatabase2D::GetOrigValues(RefTransNeigh, edgeRefData.xi1DNeigh[i],
                                             edgeRefData.eta1DNeigh[i],
                                             TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                                             &coll, (TGridCell *) neigh,
                                             edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                                             edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
                                             edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
                                             edgeRefData.xyval_refNeigh1D[i],
                                             edgeRefData.xderiv_refNeigh1D[i],
                                             edgeRefData.yderiv_refNeigh1D[i]);
              }

              // for all quadrature points
              for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
              {
                val[0] = val[1] = val[2] = 0;
                // for all basis functions
                for(int l = 0; l < N_Neigh; l++)
                {
                  // accumulate value of derivative
                  val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l];
                  // accumulate value of derivative
                  val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l];
                  // accumulate value of derivative
                  val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l];
                  //Output::print<2>(l, "  ",
                  //                 edgeRefData.xderiv_refNeigh1D[i][l], "  ",
                  //                 edgeRefData.yderiv_refNeigh1D[i][l], "  ",
                  //                 edgeRefData.FEFunctValuesNeigh[l]);
                }                                 // endfor l
                edgeRefData.xderiv_Neigh1D[i] = val[0];        // for k-th
                edgeRefData.yderiv_Neigh1D[i] = val[1];        // for k-th
                edgeRefData.xyval_Neigh1D[i] = val[2];         // for k-th
              }                                   // endfor i

              TFEDatabase2D::GetOrigFromRef(RefTransNeigh, N_QuadraturePoints1D, edgeRefData.xi1DNeigh.data(),
                                            edgeRefData.eta1DNeigh.data(),
                                            edgeRefData.X1DNeigh.data(), edgeRefData.Y1DNeigh.data(), absdet1D);
              jump = 0.0;
              // only half edge is considered
              double absdetjk1D = hE / (2.0 * N_child);
              // compute jump
              for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)
              {
                if((std::abs(edgeRefData.X1DCell[i] - edgeRefData.X1DNeigh[i])
                  + std::abs(edgeRefData.Y1DCell[i] - edgeRefData.Y1DNeigh[i])) > 1e-8)
                  Output::print(" wrong quad points 2 ",
                                edgeRefData.X1DCell[i], " , ",
                                edgeRefData.Y1DCell[i], "   ",
                                edgeRefData.X1DNeigh[i], " , ",
                                edgeRefData.Y1DNeigh[i]);
                if(check_cont)
                {
                  if (std::abs(edgeRefData.xyval_Neigh1D[i] - edgeRefData.xyval_Cell1D[i]) > 1e-8)
                  {
                    cout << "quad points b " << edgeRefData.X1DCell[i] << " , " << edgeRefData.Y1DCell[i] << endl;
                    cout << " i " << i << " valb " << edgeRefData.xyval_Cell1D[i] << " neighb " << edgeRefData.xyval_Neigh1D[i] << " " << std::abs(edgeRefData.xyval_Cell1D[i] - edgeRefData.xyval_Neigh1D[i]) << endl;
                  }
                }
                double e1 = coeff[0] * ((edgeRefData.xderiv_Cell1D[i] - edgeRefData.xderiv_Neigh1D[i]) * nx
                                        + (edgeRefData.yderiv_Cell1D[i] - edgeRefData.yderiv_Neigh1D[i]) * ny);
                //Output::print<2>(i, " jumpx ", edgeRefData.xderiv_Cell1D[i],
                //                 " ", edgeRefData.xderiv_Neigh1D[i]);
                double w = weights1D[i] * absdetjk1D;
                jump += w * e1 * e1;                   // integral on the edge
              }
              //Output::print<2>("jump ", jump);
              result += compute_estimator_weight((int)estimatorType, hE,
                                                 coeff[0], linfb, alpha, meas)
                        * jump / 2.0;
            }
          }                                       // end clipboard==-1
          else
          /*************************************************************************/
          /*  neighbour is on the finest level                                     */
          /*************************************************************************/
          {                                       // the neighbour is a member of the collection
            // find the finite element on the other side
            // find the local edge of neigh on which cell is -> l
            const int *TmpEdVerNeigh;
            int neigh_edge = 0;
            while (neigh->GetJoint(neigh_edge)->GetNeighbour(neigh) != cell) neigh_edge++;
            refdesc = neigh->GetRefDesc();
            refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerNeigh);
            ver0 = cell->GetVertex(TmpEdVer[2 * edgeIdx]);
            ver1 = cell->GetVertex(TmpEdVer[2 * edgeIdx + 1]);
            // vertices of edge
            const TVertex *ver2 = neigh->GetVertex(TmpEdVerNeigh[2 * neigh_edge]);
            const TVertex *ver3 = neigh->GetVertex(TmpEdVerNeigh[2 * neigh_edge + 1]);
            if(!(((ver0 == ver2) && (ver1 == ver3)) || ((ver0 == ver3) && (ver1 == ver2))))
            {
              cout << "wrong edge " << endl;
            }

            // compute gradient at the quadrature points on the edge of
            // the neighbour element

            int neigh_N_ = neigh->GetClipBoard();     // number of neighbour in iterator
            // finite element on neighbour
            FE2D CurrEleNeigh = fe_space->GetFE2D(neigh_N_, neigh);
            TFE2D *eleNeigh = TFEDatabase2D::GetFE2D(CurrEleNeigh);

            // basis functions on neighbout
            BaseFunct2D BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();
            int N_Neigh = eleNeigh->GetN_DOF();       // number of basis functions

            TBaseFunct2D *bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);
            // reference cell of neighbour
            BF2DRefElements bf2DrefelementsNeigh = bfNeigh->GetRefElement();

            eta_to_edge_ref<d>(bf2DrefelementsNeigh, neigh_edge,
                               N_QuadraturePoints1D, zeta, edgeRefData);
            
            //for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)
            //{
            //  Output::print<2>("xiN ", edgeRefData.xi1DNeigh[i], " etaN ",
            //                   edgeRefData.eta1DNeigh[i]);
            //}

            // compute gradients in reference cell of the neighbour
            for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)            // for all quadrature points
            {
              bfNeigh->GetDerivatives(
                D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i],
                edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
              bfNeigh->GetDerivatives(
                D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i],
                edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
              bfNeigh->GetDerivatives(
                D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i],
                edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
            }
            // reftrafo of neighbour
            RefTrans2D RefTransNeigh = eleNeigh->GetRefTransID();
            TFEDatabase2D::SetCellForRefTrans(neigh, RefTransNeigh);

            int *DOF = fe_space->GetGlobalNumbers() + fe_space->GetBeginIndex()[neigh_N_];
            for(int i = 0; i < N_Neigh; i++)
            {
              edgeRefData.FEFunctValuesNeigh[i] = Values[DOF[i]];
              //Output::print<2>(" value ", edgeRefData.FEFunctValuesNeigh[i]);
            }
            for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)            // get values and derivatives in original cell
            {
              TFEDatabase2D::GetOrigValues(RefTransNeigh, edgeRefData.xi1DNeigh[i],
                                           edgeRefData.eta1DNeigh[i],
                                           TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                                           &coll, (TGridCell *) neigh,
                                           edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                                           edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
                                           edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
                                           edgeRefData.xyval_refNeigh1D[i],
                                           edgeRefData.xderiv_refNeigh1D[i],
                                           edgeRefData.yderiv_refNeigh1D[i]);
            }

            double val[3];
            for(unsigned int i = 0; i < N_QuadraturePoints1D; i++)             // for all quadrature points
            {
              val[0] = val[1] = val[2] = 0;
              for(int l = 0; l < N_Neigh; l++)              // for all basis functions
              {
                // accumulate value of derivative
                val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l];
                // accumulate value of derivative
                val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l];
                // accumulate value of derivative
                val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l];
                //Output::print<2>(l, "  ", edgeRefData.xderiv_refNeigh1D[i][l],
                //                 "  ", edgeRefData.yderiv_refNeigh1D[i][l],
                //                 "  ", edgeRefData.FEFunctValuesNeigh[l]);
              }                                   // endfor l
              edgeRefData.xderiv_Neigh1D[i] = val[0];          // for k-th
              edgeRefData.yderiv_Neigh1D[i] = val[1];          // for k-th
              edgeRefData.xyval_Neigh1D[i] = val[2];           // for k-th
            }                                     // endfor i

            TFEDatabase2D::GetOrigFromRef(RefTransNeigh, N_QuadraturePoints1D, edgeRefData.xi1DNeigh.data(),
                                          edgeRefData.eta1DNeigh.data(),
                                          edgeRefData.X1DNeigh.data(), edgeRefData.Y1DNeigh.data(), absdet1D);

            jump = 0.0;
            double absdetjk1D = hE / 2.0;
            for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)            // compute jump
            {
              if ((std::abs(edgeData.XEdge1D[edgeIdx][i] - edgeRefData.X1DNeigh[i]) + std::abs(edgeData.YEdge1D[edgeIdx][i] - edgeRefData.Y1DNeigh[i])) > 1e-8)
              {
                cout << " wrong quad points 0 " << edgeData.XEdge1D[edgeIdx][i] << " , " << edgeData.YEdge1D[edgeIdx][i] << "   " << edgeRefData.X1DNeigh[i] << " , " << edgeRefData.Y1DNeigh[i] << endl;
              }
              if (check_cont)
              {
                //if (std::abs(edgeRefData.xyval_Neigh1D[i] - edgeData.xyval_1D[edgeIdx][i]) > 1e-8)
                //{
                //  cout << " i " << i << " valc " << edgeData.xyval_1D[edgeIdx][i] << " neighc " << edgeRefData.xyval_Neigh1D[i] << endl;
                //}
              }
              double e1 = coeff[0] * ((edgeData.xderiv_1D[edgeIdx][i] - edgeRefData.xderiv_Neigh1D[i]) * nx + (edgeData.yderiv_1D[edgeIdx][i] - edgeRefData.yderiv_Neigh1D[i]) * ny);
              //Output::print<2>(i, " jumpx ", edgeData.xderiv_1D[edgeIdx][i],
              //                 " ", edgeRefData.xderiv_Neigh1D[i]);
              double w = weights1D[i] * absdetjk1D;
              jump += w * e1 * e1;                     // integral on the edge
            }
            //Output::print<2>("jump ", jump);
            beta[0] = hE;               // weight for H^1 estimator
            beta[1] = hE * hE * hE;     // weight for L^2 estimator
            beta[2] = hE / coeff[0];    // weight for energy norm estimator
            if(TDatabase::ParamDB->INTERNAL_COERCIVITY > 0)
            {
              double w = 1 / std::sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
              if (w < beta[2])
                  beta[2] = w;
            }
            //beta[2] *= 2.0;
            if(24.0 / linfb < beta[2])
            {
              beta[4] = 24.0 / linfb;
            }
            else
            {
              beta[4] = beta[2];
            }
            /*beta[4] = 24;
            if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
            {
              linfb = std::sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY) * std::sqrt(coeff[0]);
              if (1.0/linfb < beta[4])
                beta[4] = 1.0/linfb;
            }
            linfb = hE / coeff[0];
            if (linfb < beta[4])
              beta[4] = linfb;
            */
            beta[5] = 1.0;
            beta[5] = alpha[5] * hE / (4.0 * meas);
            /*for (int i = 1; i < N_estimators; i++)
                estimated_error[i] += beta[i - 1] * jump / 2.0;*/
            result += beta[int(estimatorType) - 1] * jump / 2.0;
          }                                       // end neighbour is member of the collection
        }
      }
    }
  }
  // no delta_K, i.e., no SUPG
  if (delta_K == 4711 && (int(estimatorType) == 5 || int(estimatorType) == 6))
  {
    return 0;
  }
  return result;
}

template<int d>
bool CDErrorEstimator<d>::handleJump_BoundaryEdge(
  double *result, const Example_CD &example, const int estimatorType,
  const int N_QuadraturePoints1D, const double *const &weights1D,
  const CDErrorEstimator::JointData &edgeData,
  const double meas, const double *coeff, double linfb,
  const std::vector<double> &alpha, int edgeIdx, const TJoint *joint) const
{
  // vector holding the weights corresponding to the estimators, initialized with 0
  std::vector<double> beta(N_CD2D_ESTIMATOR_TYPES - 1, 0);
  // the edge
  const TBoundEdge *bdryEdge = (const TBoundEdge *) joint;
  // get boundary component
  const TBoundComp2D *boundComp = bdryEdge->GetBoundComp();
  // parameter interval
  double t0, t1;
  bdryEdge->GetParameters(t0, t1);
  // boundary id
  int bdryId = boundComp->GetID();
  // type of boundary condition
  BoundCond Cond0;
  example.get_bc(0)(bdryId, (t0 + t1) / 2.0, Cond0);
  double jump = 0;
  // at midpoint of boundary
  switch (Cond0)
  {
    case DIRICHLET:
    {
      // no error
      if(TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS)
      {
        return false;
      }
      break;
    }
    case NEUMANN:
    {
      double x0, x1, y0, y1;
      // coordinates at begin of parameter interval
      bdryEdge->GetXYofT(t0, x0, y0);
      // coordinates at end of parameter interval
      bdryEdge->GetXYofT(t1, x1, y1);
      // outer normal vector
      double nx = y1 - y0;
      double ny = x0 - x1;
      // length of edge
      double hE = std::sqrt(nx * nx + ny * ny);
      // normalized normal vector
      nx /= hE;
      ny /= hE;

      // compute difference to Neumann condition
      for(int i = 0; i < N_QuadraturePoints1D; i++)
      {
        x0 = edgeData.XEdge1D[edgeIdx][i];
        y0 = edgeData.YEdge1D[edgeIdx][i];
        // coordinates at quadrature points
        bdryEdge->GetTofXY(x0, y0, t0);
        // Neumann data
        double neumann_data;
        example.get_bd(0)(bdryId, t0, neumann_data);
        double e1 = coeff[0] * (edgeData.xderiv_1D[edgeIdx][i] * nx
                              + edgeData.yderiv_1D[edgeIdx][i] * ny)
                    - neumann_data;
        double w = weights1D[i] * hE / 2.0;
        // integral on the edge
        jump += w * e1 * e1;
      }
      *result += jump * compute_estimator_weight(estimatorType, hE, coeff[0],
                                                 linfb, alpha, meas);
      break;
    }
    case ROBIN:
      // TODO: Robin jump!
      break;
    case SLIP:
    case FREESURF:
    case SLIP_FRICTION_PENETRATION_RESISTANCE:
    case INTERFACE:
    case SUBDOMAIN_INTERFACE:
    case SUBDOMAIN_HALOBOUND:
    case DIRICHLET_WEAK:
    {
      ErrThrow("Only few BC implementation done ");
    }
  }
  return true;
}

template<int d>
void CDErrorEstimator<d>::info()
{
  if(this->currentCollection == nullptr)
  {
    Output::print("CDErrorEstimator<", d, ">: the method 'estimate' has not "
                  "been called yet, therefore there is no further information "
                  "available on this object.");
    return;
  }
  Output::print("CDErrorEstimator<", d, "> for ",
                this->currentCollection->GetN_Cells(), " cells");
  Output::dash("maximal_local_error: ", this->maximal_local_error);
  std::stringstream estimated_global_errors;
  for(int i = 0; i < N_CD2D_ESTIMATOR_TYPES; ++i)
  {
    estimated_global_errors << std::to_string(this->estimated_global_error[i]);
    estimated_global_errors << "  ";
  }
  Output::dash("estimated global errors: ", estimated_global_errors.str());
}


std::ostream &operator<<(std::ostream &os, CDErrorEstimatorType type)
{
  switch (type)
  {
    case CDErrorEstimatorType::GradientIndicator:
      os << "Gradient indicator";
      break;
    case CDErrorEstimatorType::H1_ResidualEstimator:
      os << "H1 residual estimator";
      break;
    case CDErrorEstimatorType::L2_ResidualEstimator:
      os << "L2 residual estimator";
      break;
    case CDErrorEstimatorType::Energy_ResidualEstimatorQuasiRobust:
      os << "Energy norm + dual norm residual estimator, Verfürth 2005";
      break;
    case CDErrorEstimatorType::Energy_ResidualEstimatorWithoutJumps:
      os << "Energy norm estimator without jumps";
      break;
    case CDErrorEstimatorType::SUPG_Upper:
      os << "SUPG estimator John/Novo, upper estimate";
      break;
    case CDErrorEstimatorType::SUPG_Lower:
      os << "SUPG estimator John/Novo, lower estimate";
      break;
    default:
      ErrThrow("unknown CDErrorEstimatorType");
  }
  os << " (database value = " << int(type) << ")";
  return os;
}


#ifdef __3D__
//template class CDErrorEstimator<3>;
#else
template class CDErrorEstimator<2>;
#endif
