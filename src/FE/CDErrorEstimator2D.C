//
// Created by Moritz Hoffmann on 05/07/15.
//

#include <AuxParam2D.h>
#include <Database.h>
#include <CD2DErrorEstimator.h>
#include <CDErrorEstimator2D.h>
#include <ConvDiff.h>
#include <FEDatabase2D.h>
#include <array>
#include <cmath>
#include <BoundEdge.h>

#define DEBUG_COMPARE_RESULTS_WITH_OLD_CODE 0

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
struct CDErrorEstimator2D::EdgeData {

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

    EdgeData() {
        // initialize structures holding values and derivatives on reference edge
        xieta_ref1D_data.resize(N_BaseFuncts2D_loc * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc * 3);
        {
            // back xietaval_ref1D, xideriv_ref1D, etaderiv_ref1D by xieta_ref1D_data
            auto *ptr = xieta_ref1D_data.data();
            xietaval_ref1D[0][0][0] = &ptr[0];
            xideriv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D_loc * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc];
            etaderiv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D_loc * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc * 2];
            for (size_t ii = 0; ii < N_BaseFuncts2D_loc; ii++) {
                for (size_t jj = 0; jj < 4; jj++) {
                    for (size_t kk = 0; kk < MaxN_QuadPoints_1D_loc; kk++) {
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
struct CDErrorEstimator2D::EdgeRefData {
private:
    std::vector<double> refNeigh1D_data;
    unsigned int refDataQuadPoints1D;
public:
    // these data structures depend on the maximal number of 2d base functions
    unsigned int max_n_base_funct_2d;
    std::array<std::array<double *, MaxN_QuadPoints_1D>, MaxN_BaseFunctions2D> xietaval_refNeigh1D;
    std::array<std::array<double *, MaxN_QuadPoints_1D>, MaxN_BaseFunctions2D> xideriv_refNeigh1D;
    std::array<std::array<double *, MaxN_QuadPoints_1D>, MaxN_BaseFunctions2D> etaderiv_refNeigh1D;

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
    std::array<double *, MaxN_QuadPoints_1D> xyval_refNeigh1D{};
    std::array<double *, MaxN_QuadPoints_1D> xderiv_refNeigh1D{};
    std::array<double *, MaxN_QuadPoints_1D> yderiv_refNeigh1D{};
    std::vector<double> FEFunctValuesNeigh;

    EdgeRefData(unsigned int max_n_base_funct_2d) : refDataQuadPoints1D{0}, max_n_base_funct_2d{max_n_base_funct_2d} {
        unsigned int k = MaxN_BaseFunctions2D * MaxN_QuadPoints_1D * max_n_base_funct_2d;
        refNeigh1D_data = std::vector<double>(3 * k);
        xietaval_refNeigh1D[0][0] = &refNeigh1D_data[0];
        xideriv_refNeigh1D[0][0] = &refNeigh1D_data[0] + k;
        etaderiv_refNeigh1D[0][0] = &refNeigh1D_data[0] + 2 * k;
        for (unsigned int i = 0; i < MaxN_BaseFunctions2D; i++) {
            unsigned int n = i * MaxN_BaseFunctions2D * max_n_base_funct_2d;
            for (unsigned int j = 0; j < MaxN_QuadPoints_1D; j++) {
                unsigned int l = n + j * max_n_base_funct_2d;
                xietaval_refNeigh1D[i][j] = xietaval_refNeigh1D[0][0] + l;
                xideriv_refNeigh1D[i][j] = xideriv_refNeigh1D[0][0] + l;
                etaderiv_refNeigh1D[i][j] = etaderiv_refNeigh1D[0][0] + l;
            }
        }
    }

    void updateQuadPointData(unsigned int refDataQuadPoints1D) {
        if (this->refDataQuadPoints1D != refDataQuadPoints1D) {
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
            for (size_t i = 0; i < refDataQuadPoints1D; i++) {
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
namespace {
    // returns the maximal number of used base functions 2D per element for a fe space
    unsigned int getMaxN_BaseFunctions2D(const TFESpace2D &fe_space) {
        unsigned int MaxN_BaseFunctions2D_loc = 0;
        // # used finite elements
        unsigned int n{(unsigned int) fe_space.GetN_UsedElements()};
        // used finite elements
        FE2D *UsedElements = fe_space.GetUsedElements();
        // for all finite elements
        for (unsigned int j = 0; j < n; j++) {
            // if number of base functions is higher, update value
            unsigned int k = (unsigned int) TFEDatabase2D::GetN_BaseFunctFromFE2D(UsedElements[j]);
            if (k > MaxN_BaseFunctions2D_loc) {
                MaxN_BaseFunctions2D_loc = k;
            }
        }
        return MaxN_BaseFunctions2D_loc;
    }

    // this returns weights for the cell residuals and calculates linfb (pass by reference)
    // alpha[0] for H^1 estimator
    // alpha[1] for L^2 estimator
    // alpha[2] for energy norm estimator
    // alpha[3] for energy norm estimator without jumps
    // alpha[4| for the SUPG norm, version 1, Novo
    // alpha[5] for the SUPG norm, version 2, Novo
    std::vector<double> getCellWeights(const double hK, const double *coeff, const double hK_tilde, double &linfb, double &delta_K) {
        std::vector<double> alpha;
        alpha.resize(6);
        alpha[0] = hK * hK;
        alpha[1] = hK * hK * hK * hK;
        alpha[2] = hK * hK / coeff[0];
        if (TDatabase::ParamDB->INTERNAL_COERCIVITY > 0) {
            if (1.0 / TDatabase::ParamDB->INTERNAL_COERCIVITY < alpha[2]) {
                // update weight for energy norm estimator
                alpha[2] = 1.0 / TDatabase::ParamDB->INTERNAL_COERCIVITY;
            }
        }
        alpha[3] = alpha[2];
        alpha[4] = alpha[3];
        linfb = fabs(coeff[1]);
        if (fabs(coeff[2]) > linfb)
            linfb = fabs(coeff[2]);
        if (TDatabase::ParamDB->DISCTYPE == SUPG) {
            // compute stabilization parameter
            delta_K = Compute_SDFEM_delta(hK, coeff[0], coeff[1], coeff[2], coeff[3], linfb);
            if (alpha[4] > 24 * delta_K) {
                alpha[4] = 24 * delta_K;
            }
            // second contribution
            alpha[4] += 24 * delta_K;
        }
        alpha[5] = hK;
        if (TDatabase::ParamDB->DISCTYPE == SUPG) {
            alpha[5] = hK * hK / (3 * sqrt(10.0) * coeff[0]);
            if (coeff[3] > 0) {
                if (1.0 / coeff[3] < alpha[5]) {
                    alpha[5] = 1.0 / coeff[3];
                }
            }
            // compute stabilization parameter
            if (hK_tilde / (sqrt(2.0) * linfb) < alpha[5]) {
                alpha[5] = hK / (sqrt(2.0) * linfb);
            }

            alpha[5] *= 2 * alpha[5];
        }
        return alpha;
    }

}

void CDErrorEstimator2D::estimate(const std::vector<MultiIndex2D> &derivatives, const TFEFunction2D &fe_function2D) {

    // remove old eta_K
    if (eta_K && eta_K != nullptr) delete[] eta_K;

    // initialization
    TCollection *coll = fe_function2D.GetFESpace2D()->GetCollection();
    // this call is obligatory,
    // otherwise the collection is unknown to the refinement strategy
    setCollection(coll);

    eta_K = new double[coll->GetN_Cells()];

    // array of pointers holding quad point -> derivatives
    std::vector<double *> derivativesPerQuadPoint;

    // actual data for quad point -> derivatives in a 1D structure
    std::vector<double> derivativesPerQuadPointData;
    derivativesPerQuadPointData.resize(MaxN_QuadPoints_2D * derivatives.size());

    // map this structure into the 2D array
    for (unsigned int j = 0; j < MaxN_QuadPoints_2D; j++) {
        derivativesPerQuadPoint.push_back(&derivativesPerQuadPointData[j * derivatives.size()]);
    }

    // vector containing quadpoint index -> coefficients
    std::vector<double *> coefficientsPerQuadPoint;

    const unsigned int stride = 20;
    std::vector<double> coefficientVectorData;
    coefficientVectorData.resize(MaxN_QuadPoints_2D * stride);
    for (unsigned int j = 0; j < MaxN_QuadPoints_2D; j++) {
        coefficientsPerQuadPoint.push_back(&coefficientVectorData[j * stride]);
    }


    // the FE basis functions
    BaseFunct2D *baseFunctions = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
    // number of basis functions
    int *n_baseFunctions = TFEDatabase2D::GetN_BaseFunctFromFE2D();

    if (DEBUG_COMPARE_RESULTS_WITH_OLD_CODE) {
        std::cout << " ---running old code--- " << std::endl;
        TAuxParam2D aux;
        int estimator = int(estimatorType);
        const TFESpace2D *fe_space = fe_function2D.GetFESpace2D();
        TCD2DErrorEstimator CDErrorEstimator(estimator, const_cast<TFEFunction2D *>(&fe_function2D), TDatabase::ParamDB->ERROR_CONTROL);
        CDErrorEstimator.GetErrorEstimate((int) derivatives.size(), const_cast<MultiIndex2D *>(derivatives.data()),
                                          example2D.get_coeffs(), const_cast<BoundCondFunct2D **>(example2D.get_bc()),
                                          const_cast<BoundValueFunct2D **>(example2D.get_bd()), &aux,
                                          1, const_cast<TFESpace2D **>(&fe_space), eta_K, &maximal_local_error,
                                          &estimated_global_error[0]);
        std::cout << " ---done running old code--- " << std::endl;
    }

    // the fe space
    const TFESpace2D *fe_space = fe_function2D.GetFESpace2D();
    // attributes of the reference cell (coordinates and determinant)
    double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];

    EdgeData edgeData{};
    // initialize and resize vectors / arrays accordingly
    const unsigned int max_n_base_funct_2d = getMaxN_BaseFunctions2D(fe_space[0]);
    EdgeRefData edgeRefData{max_n_base_funct_2d};
    {
        edgeData.xietaval_ref1D.resize(max_n_base_funct_2d);
        edgeData.xideriv_ref1D.resize(max_n_base_funct_2d);
        edgeData.etaderiv_ref1D.resize(max_n_base_funct_2d);
        edgeData.xyval_ref1D.resize(max_n_base_funct_2d);
        edgeData.xderiv_ref1D.resize(max_n_base_funct_2d);
        edgeData.yderiv_ref1D.resize(max_n_base_funct_2d);
        for (unsigned int x = 0; x < 4; x++) {
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
        TBaseCell *cell, *neigh;
        // do for all mesh cells
        for (size_t i = 0; i < coll->GetN_Cells(); i++) {
            // on the finest level
            cell = coll->GetCell((int) i);
            // for all edges
            for (int j = 0; j < cell->GetN_Edges(); j++) {
                // neighbour cell
                neigh = cell->GetJoint(j)->GetNeighbour(cell);
                // set clipboard to -1
                if (neigh) neigh->SetClipBoard(-1);
            }
            cell->SetClipBoard(-1);
        }
        // now non finest neighbours of finest cells have clipboard -1
        for (size_t i = 0; i < coll->GetN_Cells(); i++) {
            // set clipboard of cells on finest
            cell = coll->GetCell((int) i);
            cell->SetClipBoard((int) i);
        }
    }

    int N_Points;
    double *xi, *eta, *weights;
    // we do need second derivatives on the reference cell
    bool secondDerivatives[] = {true};

    // ########################################################################
    // calculate values of base functions and derivatives on ref element
    // ########################################################################

    auto n_used_elements = fe_space->GetN_UsedElements();

    // store used finite elements
    std::unique_ptr<FE2D> UsedElements{new FE2D[n_used_elements]};
    {
        std::vector<int> Used = std::vector<int>(N_FEs2D);
        auto n = fe_space->GetN_UsedElements();
        auto use = fe_space->GetUsedElements();
        for (auto j = 0; j < n; j++) {
            FE2D CurrentElement = use[j];
            Used[CurrentElement] = 1;
        }
        int N_UsedElements = 0;
        for (auto i = 0; i < N_FEs2D; i++)
            if (Used[i]) N_UsedElements++;

        auto ptr = UsedElements.get();
        size_t j = 0;
        for (auto i = 0; i < N_FEs2D; i++) {
            if (Used[i]) {
                ptr[j] = (FE2D) i;
                j++;
            }
        }
        if (n_used_elements != N_UsedElements) {
            std::cerr << "neineineineineineineinein" << std::endl;
        }
    }

    int N_Points1D;
    for (auto i = 0; i < n_used_elements; i++) {
        FE2D usedElement = UsedElements.get()[i];//fe_space->GetUsedElements()[i];
        QuadFormula1D qf = TFEDatabase2D::GetQFLineFromDegree(2 * TFEDatabase2D::GetPolynomialDegreeFromFE2D(usedElement));
        TQuadFormula1D *quadFormula = TFEDatabase2D::GetQuadFormula1D(qf);
        double *weights1D, *zeta;
        quadFormula->GetFormulaData(N_Points1D, weights1D, zeta);
        if (N_Points1D > edgeData.MaxN_QuadPoints_1D_loc) {
            Output::print<1>("CD2DErrorEstimator: too many 1D quadrature points ", N_Points1D);
            Output::print<1>("Increase  MaxN_QuadPoints_1D_loc !!!");
            exit(4711);
        }

        {
            for (auto edge = 0; edge < 4; edge++) {
                edgeData.xyval_1D[edge].resize((unsigned long) N_Points1D);
                edgeData.xderiv_1D[edge].resize((unsigned long) N_Points1D);
                edgeData.yderiv_1D[edge].resize((unsigned long) N_Points1D);
            }
        }

        BaseFunct2D baseFunct2D = baseFunctions[usedElement];
        TBaseFunct2D *bf = TFEDatabase2D::GetBaseFunct2D(baseFunct2D);// get base functions
        BF2DRefElements bf2Drefelements = bf->GetRefElement();
        auto used_ptr = UsedElements.get();
        int baseFunct_loc = 0;
        {
            unsigned int j;
            for (j = 0; j < n_used_elements; j++) {
                if ((int) baseFunct2D == (int) used_ptr[j]) {
                    break;
                }
            }
            baseFunct_loc = j;
        }

        // compute coordinates of line quadrature points in reference cell
        switch (bf2Drefelements) {
            // quadrilateral cell
            case BFUnitSquare : {
                // edge 0
                for (auto j = 0; j < N_Points1D; j++) {
                    edgeData.xi1D[baseFunct_loc][0][j] = zeta[j];
                    edgeData.eta1D[baseFunct_loc][0][j] = -1;
                    bf->GetDerivatives(D00, zeta[j], -1, edgeData.xietaval_ref1D[baseFunct_loc][0][j]);
                    bf->GetDerivatives(D10, zeta[j], -1, edgeData.xideriv_ref1D[baseFunct_loc][0][j]);
                    bf->GetDerivatives(D01, zeta[j], -1, edgeData.etaderiv_ref1D[baseFunct_loc][0][j]);
                }
                // edge 1
                for (auto j = 0; j < N_Points1D; j++) {
                    edgeData.xi1D[baseFunct_loc][1][j] = 1;
                    edgeData.eta1D[baseFunct_loc][1][j] = zeta[j];
                    bf->GetDerivatives(D00, 1, zeta[j], edgeData.xietaval_ref1D[baseFunct_loc][1][j]);
                    bf->GetDerivatives(D10, 1, zeta[j], edgeData.xideriv_ref1D[baseFunct_loc][1][j]);
                    bf->GetDerivatives(D01, 1, zeta[j], edgeData.etaderiv_ref1D[baseFunct_loc][1][j]);
                }
                // edge 2
                for (auto j = 0; j < N_Points1D; j++) {
                    edgeData.xi1D[baseFunct_loc][2][j] = -zeta[j];
                    edgeData.eta1D[baseFunct_loc][2][j] = 1;
                    bf->GetDerivatives(D00, -zeta[j], 1, edgeData.xietaval_ref1D[baseFunct_loc][2][j]);
                    bf->GetDerivatives(D10, -zeta[j], 1, edgeData.xideriv_ref1D[baseFunct_loc][2][j]);
                    bf->GetDerivatives(D01, -zeta[j], 1, edgeData.etaderiv_ref1D[baseFunct_loc][2][j]);
                }
                // edge 3
                for (auto j = 0; j < N_Points1D; j++) {
                    edgeData.xi1D[baseFunct_loc][3][j] = -1;
                    edgeData.eta1D[baseFunct_loc][3][j] = -zeta[j];
                    bf->GetDerivatives(D00, -1, -zeta[j], edgeData.xietaval_ref1D[baseFunct_loc][3][j]);
                    bf->GetDerivatives(D10, -1, -zeta[j], edgeData.xideriv_ref1D[baseFunct_loc][3][j]);
                    bf->GetDerivatives(D01, -1, -zeta[j], edgeData.etaderiv_ref1D[baseFunct_loc][3][j]);
                }
                break;
            }

                // triangular cell
            case BFUnitTriangle : {
                // edge 0
                for (auto j = 0; j < N_Points1D; j++) {
                    edgeData.xi1D[baseFunct_loc][0][j] = (zeta[j] + 1) / 2;
                    edgeData.eta1D[baseFunct_loc][0][j] = 0;
                    bf->GetDerivatives(D00, (zeta[j] + 1) / 2, 0, edgeData.xietaval_ref1D[baseFunct_loc][0][j]);
                    bf->GetDerivatives(D10, (zeta[j] + 1) / 2, 0, edgeData.xideriv_ref1D[baseFunct_loc][0][j]);
                    bf->GetDerivatives(D01, (zeta[j] + 1) / 2, 0, edgeData.etaderiv_ref1D[baseFunct_loc][0][j]);
                }
                // edge 1
                for (auto j = 0; j < N_Points1D; j++) {
                    edgeData.xi1D[baseFunct_loc][1][j] = (-zeta[j] + 1) / 2;
                    edgeData.eta1D[baseFunct_loc][1][j] = (zeta[j] + 1) / 2;
                    bf->GetDerivatives(D00, (-zeta[j] + 1) / 2, (zeta[j] + 1) / 2, edgeData.xietaval_ref1D[baseFunct_loc][1][j]);
                    bf->GetDerivatives(D10, (-zeta[j] + 1) / 2, (zeta[j] + 1) / 2, edgeData.xideriv_ref1D[baseFunct_loc][1][j]);
                    bf->GetDerivatives(D01, (-zeta[j] + 1) / 2, (zeta[j] + 1) / 2, edgeData.etaderiv_ref1D[baseFunct_loc][1][j]);
                }
                // edge 2
                for (auto j = 0; j < N_Points1D; j++) {
                    edgeData.xi1D[baseFunct_loc][2][j] = 0;
                    edgeData.eta1D[baseFunct_loc][2][j] = (-zeta[j] + 1) / 2;
                    bf->GetDerivatives(D00, 0, (-zeta[j] + 1) / 2, edgeData.xietaval_ref1D[baseFunct_loc][2][j]);
                    bf->GetDerivatives(D10, 0, (-zeta[j] + 1) / 2, edgeData.xideriv_ref1D[baseFunct_loc][2][j]);
                    bf->GetDerivatives(D01, 0, (-zeta[j] + 1) / 2, edgeData.etaderiv_ref1D[baseFunct_loc][2][j]);
                }
                break;
            }
        }
    }


    // for each cell
    for (unsigned int cellIdx = 0; cellIdx < coll->GetN_Cells(); cellIdx++) {
        // the cell
        TBaseCell *cell = coll->GetCell(cellIdx);

        // fe2d element on cell
        FE2D element = fe_space->GetFE2D(cellIdx, cell);
        // get reference transformation RefTrans2D
        RefTrans2D refTrans2D = TFEDatabase2D::GetOrig(1, &element, coll, cell, secondDerivatives, N_Points, xi, eta, weights, X, Y, AbsDetjk);

        // get fe function values for every basis function of the element
        double FEFunctValues[n_baseFunctions[element]];
        for (unsigned int l = 0; l < n_baseFunctions[element]; l++) {
            // fe values of dofs
            int *DOF = fe_space->GetGlobalNumbers() + fe_space->GetBeginIndex()[cellIdx];
            FEFunctValues[l] = fe_function2D.GetValues()[DOF[l]];
        }

        // compute values for all derivatives in all quadrature points in original mesh cell
        // for all derivatives
        for (unsigned int k = 0; k < derivatives.size(); k++) {
            // get values in original cell by dof of current mesh cell
            double **OrigFEValues = TFEDatabase2D::GetOrigElementValues(baseFunctions[element], derivatives[k]);

            // for all quadrature points
            for (unsigned long j = 0; j < N_Points; j++) {
                // value in original cell
                double *Orig = OrigFEValues[j];
                double value = 0;
                // for all basis functions
                for (unsigned int l = 0; l < n_baseFunctions[element]; l++) {
                    // accumulate value of derivative in point j
                    value += FEFunctValues[l] * Orig[l];
                }
                // for k-th derivative
                derivativesPerQuadPoint[j][k] = value;
            }
        }

        // problem's coefficients
        if (example2D.get_coeffs()) {
            example2D.get_coeffs()(N_Points, X, Y, nullptr, &coefficientsPerQuadPoint[0]);
        }

        // 1D quadrature formula
        QuadFormula1D qf = TFEDatabase2D::GetQFLineFromDegree(2 * TFEDatabase2D::GetPolynomialDegreeFromFE2D(element));
        TQuadFormula1D *quadFormula = TFEDatabase2D::GetQuadFormula1D(qf);
        int N_QuadraturePoints1D;
        double *weights1D, *zeta;
        quadFormula->GetFormulaData(N_QuadraturePoints1D, weights1D, zeta);


        // update data base
        TFEDatabase2D::GetBaseFunct2DFromFE2D(element)->MakeRefElementData(qf);


        unsigned int BaseFunct_loc = 0;
        {
            unsigned int j = 0;
            for (j = 0; j < n_used_elements; j++) {
                if ((int) baseFunctions[element] == (int) fe_space->GetUsedElements()[j]) {
                    break;
                }
            }
            BaseFunct_loc = j;
        }

        // loop over all edges of cell
        for (auto edgeIdx = 0; edgeIdx < cell->GetN_Edges(); edgeIdx++) {
            // get original coordinates of edge quad. points
            TFEDatabase2D::GetOrigFromRef(refTrans2D, N_QuadraturePoints1D,
                                          edgeData.xi1D[BaseFunct_loc][edgeIdx],
                                          edgeData.eta1D[BaseFunct_loc][edgeIdx],
                                          edgeData.XEdge1D[edgeIdx].data(),
                                          edgeData.YEdge1D[edgeIdx].data(),
                                          edgeData.AbsDetjk1D[edgeIdx].data());
            // get values and derivatives in original cell
            for (unsigned int k = 0; k < N_QuadraturePoints1D; k++) {
                TFEDatabase2D::GetOrigValues(refTrans2D, edgeData.xi1D[BaseFunct_loc][edgeIdx][k],
                                             edgeData.eta1D[BaseFunct_loc][edgeIdx][k],
                                             TFEDatabase2D::GetBaseFunct2D(baseFunctions[element]),
                                             coll, (TGridCell *) cell,
                                             edgeData.xietaval_ref1D[BaseFunct_loc][edgeIdx][k],
                                             edgeData.xideriv_ref1D[BaseFunct_loc][edgeIdx][k],
                                             edgeData.etaderiv_ref1D[BaseFunct_loc][edgeIdx][k],
                                             &edgeData.xyval_ref1D[edgeIdx][k][0],
                                             &edgeData.xderiv_ref1D[edgeIdx][k][0],
                                             &edgeData.yderiv_ref1D[edgeIdx][k][0]);
            }
            double val[3];
            // for all quadrature points
            for (auto k = 0; k < N_QuadraturePoints1D; k++) {
                val[0] = val[1] = val[2] = 0;
                // for all basis functions
                for (auto l = 0; l < n_baseFunctions[element]; l++) {
                    // accumulate value of derivative
                    val[0] += FEFunctValues[l] * edgeData.xyval_ref1D[edgeIdx][k][l];
                    // accumulate value of derivative
                    val[1] += FEFunctValues[l] * edgeData.xderiv_ref1D[edgeIdx][k][l];
                    // accumulate value of derivative
                    val[2] += FEFunctValues[l] * edgeData.yderiv_ref1D[edgeIdx][k][l];
                }
                /*std::cout << "cell=" << cellIdx << ", " << "xyval_1d[" << edgeIdx << "][" << k << "]=" << val[0] << std::endl;
                for(auto l = 0; l < n_baseFunctions[element]; l++) {
                    std::cout << "\tl="<<l<<", Fe_val="<<FEFunctValues[l] <<", " << "xyval="<<edgeData.xyval_ref1D[edgeIdx][k][l]<<std::endl;
                }*/
                edgeData.xyval_1D[edgeIdx][k] = val[0];
                edgeData.xderiv_1D[edgeIdx][k] = val[1];
                edgeData.yderiv_1D[edgeIdx][k] = val[2];
            }
        }

        if ((TDatabase::ParamDB->DISCTYPE == SDFEM) || (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM)) {
            TDatabase::ParamDB->INTERNAL_LOCAL_DOF = cellIdx;
            int cellEdges = cell->GetN_Edges();
            for (unsigned int ij = 0; ij < cellEdges; ij++) {
                TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
                TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
            }
            if (cellEdges == 3) {
                TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
            }
            TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
        }


        double result = calculateEtaK(cell, fe_function2D, *coll, derivativesPerQuadPoint, AbsDetjk, weights,
                                      coefficientsPerQuadPoint, N_Points, (const unsigned int) N_QuadraturePoints1D, weights1D, edgeData, zeta, edgeRefData);

        if (DEBUG_COMPARE_RESULTS_WITH_OLD_CODE && eta_K[cellIdx] - result > 1e-15) {
            cerr << "something wrong with the eta_k's" << endl;
        }
        eta_K[cellIdx] = result;
        estimated_global_error[int(estimatorType)] += result;
        maximal_local_error = result > maximal_local_error ? result : maximal_local_error;
    }
    estimated_global_error[int(estimatorType)] = sqrt(estimated_global_error[int(estimatorType)]);
    maximal_local_error = sqrt(maximal_local_error);
}

double CDErrorEstimator2D::calculateEtaK(TBaseCell *cell, const TFEFunction2D &fe_function, TCollection &coll, std::vector<double *> &derivativesPerQuadPoint,
                                         double (&AbsDetjk)[MaxN_QuadPoints_2D], double *(&weights),
                                         std::vector<double *> &coefficientsPerQuadPoint,
                                         int n_quadrature_points, const unsigned int N_QuadraturePoints1D, double *(&weights1D), const EdgeData &edgeData, const double *zeta, EdgeRefData &edgeRefData) const {

    // constants
    int ee_verbose = TDatabase::ParamDB->SC_VERBOSE;

    // the (current) boundary condition, where needed retrieved by example2D.get_bc()
    BoundCond Cond0;

    // SUPG parameter
    double delta_K = 0.0;

    const TFESpace2D *fe_space = fe_function.GetFESpace2D();

    if (TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS) {
        int N_Edges = cell->GetN_Edges();
        // loop over all edges of cell
        for (unsigned int j = 0; j < N_Edges; j++) {
            TVertex *ver0 = cell->GetVertex(j);
            // clip board contains information on the position of the vertex
            int w = ver0->GetClipBoard();
            // vertex not on the boundary
            if (w < -1e-8) continue;
            // component of boundary
            int comp = (int) floor(w + 1e-8);
            // parameter
            w -= comp;
            // get boundary condition
            //Output::print<1>(comp << " " << w);
            example2D.get_bc()[0](comp, w, Cond0);
            // Dirichlet
            if (Cond0 == DIRICHLET) return 0;
        }
    }

    double result = 0.0;
    switch (estimatorType) {
        case CDErrorEstimatorType::GradientIndicator: {
            // calculate the eta_K
            for (unsigned int j = 0; j < n_quadrature_points; j++) {
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
            break;
        }
        default: {
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
            if (TDatabase::ParamDB->ANSATZ_ORDER > 0)
                check_cont = 1;
            else
                check_cont = 0;

            /**
             * Calculate strong residual
             */
            const double hK = cell->GetDiameter();
            const int N_Edges = cell->GetN_Edges();
            for (unsigned int i = 0; i < N_Edges; i++) {
                TDatabase::ParamDB->INTERNAL_VERTEX_X[i] = cell->GetVertex(i)->GetX();
                TDatabase::ParamDB->INTERNAL_VERTEX_Y[i] = cell->GetVertex(i)->GetY();
            }
            if (N_Edges == 3) {
                TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
            }

            const double meas = cell->GetMeasure();
            double strong_residual = 0;
            // for all quadrature points
            for (unsigned int i = 0; i < n_quadrature_points; i++) {
                double *coeff = coefficientsPerQuadPoint[i];
                // all derivatives in quadrature points
                const double *deriv = derivativesPerQuadPoint[i];
                const double w = weights[i] * AbsDetjk[i];
                // strong residual
                double e1 = -coeff[0] * (deriv[3] + deriv[4]) + coeff[1] * deriv[0] + coeff[2] * deriv[1] + coeff[3] * deriv[2] - coeff[4];

                // time-dependent problem
                if (TDatabase::ParamDB->P15 == -4711.0) {
                    e1 += coeff[6];
                }
                // L^2 norm
                strong_residual += w * e1 * e1;
            }
            double *coeff = coefficientsPerQuadPoint[n_quadrature_points - 1];

            double hK_tilde = Mesh_size_in_convection_direction(hK, coeff[1], coeff[2]);

            // getting the cell weight in a vector [weight_0,...,weight_5], see description of the method
            double linfb;
            std::vector<double> alpha = getCellWeights(hK, coeff, hK_tilde, linfb, delta_K);
            double estimator_cell_weight = alpha[int(estimatorType) - 1];

            // this is the accordingly weighted strong residual without jumps
            result = estimator_cell_weight * strong_residual;

            /**
             * Calculate jumps
             */
            for (unsigned int edgeIdx = 0; edgeIdx < cell->GetN_Edges(); edgeIdx++) {
                TJoint *joint = cell->GetJoint(edgeIdx);

                // boundary edge
                if (joint->GetType() == BoundaryEdge || joint->GetType() == IsoBoundEdge) {
                    // This is a boundary edge
                    double boundaryJump = 0.0;
                    bool shouldNotAbortEstimation = CDErrorEstimator2D::handleJump_BoundaryEdge(&boundaryJump, example2D, int(estimatorType), N_QuadraturePoints1D, weights1D, edgeData, Cond0, meas, coeff, linfb, alpha, edgeIdx, joint);
                    if (shouldNotAbortEstimation) {
                        result += boundaryJump;
                    } else {
                        return 0;
                    }
                } else {
                    // this is an interior edge
                    TRefDesc *refdesc = cell->GetRefDesc();
                    // get refinement descriptor
                    const int *TmpEdVer;
                    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
                    // get vertices of edge edgeIdx
                    TVertex *ver0 = cell->GetVertex(TmpEdVer[2 * edgeIdx]);
                    TVertex *ver1 = cell->GetVertex(TmpEdVer[2 * edgeIdx + 1]);
                    // coordinates of edge "edgeIdx"
                    double cell_x0 = ver0->GetX();
                    double cell_y0 = ver0->GetY();
                    double cell_x1 = ver1->GetX();
                    double cell_y1 = ver1->GetY();
                    double jump = 0;
                    TBaseCell *neigh = joint->GetNeighbour(cell);
                    // compute normal
                    double nx = cell_y1 - cell_y0;
                    double ny = cell_x0 - cell_x1;
                    // length of edge
                    double hE = sqrt(nx * nx + ny * ny);
                    // normalized normal vector
                    nx /= hE;
                    ny /= hE;
                    if (ee_verbose > 1) {
                        std::cout << " A " << cell_x0 << " " << cell_y0 << std::endl;
                        std::cout << " B " << cell_x1 << " " << cell_y1 << std::endl;
                        std::cout << " n " << nx << " " << ny << std::endl;
                    }

                    /*************************************************************************/
                    /*  no neighbour, find neighbour of parent                               */
                    /*************************************************************************/
                    if (!neigh) {
                        // there is no neighbour on the same level
                        //  => finer cell in 1 regularity

                        // parent edge vertex, child edge, newEdgeOldEdge
                        const int *TmpEdVerParent, *TmpCE, *TmpoEnlE;
                        // maxlen of child edge
                        int MaxLen1;
                        // parent cell
                        TBaseCell *parent = cell->GetParent();
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

                        TJoint *parent_joint = parent->GetJoint(parent_edge);
                        // neighbour to parent
                        neigh = parent_joint->GetNeighbour(parent);
                        if (!neigh) {
                            std::cerr << "cell was no boundary cell and had no neighbour, thus its parent cell must have a neighbor" << std::endl;
                        }

                        // edge index of the parents neighbor
                        int neigh_edge = 0;
                        {
                            // iterate through the neighbours edges until the neighbor's neighbor is the parent cell itself again
                            while (neigh->GetJoint(neigh_edge)->GetNeighbour(neigh) != parent) neigh_edge++;
                        }
                        // vertices of edge
                        TVertex *ver2 = neigh->GetVertex(TmpEdVerParent[2 * neigh_edge]);
                        TVertex *ver3 = neigh->GetVertex(TmpEdVerParent[2 * neigh_edge + 1]);
                        // variable indicating if it is the first or the second part of the long edge, i.e.,
                        // comparison if ver1 (the second vertex of the cell) matches ver2 (first vertex of parent edge)
                        // or ver0 (the first vertex of the cell) matches ver3 (the second vertex of the parent edge)
                        int part = 0;
                        if (ver1 == ver2) {
                            // first part of long edge
                            part = -1;
                        } else if (ver0 == ver3) {
                            // second part of long edge
                            part = 1;
                        } else {
                            std::cerr << "\"2nd vertex of cell edge\" != \"1st vertex of parent edge\" and \"1st vertex of cell edge\" != \"2nd vertex of parent edge\"" << std::endl
                            << "this should not happen since the vertices are ordered" << std::endl;
                        }

                        // number of neighbour in iterator
                        int neigh_N_ = neigh->GetClipBoard();
                        if (neigh_N_ == -1) {
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

                        if (conform_grid) {
                            // compute coordinates of line quadrature
                            // points in reference cell
                            switch (bf2DrefelementsNeigh) {
                                case BFUnitSquare : {
                                    // edge 0
                                    if (neigh_edge == 0) {
                                        // for all quadrature points
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            edgeRefData.xi1DNeigh[i] = -zeta[i];
                                            edgeRefData.eta1DNeigh[i] = -1;
                                        }
                                    }
                                    // edge 1
                                    if (neigh_edge == 1) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            edgeRefData.xi1DNeigh[i] = 1;
                                            edgeRefData.eta1DNeigh[i] = -zeta[i];
                                        }
                                    }
                                    // edge 2
                                    if (neigh_edge == 2) {
                                        // for all quadrature points
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            edgeRefData.xi1DNeigh[i] = zeta[i];
                                            edgeRefData.eta1DNeigh[i] = 1;
                                        }
                                    }

                                    // edge 3
                                    if (neigh_edge == 3) {
                                        // for all quadrature points
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            edgeRefData.xi1DNeigh[i] = -1;
                                            edgeRefData.eta1DNeigh[i] = zeta[i];
                                        }
                                    }
                                    break;
                                }
                                case BFUnitTriangle :
                                    if (neigh_edge == 0) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            edgeRefData.xi1DNeigh[i] = (-zeta[i] + 1) / 2;
                                            edgeRefData.eta1DNeigh[i] = 0;
                                        }
                                    }
                                    if (neigh_edge == 1) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            edgeRefData.xi1DNeigh[i] = (zeta[i] + 1) / 2;
                                            edgeRefData.eta1DNeigh[i] = (-zeta[i] + 1) / 2;
                                        }
                                    }
                                    if (neigh_edge == 2) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            edgeRefData.xi1DNeigh[i] = 0;
                                            edgeRefData.eta1DNeigh[i] = (zeta[i] + 1) / 2;
                                        }
                                    }
                                    break;
                            }
                        } else {
                            // compute coordinates of line quadrature
                            // this is only for 1-regular triangulations
                            switch (bf2DrefelementsNeigh) {
                                // points in reference cell
                                case BFUnitSquare :
                                    // edge 0
                                    if (neigh_edge == 0) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            edgeRefData.xi1DNeigh[i] = (-zeta[i] + part) / 2;
                                            edgeRefData.eta1DNeigh[i] = -1;
                                        }
                                    }
                                    // edge 1
                                    if (neigh_edge == 1) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            edgeRefData.xi1DNeigh[i] = 1;
                                            edgeRefData.eta1DNeigh[i] = (-zeta[i] + part) / 2;
                                        }
                                    }
                                    // edge 2
                                    if (neigh_edge == 2) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            edgeRefData.xi1DNeigh[i] = (zeta[i] - part) / 2;
                                            edgeRefData.eta1DNeigh[i] = 1;
                                        }
                                    }
                                    // edge 3
                                    if (neigh_edge == 3) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            edgeRefData.xi1DNeigh[i] = -1;
                                            edgeRefData.eta1DNeigh[i] = (zeta[i] - part) / 2;
                                        }
                                    }
                                    break;

                                case BFUnitTriangle :
                                    if (neigh_edge == 0) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            if (part == -1) part = 0;
                                            edgeRefData.xi1DNeigh[i] = ((-zeta[i] + 1) / 2 + part) / 2;
                                            edgeRefData.eta1DNeigh[i] = 0;
                                        }
                                    }
                                    if (neigh_edge == 1) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            if (part == 1) part = 0;
                                            edgeRefData.xi1DNeigh[i] = ((zeta[i] + 1) / 2 - part) / 2;
                                            if (part == 0) part = 1;
                                            if (part == -1) part = 0;
                                            edgeRefData.eta1DNeigh[i] = ((-zeta[i] + 1) / 2 + part) / 2;
                                            if (part == 0) part = -1;
                                        }
                                    }
                                    if (neigh_edge == 2) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                            // for all quadrature points
                                            if (part == 1) part = 0;
                                            edgeRefData.xi1DNeigh[i] = 0;
                                            edgeRefData.eta1DNeigh[i] = ((zeta[i] + 1) / 2 - part) / 2;
                                        }
                                    }
                                    break;
                            }
                        }
                        if (ee_verbose > 1) {
                            for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                cout << "xiN " << edgeRefData.xi1DNeigh[i] << " etaN " << edgeRefData.eta1DNeigh[i] << endl;
                            }
                        }

                        // compute gradients in reference cell of the neighbour
                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
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
                        for (unsigned int i = 0; i < N_Neigh; i++) {
                            edgeRefData.FEFunctValuesNeigh[i] = Values[DOF[i]];
                            if (ee_verbose > 1) {
                                cout << " value " << edgeRefData.FEFunctValuesNeigh[i] << endl;
                            }
                        }

                        // get values and derivatives in original cell
                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
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
                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                            val[0] = val[1] = val[2] = 0;
                            // for all basis functions
                            for (neighIdx = 0; neighIdx < N_Neigh; neighIdx++) {
                                // accumulate value of derivative
                                val[0] += edgeRefData.FEFunctValuesNeigh[neighIdx] * edgeRefData.xderiv_refNeigh1D[i][neighIdx];
                                // accumulate value of derivative
                                val[1] += edgeRefData.FEFunctValuesNeigh[neighIdx] * edgeRefData.yderiv_refNeigh1D[i][neighIdx];
                                // accumulate value of derivative
                                val[2] += edgeRefData.FEFunctValuesNeigh[neighIdx] * edgeRefData.xyval_refNeigh1D[i][neighIdx];
                                if (ee_verbose > 1) {
                                    cout << neighIdx << "  " << edgeRefData.xderiv_refNeigh1D[i][neighIdx] << "  " << edgeRefData.yderiv_refNeigh1D[i][neighIdx] << "  " << edgeRefData.FEFunctValuesNeigh[neighIdx] << endl;
                                }
                            }
                            edgeRefData.xderiv_Neigh1D[i] = val[0];
                            edgeRefData.yderiv_Neigh1D[i] = val[1];
                            edgeRefData.xyval_Neigh1D[i] = val[2];
                        }

                        TFEDatabase2D::GetOrigFromRef(RefTransNeigh, N_QuadraturePoints1D, edgeRefData.xi1DNeigh.data(), edgeRefData.eta1DNeigh.data(), edgeRefData.X1DNeigh.data(), edgeRefData.Y1DNeigh.data(), absdet1D);
                        jump = 0.0;
                        double absdetjk1D = hE / 2.0;
                        // compute jump
                        for (size_t i = 0; i < N_QuadraturePoints1D; i++) {
                            if ((fabs(edgeData.XEdge1D[edgeIdx][i] - edgeRefData.X1DNeigh[i]) + fabs(edgeData.YEdge1D[edgeIdx][i] - edgeRefData.Y1DNeigh[i])) > 1e-8)
                                cout << " wrong quad points 1 " << edgeData.XEdge1D[edgeIdx][i] << " , " << edgeData.YEdge1D[edgeIdx][i] << "   " << edgeRefData.X1DNeigh[i] << " , " << edgeRefData.Y1DNeigh[i] << endl;
                            if (check_cont && fabs(edgeRefData.xyval_Neigh1D[i] - edgeData.xyval_1D[edgeIdx][i]) > 1e-8) {
                                cout << "quad points a " << edgeData.XEdge1D[edgeIdx][i] << " , " << edgeData.YEdge1D[edgeIdx][i] << endl;
                                cout << " i " << i << " vala " << edgeData.xyval_1D[edgeIdx][i] << " neigha " << edgeRefData.xyval_Neigh1D[i] << " " << fabs(edgeData.xyval_1D[edgeIdx][i] - edgeRefData.xyval_Neigh1D[i]) << endl;
                            }
                            double e1 = coeff[0] * ((edgeData.xderiv_1D[edgeIdx][i] - edgeRefData.xderiv_Neigh1D[i]) * nx + (edgeData.yderiv_1D[edgeIdx][i] - edgeRefData.yderiv_Neigh1D[i]) * ny);
                            if (ee_verbose > 1) {
                                cout << i << " jumpx " << edgeData.xderiv_1D[edgeIdx][i] << " " << edgeRefData.xderiv_Neigh1D[i] << endl;
                            }
                            double w = weights1D[i] * absdetjk1D;
                            jump += w * e1 * e1;                       // integral on the edge
                        }
                        if (ee_verbose > 1) {
                            cout << "jump " << jump << endl;
                        }

                        beta[0] = hE;                             // weight for H^1 estimator
                        beta[1] = hE * hE * hE;                       // weight for L^2 estimator
                        beta[2] = hE / coeff[0];                    // weight for energy norm estimator
                        if (TDatabase::ParamDB->INTERNAL_COERCIVITY > 0) {
                            double w = 1 / sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
                            if (w < beta[2]) beta[2] = w;
                        }
                        //beta[2] *= 2.0;
                        if (24.0 / linfb < beta[2]) {
                            beta[4] = 24.0 / linfb;
                        }
                        else {
                            beta[4] = beta[2];
                        }
                        /*
                        beta[4] = 24;
                        if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
                        {
                          linfb = sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY) * sqrt(coeff[0]);
                          if (1.0/linfb < beta[4])
                            beta[4] = 1.0/linfb;
                        }
                        linfb = hE / coeff[0];
                        if (linfb < beta[4])
                          beta[4] = linfb;*/
                        beta[5] = 1.0;
                        beta[5] = alpha[5] * hE / (4.0 * meas);
                        /*for (int i = 1; i < N_CD2D_ESTIMATOR_TYPES; i++)
                            estimated_error[i] += beta[i - 1] * jump / 2.0;*/
                        result += beta[int(estimatorType) - 1] * jump / 2.0;
                    }
                        /*************************************************************************/
                        /*  neighbour is not on the finest level, find children of neighbour     */
                        /*************************************************************************/
                    else                                      // there is a neighbour on the same level
                    {
                        int n = neigh->GetClipBoard();
                        if (n == -1) {
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
                            for (unsigned int r = 0; r < N_child; r++) {
                                // edge child, not general !!!
                                int edge1 = TmpoEnE[edge2neigh * MaxLen1 + r];
                                // local number of child cell
                                int chnum1 = TmpEC[edge1 * MaxLen2];
                                // child cell
                                TBaseCell *child = neigh->GetChild(chnum1);
                                // get local indices of child edge
                                refdesc->GetEdgeChildIndex(TmpECI, TmpLen3, MaxLen3);
                                // local index of child edge
                                int l_child = TmpECI[edge1 * MaxLen3];

                                // ref desc of child
                                TRefDesc *refdesc_child = child->GetRefDesc();
                                // conn. edge -> vertices
                                refdesc_child->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
                                // vertices of edge
                                TVertex *ver2 = child->GetVertex(TmpEdVer[2 * l_child]);
                                TVertex *ver3 = child->GetVertex(TmpEdVer[2 * l_child + 1]);

                                if (ee_verbose > 1) {
                                    cout << "ver 0 " << ver0->GetX() << "  " << ver0->GetY() << endl;
                                    cout << "ver 1 " << ver1->GetX() << "  " << ver1->GetY() << endl;
                                    cout << "ver 2 " << ver2->GetX() << "  " << ver2->GetY() << endl;
                                    cout << "ver 3 " << ver3->GetX() << "  " << ver3->GetY() << endl;
                                }

                                int part = 0;
                                if (ver1 == ver2) {
                                    part = 1;
                                } else if (ver0 == ver3) {
                                    part = -1;
                                } else {
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
                                if (neigh_N_ == -1) {
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
                                if (conform_grid) {
                                    switch (bf2DrefelementsNeigh)      // compute coordinates of line quadrature
                                    {                                 // points in reference cell
                                        case BFUnitSquare :             // edge 0
                                            if (neigh_edge == 0) {
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = -zeta[i];
                                                    edgeRefData.eta1DNeigh[i] = -1;
                                                }
                                            }
                                            if (neigh_edge == 1) {                             // edge 1
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = 1;
                                                    edgeRefData.eta1DNeigh[i] = -zeta[i];
                                                }
                                            }
                                            if (neigh_edge == 2) {                             // edge 2
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = zeta[i];
                                                    edgeRefData.eta1DNeigh[i] = 1;
                                                }
                                            }

                                            if (neigh_edge == 3) {                             // edge 3
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = -1;
                                                    edgeRefData.eta1DNeigh[i] = zeta[i];
                                                }
                                            }
                                            break;

                                        case BFUnitTriangle :
                                            if (neigh_edge == 0) {
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = (-zeta[i] + 1) / 2;
                                                    edgeRefData.eta1DNeigh[i] = 0;
                                                }
                                            }
                                            if (neigh_edge == 1) {
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = (zeta[i] + 1) / 2;
                                                    edgeRefData.eta1DNeigh[i] = (-zeta[i] + 1) / 2;
                                                }
                                            }
                                            if (neigh_edge == 2) {
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = 0;
                                                    edgeRefData.eta1DNeigh[i] = (zeta[i] + 1) / 2;
                                                }
                                            }
                                            break;
                                    }
                                }
                                else {
                                    switch (bf2DrefelementsNeigh)      // compute coordinates of line quadrature
                                        // this is only for 1-regular triangulations
                                    {                                 // points in reference cell
                                        case BFUnitSquare :             // edge 0
                                            if (neigh_edge == 0) {
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = (-zeta[i] + part) / 2;
                                                    edgeRefData.eta1DNeigh[i] = -1;
                                                }
                                            }
                                            if (neigh_edge == 1) {                             // edge 1
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = 1;
                                                    edgeRefData.eta1DNeigh[i] = (-zeta[i] + part) / 2;
                                                }
                                            }
                                            if (neigh_edge == 2) {                             // edge 2
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = (zeta[i] - part) / 2;
                                                    edgeRefData.eta1DNeigh[i] = 1;
                                                }
                                            }

                                            if (neigh_edge == 3) {                             // edge 3
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    edgeRefData.xi1DNeigh[i] = -1;
                                                    edgeRefData.eta1DNeigh[i] = (zeta[i] - part) / 2;
                                                }
                                            }
                                            break;

                                        case BFUnitTriangle :
                                            if (neigh_edge == 0) {
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    if (part == -1)
                                                        part = 0;
                                                    edgeRefData.xi1DNeigh[i] = ((-zeta[i] + 1) / 2 + part) / 2;
                                                    edgeRefData.eta1DNeigh[i] = 0;
                                                }
                                            }
                                            if (neigh_edge == 1) {
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    if (part == 1)
                                                        part = 0;
                                                    edgeRefData.xi1DNeigh[i] = ((zeta[i] + 1) / 2 - part) / 2;
                                                    if (part == 0)
                                                        part = 1;
                                                    if (part == -1)
                                                        part = 0;
                                                    edgeRefData.eta1DNeigh[i] = ((-zeta[i] + 1) / 2 + part) / 2;
                                                    if (part == 0)
                                                        part = -1;
                                                }
                                            }
                                            if (neigh_edge == 2) {
                                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)  // for all quadrature points
                                                {
                                                    if (part == 1)
                                                        part = 0;
                                                    edgeRefData.xi1DNeigh[i] = 0;
                                                    edgeRefData.eta1DNeigh[i] = ((zeta[i] + 1) / 2 - part) / 2;
                                                }
                                            }
                                            break;
                                    }
                                }
                                if (ee_verbose > 1) {
                                    for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                        cout << "xiN " << edgeRefData.xi1DNeigh[i] << " etaN " << edgeRefData.eta1DNeigh[i] << endl;
                                    }
                                }

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
                                for (unsigned int i = 0; i < N_Neigh; i++) {
                                    edgeRefData.FEFunctValuesNeigh[i] = Values[DOF[i]];
                                    if (ee_verbose > 1) {
                                        cout << " value " << edgeRefData.FEFunctValuesNeigh[i] << endl;
                                    }
                                }

                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)          // get values and derivatives in original cell
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
                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)           // for all quadrature points
                                {
                                    val[0] = val[1] = val[2] = 0;
                                    for (unsigned int l = 0; l < N_Neigh; l++)            // for all basis functions
                                    {
                                        // accumulate value of derivative
                                        val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l];
                                        // accumulate value of derivative
                                        val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l];
                                        // accumulate value of derivative
                                        val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l];
                                        if (ee_verbose > 1) {
                                            cout << l << "  " << edgeRefData.xderiv_refNeigh1D[i][l] << "  " << edgeRefData.yderiv_refNeigh1D[i][l] << "  " << edgeRefData.FEFunctValuesNeigh[l] << endl;
                                        }
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
                                switch (bf2DrefelementsNeigh) {
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

                                if (ee_verbose > 1) {
                                    for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                        cout << "xiN " << edgeRefData.xi1DNeigh[i] << " etaN " << edgeRefData.eta1DNeigh[i] << endl;
                                    }
                                }

                                // compute gradients in reference cell of the neighbour
                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)          // for all quadrature points
                                {
                                    bfNeigh->GetDerivatives(D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
                                    bfNeigh->GetDerivatives(D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
                                    bfNeigh->GetDerivatives(D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
                                }
                                // reftrafo of neighbour
                                RefTransNeigh = eleNeigh->GetRefTransID();
                                TFEDatabase2D::SetCellForRefTrans(child, RefTransNeigh);

                                DOF = fe_space->GetGlobalNumbers() + fe_space->GetBeginIndex()[neigh_N_];
                                for (unsigned int i = 0; i < N_Neigh; i++) {
                                    edgeRefData.FEFunctValuesNeigh[i] = Values[DOF[i]];
                                    if (ee_verbose > 1) {
                                        cout << " value " << edgeRefData.FEFunctValuesNeigh[i] << endl;
                                    }
                                }

                                // get values and derivatives in original cell
                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
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
                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                    val[0] = val[1] = val[2] = 0;
                                    // for all basis functions
                                    for (unsigned int l = 0; l < N_Neigh; l++) {
                                        // accumulate value of derivative
                                        val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l];
                                        // accumulate value of derivative
                                        val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l];
                                        // accumulate value of derivative
                                        val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l];
                                        if (ee_verbose > 1) {
                                            cout << l << "  " << edgeRefData.xderiv_refNeigh1D[i][l] << "  " << edgeRefData.yderiv_refNeigh1D[i][l] << "  " << edgeRefData.FEFunctValuesNeigh[l] << endl;
                                        }
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
                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                    if ((fabs(edgeRefData.X1DCell[i] - edgeRefData.X1DNeigh[i]) + fabs(edgeRefData.Y1DCell[i] - edgeRefData.Y1DNeigh[i])) > 1e-8)
                                        cout << " wrong quad points 2 " << edgeRefData.X1DCell[i] << " , " << edgeRefData.Y1DCell[i]
                                        << "   " << edgeRefData.X1DNeigh[i] << " , " << edgeRefData.Y1DNeigh[i] << endl;
                                    if (check_cont) if (fabs(edgeRefData.xyval_Neigh1D[i] - edgeRefData.xyval_Cell1D[i]) > 1e-8) {
                                        cout << "quad points b " << edgeRefData.X1DCell[i] << " , " << edgeRefData.Y1DCell[i] << endl;
                                        cout << " i " << i << " valb " << edgeRefData.xyval_Cell1D[i] << " neighb " << edgeRefData.xyval_Neigh1D[i] << " " << fabs(edgeRefData.xyval_Cell1D[i] - edgeRefData.xyval_Neigh1D[i]) << endl;
                                    }
                                    double e1 = coeff[0] * ((edgeRefData.xderiv_Cell1D[i] - edgeRefData.xderiv_Neigh1D[i]) * nx
                                                            + (edgeRefData.yderiv_Cell1D[i] - edgeRefData.yderiv_Neigh1D[i]) * ny);
                                    if (ee_verbose > 1) {
                                        cout << i << " jumpx " << edgeRefData.xderiv_Cell1D[i] << " " << edgeRefData.xderiv_Neigh1D[i] << endl;
                                    }
                                    double w = weights1D[i] * absdetjk1D;
                                    jump += w * e1 * e1;                   // integral on the edge
                                }
                                if (ee_verbose > 1) {
                                    cout << "jump " << jump << endl;
                                }
                                double hE2 = hE / N_child;
                                beta[0] = hE2;                        // weight for H^1 estimator
                                beta[1] = hE2 * hE2 * hE2;                // weight for L^2 estimator
                                beta[2] = hE2 / coeff[0];               // weight for energy norm estimator
                                if (TDatabase::ParamDB->INTERNAL_COERCIVITY > 0) {
                                    double w = 1 / sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
                                    if (w < beta[2])
                                        beta[2] = w;
                                }
                                //beta[2] *= 2.0;
                                if (24.0 / linfb < beta[2]) {
                                    beta[4] = 24.0 / linfb;
                                }
                                else {
                                    beta[4] = beta[2];
                                }
                                /*
                                beta[4] = 24;
                                if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
                                {
                                  linfb = sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY) * sqrt(coeff[0]);
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
                            TVertex *ver2 = neigh->GetVertex(TmpEdVerNeigh[2 * neigh_edge]);
                            TVertex *ver3 = neigh->GetVertex(TmpEdVerNeigh[2 * neigh_edge + 1]);
                            if (!(((ver0 == ver2) && (ver1 == ver3)) || ((ver0 == ver3) && (ver1 == ver2)))) {
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

                            switch (bf2DrefelementsNeigh)          // compute coordinates of line quadrature
                            {                                     // points in reference cell
                                case BFUnitSquare :                 // edge 0
                                    if (neigh_edge == 0) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)      // for all quadrature points
                                        {
                                            edgeRefData.xi1DNeigh[i] = -zeta[i];
                                            edgeRefData.eta1DNeigh[i] = -1;
                                        }
                                    }
                                    if (neigh_edge == 1) {                                 // edge 1
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)      // for all quadrature points
                                        {
                                            edgeRefData.xi1DNeigh[i] = 1;
                                            edgeRefData.eta1DNeigh[i] = -zeta[i];
                                        }
                                    }
                                    if (neigh_edge == 2) {                                 // edge 2
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)      // for all quadrature points
                                        {
                                            edgeRefData.xi1DNeigh[i] = zeta[i];
                                            edgeRefData.eta1DNeigh[i] = 1;
                                        }
                                    }

                                    if (neigh_edge == 3) {                                 // edge 3
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)      // for all quadrature points
                                        {
                                            edgeRefData.xi1DNeigh[i] = -1;
                                            edgeRefData.eta1DNeigh[i] = zeta[i];
                                        }
                                    }
                                    break;

                                case BFUnitTriangle :
                                    if (neigh_edge == 0) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)      // for all quadrature points
                                        {
                                            edgeRefData.xi1DNeigh[i] = (-zeta[i] + 1) / 2;
                                            edgeRefData.eta1DNeigh[i] = 0;
                                        }
                                    }
                                    if (neigh_edge == 1) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)      // for all quadrature points
                                        {
                                            edgeRefData.xi1DNeigh[i] = (zeta[i] + 1) / 2;
                                            edgeRefData.eta1DNeigh[i] = (-zeta[i] + 1) / 2;
                                        }
                                    }
                                    if (neigh_edge == 2) {
                                        for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)      // for all quadrature points
                                        {
                                            edgeRefData.xi1DNeigh[i] = 0;
                                            edgeRefData.eta1DNeigh[i] = (zeta[i] + 1) / 2;
                                        }
                                    }
                                    break;
                            }

                            if (ee_verbose > 1) {
                                for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                                    cout << "xiN " << edgeRefData.xi1DNeigh[i] << " etaN " << edgeRefData.eta1DNeigh[i] << endl;
                                }
                            }

                            // compute gradients in reference cell of the neighbour
                            for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)            // for all quadrature points
                            {
                                bfNeigh->GetDerivatives(D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
                                bfNeigh->GetDerivatives(D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
                                bfNeigh->GetDerivatives(D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
                            }
                            // reftrafo of neighbour
                            RefTrans2D RefTransNeigh = eleNeigh->GetRefTransID();
                            TFEDatabase2D::SetCellForRefTrans(neigh, RefTransNeigh);

                            int *DOF = fe_space->GetGlobalNumbers() + fe_space->GetBeginIndex()[neigh_N_];
                            for (unsigned int i = 0; i < N_Neigh; i++) {
                                edgeRefData.FEFunctValuesNeigh[i] = Values[DOF[i]];
                                if (ee_verbose > 1) {
                                    cout << " value " << edgeRefData.FEFunctValuesNeigh[i] << endl;
                                }
                            }
                            for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)            // get values and derivatives in original cell
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
                            for (unsigned int i = 0; i < N_QuadraturePoints1D; i++)             // for all quadrature points
                            {
                                val[0] = val[1] = val[2] = 0;
                                for (unsigned int l = 0; l < N_Neigh; l++)              // for all basis functions
                                {
                                    // accumulate value of derivative
                                    val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l];
                                    // accumulate value of derivative
                                    val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l];
                                    // accumulate value of derivative
                                    val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l];
                                    if (ee_verbose > 1) {
                                        cout << l << "  " << edgeRefData.xderiv_refNeigh1D[i][l] << "  " << edgeRefData.yderiv_refNeigh1D[i][l] << "  " << edgeRefData.FEFunctValuesNeigh[l] << endl;
                                    }
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
                                if ((fabs(edgeData.XEdge1D[edgeIdx][i] - edgeRefData.X1DNeigh[i]) + fabs(edgeData.YEdge1D[edgeIdx][i] - edgeRefData.Y1DNeigh[i])) > 1e-8) {
                                    cout << " wrong quad points 0 " << edgeData.XEdge1D[edgeIdx][i] << " , " << edgeData.YEdge1D[edgeIdx][i] << "   " << edgeRefData.X1DNeigh[i] << " , " << edgeRefData.Y1DNeigh[i] << endl;
                                }
                                if (check_cont) {
                                    if (fabs(edgeRefData.xyval_Neigh1D[i] - edgeData.xyval_1D[edgeIdx][i]) > 1e-8) {
                                        cout << " i " << i << " valc " << edgeData.xyval_1D[edgeIdx][i] << " neighc " << edgeRefData.xyval_Neigh1D[i] << endl;
                                    }
                                }
                                double e1 = coeff[0] * ((edgeData.xderiv_1D[edgeIdx][i] - edgeRefData.xderiv_Neigh1D[i]) * nx + (edgeData.yderiv_1D[edgeIdx][i] - edgeRefData.yderiv_Neigh1D[i]) * ny);
                                if (ee_verbose > 1) {
                                    cout << i << " jumpx " << edgeData.xderiv_1D[edgeIdx][i] << " " << edgeRefData.xderiv_Neigh1D[i] << endl;
                                }
                                double w = weights1D[i] * absdetjk1D;
                                jump += w * e1 * e1;                     // integral on the edge
                            }
                            if (ee_verbose > 1) {
                                cout << "jump " << jump << endl;
                            }
                            beta[0] = hE;                           // weight for H^1 estimator
                            beta[1] = hE * hE * hE;                     // weight for L^2 estimator
                            beta[2] = hE / coeff[0];                  // weight for energy norm estimator
                            if (TDatabase::ParamDB->INTERNAL_COERCIVITY > 0) {
                                double w = 1 / sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
                                if (w < beta[2])
                                    beta[2] = w;
                            }
                            //beta[2] *= 2.0;
                            if (24.0 / linfb < beta[2]) {
                                beta[4] = 24.0 / linfb;
                            }
                            else {
                                beta[4] = beta[2];
                            }
                            /*beta[4] = 24;
                            if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
                            {
                              linfb = sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY) * sqrt(coeff[0]);
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
    }
    // no delta_K, i.e., no SUPG
    if (delta_K == 4711 && (int(estimatorType) == 5 || int(estimatorType) == 6)) {
        return 0;
    }
    return result;
}

bool CDErrorEstimator2D::handleJump_BoundaryEdge(double *result, Example2D &example2D, const int estimatorType, const int N_QuadraturePoints1D, double *const &weights1D, const CDErrorEstimator2D::EdgeData &edgeData, BoundCond &Cond0,
                                                 const double meas, const double *coeff, double linfb, const std::vector<double> &alpha, int edgeIdx, const TJoint *joint) const {
    // vector holding the weights corresponding to the estimators, initialized with 0
    std::vector<double> beta(N_CD2D_ESTIMATOR_TYPES - 1, 0);
    // the edge
    TBoundEdge *bdryEdge = (TBoundEdge *) joint;
    // get boundary component
    TBoundComp2D *boundComp = bdryEdge->GetBoundComp();
    // parameter interval
    double t0, t1;
    bdryEdge->GetParameters(t0, t1);
    // boundary id
    int bdryId = boundComp->GetID();
    // type of boundary condition
    example2D.get_bc()[0](bdryId, (t0 + t1) / 2.0, Cond0);
    double jump = 0;
    // at midpoint of boundary
    switch (Cond0) {
        case DIRICHLET: {
            // no error
            if (TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS) {
                return false;
            }
            break;
        }
        case NEUMANN: {
            double x0, x1, y0, y1;
            // coordinates at begin of parameter interval
            bdryEdge->GetXYofT(t0, x0, y0);
            // coordinates at end of parameter interval
            bdryEdge->GetXYofT(t1, x1, y1);
            // outer normal vector
            double nx = y1 - y0;
            double ny = x0 - x1;
            // length of edge
            double hE = sqrt(nx * nx + ny * ny);
            // normalized normal vector
            nx /= hE;
            ny /= hE;

            // compute difference to Neumann condition
            for (unsigned int i = 0; i < N_QuadraturePoints1D; i++) {
                x0 = edgeData.XEdge1D[edgeIdx][i];
                y0 = edgeData.YEdge1D[edgeIdx][i];
                // coordinates at quadrature points
                bdryEdge->GetXYofT(t0, x0, y0);
                // Neumann data
                double neumann_data;
                example2D.boundary_data[0](bdryId, t0, neumann_data);
                double e1 = coeff[0] * (edgeData.xderiv_1D[edgeIdx][i] * nx + edgeData.yderiv_1D[edgeIdx][i] * ny) - neumann_data;
                double w = weights1D[i] * hE / 2.0;
                // integral on the edge
                jump += w * e1 * e1;
            }
            // weight for H^1 estimator
            beta[0] = hE;
            // weight for L^2 estimator
            beta[1] = hE * hE * hE;
            // weight for energy norm estimator
            beta[2] = hE / coeff[0];
            if (TDatabase::ParamDB->INTERNAL_COERCIVITY > 0) {
                double w = 1 / sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
                if (w < beta[2])
                    beta[2] = w;
            }
            if (24.0 / linfb < beta[2]) {
                beta[4] = 24.0 / linfb;
            }
            else {
                beta[4] = beta[2];
            }

            beta[5] = 1.0;
            beta[5] = alpha[5] * hE / (4.0 * meas);

            // fall through
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
        case DIRICHLET_WEAK: {
            cerr << "Only few BC implementation done " << endl;
            exit(-1);
        }
    }
    *result += beta[estimatorType - 1] * jump;
    return true;
}

std::ostream &operator<<(std::ostream &os, CDErrorEstimatorType &type) {
    switch (type) {
        case CDErrorEstimatorType::GradientIndicator: {
            os << "Gradient indicator";
            break;
        }
        case CDErrorEstimatorType::H1_ResidualEstimator: {
            os << "H1 residual estimator";
            break;
        }
        case CDErrorEstimatorType::L2_ResidualEstimator: {
            os << "L2 residual estimator";
            break;
        }
        case CDErrorEstimatorType::Energy_ResidualEstimatorQuasiRobust: {
            os << "Energy norm + dual norm residual estimator, Verfürth 2005";
            break;
        }
        case CDErrorEstimatorType::Energy_ResidualEstimatorWithoutJumps: {
            os << "Energy norm estimator without jumps";
            break;
        }
        case CDErrorEstimatorType::SUPG_Upper: {
            os << "SUPG estimator John/Novo, upper estimate";
            break;
        }
        case CDErrorEstimatorType::SUPG_Lower: {
            os << "SUPG estimator John/Novo, lower estimate";
            break;
        }
    }
    os << " (database value = " << int(type) << ")";
    return os;
}

CDErrorEstimator2D::CDErrorEstimator2D(Example2D &ex, int type) :
        ErrorEstimator2D(ex), estimatorType{CDErrorEstimatorType(type)} {
    estimated_global_error.resize(N_CD2D_ESTIMATOR_TYPES);
    conform_grid = TDatabase::ParamDB->GRID_TYPE;
}

int CDErrorEstimator2D::isConformGrid() const {
    return conform_grid;
}

void CDErrorEstimator2D::setConformGrid(int conform_grid) {
    this->conform_grid = conform_grid;
}