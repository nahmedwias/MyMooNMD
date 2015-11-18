#include <Enumerations.h>
#include <FEFunction2D.h>
#include <FEVectFunct2D.h>
#include <Database.h>
#include <NSEErrorEstimator2D.h>
#include <array>
#include <FEDatabase2D.h>

#define DEBUG_COMPARE_RESULTS_WITH_OLD_CODE 0

// TODO: Remove.
#include <NSE2DErrorEstimator.h>

#if DEBUG_COMPARE_RESULTS_WITH_OLD_CODE == 1
#include <NSE2DErrorEstimator.h>
#endif

namespace {
    unsigned int get_max_n_base_functions(const TFESpace2D &fe_space) {
        unsigned int max_n_base_functions_loc = 0;
        // # used finite elements
        unsigned int n{(unsigned int) fe_space.GetN_UsedElements()};
        // used finite elements
        FE2D *UsedElements = fe_space.GetUsedElements();
        // for all finite elements
        for (unsigned int j = 0; j < n; j++) {
            // if number of base functions is higher, update value
            unsigned int k = (unsigned int) TFEDatabase2D::GetN_BaseFunctFromFE2D(UsedElements[j]);
            if (k > max_n_base_functions_loc) {
                max_n_base_functions_loc = k;
            }
        }
        return max_n_base_functions_loc;
    }
}

/**
 * Data structure holding relevant data of the edges.
 *
 * XEdge1D, YEdge1D: coordinates of the edges
 * xi1D, eta1D: arrays holding the line quadrature points on reference cell of the form arr[BaseFunction][edgeIdx][quadraturePoint] = coordinate
 * xietaval, xideriv, etaderiv: arrays holding values of the derivatives at the quadrature points, same structure as xi1D and eta1D
 * AbsDetjk1D: determinant of the affine mapping for each edge of the form arr[edgeIdx][quadPoint] = value
 * xyval_ref1D: values in reference cell
 * xderiv_ref1D, yderiv_ref1D: derivative values in reference cell
 * xyval_1D: values in original cell
 * xderiv_1D, yderiv_1D: values of derivatives in original cell
 */
struct NSEErrorEstimator2D::EdgeData {

protected:
    // pointer to backing array for xietaval_ref1D, xideriv_ref1D, etaderiv_ref1D
    std::shared_ptr<double> xieta_ref1D_data;
    // pointer to backing array for xi1D, eta1D
    std::shared_ptr<double> xi_eta_1D_data;
    // pointer to backing array for gradient and value of function on edge
    std::shared_ptr<double> xyval_xyderiv_1D_data;
    // pointer to backing array for gradient and value of function on reference edge
    std::shared_ptr<double> xyval_xyderiv_ref_1D_data;
public:
    // coordinates of the edges
    std::vector<std::vector<double>> XEdge1D, YEdge1D;

    // determinant of the affine mapping for each edge
    std::vector<std::vector<double>> AbsDetjk1D{4, std::vector<double>(MaxN_QuadPoints_2D)};
    // values and derivative values in reference cell
    std::vector<std::array<double *, MaxN_QuadPoints_1D>> xyval_ref1D{4};
    std::vector<std::array<double *, MaxN_QuadPoints_1D>> xderiv_ref1D{4};
    std::vector<std::array<double *, MaxN_QuadPoints_1D>> yderiv_ref1D{4};
    // values and derivative values in original cell
    std::array<double *, 4> xyval_1D{};
    std::array<double *, 4> xderiv_1D{};
    std::array<double *, 4> yderiv_1D{};

    // mapping (base_function, edge, quadPoint) -> value
    std::vector<std::array<double *, 4>> xi1D{N_BaseFuncts2D}, eta1D{N_BaseFuncts2D};

    // mapping (base_function2D, edge, quadPoint, BaseFunction) -> value
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D>, 4>> xietaval_ref1D{N_BaseFuncts2D};
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D>, 4>> xideriv_ref1D{N_BaseFuncts2D};
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D>, 4>> etaderiv_ref1D{N_BaseFuncts2D};

    double *weights1D;
    double *zeta;
    unsigned int n_points_1D;
    // maximal number of base functions 2d
    unsigned int max_n_base_funct_2d;

    EdgeData(unsigned int max_n_base_functions_2d) {
        this->max_n_base_funct_2d = max_n_base_functions_2d;

        // initialize structures holding function values and derivatives on reference edge
        xyval_xyderiv_ref_1D_data = std::shared_ptr<double>(new double[3 * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d], std::default_delete<double[]>());
        {
            double *ptr = xyval_xyderiv_ref_1D_data.get();
            xyval_ref1D[0][0] = &ptr[0];
            xderiv_ref1D[0][0] = &ptr[4 * MaxN_QuadPoints_1D * max_n_base_funct_2d];
            yderiv_ref1D[0][0] = &ptr[2 * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d];
            for (size_t ii = 0; ii < 4; ii++) {
                for (size_t jj = 0; jj < MaxN_QuadPoints_1D; jj++) {
                    xyval_ref1D[ii][jj] = xyval_ref1D[0][0] + (ii * MaxN_QuadPoints_1D * max_n_base_funct_2d + jj * max_n_base_funct_2d);
                    xderiv_ref1D[ii][jj] = xderiv_ref1D[0][0] + (ii * MaxN_QuadPoints_1D * max_n_base_funct_2d + jj * max_n_base_funct_2d);
                    yderiv_ref1D[ii][jj] = yderiv_ref1D[0][0] + (ii * MaxN_QuadPoints_1D * max_n_base_funct_2d + jj * max_n_base_funct_2d);
                }
            }
        }

        // initialize structures holding (xi, eta) values on edge
        xi_eta_1D_data = std::shared_ptr<double>(new double[N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D * 2], std::default_delete<double[]>());
        {
            // back xi1D and eta1D by xi_eta_1D_data
            double *ptr = xi_eta_1D_data.get();
            xi1D[0][0] = &ptr[0];
            eta1D[0][0] = &ptr[N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D];
            for (size_t ii = 0; ii < N_BaseFuncts2D; ii++) {
                for (size_t jj = 0; jj < 4; jj++) {
                    xi1D[ii][jj] = xi1D[0][0] + (ii * 4 * MaxN_QuadPoints_1D + jj * MaxN_QuadPoints_1D);
                    eta1D[ii][jj] = eta1D[0][0] + (ii * 4 * MaxN_QuadPoints_1D + jj * MaxN_QuadPoints_1D);
                }
            }
        }
        // initialize structures holding values and derivatives on reference edge
        xieta_ref1D_data = std::shared_ptr<double>(new double[N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d * 3], std::default_delete<double[]>());
        {
            // back xietaval_ref1D, xideriv_ref1D, etaderiv_ref1D by xieta_ref1D_data
            double *ptr = xieta_ref1D_data.get();
            xietaval_ref1D[0][0][0] = &ptr[0];
            xideriv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d];
            etaderiv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d * 2];
            for (size_t ii = 0; ii < N_BaseFuncts2D; ii++) {
                for (size_t jj = 0; jj < 4; jj++) {
                    for (size_t kk = 0; kk < MaxN_QuadPoints_1D; kk++) {
                        xietaval_ref1D[ii][jj][kk] = xietaval_ref1D[0][0][0] + (kk * max_n_base_funct_2d + jj * MaxN_QuadPoints_1D * max_n_base_funct_2d + ii * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d);
                        xideriv_ref1D[ii][jj][kk] = xideriv_ref1D[0][0][0] + (kk * max_n_base_funct_2d + jj * MaxN_QuadPoints_1D * max_n_base_funct_2d + ii * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d);
                        etaderiv_ref1D[ii][jj][kk] = etaderiv_ref1D[0][0][0] + (kk * max_n_base_funct_2d + jj * MaxN_QuadPoints_1D * max_n_base_funct_2d + ii * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d);
                    }
                }
            }
        }
    }

    void setPoints_1D(unsigned int points_1D) {
        n_points_1D = points_1D;

        // initialize structures holding function values and derivatives on edge
        xyval_xyderiv_1D_data = std::shared_ptr<double>(new double[3 * 4 * 3 * n_points_1D]);
        {
            double *ptr = xyval_xyderiv_1D_data.get();
            for (size_t i = 0; i < 4; i++) {
                xyval_1D[i] = &ptr[i * 3 * n_points_1D];
                xderiv_1D[i] = &ptr[1 * 4 * 3 * n_points_1D + i * 3 * n_points_1D];
                yderiv_1D[i] = &ptr[2 * 4 * 3 * n_points_1D + i * 3 * n_points_1D];
            }
        }
    }

    void setQuadratureFormula(TQuadFormula1D &quadFormula) {
        quadFormula.GetFormulaData((int &) n_points_1D, weights1D, zeta);
        XEdge1D = std::vector<std::vector<double>>(4, std::vector<double>(n_points_1D));
        YEdge1D = std::vector<std::vector<double>>(4, std::vector<double>(n_points_1D));
    }
};

void NSEErrorEstimator2D::estimate(TFEVectFunct2D &fe_function2D_u, TFEFunction2D &fe_function2D_p, TAuxParam2D &Aux) {
    // remove old eta_K
    if (eta_K && eta_K == nullptr) delete[] eta_K;

    // initialization
    TCollection *coll = fe_function2D_u.GetFESpace2D()->GetCollection();
    eta_K = new double[coll->GetN_Cells()];

    // spaces u and p
    constexpr int n_fespaces = 2;

    TFESpace2D const *fe_spaces[2] = {fe_function2D_u.GetFESpace2D(), fe_function2D_p.GetFESpace2D()};

    /**
    * Maximal number of base functions for u and p, respectively
    */
    unsigned int max_n_base_functions_u = get_max_n_base_functions(fe_function2D_u.GetFESpace2D()[0]);
    unsigned int max_n_base_functions_p = get_max_n_base_functions(fe_function2D_p.GetFESpace2D()[0]);
    unsigned int max_n_base_functions = max_n_base_functions_u + max_n_base_functions_p;

    int NUsedElements, N_LocalUsedElements;
    int N_Points, N_Parameters, N_Edges, N_;
    std::vector<int> Used(N_FEs2D);
    int *N_BaseFunct;
    BaseFunct2D *BaseFuncts;
    std::vector<FE2D> LocalUsedElements(N_FEs2D);
    FE2D CurrentElement;
    QuadFormula1D LineQuadFormula;
    TQuadFormula1D *qf1D;
    BaseFunct2D BaseFunct, BaseFunctP;
    TBaseFunct2D *bf;
    BF2DRefElements bf2Drefelements;
    double *weights, *xi, *eta;

    double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
    double AbsDetjk[MaxN_QuadPoints_2D];
    RefTrans2D RefTrans;
    double *Param[MaxN_QuadPoints_2D];
    double *Derivatives[3 * MaxN_QuadPoints_2D];
    double *AuxArray[MaxN_QuadPoints_2D];
    std::vector<double> FEFunctValues(3 * max_n_base_functions);
    double max_loc_err;
    double estimated_global_errors[N_NSE2D_ESTIMATOR_TYPES], estimated_local_errors[N_NSE2D_ESTIMATOR_TYPES];
    int LocN_BF[N_BaseFuncts2D];
    BaseFunct2D LocBF[N_BaseFuncts2D];

    int *DOF_u, *DOF_p;

    int ee_verbose = 2;                             // verbosity

    {
        // ########################################################################
        // store information in local arrays
        // ########################################################################
        BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
        N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
        // note that char is used instead of bool to prevent specialization of the vector
        std::vector<char> NeedsSecondDer(n_fespaces);


        FE2D *FEUsedElements;
        NUsedElements = 0;
        for (TFESpace2D const *fespace : fe_spaces) {
            FEUsedElements = fespace->GetUsedElements();
            for (size_t j = 0; j < fespace->GetN_UsedElements(); j++) {
                Used[FEUsedElements[j]] = 1;
            }
        }

        // compute number of used elements
        std::vector<FE2D> UsedElements;
        {
            for (size_t i = 0; i < N_FEs2D; i++) {
                if (Used[i]) {
                    NUsedElements++;
                    UsedElements.push_back((FE2D) i);
                }
            }

            if (ee_verbose > 1) {
                cout << "estimator number of used elements: " << NUsedElements << endl;
                for (size_t i = 0; i < NUsedElements; i++)
                    cout << "UsedElements[" << i << "]: " << UsedElements[i] << endl;
            }
        }


        // ########################################################################
        // prepare error estimates
        // ########################################################################

        EdgeData edgeData{max_n_base_functions};

        double *values_u = fe_function2D_u.GetValues();                      // values of fe function
        double *values_p = fe_function2D_p.GetValues();

        {
            /**
             * prepare clip board to mark cells on the finest level with integers > -1 and everything
             * on non finest level with -1
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

        // initialize some quantities
        for (size_t i = 0; i < N_NSE2D_ESTIMATOR_TYPES; i++) {
            estimated_global_errors[i] = 0.0;
        }
        max_loc_err = 0;

        // ########################################################################
        // calculate values of base functions and derivatives on ref element
        // ########################################################################
        for (int i = 0; i < NUsedElements; i++) {
            CurrentElement = UsedElements[i];
            int l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement); // get 1d quad fromula
            LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2 * l);
            qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
            int N_Points1D = 0;
            double *weights1D, *zeta;
            qf1D->GetFormulaData((int &) N_Points1D, weights1D, zeta);

            edgeData.setPoints_1D((unsigned int) N_Points1D);
            edgeData.setQuadratureFormula(qf1D[0]);

            BaseFunct = BaseFuncts[CurrentElement];

            bf = TFEDatabase2D::GetBaseFunct2D(BaseFunct); // get base functions
            bf2Drefelements = bf->GetRefElement();

            switch (bf2Drefelements)                      // compute coordinates of line quadrature
            {                                            // points in reference cell
                // quadrilateral cell
                case BFUnitSquare :                          // edge 0
                    for (size_t j = 0; j < N_Points1D; j++)                 // for all quadrature points
                    {
                        edgeData.xi1D[BaseFunct][0][j] = edgeData.zeta[j];
                        edgeData.eta1D[BaseFunct][0][j] = -1;
                        bf->GetDerivatives(D00, edgeData.zeta[j], -1, edgeData.xietaval_ref1D[BaseFunct][0][j]);
                        bf->GetDerivatives(D10, edgeData.zeta[j], -1, edgeData.xideriv_ref1D[BaseFunct][0][j]);
                        bf->GetDerivatives(D01, edgeData.zeta[j], -1, edgeData.etaderiv_ref1D[BaseFunct][0][j]);
                    }                                        // edge 1
                    for (size_t j = 0; j < N_Points1D; j++)               // for all quadrature points
                    {
                        edgeData.xi1D[BaseFunct][1][j] = 1;
                        edgeData.eta1D[BaseFunct][1][j] = edgeData.zeta[j];
                        bf->GetDerivatives(D00, 1, edgeData.zeta[j], edgeData.xietaval_ref1D[BaseFunct][1][j]);
                        bf->GetDerivatives(D10, 1, edgeData.zeta[j], edgeData.xideriv_ref1D[BaseFunct][1][j]);
                        bf->GetDerivatives(D01, 1, edgeData.zeta[j], edgeData.etaderiv_ref1D[BaseFunct][1][j]);
                    }                                        // edge 2
                    for (size_t j = 0; j < N_Points1D; j++)                 // for all quadrature points
                    {
                        edgeData.xi1D[BaseFunct][2][j] = -edgeData.zeta[j];
                        edgeData.eta1D[BaseFunct][2][j] = 1;
                        bf->GetDerivatives(D00, -edgeData.zeta[j], 1, edgeData.xietaval_ref1D[BaseFunct][2][j]);
                        bf->GetDerivatives(D10, -edgeData.zeta[j], 1, edgeData.xideriv_ref1D[BaseFunct][2][j]);
                        bf->GetDerivatives(D01, -edgeData.zeta[j], 1, edgeData.etaderiv_ref1D[BaseFunct][2][j]);
                    }                                         // edge 3
                    for (size_t j = 0; j < N_Points1D; j++)                  // for all quadrature points
                    {
                        edgeData.xi1D[BaseFunct][3][j] = -1;
                        edgeData.eta1D[BaseFunct][3][j] = -edgeData.zeta[j];
                        bf->GetDerivatives(D00, -1, -edgeData.zeta[j], edgeData.xietaval_ref1D[BaseFunct][3][j]);
                        bf->GetDerivatives(D10, -1, -edgeData.zeta[j], edgeData.xideriv_ref1D[BaseFunct][3][j]);
                        bf->GetDerivatives(D01, -1, -edgeData.zeta[j], edgeData.etaderiv_ref1D[BaseFunct][3][j]);
                    }
                    break;

                case BFUnitTriangle :                        // triangular cell
                    for (size_t j = 0; j < N_Points1D; j++)                 // for all quadrature points
                    {
                        edgeData.xi1D[BaseFunct][0][j] = (edgeData.zeta[j] + 1) / 2;
                        edgeData.eta1D[BaseFunct][0][j] = 0;
                        bf->GetDerivatives(D00, (edgeData.zeta[j] + 1) / 2, 0, edgeData.xietaval_ref1D[BaseFunct][0][j]);
                        bf->GetDerivatives(D10, (edgeData.zeta[j] + 1) / 2, 0, edgeData.xideriv_ref1D[BaseFunct][0][j]);
                        bf->GetDerivatives(D01, (edgeData.zeta[j] + 1) / 2, 0, edgeData.etaderiv_ref1D[BaseFunct][0][j]);
                    }                                       // edge 1
                    for (size_t j = 0; j < N_Points1D; j++)                // for all quadrature points
                    {
                        edgeData.xi1D[BaseFunct][1][j] = (-edgeData.zeta[j] + 1) / 2;
                        edgeData.eta1D[BaseFunct][1][j] = (edgeData.zeta[j] + 1) / 2;
                        bf->GetDerivatives(D00, (-edgeData.zeta[j] + 1) / 2, (edgeData.zeta[j] + 1) / 2, edgeData.xietaval_ref1D[BaseFunct][1][j]);
                        bf->GetDerivatives(D10, (-edgeData.zeta[j] + 1) / 2, (edgeData.zeta[j] + 1) / 2, edgeData.xideriv_ref1D[BaseFunct][1][j]);
                        bf->GetDerivatives(D01, (-edgeData.zeta[j] + 1) / 2, (edgeData.zeta[j] + 1) / 2, edgeData.etaderiv_ref1D[BaseFunct][1][j]);
                    }                                       // edge 2
                    for (size_t j = 0; j < N_Points1D; j++)                // for all quadrature points
                    {
                        edgeData.xi1D[BaseFunct][2][j] = 0;
                        edgeData.eta1D[BaseFunct][2][j] = (-edgeData.zeta[j] + 1) / 2;
                        bf->GetDerivatives(D00, 0, (-edgeData.zeta[j] + 1) / 2, edgeData.xietaval_ref1D[BaseFunct][2][j]);
                        bf->GetDerivatives(D10, 0, (-edgeData.zeta[j] + 1) / 2, edgeData.xideriv_ref1D[BaseFunct][2][j]);
                        bf->GetDerivatives(D01, 0, (-edgeData.zeta[j] + 1) / 2, edgeData.etaderiv_ref1D[BaseFunct][2][j]);
                    }
                    break;
            }
        }

        // ########################################################################
        // loop over all cells
        // ########################################################################
        // for all cells in the collection
        for (size_t cellIdx = 0; cellIdx < coll->GetN_Cells(); cellIdx++) {
            // next cell
            TBaseCell *cell = coll->GetCell((int) cellIdx);

            // initialize local estimate
            eta_K[cellIdx] = 0.0;

            // get number of parameters of equation
            N_Parameters = Aux.GetN_Parameters();
            std::shared_ptr<double> ParamData(new double[MaxN_QuadPoints_2D * N_Parameters + 3 * MaxN_QuadPoints_2D * derivatives_u.size() + MaxN_QuadPoints_2D * 20], std::default_delete<double[]>());
            double *ParamDataPtr = ParamData.get();
            // set pointers
            for (size_t j = 0; j < MaxN_QuadPoints_2D; j++) {
                Param[j] = &ParamDataPtr[0] + j * N_Parameters;
                AuxArray[j] = &ParamDataPtr[MaxN_QuadPoints_2D * N_Parameters + 3 * MaxN_QuadPoints_2D * derivatives_u.size()] + j * 20;
            }
            for (size_t j = 0; j < 3 * MaxN_QuadPoints_2D; j++) {
                Derivatives[j] = &ParamDataPtr[MaxN_QuadPoints_2D * N_Parameters] + j * derivatives_u.size();
            }

            const TFESpace2D *fe_space_u = fe_function2D_u.GetFESpace2D();
            int *global_numbers_u = fe_space_u->GetGlobalNumbers();
            int *begin_index_u = fe_space_u->GetBeginIndex();

            const TFESpace2D *fe_space_p = fe_function2D_p.GetFESpace2D();
            int *global_numbers_p = fe_space_p->GetGlobalNumbers();
            int *begin_index_p = fe_space_p->GetBeginIndex();

            // ####################################################################
            // find local used elements on this cell
            // ####################################################################
            for (size_t j = 0; j < n_fespaces; j++) {
                // j == 0: velocity
                // j == 1: pressure
                CurrentElement = fe_spaces[j]->GetFE2D((int) cellIdx, cell);
                LocalUsedElements[j] = CurrentElement;
                LocN_BF[j] = N_BaseFunct[CurrentElement]; // local basis functions
                LocBF[j] = BaseFuncts[CurrentElement];
                // velocity needs 2nd derivative
                NeedsSecondDer[j] = j == 0;
            }
            N_LocalUsedElements = n_fespaces;

            // ####################################################################
            // calculate values on original element
            // ####################################################################

            // get reference transformation
            RefTrans = TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements.data(),
                                              coll, cell, (bool *) NeedsSecondDer.data(),
                                              N_Points, xi, eta, weights, X, Y, AbsDetjk);
            // get parameters of equ.
            if (N_Parameters > 0) {
                Aux.GetParameters(N_Points, coll, cell, (int) cellIdx, xi, eta, X, Y, Param);
            }

            // velocity first
            {
                // calculate all needed derivatives of this FE function
                CurrentElement = fe_space_u->GetFE2D((int) cellIdx, cell);
                BaseFunct = BaseFuncts[CurrentElement];
                N_ = N_BaseFunct[CurrentElement];
                // dof of current mesh cell
                DOF_u = global_numbers_u + begin_index_u[cellIdx];
                for (size_t l = 0; l < N_; l++) {
                    // u fe values of dofs
                    FEFunctValues[l] = values_u[DOF_u[l]];
                    // v fe values of dofs
                    FEFunctValues[l + max_n_base_functions] = values_u[DOF_u[l] + fe_function2D_u.GetLength()];
                }

                // compute values for all derivatives
                // in all quadrature points
                // in original mesh cell
                {
                    double **OrigFEValues;
                    double value[2];
                    // for all derivatives
                    for (size_t k = 0; k < derivatives_u.size(); k++) {
                        // get values in original cell
                        OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, derivatives_u[k]);
                        // for all quadrature points
                        for (size_t j = 0; j < N_Points; j++) {
                            // value in original cell
                            double *Orig = OrigFEValues[j];
                            value[0] = value[1] = 0;
                            // for all basis functions
                            for (size_t l = 0; l < N_; l++) {
                                // accumulate u value of derivative in point j
                                value[0] += FEFunctValues[l] * Orig[l];
                                // accumulate v value of derivative in point j
                                value[1] += FEFunctValues[l + max_n_base_functions] * Orig[l];
                            }
                            // k-th u derivative
                            Derivatives[j][k] = value[0];
                            // k-th v derivative
                            Derivatives[j + MaxN_QuadPoints_2D][k] = value[1];
                        }
                    }
                }

                // get coefficients of pde
                if (example2D.get_coeffs()) {
                    // 0 - eps, 1,2 - rhs
                    (example2D.get_coeffs())(N_Points, X, Y, Param, AuxArray);
                }
            }
            // prepare 1D quadrature formula
            // pressure first
            int N_P;
            {
                FE2D CurrentElementP = fe_space_p->GetFE2D((int) cellIdx, cell);  // finite element on cell
                BaseFunctP = BaseFuncts[CurrentElementP];       // basis functions
                N_P = N_BaseFunct[CurrentElementP];             // # basis functions
                DOF_p = global_numbers_p + begin_index_p[cellIdx];    // dof of current mesh cell
                //cout << "cell no " << i << " begin " << *DOF_p <<
                // " " << *GlobalNumbersP << " " << BeginIndexP[i] << endl;
                for (size_t l = 0; l < N_P; l++) {
                    FEFunctValues[l + 2 * max_n_base_functions] = values_p[DOF_p[l]];    // p fe values of dofs
                    //cout <<  FEFunctValues[l+2*max_n_base_functions] << endl;
                }

                {
                    double value;
                    double **OrigFEValues;
                    for (size_t k = 0; k < derivatives_p.size(); k++)            // for all derivatives
                    {                                       // get values in original cell
                        OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunctP, derivatives_p[k]);
                        for (size_t j = 0; j < N_Points; j++)               // for all quadrature points
                        {
                            double *Orig = OrigFEValues[j];             // value in original cell
                            value = 0;
                            for (size_t l = 0; l < N_P; l++)                   // for all basis functions
                            {
                                value += FEFunctValues[l + 2 * max_n_base_functions] * Orig[l]; // accumulate p value of derivative in point j
                            }
                            Derivatives[j + 2 * MaxN_QuadPoints_2D][k] = value;// for k-th v derivative
                        }                                     // endfor j
                    } // endfor k
                }
            }

            unsigned int n_quadrature_points_1d;
            {
                unsigned int l = (unsigned int) TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2 * l);
                qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                {
                    int n_points1d_val;
                    qf1D->GetFormulaData(n_points1d_val, edgeData.weights1D, edgeData.zeta);
                    n_quadrature_points_1d = (unsigned int) n_points1d_val;
                }
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)->MakeRefElementData(LineQuadFormula);
            }

            // # edges
            N_Edges = cell->GetN_Edges();
            // pressure second
            {
                // loop over all edges of cell
                for (size_t j = 0; j < N_Edges; j++) {
                    // get original coordinates of edge quad. points
                    TFEDatabase2D::GetOrigFromRef(RefTrans, (int) n_quadrature_points_1d, edgeData.xi1D[BaseFunctP][j], edgeData.eta1D[BaseFunctP][j], edgeData.XEdge1D[j].data(), edgeData.YEdge1D[j].data(), edgeData.AbsDetjk1D[j].data());
                    for (size_t k = 0; k < n_quadrature_points_1d; k++)  // get values and derivatives in original cell
                    {
                        TFEDatabase2D::GetOrigValues(RefTrans, edgeData.xi1D[BaseFunctP][j][k],
                                                     edgeData.eta1D[BaseFunctP][j][k],
                                                     TFEDatabase2D::GetBaseFunct2D(BaseFunctP),
                                                     coll, (TGridCell *) cell,
                                                     edgeData.xietaval_ref1D[BaseFunctP][j][k],
                                                     edgeData.xideriv_ref1D[BaseFunctP][j][k],
                                                     edgeData.etaderiv_ref1D[BaseFunctP][j][k],
                                                     edgeData.xyval_ref1D[j][k],
                                                     edgeData.xderiv_ref1D[j][k],
                                                     edgeData.yderiv_ref1D[j][k]);
                    }

                    {
                        double value;
                        size_t m;
                        // for all quadrature points
                        for (size_t k = 0; k < n_quadrature_points_1d; k++) {
                            value = 0;
                            // for all basis functions
                            for (size_t l = 0; l < N_P; l++) {
                                m = (size_t) (l + 2 * max_n_base_functions);
                                value += FEFunctValues[m] * edgeData.xyval_ref1D[j][k][l]; // accumulate value of derivative
                            } // endfor l
                            m = k + 2 * n_quadrature_points_1d;
                            edgeData.xyval_1D[j][m] = value; // for k-th
                        } // endfor k
                    }
                }     // endfor j
            }
            // velocity second
            {
                unsigned int m;
                double val[6];
                for (size_t j = 0; j < N_Edges; j++)                  // loop over all edges of cell
                {                                     // get original coordinates of edge quad. points
                    TFEDatabase2D::GetOrigFromRef(RefTrans, (int) n_quadrature_points_1d, edgeData.xi1D[BaseFunct][j],
                                                  edgeData.eta1D[BaseFunct][j],
                                                  edgeData.XEdge1D[j].data(), edgeData.YEdge1D[j].data(), edgeData.AbsDetjk1D[j].data());
                    for (size_t k = 0; k < n_quadrature_points_1d; k++)  // get values and derivatives in original cell
                    {
                        TFEDatabase2D::GetOrigValues(RefTrans, edgeData.xi1D[BaseFunct][j][k], edgeData.eta1D[BaseFunct][j][k],
                                                     TFEDatabase2D::GetBaseFunct2D(BaseFunct), coll, (TGridCell *) cell,
                                                     edgeData.xietaval_ref1D[BaseFunct][j][k], edgeData.xideriv_ref1D[BaseFunct][j][k],
                                                     edgeData.etaderiv_ref1D[BaseFunct][j][k], edgeData.xyval_ref1D[j][k],
                                                     edgeData.xderiv_ref1D[j][k], edgeData.yderiv_ref1D[j][k]);
                    }

                    for (size_t k = 0; k < n_quadrature_points_1d; k++)     // for all quadrature points
                    {
                        val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = 0;
                        for (size_t l = 0; l < N_; l++)       // for all basis functions
                        {
                            m = (int) (l + max_n_base_functions);
                            val[0] += FEFunctValues[l] * edgeData.xyval_ref1D[j][k][l]; // accumulate value of derivative
                            val[1] += FEFunctValues[l] * edgeData.xderiv_ref1D[j][k][l]; // accumulate value of derivative
                            val[2] += FEFunctValues[l] * edgeData.yderiv_ref1D[j][k][l]; // accumulate value of derivative
                            val[3] += FEFunctValues[m] * edgeData.xyval_ref1D[j][k][l]; // accumulate value of derivative
                            val[4] += FEFunctValues[m] * edgeData.xderiv_ref1D[j][k][l]; // accumulate value of derivative
                            val[5] += FEFunctValues[m] * edgeData.yderiv_ref1D[j][k][l]; // accumulate value of derivative
                        } // endfor l
                        m = (int) (k + n_quadrature_points_1d);
                        edgeData.xyval_1D[j][k] = val[0]; // for k-th
                        edgeData.xderiv_1D[j][k] = val[1]; // for k-th
                        edgeData.yderiv_1D[j][k] = val[2]; // for k-th
                        edgeData.xyval_1D[j][m] = val[3]; // for k-th
                        edgeData.xderiv_1D[j][m] = val[4]; // for k-th
                        edgeData.yderiv_1D[j][m] = val[5]; // for k-th
                    } // endfor k
                }     // endfor j
            }

            {
                // estimate local errors
                calculateEtaK(fe_space_u[0], fe_space_p[0],
                              cell, N_Points, X, Y, AbsDetjk, weights, Derivatives, AuxArray,
                              example2D, edgeData,
                              global_numbers_u, begin_index_u, DOF_u, values_u,
                              global_numbers_p, begin_index_p, DOF_p, values_p,
                              &estimated_local_errors[0]);
            }


            for (size_t k = 0; k < N_NSE2D_ESTIMATOR_TYPES; k++) {
                // update global error estimates
                estimated_global_errors[k] += estimated_local_errors[k];
            }
            // update maximal local error estimate
            if (estimated_local_errors[int(estimatorType)] > max_loc_err) {
                max_loc_err = estimated_local_errors[int(estimatorType)];
            }

            eta_K[cellIdx] = estimated_local_errors[int(estimatorType)];
        }

        // compute global error estimates
        for (size_t i = 1; i < N_NSE2D_ESTIMATOR_TYPES; i++) {
            estimated_global_error[i] = sqrt(estimated_global_errors[i]);
        }
        estimated_global_error[0] = estimated_global_errors[0];
        // compute maximal local error estimate
        maximal_local_error = sqrt(max_loc_err);
    }

#if DEBUG_COMPARE_RESULTS_WITH_OLD_CODE == 0
    int type = int(estimatorType);
    int errorControl = TDatabase::ParamDB->ERROR_CONTROL;
    int is_nse = TDatabase::ParamDB->PROBLEM_TYPE != 3;
    double *debug_eta_K = new double[coll->GetN_Cells()];
    double debug_maximal_local_error;
    std::vector<double> debug_estimated_global_error (N_NSE2D_ESTIMATOR_TYPES);
    TNS2DErrorEstimator debugimator(type, &fe_function2D_u, &fe_function2D_p, errorControl, is_nse);
    TAuxParam2D aux;
    debugimator.GetErrorEstimate((int) derivatives_u.size(), const_cast<MultiIndex2D *>(derivatives_u.data()),
                                 (int) derivatives_p.size(), const_cast<MultiIndex2D *>(derivatives_p.data()),
                                 example2D.get_coeffs(), example2D.get_bc(), example2D.get_bd(),
                                 &aux, 2, const_cast<TFESpace2D **>(fe_spaces), &debug_eta_K[0], &debug_maximal_local_error, &debug_estimated_global_error[0]);

    for(int i = 0; i < coll->GetN_Cells(); i++) {
        if(fabs(debug_eta_K[i] - eta_K[i])> 1e-15) {
            std::cerr << "something wrong with the eta_k's" << std::endl;
        }
    }
    delete[] debug_eta_K;
#endif
}

void NSEErrorEstimator2D::calculateEtaK(const TFESpace2D &fe_space_u, const TFESpace2D &fe_space_p, TBaseCell *cell,
                                        int N_Points, double *X, double *Y, double *AbsDetjk, double *weights, double **Derivatives, double **coeffs,
                                        Example2D &example,
                                        EdgeData &edgeData,
                                        int *global_numbers_u, int *begin_index_u, int *DOF_u, double *values_u,
                                        int *global_numbers_p, int *begin_index_p, int *DOF_p, double *values_p,
                                        double *estimated_local_error) {

    double estimated_error[N_NSE2D_ESTIMATOR_TYPES];

    // initalize local error estimates
    for (size_t i = 0; i < N_NSE2D_ESTIMATOR_TYPES; i++)
        estimated_error[i] = 0.0;

    if (estimatorType == NSE2DErrorEstimatorType::gradient_indicator) {
        double w;
        double e1, e2, e3, e4;
        // for all quadrature points
        for (int i = 0; i < N_Points; i++) {
            w = weights[i] * AbsDetjk[i];
            // all u derivatives in quadrature points
            double *deriv = Derivatives[i];
            // x derivative
            e1 = deriv[0];
            // y derivative
            e2 = deriv[1];
            deriv = Derivatives[i + MaxN_QuadPoints_2D];
            // all v derivatives in quadrature points
            e3 = deriv[0];
            e4 = deriv[1];
            estimated_error[int(NSE2DErrorEstimatorType::gradient_indicator)] += w * (e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4);
        } // endfor i
    }

    for (int i = 0; i < N_NSE2D_ESTIMATOR_TYPES; i++) {
        estimated_local_error[i] = estimated_error[i];
    }

}

NSEErrorEstimator2D::NSEErrorEstimator2D(Example2D &ex, TDomain &domain, int type) : ErrorEstimator2D(ex, domain), estimatorType{NSE2DErrorEstimatorType(type)} {
    estimated_global_error.resize(N_NSE2D_ESTIMATOR_TYPES);
    conform_grid = TDatabase::ParamDB->GRID_TYPE;
}

std::ostream &operator<<(std::ostream &os, NSE2DErrorEstimatorType &type) {
    switch (type) {
        case NSE2DErrorEstimatorType::gradient_indicator: {
            os << "Gradient indicator";
            break;
        }
        case NSE2DErrorEstimatorType::residual_estimator_h1: {
            os << "H1 residual estimator";
            break;
        }
        case NSE2DErrorEstimatorType::residual_estimator_l2: {
            os << "L2 residual estimator";
            break;
        }
        case NSE2DErrorEstimatorType::residual_estimator_energy_quasi_robust: {
            os << "Energy residual estimator quasi robust";
            break;
        }
        case NSE2DErrorEstimatorType::gradient_recovery: {
            os << "Gradient recovery";
            break;
        }
        case NSE2DErrorEstimatorType::implicit_estimator_neumann: {
            os << "Implicit neumann estimator";
            break;
        }
    }
    os << " (database value = " << int(type) << ")";
    return os;
}

