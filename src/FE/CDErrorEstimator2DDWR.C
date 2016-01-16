#include <CDErrorEstimator2DDWR.h>
#include <MainUtilities.h>
#include <Output2D.h>
#include <Solver.h>
#include <LinAlg.h>
#include <FEDatabase2D.h>

namespace {

    struct OrigValues {
        int n_points;
        double *xi;
        double *eta;
        double *weights;
        std::vector<double> X = std::vector<double>(MaxN_QuadPoints_2D, 0);
        std::vector<double> Y = std::vector<double>(MaxN_QuadPoints_2D, 0);
        std::vector<double> AbsDetjk = std::vector<double>(MaxN_QuadPoints_2D, 0);
    };

    // calculates the l2 norm of the laplace in the cell
    double calcL2LaplaceInCell(const TFEFunction2D* dual_sol, TBaseCell &cell, OrigValues &origValues, BaseFunct2D *baseFunctions, int *n_baseFunctions) {
        auto cellIdx = cell.GetCellIndex();
        auto fe_space = dual_sol->GetFESpace2D();
        // fe2d element on cell
        FE2D element = fe_space->GetFE2D(cellIdx, &cell);
        // we do require second derivatives
        bool secondDerivatives[] {true};
        // get reference transformation RefTrans2D
        TFEDatabase2D::GetOrig(1, &element, fe_space->GetCollection(), &cell, secondDerivatives,
                                                       origValues.n_points, origValues.xi, origValues.eta,
                                                       origValues.weights, origValues.X.data(),
                                                       origValues.Y.data(), origValues.AbsDetjk.data());
        // get fe function values for every basis function of the element
        double FEFunctValues[n_baseFunctions[element]];
        for (unsigned int l = 0; l < n_baseFunctions[element]; l++) {
            // fe values of dofs
            int *DOF = fe_space->GetGlobalNumbers() + fe_space->GetBeginIndex()[cellIdx];
            FEFunctValues[l] = dual_sol->GetValues()[DOF[l]];
        }
        // get values of 2nd derivatives
        double **orig_values_D20 = TFEDatabase2D::GetOrigElementValues(baseFunctions[element], D20);
        double **orig_values_D02 = TFEDatabase2D::GetOrigElementValues(baseFunctions[element], D02);
        double result = 0;
        for(auto quad_point = 0; quad_point < origValues.n_points; quad_point++) {
            double *orig_quad_point_D20 = orig_values_D20[quad_point];
            double *orig_quad_point_D02 = orig_values_D02[quad_point];
            double D_20 = 0;
            double D_02 = 0;
            for (auto base_function_idx = 0; base_function_idx < n_baseFunctions[element]; base_function_idx++) {
                D_20 += FEFunctValues[base_function_idx] * orig_quad_point_D20[base_function_idx];
                D_02 += FEFunctValues[base_function_idx] * orig_quad_point_D02[base_function_idx];
            }
            // weight of the gradient obtained by transformation of the element
            const double w = origValues.weights[quad_point] * origValues.AbsDetjk[quad_point];
            const double e1 = D_20 + D_02;
            // L^2 norm squared
            result += w * e1 * e1;
        }
        return result;
    }
}

void CDErrorEstimator2DDWR::estimate(const std::vector<MultiIndex2D> &derivatives, const TFEFunction2D &fe_function2D) {
    // for copy semantics
    std::vector<const TFEFunction2D*> sol_vector = {&fe_function2D};

    /**
     * Assemble and solve dual problem
     */

    // const reference to manipulate function for rhs
    const ManipulateAssemblingFct2D &manipulateFct = std::bind(&CDErrorEstimator2DDWR::manipulate, this,
                                                               sol_vector[0],
                                                               std::placeholders::_1,
                                                               std::placeholders::_2,
                                                               std::placeholders::_3,
                                                               std::placeholders::_4,
                                                               std::placeholders::_5,
                                                               std::placeholders::_6);
    std::vector<BoundCondFunct2D *> bc(1, example2D.get_bc()[0]);
    std::vector<BoundValueFunct2D *> bd(1, example2D.get_bd()[0]);
    // set up dual example as copy of original example
    Example_CD2D dualExample(example2D.get_exact(), bc, bd, example2D.get_coeffs());

    // dual problem class instance
    CD2DDual cd2dDual(domain, dualExample);
    // manipulated assemble
    cd2dDual.assemble(manipulateFct);
    // solve
    cd2dDual.solve_transposed();
    // solution
    const TFEFunction2D &z_h = cd2dDual.get_function();


    /**
     * Calculate resulting weights
     */

    // the fe space
    auto fe_space = fe_function2D.GetFESpace2D();
    // collection of fe space
    auto coll = fe_space->GetCollection();
    // array holding the dual weights for eta_K
    std::unique_ptr<double> dual_weights (new double[coll->GetN_Cells()]);
    double* dual_weights_ptr = dual_weights.get();
    // store values from ref transformation
    OrigValues origValues;
    // calculate 2nd derivatives, get l2 norm of laplace per cell and store resulting weights
    {
        // get base function set from fe space
        BaseFunct2D *baseFunctions = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
        // number of basis functions
        int *n_baseFunctions = TFEDatabase2D::GetN_BaseFunctFromFE2D();
        // for each cell
        for(auto cellIdx = 0; cellIdx < coll->GetN_Cells(); cellIdx++) {
            // current cell
            TBaseCell *cell = coll->GetCell(cellIdx);
            // diameter
            const double h_K = cell->GetDiameter();
            // calculate L^2 squared norm of laplace
            const double laplaceInCell = calcL2LaplaceInCell(&z_h, cell[0], origValues, baseFunctions, n_baseFunctions);
            // save weight (additional h_K per cell)
            dual_weights_ptr[cellIdx] = laplaceInCell;
        }
    }

    /**
     * Apply weights
     */
    // call super estimate
    CDErrorEstimator2D::estimate(derivatives, fe_function2D);
    // post processing with weights
    {
        for(auto x = 0; x < N_CD2D_ESTIMATOR_TYPES; x++) {
            estimated_global_error[x] = 0.0;
        }
        maximal_local_error = 0;
        for (auto cellIdx = 0; cellIdx < coll->GetN_Cells(); cellIdx++) {
            eta_K[cellIdx] *= dual_weights_ptr[cellIdx];
            for(auto x = 0; x < N_CD2D_ESTIMATOR_TYPES; x++) {
                estimated_global_error[x] += eta_K[cellIdx];
            }
            maximal_local_error = maximal_local_error < eta_K[cellIdx] ? eta_K[cellIdx] : maximal_local_error;
        }
        for(auto x = 0; x < N_CD2D_ESTIMATOR_TYPES; x++) {
            estimated_global_error[x] = sqrt(estimated_global_error[int(estimatorType)]);
        }
        maximal_local_error = sqrt(maximal_local_error);
    }
}

void CDErrorEstimator2DDWR::manipulate(const TFEFunction2D *sol, int N_Points, double* X, double* Y, double **Coeff, double **parameters, TBaseCell *cell) {
    for(size_t i = 0; i < N_Points; i++) {
        double x = X[i];
        double y = Y[i];
        Coeff[i][4] = dwrFunctional(sol, x, y, cell[0]);
    }
}

void CD2DDual::assemble(const ManipulateAssemblingFct2D &manipulate) {
    LocalAssembling2D_type t = LocalAssembling2D_type::ConvDiff;

    // this loop has more than one iteration only in case of multigrid
    for(auto & s : this->systems)
    {
        TFEFunction2D * pointer_to_function = &s.fe_function;
        // create a local assembling object which is needed to assemble the matrix
        LocalAssembling2D la(t, &pointer_to_function, example.get_coeffs());
        // set manipulate
        la.SetManipulateAssembling(const_cast<ManipulateAssemblingFct2D *>(&manipulate));
        // assemble the system matrix with given local assembling, solution and rhs
        s.matrix.Assemble(la, s.solution, s.rhs);
    }

    // when using afc, do it now
    if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION > 0)
    {
        this->performAlgebraicFluxCorrection();
    }
}

CD2DDual::CD2DDual(const TDomain& domain, Example_CD2D example, int reference_id)
        : CD2D(domain, example, reference_id) {
}

void CD2DDual::solve_transposed() {
    double t = GetTime();
    System_per_grid& s = this->systems.front();
    TSquareMatrix2D *SqMat[1] = {(TSquareMatrix2D *) s.matrix.get_matrix()->GetTransposed()};
    Solver((TSquareMatrix **)SqMat, NULL, s.rhs.get_entries(),
           s.solution.get_entries(), MatVect_Scalar, Defect_Scalar,
           this->multigrid.get(), this->get_size(), 0);

    t = GetTime() - t;
    delete SqMat[0];
    Output::print<2>(" solving of a dual CD2D problem done in ", t, " seconds");
}

//TODO remove me
static int llevel = 0;
void CD2DDual::output_dual(size_t level) {
    llevel++;
    level = (size_t) llevel;
    if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
        return;

    // print the value of the largest and smallest entry in the finite element
    // vector
    TFEFunction2D & fe_function = this->systems.front().fe_function;

    // write solution to a vtk file
    if(TDatabase::ParamDB->WRITE_VTK)
    {
        // last argument in the following is domain, but is never used in this class
        TOutput2D Output(1, 1, 0, 0, NULL);
        Output.AddFEFunction(&fe_function);
        std::string filename(TDatabase::ParamDB->OUTPUTDIR);
        filename += "/" + std::string(TDatabase::ParamDB->BASENAME) + "_dual";
        if(level >= 0)
            filename += "_" + std::to_string(level);
        filename += ".vtk";
        Output.WriteVtk(filename.c_str());
    }
}

CDErrorEstimator2DDWR::CDErrorEstimator2DDWR(Example2D &ex, CD2DDwrFunctional &functional, TDomain &domain)
        : CDErrorEstimator2D(ex, int(CDErrorEstimatorType::L2_ResidualEstimator)) {
    this->dwrFunctional = functional;
    this->domain = domain;
}