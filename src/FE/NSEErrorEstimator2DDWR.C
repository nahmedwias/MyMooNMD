#include <NSEErrorEstimator2DDWR.h>
#include <NSEErrorEstimator2D.h>
#include <FEDatabase2D.h>
#include <Solver.h>
#include <LinAlg.h>

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

void NSEErrorEstimator2DDWR::estimate(const TFEVectFunct2D &u, const TFEFunction2D &p, TAuxParam2D &Aux) {
    // for copy semantics
    std::vector<const TFEVectFunct2D*> sol_vector_u = {&u};
    std::vector<const TFEFunction2D*> sol_vector_p = {&p};

    /**
     * Assemble and solve dual problem
     */

    // const reference to manipulate function for rhs
    const ManipulateAssemblingFct2D &manipulateFct = std::bind(&NSEErrorEstimator2DDWR::manipulate, this,
                                                               sol_vector_u[0],
                                                               sol_vector_p[0],
                                                               std::placeholders::_1,
                                                               std::placeholders::_2,
                                                               std::placeholders::_3,
                                                               std::placeholders::_4,
                                                               std::placeholders::_5,
                                                               std::placeholders::_6);
    std::vector<BoundCondFunct2D *> bc(example2D.boundary_conditions);
    std::vector<BoundValueFunct2D *> bd(example2D.boundary_data);
    // set up dual example as copy of original example
    Example_NSE2D dualExample(example2D.get_exact(), bc, bd, example2D.get_coeffs());

    // dual problem class instance
    NSE2DDual nse2dDual(domain, dualExample);
    // manipulated assemble
    nse2dDual.assemble(manipulateFct);
    // solve
    nse2dDual.solve_transposed();
    // solution
    const TFEVectFunct2D &z_h_u = nse2dDual.get_velocity();
    const TFEFunction2D &z_h_p = nse2dDual.get_pressure();


    /**
     * Calculate resulting weights
     */

    // the fe space
    auto fe_space = u.GetFESpace2D();
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
            const double laplaceInCell_u = calcL2LaplaceInCell(&z_h_u, cell[0], origValues, baseFunctions, n_baseFunctions);
            const double laplaceInCell_p = calcL2LaplaceInCell(&z_h_p, cell[0], origValues, baseFunctions, n_baseFunctions);
            // save weight (additional h_K per cell)
            dual_weights_ptr[cellIdx] = h_K * (laplaceInCell_u * laplaceInCell_u + laplaceInCell_p * laplaceInCell_p);
        }
    }

    /**
     * Apply weights
     */
    // call super estimate
    NSEErrorEstimator2D::estimate(u, p, Aux);
    // post processing with weights
    {
        for(auto x = 0; x < N_NSE2D_ESTIMATOR_TYPES; x++) {
            estimated_global_error[x] = 0.0;
        }
        maximal_local_error = 0;
        for (auto cellIdx = 0; cellIdx < coll->GetN_Cells(); cellIdx++) {
            eta_K[cellIdx] *= dual_weights_ptr[cellIdx];
            for(auto x = 0; x < N_NSE2D_ESTIMATOR_TYPES; x++) {
                estimated_global_error[x] += eta_K[cellIdx];
            }
            maximal_local_error = maximal_local_error < eta_K[cellIdx] ? eta_K[cellIdx] : maximal_local_error;
        }
        for(auto x = 0; x < N_NSE2D_ESTIMATOR_TYPES; x++) {
            estimated_global_error[x] = sqrt(estimated_global_error[int(estimatorType)]);
        }
        maximal_local_error = sqrt(maximal_local_error);
    }

}

void NSEErrorEstimator2DDWR::manipulate(const TFEVectFunct2D *sol_u , const TFEFunction2D *sol_p, int N_Points, double* X, double* Y, double **Coeff, double **parameters, TBaseCell *cell) {
    for(size_t i = 0; i < N_Points; i++) {
        double x = X[i];
        double y = Y[i];
        std::vector<double> functVals = dwrFunctional(sol_u, sol_p, x, y, cell[0]);
        Coeff[i][1] = functVals[0];
        Coeff[i][2] = functVals[1];
        Coeff[i][3] = functVals[2];
        Coeff[i][4] = functVals[3];
        Coeff[i][5] = functVals[4];
        Coeff[i][6] = functVals[5];
    }
}

NSEErrorEstimator2DDWR::NSEErrorEstimator2DDWR(Example2D &ex, NSE2DDwrFunctional &funct, TDomain &domain) : NSEErrorEstimator2D(ex, int(NSE2DErrorEstimatorType::residual_estimator_h1), false) {
    this->dwrFunctional = funct;
    this->domain = domain;
}

NSE2DDual::NSE2DDual(const TDomain &domain, Example_NSE2D example, unsigned int reference_id) : NSE2D(domain, example, reference_id) { }

void NSE2DDual::assemble(const ManipulateAssemblingFct2D &manipulate) {
    // the class LocalAssembling2D which we will need next, requires an array of
    // pointers to finite element functions, i.e. TFEFunction2D **.
    for(System_per_grid& s : this->systems)
    {
        TFEFunction2D *fe_functions[3] =
                { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
        LocalAssembling2D la(NSE2D_Galerkin, fe_functions,
                             this->example.get_coeffs());
        la.SetManipulateAssembling(const_cast<ManipulateAssemblingFct2D *>(&manipulate));
        s.matrix.Assemble(la, s.rhs);
        // set rhs for Dirichlet nodes
        s.solution.copy_nonactive(s.rhs);
        delete fe_functions[0];
        delete fe_functions[1];
    }
}

void NSE2DDual::solve_transposed() {
    System_per_grid& s = this->systems.front();
    s.matrix.get_BT_block(0)->scale(-1.);
    s.matrix.get_BT_block(1)->scale(-1.);
    NSE2D::solve();
    s.matrix.get_BT_block(0)->scale(-1.);
    s.matrix.get_BT_block(1)->scale(-1.);
}