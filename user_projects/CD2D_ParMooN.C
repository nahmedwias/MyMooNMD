// =============================================================================
//
// Purpose:  main program for solving a stationary scalar equation using ParMooN
//
// Author:   Sashikumaar Ganesan
//
// History:  Implementation started on 22.08.2014
//           rewrite using CD2D class, Ulrich, March, 2015
// =============================================================================
#include <Domain.h>
#include <FEDatabase2D.h>
#include <MooNMD_Io.h>

#include <sys/stat.h>

#include <Output2D.h>
#include <CD2D.h>
#include <CDErrorEstimator2D.h>
#include <CDErrorEstimator2DDWR.h>
#include <functional>
#include <MainUtilities.h>

# define COMPARE_TO_VERY_FINE 0

// =============================================================================
// main program
// =============================================================================
int main(int argc, char *argv[]) {
    //  declaration of database, you need this in every program
    TDatabase Database;
    TFEDatabase2D FEDatabase;

#if COMPARE_TO_VERY_FINE
    Output::print<1>("solver type = ", TDatabase::ParamDB->SOLVER_TYPE );
    Output::print<1>("----------------solving on very fine grid----------------");
    TDomain domain2(argv[1]);
    domain2.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); //"hemker_fine.01.GEO
    Output::print<1>("Refining regularly:");
    for(int j = 0; j < 5; j++) {
        domain2.RegRefineAll();
    }
    std::cout << "done";
    Example_CD2D example2;
    CD2D cd2d_2(domain2, example2);
    std::cout << " assemble .. ";
    cd2d_2.assemble();
    std::cout << " solve ..";
    cd2d_2.solve();
    std::cout << "done" << std::endl;
    {
        TOutput2D Output(1, 1, 0, 0, NULL);
        Output.AddFEFunction(&cd2d_2.get_function());
        std::string filename(TDatabase::ParamDB->OUTPUTDIR);
        filename += "/uniform.vtk";
        Output.WriteVtk(filename.c_str());
    }
#endif

    /** set variables' value in TDatabase using argv[1] (*.dat file) */
    TDomain Domain(argv[1]);

    //set PROBLEM_TYPE to CD (convection diffusion) if not yet set
    if (TDatabase::ParamDB->PROBLEM_TYPE == 0)
        TDatabase::ParamDB->PROBLEM_TYPE = 1;
    //open OUTFILE, where all output is written to (in addition to console)
    Output::set_outfile(TDatabase::ParamDB->OUTFILE);

    // write all Parameters to the OUTFILE (not to console) for later reference
    Database.WriteParamDB(argv[0]);

    /* include the mesh from a mesh generator, for a standard mesh use the
     * build-in function. The GEOFILE describes the boundary of the domain. */
    Domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE);

    // refine grid up to the coarsest level
    for (int i = 0; i < TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR; i++)
        Domain.RegRefineAll();

    // write grid into an Postscript file
    if (TDatabase::ParamDB->WRITE_PS)
        Domain.PS("Domain.ps", It_Finest, 0);

    // create output directory, if not already existing
    if (TDatabase::ParamDB->WRITE_VTK)
        mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

    // for adaptive grid refinement:
    int current_estimator = TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION;
    // flag for either conforming closures or hanging nodes
    bool conforming = TDatabase::ParamDB->GRID_TYPE != 0;

    // choose example according to the value of TDatabase::ParamDB->EXAMPLE
    Example_CD2D example;

    // initialize the error estimator
    // this functional evaluates the integral of the solution in a ball around the center of radius 0.2
    std::function<double(const TFEFunction2D *, double, double, TBaseCell &)> dwrFunctional = [](const TFEFunction2D *sol, double x, double y, TBaseCell &cell) {
        if ((x - 0.5) * (x - 0.5) + y * y <= 0.2 * 0.2) {
            double integral = 0;
            FE2D feID = sol->GetFESpace2D()->GetFE2D(cell.GetCellIndex(), &cell); // id of finite element
            // calculate values on original element (i.e. prepare reference transformation)
            bool SecondDer = false; // defined in include/General/Constants.h
            double *weights, *xi, *eta;//quadrature weights and points in reference cell
            double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D]; // quadrature points
            double AbsDetjk[MaxN_QuadPoints_2D]; // determinant of transformation
            int n_points = 0;
            TFEDatabase2D::GetOrig(1, &feID, sol->GetFESpace2D()->GetCollection(), &cell, &SecondDer,
                                   n_points, xi, eta, weights, X, Y, AbsDetjk);
            // finite element on the current cell
            TFE2D *fe = TFEDatabase2D::GetFE2D(feID);
            const int n_loc_dof = fe->GetN_DOF(); // number of local dofs
            int *DOF = sol->GetFESpace2D()->GetGlobalDOF(cell.GetCellIndex());
            // id of the local basis functions
            BaseFunct2D base_fc_id = fe->GetBaseFunct2D()->GetID();
            // transformed values of basis functions
            double **orig_values = TFEDatabase2D::GetOrigElementValues(base_fc_id, D00);
            // local integration (loop over all quadrature points)
            for (auto j = 0; j < n_points; j++) {
                // local transformed values on this quadrature point
                double *orig = orig_values[j];
                double value = 0; // value of this TFEFunction2D at this quadrature point
                for (int l = 0; l < n_loc_dof; l++) {
                    // entry in the vector of this TFEFunction2D times basis function
                    value += sol->GetValues()[DOF[l]] * orig[l];
                }

                const double w = weights[j] * AbsDetjk[j];
                integral += w * value;
            }
            return integral; //* pow(10, 8)
        }
        return 0.0;
    };
    std::function<double(const TFEFunction2D *, double, double, TBaseCell &)> dwrFunctional2 = [](const TFEFunction2D *sol, double x, double y, TBaseCell &cell) {
        return 1.0;
    };
    //CDErrorEstimator2DDWR estimator{example, dwrFunctional2, Domain};
    CDErrorEstimator2D estimator {example, TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION };

    // refinement
    RefinementStrategy refinementStrategy;

    int LEVELS = TDatabase::ParamDB->LEVELS;
    for (int curr_level = 0; curr_level < LEVELS; curr_level++) {
        double t_start_level = GetTime(); // time for level
        {
            std::stringstream out;
            out << std::endl << "***********************************************************" << std::endl;
            out << "GEOMETRY  LEVEL " << curr_level;
            Output::print<1>(out.str());
        }


        // refine grid if level is greater than 0
        if ((curr_level > 0)) {
            // regular refinement if
            // adaptive procedure should not yet start
            // or no error estimation
            if ((curr_level <= TDatabase::ParamDB->UNIFORM_STEPS) || (!TDatabase::ParamDB->ESTIMATE_ERRORS)) {
                Domain.RegRefineAll();
            } else {
                Domain.RefineByRefinementStrategy(&refinementStrategy, conforming);
            }
        }

        // run cd2d on current grid
        CD2D cd2d(Domain, example);
        cd2d.assemble();
        cd2d.solve();
#if !COMPARE_TO_VERY_FINE
        cd2d.output(curr_level);
#endif

        if(TDatabase::ParamDB->ESTIMATE_ERRORS){
            estimator.estimate(cd2d.get_function());
            {
                if(TDatabase::ParamDB->WRITE_VTK) {
                    // the FE basis functions
                    BaseFunct2D *baseFunctions = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
                    // number of basis functions
                    int *n_baseFunctions = TFEDatabase2D::GetN_BaseFunctFromFE2D();
                    std::vector<double> eta_fct_values;
                    auto coll = cd2d.get_function().GetFESpace2D()->GetCollection();
                    eta_fct_values.resize((unsigned long) cd2d.get_function().GetLength());
                    for (auto cellIdx = 0; cellIdx < coll->GetN_Cells(); cellIdx++) {
                        // the cell
                        TBaseCell *cell = coll->GetCell(cellIdx);
                        // fe2d element on cell
                        FE2D element = cd2d.get_function().GetFESpace2D()->GetFE2D(cellIdx, cell);
                        for (unsigned int l = 0; l < n_baseFunctions[element]; l++) {
                            // fe values of dofs
                            auto *DOF = cd2d.get_function().GetFESpace2D()->GetGlobalNumbers() + cd2d.get_function().GetFESpace2D()->GetBeginIndex()[cellIdx];
                            eta_fct_values[DOF[l]] = sqrt(estimator.GetEta_K()[cellIdx]) / estimator.GetMaximalLocalError();
                        }
                    }

                    TFEFunction2D eta_fct = TFEFunction2D(cd2d.get_function().GetFESpace2D(), (char *) "etas", (char *) "etas", eta_fct_values.data(), (int) eta_fct_values.size());
                    TOutput2D Output(1, 1, 0, 0, NULL);
                    Output.AddFEFunction(&eta_fct);
                    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
                    filename += "/eta_" + std::to_string(curr_level);
                    filename += ".vtk";
                    Output.WriteVtk(filename.c_str());
                }
            }
#if COMPARE_TO_VERY_FINE
            {
                // calc errors to very fine sol'n
                TFEFunction2D *foo = const_cast<TFEFunction2D*>(&cd2d.get_function());
                double vals [foo->GetLength()];
                TFEFunction2D foutzz = TFEFunction2D(&cd2d.get_space(), "name", "dsc", vals, foo->GetLength());
                TFEFunction2D *foo2 = const_cast<TFEFunction2D*>(&cd2d_2.get_function());
                foutzz.Interpolate(foo2);
                foo->operator*=(-1.);
                foo->operator+=(foutzz);
                {
                    double errors[5];
                    TAuxParam2D aux;
                    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
                    const TFESpace2D* space = foo->GetFESpace2D();

                    foo->GetErrors(example.get_exact(0), 3, AllDerivatives, 4,
                                          SDFEMErrors, example.get_coeffs(), &aux, 1,
                                          &space, errors);

                    Output::print<1>("exact L2     : ", errors[0]);
                    Output::print<1>("exact H1-semi: ", errors[1]);
                    Output::print<1>("exact SD     : ", errors[2]);
                    Output::print<1>("exact L_inf  : ", errors[3]);
                }
            }
#endif
            refinementStrategy.applyEstimator(estimator);
            std::stringstream out;
            out << "estimated global error " << setw(10) << estimator.GetEstimatedGlobalError()[current_estimator];
            Output::print<1>(out.str());
        }
        {
            std::stringstream out;
            out << "Time for level " << curr_level << ": " << GetTime() - t_start_level;
            Output::print<1>(out.str());
        }
    }

    Output::close_file();
    return 0;
} // end main
