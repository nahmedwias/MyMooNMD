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

static double hmin = 0;
static double hmax = 0;

// =============================================================================
// main program
// =============================================================================
int main(int argc, char *argv[]) {
    //  declaration of database, you need this in every program
    TDatabase Database;
    TFEDatabase2D FEDatabase;

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

    std::function<double(const TFEFunction2D *, double, double, TBaseCell &)> dwrFunctional2 = [](const TFEFunction2D *sol, double x, double y, TBaseCell &cell) {
        if (hmin == 0) {
            std::cout << " --------- <<<<< assigning hmin >>>>> -----------" << std::endl;
            //double _hmin, _hmax;
            //sol->GetFESpace2D()->GetCollection()->GetHminHmax(&_hmin, &_hmax);
            hmin = 1.;
        }
        double r = .01; ////hmin;
        if ((y - .3) * (y - .3) + (x - .7) * (x - .7) <= r * r) {
            return 1.0 / (Pi * r * r);
        } else {
            return 0.;
        }
    };
    CDErrorEstimator2DDWR estimator{example, dwrFunctional2, Domain};
    //CDErrorEstimator2D estimator {example, TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION };

    // refinement
    RefinementStrategy refinementStrategy;

    int LEVELS = TDatabase::ParamDB->LEVELS;
    for (int curr_level = 0; curr_level < LEVELS; curr_level++) {
        hmin = 0;
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
        cd2d.output(curr_level);

        if (TDatabase::ParamDB->ESTIMATE_ERRORS) {
            estimator.estimate(cd2d.get_function());
            {
                if (TDatabase::ParamDB->WRITE_VTK) {
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
            {
                double estimated = estimator.GetEstimatedGlobalError()[current_estimator];
                double trueErr = 0;
                double FEFunctValues[MaxN_BaseFunctions2D];
                std::vector<double> values = std::vector<double>((unsigned long) cd2d.get_function().GetLength());
                // calculate efficiency index
                TFEFunction2D err_fct = TFEFunction2D(cd2d.get_function().GetFESpace2D(), (char *) "errf", (char *) "errf", values.data(), (int) values.size());
                err_fct.Interpolate(example.exact_solution[0]);
                err_fct *= -1.;
                err_fct += cd2d.get_function();

                double *weights, *xi, *eta;
                double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
                double AbsDetjk[MaxN_QuadPoints_2D];
                double *righthand;
                auto space = cd2d.get_function().GetFESpace2D();

                auto BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
                auto N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

                auto GlobalNumbers = space->GetGlobalNumbers();
                auto BeginIndex = space->GetBeginIndex();
                for (auto cellIdx = 0; cellIdx < space->GetN_Cells(); cellIdx++) {
                    auto cell = space->GetCollection()->GetCell(cellIdx);
                    constexpr int n_fespaces = 1;

                    auto DOF = GlobalNumbers + BeginIndex[cellIdx];

                    std::vector<FE2D> LocalUsedElements(n_fespaces);
                    std::vector<int> LocN_BF(n_fespaces);
                    std::vector<BaseFunct2D> LocBF(n_fespaces);
                    int iSpace = 0;

                    FE2D CurrentElement = space->GetFE2D(cellIdx, cell);
                    TFE2D *fe = TFEDatabase2D::GetFE2D(CurrentElement);
                    LocalUsedElements[iSpace] = CurrentElement;
                    LocN_BF[iSpace] = fe->GetSize();
                    LocBF[iSpace] = fe->GetBaseFunct2D_ID();

                    for (auto l = 0; l < N_BaseFunct[CurrentElement]; l++) {
                        FEFunctValues[l] = err_fct.GetValues()[DOF[l]];
                    }

                    int N_LocalUsedElements = n_fespaces;

                    int N_Points; // number of quadrature points
                    bool secondDer[]{false};
                    TFEDatabase2D::GetOrig(N_LocalUsedElements, &LocalUsedElements[0], cd2d.get_space().GetCollection(),
                                           cell, secondDer, N_Points, xi, eta, weights, X, Y,
                                           AbsDetjk);
                    double **elementVals = TFEDatabase2D::GetOrigElementValues(LocBF[0], D00);
                    double *orig = elementVals[0];
                    for (int i = 0; i < N_Points; i++) {
                        double x = X[i];
                        double y = Y[i];
                        double Mult = weights[i] * AbsDetjk[i];
                        double functVal = 0;
                        functVal = dwrFunctional2(&err_fct, x, y, cell[0]);
                        trueErr += Mult * (functVal * FEFunctValues[i] * orig[i]);
                    }

                }
                {
                    std::stringstream out;
                    out << "efficiency index " << setw(10) << (fabs(estimated) / fabs(trueErr)) << ", trueErr " << fabs(trueErr);
                    Output::print<1>(out.str());
                }
            }
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
