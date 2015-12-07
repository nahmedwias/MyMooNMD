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

    // initialize the error estimator
    //CDErrorEstimator2D estimator {example, TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION };
    // this functional evaluates the integral of the solution in a ball around the center of radius 0.2
    std::function<double(const TFEFunction2D*, double, double, TBaseCell&)> dwrFunctional = [](const TFEFunction2D* sol, double x, double y, TBaseCell& cell) {
        if( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) <= 0.2 * 0.2 ) {
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
            //BaseFunct2D base_fc_id = fe->GetBaseFunct2D()->GetID();
            // transformed values of basis functions
            //double **orig_values = TFEDatabase2D::GetOrigElementValues(base_fc_id, D00);
            // local integration (loop over all quadrature points)
            for (auto j = 0; j < n_points; j++) {
                // local transformed values on this quadrature point
                //double *orig = orig_values[j];
                double value = 0; // value of this TFEFunction2D at this quadrature point
                for (int l = 0; l < n_loc_dof; l++) {
                    value += sol->GetValues()[DOF[l]];
                }

                integral += value / n_loc_dof;
            }
            return integral / n_points;
        }
        return 0.0;
    };
    CDErrorEstimator2DDWR estimator {example, dwrFunctional, Domain};

    // refinement
    RefinementStrategy refinementStrategy;

    int LEVELS = TDatabase::ParamDB->LEVELS;
    for (int curr_level = 0; curr_level < LEVELS; curr_level++) {
        double t_start_level = GetTime(); // time for level
        {
            std::stringstream out;
            out << std::endl << "***********************************************************"<<std::endl;
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

        {
            estimator.estimate(cd2d.get_function());
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
