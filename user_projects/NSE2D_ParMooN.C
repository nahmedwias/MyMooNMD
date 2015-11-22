// =======================================================================
//
// Purpose:     main program for solving a stationary NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 23.08.2014

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Output2D.h>
#include <NSE2D.h>

#include <sys/stat.h>
#include <RefinementStrategy.h>
#include "NSEErrorEstimator2D.h"

// =======================================================================
// main program
// =======================================================================
int main(int argc, char *argv[]) {
    //  declaration of database, you need this in every program
    TDatabase Database;
    TFEDatabase2D FEDatabase;

    /** set variables' value in TDatabase using argv[1] (*.dat file) */
    TDomain Domain(argv[1]);

    //set PROBLEM_TYPE to NSE if not yet set (3 means Stokes, 5 Naver-Stokes)
    if (TDatabase::ParamDB->PROBLEM_TYPE != 3 && TDatabase::ParamDB->PROBLEM_TYPE != 5)
        TDatabase::ParamDB->PROBLEM_TYPE = 5;
    //open OUTFILE, this is where all output is written to (addionally to console)
    OpenFiles();

    // possibly change parameters in the database, if they are not meaningful now
    Database.CheckParameterConsistencyNSE();
    // write all Parameters to the OUTFILE (not to console) for later reference
    Database.WriteParamDB(argv[0]);

    /* include the mesh from a mesh generator, for a standard mesh use the
     * build-in function. The GEOFILE describes the boundary of the domain. */
    Domain.Init(NULL, TDatabase::ParamDB->GEOFILE); // call mesh generator

    // refine grid up to the coarsest level
    for (int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain.RegRefineAll();

    // write grid into an Postscript file
    if (TDatabase::ParamDB->WRITE_PS)
        Domain.PS("Domain.ps", It_Finest, 0);

    // create output directory, if not already existing
    if (TDatabase::ParamDB->WRITE_VTK)
        mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

    //======================================================================
    // create an example (depending on TDatabase::ParamDB->EXAMPLE)
    Example_NSE2D example;

    // for adaptive grid refinement:
    int current_estimator = TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION;
    // flag for either conforming closures or hanging nodes
    bool conforming = TDatabase::ParamDB->GRID_TYPE != 0;

    NSEErrorEstimator2D estimator {example, Domain, TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION, TDatabase::ParamDB->PROBLEM_TYPE != 3 };
    RefinementStrategy refinementStrategy;
    TAuxParam2D aux;

    int LEVELS = TDatabase::ParamDB->LEVELS;
    for (int curr_level = 0; curr_level < LEVELS; curr_level++) {
        double t_start_level = GetTime(); // time for level
        OutPut(endl << endl << "***********************************************************" << endl);
        OutPut("GEOMETRY  LEVEL " << curr_level << endl);


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

        NSE2D ns(Domain, example);
        ns.assemble();
        // if solution was not zero up to here, you should call
        //ns.assemble_nonlinear_term();

        OutPut("Nonlinear iteration step   0\t");
        ns.stopIt(0);

        //======================================================================
        // nonlinear loop
        // in function 'stopIt' termination condition is checked
        for (unsigned int k = 1; ; k++) {
            ns.solve();

            //no nonlinear iteration for Stokes problem
            if (TDatabase::ParamDB->PROBLEM_TYPE == 3)
                break;

            ns.assemble_nonlinear_term();

            OutPut("nonlinear iteration step " << setw(3) << k << "\t");
            if (ns.stopIt(k))
                break;
        } // end for k

        if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
            ns.get_pressure().project_into_L20(0.0);

        ns.output(curr_level);

        {
            estimator.estimate(ns.get_velocity(), ns.get_pressure(), aux);
            refinementStrategy.applyEstimator(estimator);
            OutPut("estimated global error " << setw(10) << estimator.GetEstimatedGlobalError()[current_estimator] << endl);
        }
        OutPut("Time for level " << curr_level << ": " << GetTime() - t_start_level << endl);
    }

    CloseFiles();
    return 0;
} // end main
