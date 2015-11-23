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

#include <sys/stat.h>

#include <CD2D.h>
#include <CDErrorEstimator2D.h>

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
    OpenFiles(); // call CloseFiles() at the end of the program

    // write all Parameters to the OUTFILE (not to console) for later reference
    Database.WriteParamDB(argv[0]);

    /* include the mesh from a mesh generator, for a standard mesh use the
     * build-in function. The GEOFILE describes the boundary of the domain. */
    Domain.Init(NULL, TDatabase::ParamDB->GEOFILE); // call mesh generator

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
    CDErrorEstimator2D estimator {example, Domain.GetCollection(It_Finest, 0)[0], TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION };

    // refinement
    RefinementStrategy refinementStrategy;

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

        // run cd2d on current grid
        CD2D cd2d(Domain, example);
        cd2d.assemble();
        cd2d.solve();
        cd2d.output(curr_level);

        {
            estimator.estimate(cd2d.get_function());
            refinementStrategy.applyEstimator(estimator);
            OutPut("estimated global error " << setw(10) << estimator.GetEstimatedGlobalError()[current_estimator] << endl);
        }
        OutPut("Time for level " << curr_level << ": " << GetTime() - t_start_level << endl);
    }

    CloseFiles();
    return 0;
} // end main
