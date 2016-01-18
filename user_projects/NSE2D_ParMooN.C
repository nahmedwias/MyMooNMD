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
#include <NSEErrorEstimator2DDWR.h>
#include "NSEErrorEstimator2D.h"

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{

    //  declaration of database, you need this in every program
    TDatabase Database;
    TFEDatabase2D FEDatabase;

    /** set variables' value in TDatabase using argv[1] (*.dat file) */
    TDomain Domain(argv[1]);

    //set PROBLEM_TYPE to NSE if not yet set (3 means Stokes, 5 Naver-Stokes)
    if(TDatabase::ParamDB->PROBLEM_TYPE!=3 && TDatabase::ParamDB->PROBLEM_TYPE!=5)
        TDatabase::ParamDB->PROBLEM_TYPE = 5;
    //open OUTFILE, this is where all output is written to (addionally to console)
    Output::set_outfile(TDatabase::ParamDB->OUTFILE);

    // possibly change parameters in the database, if they are not meaningful now
    Database.CheckParameterConsistencyNSE();
    // write all Parameters to the OUTFILE (not to console) for later reference
    Database.WriteParamDB(argv[0]);

    /* include the mesh from a mesh generator, for a standard mesh use the
     * build-in function. The GEOFILE describes the boundary of the domain. */
    Domain.Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); // call mesh generator

    // refine grid up to the coarsest level
    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
        Domain.RegRefineAll();

    // write grid into an Postscript file
    if(TDatabase::ParamDB->WRITE_PS)
        Domain.PS("Domain.ps",It_Finest,0);

    // create output directory, if not already existing
    if(TDatabase::ParamDB->WRITE_VTK)
        mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

    Example_NSE2D example;

    // flag for either conforming closures or hanging nodes
    bool conforming = TDatabase::ParamDB->GRID_TYPE != 0;

    std::function<std::vector<double>(const TFEVectFunct2D*, const TFEFunction2D*, double, double, TBaseCell&)> functional = [](const TFEVectFunct2D *u, const TFEFunction2D *p, double x, double y, TBaseCell &cell) {
        std::vector<double> result = std::vector<double>(6);

        result[0] = x*x + y*y <= 3 ? 5 : 0; // f1
        result[1] = x*x + y*y <= 3 ? 5 : 0; // f2
        result[2] = 0; // grad f1
        result[3] = 0; // grad f1
        result[4] = 0; // grad f2
        result[5] = 0; // grad f2

        return result;
    };

    //NSEErrorEstimator2DDWR estimator = NSEErrorEstimator2DDWR(example, functional, Domain);
    NSEErrorEstimator2D estimator {example};
    // for adaptive grid refinement
    auto current_estimator = int(estimator.GetEstimatorType());
    RefinementStrategy refinementStrategy;
    TAuxParam2D aux;

    {
        int LEVELS = TDatabase::ParamDB->LEVELS;
        for (int curr_level = 0; curr_level < LEVELS; curr_level++) {
            double t_start_level = GetTime(); // time for level
            {
                std::stringstream out;
                out << std::endl << std::endl << "***********************************************************" << std::endl;
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

            // create an object of the Navier-Stokes class
            NSE2D ns(Domain, example);
            ns.assemble();
            // if solution was not zero up to here, you should call
            //ns.assemble_nonlinear_term();

            ns.stopIt(0);
            Output::print<1>("Nonlinear iteration step   0\t", ns.getResiduals());

            //======================================================================
            // nonlinear loop
            // in function 'stopIt' termination condition is checked
            for(unsigned int k = 1;; k++)
            {
                ns.solve();

                //no nonlinear iteration for Stokes problem
                if(TDatabase::ParamDB->PROBLEM_TYPE == 3)
                    break;

                ns.assemble_nonlinear_term();

                Output::print<1>("nonlinear iteration step ", setw(3), k, "\t",
                                 ns.getResiduals());
                if(ns.stopIt(k))
                    break;
            } // end for k

            if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
                ns.get_pressure().project_into_L20(0.0);

            ns.output(curr_level);

            double pmin, pmax;
            ns.get_pressure().MinMax(pmin, pmax);
            Output::print<1>("pmin ", pmin, " pmax ", pmax);

            {
                estimator.estimate(ns.get_velocity(), ns.get_pressure(), aux);
                refinementStrategy.applyEstimator(estimator);
                Output::print<1>("estimated global error ", setw(10), estimator.GetEstimatedGlobalError()[current_estimator]);
            }
            Output::print<1>("Time for level ", curr_level, ": ", GetTime() - t_start_level);
        }
    }




    Output::close_file();
    return 0;
} // end main