/**
 * @brief A test program for the solving of CD2D problems.
 *
 * This serves as a test for the solving of CD2D problems. It is intended to
 * perform CD2D calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
 * So far only one such test is implemented.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes wrong)
 * the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then this tests
 * shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 *
 * @date 2015/11/11
 * @author Clemens Bartsch
 *
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <CD2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_CD2D.h>

#include <MainUtilities.h> //for error measuring

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  /** Program 1
   *  This program tests direct solve with galerkin discretization and
   *  fem-tvd-type algebraic flux correction.
   */
  {

    //  declaration of databases
    TDatabase Database;
    TFEDatabase2D FEDatabase;

    // default construct a domain object
    TDomain domain;
    
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db["problem_type"] = 1;
    db.add("solver_type", std::string("direct"), "");

    // Set Database values (this is what is usually done by the input-file)
    TDatabase::ParamDB->PROBLEM_TYPE = 1; //CDR Problem
    TDatabase::ParamDB->EXAMPLE = 3; //Sharp Boundary Layer Example
    TDatabase::ParamDB->UNIFORM_STEPS = 1;
    TDatabase::ParamDB->LEVELS = 2;
    TDatabase::ParamDB->ANSATZ_ORDER = 1; //P1 elements
    TDatabase::ParamDB->DISCTYPE = 1; //Galerkin Desicreitzation
    TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION = 1; //FEM-TVD type afc
    TDatabase::ParamDB->MEASURE_ERRORS = 1;
    TDatabase::ParamDB->DELTA0 = 0.3;
    TDatabase::ParamDB->DELTA1 = 0.;
    TDatabase::ParamDB->SDFEM_TYPE = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT = 1;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0.5;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;
    TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = 1;
    TDatabase::ParamDB->SC_VERBOSE = 0; // supress solver output
    TDatabase::ParamDB->SOLVER_TYPE = 2; // use direct solver
    
    // the domain is initialised with default description and default
    // initial mesh
	  domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");

    // refine grid up to the coarsest level
    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS + TDatabase::ParamDB->LEVELS; i++)
    {
      domain.RegRefineAll();
    }

    //Here the actual computations take place
    //=========================================================================
    CD2D cd2d(domain, db);
    cd2d.assemble();
    cd2d.solve();
    //=========================================================================

    // instead of calling CD2D::output() - compare errors against reference errors
    // TODO (when there is more than one programming running, this will be encapsulated
    // and out-sourced)

    double errors[5];
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    const TFEFunction2D& function = cd2d.get_function();
    const TFESpace2D* space = function.GetFESpace2D();

    function.GetErrors(cd2d.get_example().get_exact(0), 3, AllDerivatives, 4,
                       SDFEMErrors, cd2d.get_example().get_coeffs(), &aux, 1,
                       &space, errors);
    // check L2 error
    if( fabs(errors[0] - 0.500273) > 1e-6 )
    {
      ErrThrow("Program 1: L2 norm not correct.");
    }
    // check H1-semi error
    if( fabs(errors[1] - 1.90143) > 1e-5 )
    {
      ErrThrow("Program 1: H1-semi norm not correct.");
    }
    // check SD error
    if( fabs(errors[2] - 0.275256) > 1e-6)
    {
      ErrThrow("Program 1: SD norm not correct.");
    }
    // check L_inf error
    if(fabs( errors[3] - 0.999909) > 1e-6 )
    {
      ErrThrow("Program 1: L_inf norm not correct.");
    }

  } // end program 1

  return 0;
}
