/**
 * @brief A test program for the solving of NSE2D problems.
 *
 * This serves as a test for the solving of NSE2D problems. It is intended to
 * perform NSE2D calculations with different examples in different setups to test
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
 * @author Naveed, Ulrich, Clemens
 *
 */

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <NSE2D.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_NSE2D.h>

#include <MainUtilities.h> //for error measuring

void compare(const NSE2D& nse2d, std::array<double, int(6)> errors)
{
  std::array<double, int(6)> computed_errors;
  computed_errors = nse2d.get_errors();
  
  // check the L2-error of the velcoity
  if( fabs(computed_errors[0]-errors[0]) > 1e-6 )
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the H1-error of the velcoity
  if( fabs(computed_errors[1] - errors[1]) > 1e-6 )
  {
    ErrThrow("H1 norm of velocity: ", computed_errors[1], "  ", errors[1]);
  }    
  // check the L2-error of the pressure
  if( fabs(computed_errors[2] - errors[2]) > 1e-6)
  {
    ErrThrow("L2 norm of pressure: ", computed_errors[2], "  ", errors[2]);
  }
  // check the H1-error of the pressure
  if(fabs(computed_errors[3] - errors[3]) > 1e-6 )
  {
    ErrThrow("H1 norm of pressure: ", computed_errors[3], "  ", errors[3]);
  }  
}
void check(TDomain &domain, int velocity_order, int nstype, int laplace_type,
           int nonlinear_form,
           std::array<double, int(6)> errors)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->NSTYPE = nstype;
  TDatabase::ParamDB->LAPLACETYPE = laplace_type;
  TDatabase::ParamDB->NSE_NONLINEAR_FORM = nonlinear_form;
  
  NSE2D nse2d(domain);
  nse2d.assemble();
  
  // check stopping criterion
  nse2d.stopIt(0);
  for(unsigned int k=1;; k++)
  {
    Output::print<1>("nonlinear step " , setw(3), k-1, "\t",
                     nse2d.getResiduals());
    nse2d.solve();
    
    // checking the first nonlinear iteration    
    nse2d.assemble_nonlinear_term();;
    if(nse2d.stopIt(k))
      break;
  }
  nse2d.output();  
  // compare now the errors
  compare(nse2d, errors);
}

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  /** Program 1
   *  This program tests direct solve with galerkin discretization
   * direct solver; Tests for Triangles
   * The VELOCITY_SPACE is fixed and tests for different NSTYPE's
   */
  { //  declaration of databases
    TDatabase Database;
    TFEDatabase2D FEDatabase;
    // default construct a domain object
    TDomain domain;

    TDatabase::ParamDB->PROBLEM_TYPE = 5; //NSE Problem
    TDatabase::ParamDB->EXAMPLE = 2; 
    TDatabase::ParamDB->UNIFORM_STEPS = 3;
    TDatabase::ParamDB->RE_NR=1;
    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LEVELS =1;
    TDatabase::ParamDB->SOLVER_TYPE = 2;
    TDatabase::ParamDB->LAPLACETYPE = 0;
    TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE= 1e-10;
    TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE=50;
    TDatabase::ParamDB->MEASURE_ERRORS = 1;
    
    // possibly parameters in the database
    Database.CheckParameterConsistencyNSE();
    // the domain is initialised with default description and default
    // initial mesh
    domain.Init((char*)"Default_UnitSquare", (char*)"TwoTriangles");
    // refine grid up to the coarsest level
    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    {
      domain.RegRefineAll();
    }
    //===========================================================
    Output::print<1>("Testing the P2/P1 elements");
    //===========================================================
    std::array<double, int(6)> errors;
    errors = {{0.000610487, 0.0389713, 0.0107332,  0.512621}};
    // VELOCITY_SPACE = 2 and the pressure space is chosen in the class NSE2D
    // NSTYPE = 1
    check(domain, 2, 1, 0, 0, errors);
    // NSTYPE = 2
    check(domain, 2, 2, 0, 0, errors);
    // NSTYPE = 3
    check(domain, 2, 3, 0, 0, errors);
    // NSTYPE is 4
    check(domain, 2, 4, 0, 0, errors);
    //===========================================================
    Output::print<1>("Testing the P3/P2 elements");
    //===========================================================
    errors = {{1.84056e-05, 0.00142784, 0.000556315, 0.0430551}};
    // NSTYPE = 1
    check(domain, 3, 1, 0, 0, errors);
    // NSTYPE = 2
    check(domain, 3, 2, 0, 0, errors);
    // NSTYPE = 3
    check(domain, 3, 3, 0, 0, errors);
    // NSTYPE = 4
    check(domain, 3, 4, 0, 0, errors);
    
  } // end program 1
//================================================================================  
  /** Program 2
   *  This program tests direct solve with galerkin discretization
   * direct solver; Test for Quad's
   * The VELOCITY_SPACE is fixed and tests for different NSTYPE's
   */
  {
    //  declaration of databases
    TDatabase Database;
    TFEDatabase2D FEDatabase;
    // default construct a domain object
    TDomain domain;
    // parameters used for this test
    TDatabase::ParamDB->PROBLEM_TYPE = 5; //NSE Problem
    TDatabase::ParamDB->EXAMPLE = 2; 
    TDatabase::ParamDB->UNIFORM_STEPS = 4;
    TDatabase::ParamDB->RE_NR=1;
    TDatabase::ParamDB->DISCTYPE=1;
    TDatabase::ParamDB->LEVELS =1;
    TDatabase::ParamDB->SOLVER_TYPE = 2;
    TDatabase::ParamDB->LAPLACETYPE = 0;
    TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE= 1e-10;
    TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE=50;
    TDatabase::ParamDB->MEASURE_ERRORS = 1;
        
    // possibly parameters in the database
    Database.CheckParameterConsistencyNSE();

    // the domain is initialised with default description and default
    // initial mesh
    domain.Init((char*)"Default_UnitSquare", (char*)"UnitSquare");
    // refine grid up to the coarsest level
    for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    {
      domain.RegRefineAll();
    }
    //===========================================================
    Output::print<1>("Testing the Q2/P1-disc elements");
    //===========================================================
    std::array<double, int(6)> errors;
    errors = {{6.37755e-05, 0.00661116, 0.00190367,  0.17806}};
    // VELOCITY_SPACE  = 12
    // NSTYPE = 1
    check(domain, 12, 1, 0, 0, errors);
    // NSTYPE = 2
    check(domain, 12, 2, 0, 0, errors);
    // NSTYPE = 3
    check(domain, 12, 3, 0, 0, errors);
    // NSTYPE 4 
    check(domain, 12, 4, 0, 0, errors);
    //===========================================================
    Output::print<1>("Testing the Q3/P2-disc elements");
    //===========================================================    
    errors = {{8.44484e-07, 0.000121079, 6.07644e-05,  0.00859989}};
    // VELOCITY_SPACE  = 12
    // NSTYPE = 1
    check(domain, 13, 1, 0, 0, errors);
    // NSTYPE = 2
    check(domain, 13, 2, 0, 0, errors);
    // NSTYPE = 3
    check(domain, 13, 3, 0, 0, errors);
    // NSTYPE 4 
    check(domain, 13, 4, 0, 0, errors);
    //===========================================================
    Output::print<1>("Testing the Q4/P3-disc elements");
    //===========================================================
    errors = {{1.02856e-08, 1.81526e-06, 1.45146e-06, 0.000287603}};
    // VELOCITY_SPACE  = 12
    // NSTYPE = 1
    check(domain, 14, 1, 0, 0, errors);
    // NSTYPE = 2
    check(domain, 14, 2, 0, 0, errors);
    // NSTYPE = 3
    check(domain, 14, 3, 0, 0, errors);
    // NSTYPE 4 
    check(domain, 14, 4, 0, 0, errors);
    //===========================================================
    Output::print<1>("Testing the Q5/P4-disc elements");
    //===========================================================    
    errors = {{1.41762e-10, 2.90308e-08, 2.7216e-08,  7.02901e-06}};
    // VELOCITY_SPACE  = 12
    // NSTYPE = 1
    check(domain, 15, 1, 0, 0, errors);
    // NSTYPE = 2
    check(domain, 15, 2, 0, 0, errors);
    // NSTYPE = 3
    check(domain, 15, 3, 0, 0, errors);
    // NSTYPE 4 
    check(domain, 15, 4, 0, 0, errors);
  }
  
  return 0;
}
