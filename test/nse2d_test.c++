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
#include <ParameterDatabase.h>
#include <FEDatabase2D.h>
#include <NSE2D.h>
#include <Example_NSE2D.h>
#include <Multigrid.h>
#include <Chrono.h>
#include <algorithm>
#include "LocalAssembling.h"

double accuracy = 1e-6;

void compare(const NSE2D& nse2d, std::array<double, int(5)> errors)
{
  std::array<double, int(5)> computed_errors = nse2d.get_errors();
  
  // check the L2-error of the velcoity
  if( fabs(computed_errors[0]-errors[0]) > accuracy )
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the L2-error of the divergence of the velocity
  if( fabs(computed_errors[1] - errors[1]) > accuracy )
  {
    ErrThrow("L2 norm of divergence: ", std::setprecision(14), computed_errors[1], "  ", errors[1]);
  }
  // check the H1-error of the velcoity
  if( fabs(computed_errors[2] - errors[2]) > accuracy )
  {
    ErrThrow("H1 norm of velocity: ", computed_errors[2], "  ", errors[2]);
  }
  // check the L2-error of the pressure
  if( fabs(computed_errors[3] - errors[3]) > accuracy)
  {
    ErrThrow("L2 norm of pressure: ", computed_errors[3], "  ", errors[3]);
  }
  // check the H1-error of the pressure
  if(fabs(computed_errors[4] - errors[4]) > accuracy )
  {
    ErrThrow("H1 norm of pressure: ", computed_errors[4], "  ", errors[4]);
  }
}

void compute(TDomain &domain, ParameterDatabase& db,
             std::array<double, int(5)> errors)
{
  NSE2D nse2d(domain, db);
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

void check(TDomain &domain, ParameterDatabase db,
           int velocity_order, int nstype, int laplace_type,
           int nonlinear_form,
           std::array<double, int(5)> errors)
{
  Output::print("\n\nCalling check with velocity_order=", velocity_order,
                ", nstype=", nstype, ", laplace_type=", laplace_type, 
                ", and nonlinear_form=", nonlinear_form);
  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(Multigrid::default_multigrid_database());
  db["problem_type"] = 5;
  db["solver_type"] = "direct";
  db["iterative_solver_type"] = "fgmres";
  db["residual_tolerance"] = 1.e-12;
  db["preconditioner"] = "least_squares_commutator";
  
  db["nonlinloop_maxit"] = 50;
  db["nonlinloop_epsilon"] = 1e-10;

  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->NSTYPE = nstype;
  TDatabase::ParamDB->LAPLACETYPE = laplace_type;
  TDatabase::ParamDB->NSE_NONLINEAR_FORM = nonlinear_form;
  
  Chrono timer;
  compute(domain, db, errors);
  timer.restart_and_print("nse2d direct solver,                    velocity "
                          + std::to_string(velocity_order) + ", nstype "
                          + std::to_string(nstype));
  
  // we have to reset the space codes because they are changed in nse2d
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  
  db["solver_type"] = "iterative";
  compute(domain, db, errors);
  timer.restart_and_print("nse2d fgmres(lsc preconditioner),       velocity "
                          + std::to_string(velocity_order) + ", nstype "
                          + std::to_string(nstype));
  
  // we have to reset the space codes because they are changed in nse2d
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  

  db["preconditioner"] = "multigrid";
  db["multigrid_n_levels"] = db["refinement_n_initial_steps"].get<size_t>();
  //choose smoother on fine grid according to element
  std::vector<int> disc_p = {12,13,14,15,22,23,24};
  if(std::find(disc_p.begin(), disc_p.end(), velocity_order) != disc_p.end())
    db["multigrid_smoother"] = "cell_vanka_store";
  else
    db["multigrid_smoother"] = "patch_vanka_store";

  db["multigrid_type"] = "standard";
  db["multigrid_smoother_coarse"] = "direct_solve";
  db["multigrid_n_pre_smooth"] = 0;
  db["multigrid_n_post_smooth"] = 1;
  db["multigrid_correction_damp_factor"] = 1.0;
  db["multigrid_vanka_damp_factor"] = 1.0;
  compute(domain, db, errors);
  timer.restart_and_print("nse2d fgmres(multigrid preconditioner), velocity "
                          + std::to_string(velocity_order) + ", nstype "
                          + std::to_string(nstype));
}

template <int laplace_type>
void check_one_element(TDomain& domain, ParameterDatabase db,
                       int velocity_order,
                       std::array<double, int(5)> errors)
{
  int nonlinear_form = 0;
  if(laplace_type == 0)
  {
    // NSTYPE = 1
    check(domain, db, velocity_order, 1, laplace_type, nonlinear_form, errors);
    // NSTYPE = 2
    check(domain, db, velocity_order, 2, laplace_type, nonlinear_form, errors);
  }
  // NSTYPE = 3
  check(domain, db, velocity_order, 3, laplace_type, nonlinear_form, errors);
  // NSTYPE is 4
  check(domain, db, velocity_order, 4, laplace_type, nonlinear_form, errors);
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

    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(ParameterDatabase::default_nonlinit_database());
    db.merge(ParameterDatabase::default_output_database());
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());

    db["problem_type"].set<size_t>(5);
    db["example"] = 2;

    db.add("refinement_n_initial_steps", (size_t) 2,"");
    
    db["nonlinloop_maxit"] = 100;
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_slowfactor"] = 1.;

    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "TwoTriangles", "", {"UnitSquare", "TwoTriangles"});
    
    // default construct a domain object
    TDomain domain(db);

    db["reynolds_number"] = 1;
    TDatabase::ParamDB->FLOW_PROBLEM_TYPE=5;
    db["space_discretization_type"] = "galerkin";
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE = 0;
    
    // possibly parameters in the database
    check_parameters_consistency_NSE(db);
    // refine grid
    size_t n_ref = domain.get_n_initial_refinement_steps();
    for(unsigned int i=0; i < n_ref; i++)
    {
      domain.RegRefineAll();
    }
    std::array<double, int(5)> errors;
    
    //=========================================================================
    Output::print<1>("\nTesting the P2/P1 elements");
    errors = {{ 0.005005607397208, 0.085596036177161, 0.15666212408257,
                0.071089608676709, 1.3407900222228 }};
    // VELOCITY_SPACE = 2 and the pressure space is chosen in the class NSE2D
    check_one_element<0>(domain, db, 2, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the P3/P2 elements");
    errors = {{ 0.00028020829779642, 0.0061977013305504, 0.011391210186952,
                0.0053396813413967, 0.20779333160236 }};
    // VELOCITY_SPACE = 3 and the pressure space is chosen in the class NSE2D
    check_one_element<0>(domain, db, 3, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the P4/P3 elements");
    errors = {{ 1.1817023010728e-05, 0.00029575269112296, 0.0006435418450572, 
                0.00050496270735108, 0.026998702772064 }};
    // VELOCITY_SPACE = 4 and the pressure space is chosen in the class NSE2D
    check_one_element<0>(domain, db, 4, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the P5/P4 elements");
    errors = {{ 4.3466391168252e-07, 1.4197002168949e-05, 2.793323812439e-05,
                2.1824211773585e-05, 0.0016936362911126 }};
    // VELOCITY_SPACE = 5 and the pressure space is chosen in the class NSE2D
    check_one_element<0>(domain, db, 5, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the P2-bubble/P1-disc elements");
    errors = {{ 0.0071886299046824, 0.081842610256743, 0.21185558654462,
                0.36754876295023, 5.3058522557418 }};
    // VELOCITY_SPACE = 22 and the pressure space is chosen in the class NSE2D
    check_one_element<0>(domain, db, 22, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the P3-bubble/P2-disc elements");
    errors = {{ 0.00026037876329326, 0.0057269486642841, 0.010856952083039,
                0.013201231959059, 0.39888576555041 }};
    // VELOCITY_SPACE = 23 and the pressure space is chosen in the class NSE2D
    check_one_element<0>(domain, db, 23, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the P4-bubble/P3-disc elements");
    errors = {{ 1.2032722339771e-05, 0.00024376809624887, 0.00055164963287203,
                0.00063706731983293, 0.027783948983068 }};
    // VELOCITY_SPACE = 24 and the pressure space is chosen in the class NSE2D
    check_one_element<0>(domain, db, 24, errors);
  } // end program 1
  //=========================================================================
  /** Program 2
   *  This program tests direct solve with galerkin discretization
   * direct solver; Test for Quad's
   * The VELOCITY_SPACE is fixed and tests for different NSTYPE's
   */
  {
    //  declaration of databases
    TDatabase Database;
    TFEDatabase2D FEDatabase;

    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(ParameterDatabase::default_nonlinit_database());
    db.merge(ParameterDatabase::default_output_database());
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());
    db["problem_type"].set<size_t>(5);
    db["example"] = 2;

    db.add("refinement_n_initial_steps", (size_t) 2,"");
    
    db["nonlinloop_maxit"] = 100;
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_slowfactor"] = 1.;

    // default construct a domain object
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    TDomain domain(db);

    // parameters used for this test
    db["reynolds_number"] = 1;
    TDatabase::ParamDB->FLOW_PROBLEM_TYPE=5;
    db["space_discretization_type"] = "galerkin";
    TDatabase::ParamDB->LAPLACETYPE = 0;
        
    // possibly parameters in the database
    check_parameters_consistency_NSE(db);

    // refine grid
    size_t n_ref = domain.get_n_initial_refinement_steps();
    for(unsigned int i=0; i < n_ref; i++)
    {
      domain.RegRefineAll();
    }
    std::array<double, int(5)> errors;
    
    //=========================================================================
    Output::print<1>("\nTesting the Q2/Q1 elements");
    errors = {{ 0.004083204524442, 0.05101317284072, 0.10522635824261,
                0.017686667902813, 0.51182308944019 }};
    // VELOCITY_SPACE  = 2
    check_one_element<0>(domain, db, 2, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the Q3/Q2 elements");
    errors = {{ 0.00019319433716041, 0.003563918945437, 0.0071078507849009,
                0.0018446328461379, 0.057123632497266 }};
    // VELOCITY_SPACE  = 3
    check_one_element<0>(domain, db, 3, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the Q4/Q3 elements");
    errors = {{ 6.9434747253041e-06, 0.0001677182257252, 0.00035212311261646,
                9.4703269177756e-05, 0.0048368160352994 }};
    // VELOCITY_SPACE  = 4
    check_one_element<0>(domain, db, 4, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the Q5/Q4 elements");
    errors = {{ 2.3951237974726e-07, 6.8246180893733e-06, 1.4163394749406e-05,
                4.9000676557526e-06, 0.00030469183949993 }};
    // VELOCITY_SPACE  = 5
    check_one_element<0>(domain, db, 5, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the Q2/P1-disc elements");
    errors = {{ 0.0040557267369371, 0.05094007770388, 0.10564123627325,
                0.030975452768144, 0.70842776239597 }};
    // VELOCITY_SPACE  = 22
    check_one_element<0>(domain, db, 22, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the Q3/P2-disc elements");
    errors = {{ 0.00020642736694367, 0.0042373143298489, 0.0075794259144329,
                0.0039793392063018, 0.13655739379564 }};
    // VELOCITY_SPACE  = 23
    check_one_element<0>(domain, db, 23, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the Q4/P3-disc elements");
    errors = {{ 9.9544541389566e-06, 0.00032338687271106, 0.00045888380841742,
                0.00037625559674667, 0.01835259073302 }};
    // VELOCITY_SPACE  = 24
    check_one_element<0>(domain, db, 24, errors);
    
    //=========================================================================
    Output::print<1>("\nTesting the Q5/P4-disc elements");
    errors = {{ 5.3747286967048e-07, 2.4122655352531e-05, 2.7997005479668e-05,
                2.81141239001e-05, 0.0018176331985532 }};
    // VELOCITY_SPACE  = 25
    check_one_element<0>(domain, db, 25, errors);
  }
  
  return 0;
}
