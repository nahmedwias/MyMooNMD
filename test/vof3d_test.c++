/**
 * @brief A test program for the solving of VOF3D problems.
 *

 *
 * So far the test is adapted to:
 *  - examples number ......
 *  - the finite elements ......
 *  - the nstype .....
 *  - the time stepping ....
 *  - the discretization type ....
 *  - the solver .....
 *  - the geometry (quads,trias)....
 *  - the laplace type ....
 *  - the couplings ...
 *
 * @author Najib Alia
 *
 * @date 2017/04/26
 */

#include <VOF_TwoPhase3D.h>
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <TimeDiscRout.h>
#include <AlgebraicFluxCorrection.h>

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

void check_errors(const VOF_TwoPhase3D& vof3d,
                  std::array<double, int(6)> errors,
                  double tol)
{
  std::array<double, int(6)> computed_errors_tnse ;
  computed_errors_tnse = vof3d.tnse3d_.get_errors();
  // check the Errors of velocity and pressure
  if( fabs(computed_errors_tnse[0]-errors[0]) > tol ||
      computed_errors_tnse[0] != computed_errors_tnse[0]) //check for nan!
    ErrThrow("L2 norm of velocity: ", computed_errors_tnse[0], "  ", errors[0]);

  if( fabs(computed_errors_tnse[1] - errors[1]) > tol )
    ErrThrow("H1 norm of velocity: ", computed_errors_tnse[1], "  ", errors[1]);

  if( fabs(computed_errors_tnse[2] - errors[2]) > tol)
    ErrThrow("L2 norm of pressure: ", computed_errors_tnse[2], "  ", errors[2]);

  if(fabs(computed_errors_tnse[3] - errors[3]) > tol )
    ErrThrow("H1 norm of pressure: ", computed_errors_tnse[3], "  ", errors[3]);

  std::array<double, int(3)> computed_errors_tcd;
  computed_errors_tcd = vof3d.phaseconvection3d_.get_errors();
  // check the Errors of concentration
  if( fabs(computed_errors_tcd[0]-errors[4]) > tol )
    ErrThrow("L2t norm of concentration: ", computed_errors_tcd[0], "  ", errors[4]);
  if( fabs(computed_errors_tcd[1]-errors[5]) > tol )
    ErrThrow("L2H1t norm of concentration: ", computed_errors_tcd[1], "  ", errors[5]);
}


void compute(ParameterDatabase& db,
             std::array<double, int(6)> errors,
             double tol)
{
#ifdef _MPI
  int size, my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  db["verbosity"]  = 5;
  ParameterDatabase tnse_db("Navier-Stokes Database");
  ParameterDatabase tcd_db("Convection Diffusion Database");
  tnse_db.merge(db,true);
  tcd_db.merge(db,true);
  tcd_db["problem_type"]     = 2;
  tnse_db["problem_type"]    = 6;

  check_parameters_consistency_NSE(tnse_db);
  //  declaration of databases
  TDomain domain(db);
  // Initial refinement and grid collection
#ifdef _MPI
   int maxSubDomainPerDof = 0;
#endif
   std::list<TCollection*> gridCollections
      = domain.refine_and_get_hierarchy_of_collections(
        db
#ifdef _MPI
        , maxSubDomainPerDof
#endif
     );

  /********************************************************************
   * Creating VOF object, which contains both TimeNSE2D and TimeCD2D
   ********************************************************************/
  SetTimeDiscParameters(0);                      // Initialize parameters for time discretization
#ifdef _MPI
  VOF_TwoPhase3D vof(gridCollections,tnse_db,tcd_db, maxSubDomainPerDof);
#else
  VOF_TwoPhase3D vof(gridCollections,tnse_db,tcd_db);
#endif

  vof.manage_example_parameters();
  vof.update_field_vectors();

  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=0;
  vof.tnse3d_.assemble_initial_time();      // assemble linear term

  if (!tcd_db["algebraic_flux_correction"].is("none"))
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=1;
  if (vof.solve_convection_ == true)
  {
    if (vof.nse2cd_coupling_ == true)
    {
      vof.phaseconvection3d_.assemble_initial_time_with_convection(&vof.tnse3d_.get_velocity());
    }
    else
      vof.phaseconvection3d_.assemble_initial_time();
  }

   double end_time = TDatabase::TimeDB->ENDTIME;
   int step = 0, image = 0;
   int n_substeps = GetN_SubSteps();
   vof.tnse3d_.current_step_ = 0;
   TDatabase::TimeDB->CURRENTTIME = 0.0;
  /********************************************************************
   * TIME ITERATION LOOP
   ********************************************************************/
   while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
   {
     step++;
     vof.tnse3d_.current_step_++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
     for(int j=0; j < n_substeps; ++j)
     {
       SetTimeDiscParameters(1);            // setting the time disc parameters
       double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
       TDatabase::TimeDB->CURRENTTIME += tau;
       TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=0;
       vof.tnse3d_.assemble_rhs();
       if (vof.tnse_variable_fluid_ == true)
       {
         ErrThrow("NOT IMPLEMENTED YET!");
//        vof.tnse2d_.assemble_massmatrix_withfields(&vof.rho_fefunction_);
//        if( TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION == 1 )
//            vof.tnse2d_.apply_slip_penetration_bc(true,true);
       }

       vof.tnse3d_.assemble_nonlinear_term();
       vof.tnse3d_.assemble_system();

       /********************************************************************
        * NON LINEAR LOOP
        ********************************************************************/
       for(unsigned int k = 0;; k++)
       {
         vof.tnse3d_.compute_residuals();

//         if (my_rank==0) // some outputs
//         {
//           Output::print<1>("\nNONLINEAR ITERATION :", setw(3), k);
//           Output::print<1>("Residuals :", vof.tnse3d_.get_residuals());
//         }

         if(vof.tnse3d_.stop_it(k))
           break;
         vof.tnse3d_.solve();
         if(vof.tnse3d_.imex_scheme(1))
           continue; // this interrupts the NL-Loop
         if (vof.tnse_variable_fluid_ == true)
         {
           ErrThrow("Slip BC Not Implemented yet!");
//        if( TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION == 1 )
//          vof.tnse3d_.apply_slip_penetration_bc(false,false);
         }
         vof.tnse3d_.assemble_nonlinear_term();
         vof.tnse3d_.assemble_system();
       } // end for k, non linear loop

       /********************************************************************
        * SOLVING CD3D WITH NSE3D SOLUTION
        ********************************************************************/
       if (!tcd_db["algebraic_flux_correction"].is("none"))
         TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=1;
       if (vof.solve_convection_ == true )
       {
         if (vof.nse2cd_coupling_ == true)
         {
           ErrThrow("Not implemented yet!");
//        tcd2d.assemble_rhs_vector(&tnse2d.get_velocity()); // once per time step
//        tcd2d.assemble_stiffness_matrix_alone_with_convection(&tnse2d.get_velocity());
//        tcd2d.scale_stiffness_matrix();
        vof.phaseconvection3d_.assemble_with_convection(&vof.tnse3d_.get_velocity());
         }
         else  // if we solve TCD2D standard, without any coupling
         { vof.phaseconvection3d_.assemble();}

         int ipp_was_1 = 0;
         if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE==1)
           {ipp_was_1 = 1;
           TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;}
         vof.phaseconvection3d_.solve();
         if (ipp_was_1)
           {TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;}
         vof.phaseconvection3d_.descale_stiffness(); //needed once per time loop

         /********************************************************************
          * UPDATING VELOCITY VECTOR WITH CD3D SOLUTION
          ********************************************************************/
         if (vof.cd2nse_coupling_ == true )
           vof.update_field_vectors();
       }
       vof.tnse3d_.output(step,image);
       vof.phaseconvection3d_.output(step,image);
     }
#ifdef _MPI
     if (my_rank==0)
     {
     if(TDatabase::TimeDB->CURRENTTIME >= end_time - 1e-10)
       check_errors(vof, errors,tol);
     }
#endif
   } // end for step, time loop
}


void set_errors(int example,
                std::array<double, int(6)>& errors)
{
  switch(example)
  {
    case 10:
      errors= {0.0001551214835, /*L2(u)*/ 0.004173239659, /*H1(u)*/
               0.2130611644, /*L2(p)*/ 2.798457879, /*H1(p)*/
               5.21837e-06, /*L2t(c)*/ 0.000161713  /*L2H1t(c)*/ };
      break;
    default:
      ErrThrow("Unknown example number!");
      break;
  }
}


int main(int argc, char* argv[])
{
  //declaration of databases
  TDatabase Database;

#ifdef _MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  TDatabase::ParamDB->Comm = comm;

  int size;//my_rank,
//  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  TFEDatabase3D FEDatabase;
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(ParameterDatabase::default_output_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(Solver<>::default_solver_database());
  db.merge(Example3D::default_example_database());
  db.merge(AlgebraicFluxCorrection::default_afc_database());

  db["output_compute_errors"] = true;
  db["verbosity"]             = 5;
  db.add("refinement_n_initial_steps", (size_t) 1,"", (size_t) 0, (size_t) 10);;
  db.add("boundary_file", "Default_UnitCube", "", {"Default_UnitCube"});
  db.add("geo_file", "Default_UnitCube_Hexa", "", {"Default_UnitCube_Hexa"});

  std::array<double, int(6)> errors;
  double tol = 1e-6;
  //=======================================================================
  /* ===== EXAMPLE 10   =========
   *
   *     ================= */
  // ======================================================================
  {
    db["example"] = 10;
    db["dimensional_nse"] = false;
    db["solve_cd"] = true;
    db["refinement_n_initial_steps"] = 3;
    TDatabase::TimeDB->STARTTIME     = 0.;
    TDatabase::TimeDB->ENDTIME       = 0.1;
    TDatabase::TimeDB->TIMESTEPLENGTH= 0.1;
    db["fluid_density"]           = 1;  // this is then compared to
    db["fluid_dynamic_viscosity"] = 1;    // TNSE2D with Re=200
    db["reynolds_number"]         = 1;
    TDatabase::ParamDB->RE_NR     = 1;
    db["solver_type"]             = "direct";
    db["direct_solver_type"]      = std::string(argv[1]); //umfpack or mumps
    db["nonlinloop_epsilon"]      = 1.e-7;
    TDatabase::ParamDB->VELOCITY_SPACE = 12;
    TDatabase::ParamDB->PRESSURE_SPACE = -4711;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE = 0;
    TDatabase::ParamDB->ANSATZ_ORDER= 1;
    db["algebraic_flux_correction"] = "none";
    set_errors(10,errors);
    compute(db, errors, tol);
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  }

#ifdef _MPI
  MPI_Finalize();
#endif

  return 0;
}
