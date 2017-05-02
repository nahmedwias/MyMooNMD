// =======================================================================
// Purpose:
//
// Author:      NA
//
// History:     Start 25.04.2017
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <Time_NSE3D.h>
#include <Time_CD3D.h>
#include <Example_TimeNSE3D.h>
#include <Chrono.h>
#include <LoopInfo.h>
#include <ParameterDatabase.h>
#include <TimeDiscRout.h>

#include <VOF_TwoPhase3D.h>

#include <Output3D.h>
#include <MeshPartition.h>
#include <sys/stat.h>

using namespace std;

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif



// ***** MAIN PROGRAM ***** //
int main(int argc, char* argv[])
{
#ifdef _MPI
  //Construct and initialise the default MPI communicator.
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  // Hold mpi rank and size ready, check whether the current processor
  // is responsible for output (usually root, 0).
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(my_rank==0)
  {
    Output::print("<<<<< Running ParMooN: TimeMultiphase3D Main Program >>>>>");
    Output::info("Time_Multiphase3D", "MPI, using ", size, " processes");
  }
#else
  int my_rank = 0;
  Output::print("<<<<< Running ParMooN: TimeMultiphase3D Main Program >>>>>");
  Output::info("Time_Multiphase3D", "SEQUENTIAL (or OMP...)");
#endif
  double t_start = GetTime();

  Chrono  stopwatch;   // Start a stopwatch for time measurement during execution

  TDatabase     Database;        // Initialize User Input Databases.
  TFEDatabase3D FEDatabase;      // Initialize FE Database.

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  ParameterDatabase tnse_db("Navier-Stokes Database");
  ParameterDatabase tcd_db("Convection Diffusion Database");

  std::ifstream fs(argv[1]); parmoon_db.read(fs);  fs.close();
  tnse_db.merge(parmoon_db,true);
  tcd_db.merge(parmoon_db,true);

//  tcd_db["example"]          = -1;
  tcd_db["problem_type"]     = 2;
  tcd_db["output_basename"]  = "multiphase_tconvection_output";
  //  tcd_db["space_discretization_type"] = "galerkin";

  //  tnse_db["example"]         = 18;
  tnse_db["problem_type"]    = 6;
  tnse_db["output_basename"] = "multiphase_tnse_output";


#ifdef _MPI
  TDatabase::ParamDB->Comm = comm;
#endif

  /* =====================================================================
   * set the database values and generate mesh
   * =====================================================================*/
  // Construct domain, thereby read in controls from the input file.
  TDomain domain(argv[1], parmoon_db);      // Initialize geometry

  /********************************************************************
   * WRITE PARAMETERS TO OUTFILE
   ********************************************************************/
  check_parameters_consistency_NSE(tnse_db);
  Output::setVerbosity(parmoon_db["verbosity"]);
  if(my_rank==0)
  {
    Output::set_outfile(parmoon_db["outfile"]);
    parmoon_db.write(Output::get_outfile());
    Database.WriteParamDB(argv[0]);
    Database.WriteTimeDB();
  }

  // Initial refinement and grid collection
#ifdef _MPI
  int maxSubDomainPerDof = 0;
#endif
  std::list<TCollection*> gridCollections
     = domain.refine_and_get_hierarchy_of_collections(
       parmoon_db
#ifdef _MPI
       , maxSubDomainPerDof
#endif
    );
  TCollection* coll = gridCollections.front();
  TOutput3D output(0,0,0,0,std::addressof(domain),coll);

  if (my_rank==0)
    output.WriteVtk("mesh.vtk"); // run SEQ to see the generated mesh
  domain.print_info("Multiphase3D domain");      // Output domain info

  /********************************************************************
   * Creating VOF object, which contains both TimeNSE3D and TimeCD3D
   ********************************************************************/
  SetTimeDiscParameters(0);     // Initialize parameters for time discretization
#ifdef _MPI
  VOF_TwoPhase3D vof(gridCollections,tnse_db,tcd_db, maxSubDomainPerDof);
#else
  VOF_TwoPhase3D vof(gridCollections,tnse_db,tcd_db);
#endif

  vof.manage_example_parameters();
  /* This calculates rho and mu vectors depending on example number
   * Check that the vectors are as expected using "output_vectors(..)" */
  vof.update_field_vectors(); // should be done by each process
  /* SOME OUTPUT AND INFORMATION SET */
  if (my_rank==0)
  {
    vof.output_initial_info();
//    vof.output_vectors("vector_phi_init","vector_rho_init","vector_mu_init");
  }

  /********************************************************************
   * START ASSEMBLING TimeNSE3D WITH GIVEN FIELDS
   ********************************************************************/
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
   * SOME OUTPUT AND INFORMATION SET FOR THE LOOP
   ********************************************************************/
  LoopInfo  loop_info("nonlinear");
  loop_info.print_time_every_step = true;
  loop_info.verbosity_threshold   = 1;            // full verbosity
//  loop_info.print(0, vof.tnse3d_.get_full_residual());

  stopwatch.print_total_time("setting up spaces, matrices, linear assemble");
  stopwatch.reset();
  stopwatch.start();

  Chrono nse_nl_stopwatch;
  Chrono nse_timeit_stopwatch;



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
      if (my_rank ==0)
        Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=0;
      vof.tnse3d_.assemble_rhs();
      if (vof.tnse_variable_fluid_ == true)
      {
        vof.tnse3d_.assemble_massmatrix_withfields();
//        if( TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION == 1 )
//            vof.tnse2d_.apply_slip_penetration_bc(true,true);
      }

      vof.tnse3d_.assemble_nonlinear_term();
      vof.tnse3d_.assemble_system();


    /********************************************************************
     * NON LINEAR LOOP
     ********************************************************************/
     nse_timeit_stopwatch.restart_and_print("preparation of NSE iterations");
    for(unsigned int k = 0;; k++)
    {
      vof.tnse3d_.compute_residuals();

      if (my_rank==0) // some outputs
      {
        Output::print<1>("\nNONLINEAR ITERATION :", setw(3), k);
        Output::print<1>("Residuals :", vof.tnse3d_.get_residuals());
      }

      if(vof.tnse3d_.stop_it(k))
        break;

      vof.tnse3d_.solve();

      if(vof.tnse3d_.imex_scheme(1))
        continue; // this interrupts the NL-Loop

      if (vof.tnse_variable_fluid_ == true)
      {
        if( TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION == 1 )
        {ErrThrow("Slip BC Not Implemented yet!");
//          vof.tnse3d_.apply_slip_penetration_bc(false,false);
        }
      }

      vof.tnse3d_.assemble_nonlinear_term();
      vof.tnse3d_.assemble_system();
      nse_nl_stopwatch.restart_and_print("solving and reassembling NL iter. "
                                          + std::to_string(k));
    } // end for k, non linear loop

    nse_timeit_stopwatch.restart_and_print("total NSE time iter. "
                                  +std::to_string(TDatabase::TimeDB->CURRENTTIME));


    /********************************************************************
     * SOLVING CD3D WITH NSE3D SOLUTION
     ********************************************************************/
    if (!tcd_db["algebraic_flux_correction"].is("none"))
      TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE=1;
    if (vof.solve_convection_ == true )
    {
      if (my_rank==0)
        Output::print<1>("<<<<<<<<<<<<<<<<<< NOW SOLVING CONVECTION  >>>>>>>>>>>>>");
      if (vof.nse2cd_coupling_ == true)
      {
//        tcd2d.assemble_rhs_vector(&tnse2d.get_velocity()); // once per time step
//        tcd2d.assemble_stiffness_matrix_alone_with_convection(&tnse2d.get_velocity());
//        tcd2d.scale_stiffness_matrix();
        vof.phaseconvection3d_.assemble_with_convection(&vof.tnse3d_.get_velocity());
      }
      else  // if we solve TCD3D standard, without any coupling
      { vof.phaseconvection3d_.assemble();}

      /* this is a quick workaround needed because of the
       * global parameter INTERNAL_PROJECT_PRESSURE, which disturbs
       * normal use of Mumps in TCD3D
       */
      int ipp_was_1 = 0;
      if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE==1)
      {
        ipp_was_1 = 1;
        TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;
      }
      vof.phaseconvection3d_.solve();
      if (ipp_was_1)
      {
        TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;
      }
      vof.phaseconvection3d_.descale_stiffness(); //needed once per time loop

      if (my_rank==0)
        Output::print<1>("<<<<<<<<<<<<<<<<<< END SOLVING CONVECTION >>>>>>>>>>>>>>");


      /********************************************************************
       * UPDATING VELOCITY VECTOR WITH CD3D SOLUTION
       ********************************************************************/
      if (vof.cd2nse_coupling_ == true )
      {
        vof.update_field_vectors();
//        vof.output_vectors("vector_phi_updated","vector_rho_updated","vector_mu_updated");
      }
    }

    stopwatch.restart_and_print("total whole iter. " +
                                std::to_string(TDatabase::TimeDB->CURRENTTIME));

    vof.tnse3d_.output(step,image);

    if(vof.solve_convection_ == true)
    {
      if((step-1) % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
        vof.phaseconvection3d_.output(step,image,&vof.tnse3d_.get_velocity());
    }
//        vof.tnse3d_.get_solution().write("solution_velocity");
    }
  } // end for step, time loop

  stopwatch.print_total_time("total solving duration: ");
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();

#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;
}
// end main
