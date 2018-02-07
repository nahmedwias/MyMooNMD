/**
 * Main program for solving a 3D coupled problem of a time-dependent flow,
 * a transported concentration/temperature and a transported particle population.
 * The file is, due to set-up and geometry, adapted to a special example: the
 * 3d simulation of the Kalialaun batch crystallizer by V. Wiedmeyer et al.
 *
 * The program should be kept compileable in both SEQ and MPI, because it is
 * planned to compute the flow in parallel sooner or later, while keeping the
 * interface to Brush sequential, i.e., on root process for the moment.
 *
 * @author Clemens Bartsch
 *
 * Implementation started on 2018/1/15.
 *
 */
#include <BrushWrapper.h>
#include <Domain.h>
#include <FEDatabase3D.h>
#include <LoopInfo.h>
#include <Output3D.h>
#include <Time_NSE3D.h>
#include <TimeDiscRout.h>

#include <math.h>

//CB DEBUG
#include <GridTransferTool.h>
void DirichletBoundaryConditions_TEST(int BdComp, double t, BoundCond &cond)
{
      cond = DIRICHLET;
}
//END DEBUG

// CB EXAMPLE
void transform_to_crystallizer_geometry(TCollection *coll, double outflow_stretch, bool cut_off_entry);
// END EXAMPLE

// main program
// =======================================================================
int main(int argc, char* argv[])
{
  int my_rank=0;
#ifdef _MPI
  //Construct and initialise the default MPI communicator.
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  // Hold mpi rank and size ready, check whether the current processor
  // is responsible for output (usually root, 0).
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  {
    // Construct the ParMooN Databases.
    // FIXME As long as it is only possible to read in one level of nested databases,
    // we are forced to treat the flow database as the general database, and
    // store particle database as one of its nested databases.
    std::ifstream fs(argv[1]);
    ParameterDatabase flow_database("Flow Database");
    flow_database.read(fs);
    fs.close();
    ParameterDatabase particle_database = flow_database.get_nested_database("Particle Database");
    //enrich the particle database by some general parameters from the flow database
    particle_database.add(flow_database["output_vtk_directory"]);

    //This has to do with the old database
    TDatabase Database;
    TFEDatabase3D feDatabase;
#ifdef _MPI
    TDatabase::ParamDB->Comm = comm;
#endif

    //open OUTFILE, this is where all output is written to (additionally to console)
    if(my_rank==0)
      Output::set_outfile(flow_database["outfile"]);
    Output::setVerbosity(flow_database["verbosity"]);

    // Here starts the actual program
    Output::root_info("<<<<< Running ParMooN","Population Balance System Main Program >>>>>");
#ifdef _MPI
    Output::root_info("RUN with","MPI, using ", size, " processes");
#else
    Output::info("RUN with", "SEQUENTIAL");
#endif
    //start a stopwatch which measures time spent in program parts
    Chrono timer;

    if(my_rank==0) //Only one process should do that.
    {
      flow_database.write(Output::get_outfile());
      particle_database.write(Output::get_outfile());
      //Database.WriteParamDB(argv[0]);
      //Database.WriteTimeDB();
    }

    //CB EXAMPLE
    //This code is specific to the WiedmeyerBatchCrystallizer Example.
    if(flow_database["sandwich_grid"])
    {//prepare parameters for an example specific sandwich grid
      ParameterDatabase& sw_db = flow_database.get_nested_database("Sandwich Grid Database");
      sw_db.add("lambda", {0.0,1.0}, "Check default_sandwich_grid_parameters for description.");

      int n_layers_inflow = sw_db["n_layers_inflow"];
      int n_layers_cone = sw_db["n_layers_cone"];
      int n_layers_outflow = sw_db["n_layers_outflow"];
      double inflow_part = 1.0/10;
      double cone_part = 6.0/10;
      double outflow_part = 3.0/10;
      if(flow_database["cut_off_entry"])
      {
        n_layers_inflow = 0;
        inflow_part = 0;
        cone_part = 2.0/3;
        outflow_part = 1.0/3;
      }

      std::vector<double> lambda(n_layers_inflow + n_layers_cone + n_layers_outflow + 1);
      for(int i=0; i<(int)lambda.size(); ++i)
      {//fill lambda
        if(i < n_layers_inflow)
          lambda[i] = inflow_part * i * (1.0/n_layers_inflow);
        else if(i < n_layers_inflow + n_layers_cone)
          lambda[i] = inflow_part + cone_part * (i - n_layers_inflow) * (1.0/n_layers_cone);
        else
          lambda[i] = inflow_part + cone_part
          + outflow_part * (i - n_layers_inflow - n_layers_cone) * (1.0/n_layers_outflow);
      }
      //put the new lambda into the database
      sw_db["lambda"] = lambda;

      if(flow_database["cut_off_entry"])
      {
        sw_db["drift_z"] = 45;
        Output::root_info("SANDWICH GRID", "drift_z was set to 45 due to 'cut_off_entry'");
      }
    }
    //END EXAMPLE

    // Construct domain, thereby read in the old parameter controls from the input file.
    Output::decreaseVerbosity(1);
    TDomain domain(flow_database, argv[1]);
    Output::increaseVerbosity(1);

    // Initial refinement and grid collection
#ifdef _MPI
    int maxSubDomainPerDof = 0;
    std::list<TCollection* > flow_grids =
        domain.refine_and_get_hierarchy_of_collections(flow_database, maxSubDomainPerDof );
#else
    std::list<TCollection* > flow_grids =
            domain.refine_and_get_hierarchy_of_collections(flow_database);
#endif

    std::string brush_grid_control = particle_database["brush_grid_control"];
    TCollection* brush_grid;
    if(brush_grid_control == "coarser")
    {
      Output::root_info("GRID", "Picking one level coarser grid for Brush.");
      //The iterator It_OCAF picks cells which lie one level under the finest one -
      // exactly what I want here!
      // TODO this will not work in MPI case - but anyway, in MPI case the brush grid will be only on process 0
      brush_grid = domain.GetCollection(It_OCAF, 0); //the number 0 is unused in OCAF iterator
    }
    else if(brush_grid_control == "same")
    {
     Output::info("GRID", "Picking same grid for Brush.");
     brush_grid = domain.GetCollection(It_Finest, 0);
    }
    else if(brush_grid_control == "finer")
    {
     Output::info("GRID", "Finer grid for Brush required, performing additional refinement step.");
     domain.RegRefineAll();
     brush_grid = domain.GetCollection(It_Finest, 0);
    }
    else
      ErrThrow(brush_grid_control, " is an unknown control parameter for the Brush grid.");

    //CB EXAMPLE
    // This code is specific to the WiedmeyerBatchCrystallizer Example.
    // This will be done for only the finest grid -
    // the other grids share its vertices!
    if(flow_database["sandwich_grid"])
    {
      TCollection* finest_collection = domain.GetCollection(It_Finest, 0);
      double outflow_stretch = 7.5;
      if(flow_database.contains("outflow_stretch"))
        outflow_stretch = flow_database["outflow_stretch"];
      transform_to_crystallizer_geometry(finest_collection,
                                         outflow_stretch,
                                         flow_database["cut_off_entry"]);
      delete finest_collection;
    }
    //END EXAMPLE

    //print information on the mesh partition on the finest grid
    domain.print_info("PBS domain");

    // set some parameters for time stepping
    SetTimeDiscParameters(0);

    // Choose and construct the flow example.
    Example_TimeNSE3D flow_example(flow_database);
    check_parameters_consistency_NSE(flow_database);

    // Construct an object of the Time_NSE3D-problem type.
#ifdef _MPI
    Time_NSE3D flow_object(flow_grids, flow_database, flow_example, maxSubDomainPerDof);
#else
    Time_NSE3D flow_object(flow_grids, flow_database, flow_example);
#endif

    // the particles object which wraps up Brush
    BrushWrapper part_object(brush_grid, domain.GetCollection(It_Finest, 0), particle_database);

    // PARTS: SET UP INITIAL STATES /////////////////////////////////////////////
    Output::info("PROGRAM PART", "Setting up initial states.");

    //first assembling of the flow object
    if(!flow_database["force_stationary"])
      flow_object.assemble_initial_time();

    //set fluid phase in the particle object
    const TFEFunction3D& ux = *flow_object.get_velocity_component(0);
    const TFEFunction3D& uy = *flow_object.get_velocity_component(1);
    const TFEFunction3D& uz = *flow_object.get_velocity_component(2);
    const TFEFunction3D& pressure = flow_object.get_pressure();
    TFEFunction3D* dummy_T    = &flow_object.get_pressure();    //TODO Currently the pressure will be our replacement for temperature and concentration
    TFEFunction3D* dummy_conc = &flow_object.get_pressure();    //TODO Currently the pressure will be our replacement for temperature and concentration
    std::vector<TFEFunction3D*> fcts = {dummy_T, dummy_conc};
    part_object.reset_fluid_phase(ux, uy, uz, pressure, fcts);

    timer.restart_and_print("setting up spaces, matrices and initial assembling");

    // PART: SOLVE THE SYSTEM IN A TIME LOOP ///////////////////////////////////
    Output::info("PROGRAM PART", "Solving the coupled system.");

    double end_time = TDatabase::TimeDB->ENDTIME;
    flow_object.current_step_ = 0;

    // Those two variables are used for the .VTK output filename control
    int flow_image = 0;
    int part_image = 0;

    if(flow_database["force_stationary"]) //put it out just once
      flow_object.output(flow_object.current_step_,flow_image);

    LoopInfo loop_info_time("time loop", true, true, 1);
    int linear_iterations = 0;
    TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;

    int output_steps_parts = particle_database["output_all_k_steps"];
    int output_steps_flow = flow_database["output_all_k_steps"];
    int step = 0;

    //======================================================================
    // time iteration
    //======================================================================
    while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
    {
      // time measuring during every time iteration
      Chrono timer_timeit;

      step++;

      //////////////////// THE PARTICLES //////////////////////
      Output::print(" BRUSH - UPDATE AND SOLVE");
      part_object.reset_fluid_phase(ux,uy,uz, pressure, fcts);
      Output::print(" resetting complete");
      part_object.solve(TDatabase::TimeDB->CURRENTTIME, TDatabase::TimeDB->CURRENTTIME + TDatabase::TimeDB->CURRENTTIMESTEPLENGTH);
      Output::print(" solving complete");

      if(!flow_database["force_stationary"])
        flow_object.current_step_++;

      // setting the time discretization parameters
      SetTimeDiscParameters(1);

      // tau may change depending on the time discretization (adaptive time)
      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;

      Output::root_info("CURRENT TIME", TDatabase::TimeDB->CURRENTTIME);

      if(!flow_database["force_stationary"])
      {
        ////////////////////   THE FLOW   //////////////////////
        Output::print(" PARMOON - UPDATE AND SOLVE");
        // prepare the right hand side vector - needed only once per time step
        flow_object.assemble_rhs();

        // assemble the nonlinear matrices
        flow_object.assemble_nonlinear_term();

        // prepare the matrices for defect computations and solvers
        flow_object.assemble_system();

        //LoopInfo for the nonlinear loop, printing full verbose info every step.
        LoopInfo loop_info("nonlinear", true, true, 1);

        // nonlinear iteration
        for(unsigned int k=0; ; k++)
        {
          flow_object.compute_residuals();

          Output::root_info<1>("NONLINEAR ITERATION ", setw(3), k);
          Output::root_info<1>("Residuals ", flow_object.get_residuals());

          // checking residuals and stop conditions
          if(flow_object.stop_it(k))
          {
            loop_info.finish(k, flow_object.get_full_residual());
            linear_iterations+=k;
            loop_info_time.print(linear_iterations, flow_object.get_full_residual());
            break;
          }
          else
            loop_info.print(k, flow_object.get_full_residual());

          flow_object.solve();

          if(flow_object.imex_scheme(1))
            continue;

          flow_object.assemble_nonlinear_term();

          flow_object.assemble_system();

          timer_timeit.restart_and_print("solving and reassembling in the "
              "nonlinear iteration " +
              std::to_string(k));
        }  // end of nonlinear loop
      }
      timer_timeit.restart_and_print(
          "solving the time iteration " +
          std::to_string(TDatabase::TimeDB->CURRENTTIME));

      if(step %  output_steps_parts == 0)
        part_object.output(part_image,TDatabase::TimeDB->CURRENTTIME);
      if(step % output_steps_flow == 0 && !flow_database["force_stationary"])
        flow_object.output(flow_object.current_step_,flow_image);

      timer_timeit.print_total_time(
          "time step " + std::to_string(TDatabase::TimeDB->CURRENTTIME));
    } // end of time loop
    loop_info_time.finish(linear_iterations, flow_object.get_full_residual());

    timer.print_total_time("whole solving procedure ");

    Output::close_file();

  }
#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;
}


//CB EXAMPLE
//This code is specific to the WiedmeyerBatchCrystallizer Example.
void compute_position_in_crystallizer_geometry(
    double x, double y, double z,
    double& x_trans, double& y_trans, double& z_trans,
    double outflow_stretch, bool cut_off_entry
)
{// We assume that the input geometry is a cylinder with radius 1 (cm)
  // and height 50 cm OR 45 cm (if cut_off_entry) , z being the height direction.
  // The cylinder is further assumed to 'fit' the inflow of
  // the crystallizer geometry, i.e., the conical
  // part and the outflow part are gained by a stretching.

  double tol = 1e-6;
  double inflow_end = 5;
  double cone_end = 35;
  double outflow_end = 50;
  if(cut_off_entry)
  {
    inflow_end = 0;
    cone_end = 30;
    outflow_end = 45;
  }

  if(z < inflow_end + tol)//inflow piece,keep it unchanged.
  {
    x_trans = x;
    y_trans = y;
    z_trans = z;
  }
  else if (z < cone_end)//conical piece, stretch it linearly
  {
    double lincomb = (1 - (z - inflow_end)/(cone_end - inflow_end)) * 1.0 +
        (z - inflow_end)/(cone_end - inflow_end) * outflow_stretch;
    x_trans = x * lincomb;
    y_trans = y * lincomb;
    //with those two lines enforce convexity in the conical section
    //x_trans *= 1+std::abs(z - 30)/30;
    //y_trans *= 1+std::abs(z - 30)/30;
    z_trans = z;
  }
  else if (z < outflow_end + tol)//outflow piece, stretch it constantly
  {

    x_trans = x * outflow_stretch;
    y_trans = y * outflow_stretch;
    z_trans = z;
  }
  else
  {
    ErrThrow("z is out of range!");
  }

  // Finally, divide everything by 100 - input geometry is expected
  // in [cm], but internally we use [m] currently
  //TODO this could be changed to cm eventually!
  x_trans/=100;
  y_trans/=100;
  z_trans/=100;

}

/**
 * This method transforms a cylindrical input grid to
 * the 'batch crystallizer geometry'.
 * No checks whatsoever are performed, the method it is
 * extremely specific. You must know what you are doing.
 */
void transform_to_crystallizer_geometry(
    TCollection *coll, double outflow_stretch, bool cut_off_entry)
{
  int N_Cells = coll->GetN_Cells();

  // initialise ClipBoard
  for(int i=0 ; i<N_Cells ; i++)
  {
    TBaseCell* cell = coll->GetCell(i);
    int n_verts = cell->GetN_Vertices();

    for (int j=0 ; j<n_verts ; j++)
    {
      TVertex* vertex = cell->GetVertex(j);
      vertex->SetClipBoard(0);
    }
  }

  // visit each vertex, reset its coordinates to its new position
  for(int i=0 ; i < N_Cells ; i++)
  {
    TBaseCell* cell = coll->GetCell(i);
    int n_verts = cell->GetN_Vertices();

    for (int j=0 ; j < n_verts ; j++)
    {
      TVertex* vertex = cell->GetVertex(j);
      if (!vertex->GetClipBoard())
      {
        double x, y, z;
        double x_transf, y_transf, z_transf;

        vertex->GetCoords(x, y, z);
        compute_position_in_crystallizer_geometry(x, y, z, x_transf, y_transf, z_transf,
                                                  outflow_stretch, cut_off_entry);
        vertex->SetCoords(x_transf, y_transf, z_transf);

        //mark this vertex as treated
        vertex->SetClipBoard(1);
      }
    }
  }
}
//END EXAMPLE
