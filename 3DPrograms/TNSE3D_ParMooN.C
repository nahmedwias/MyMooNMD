/**
 * @brief Main program for solving a 3D time-dependent Navier Stokes equation using ParMooN.
 * @author Najib Alia
 *
 * Implementation started on 2016/04/15.
 *
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <Time_NSE3D.h>
#include <MeshPartition.h>
#include <Chrono.h>
#include <TetGenMeshLoader.h>
#include <Output3D.h>
#include <LoopInfo.h>

#include <sys/stat.h>

#include <TimeDiscRout.h>

using namespace std;

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

// CB EXAMPLE
void transform_to_crystallizer_geometry(TCollection *coll, double outflow_stretch);
// END EXAMPLE

// main program
// =======================================================================
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
    Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
    Output::info("Time_NSE3D", "MPI, using ", size, " processes");
  }
#else
  int my_rank = 0;
  Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
  Output::info("Time_NSE3D", "SEQUENTIAL (or OMP...)");
#endif
  double t_start = GetTime();
  //start a stopwatch which measures time spent in program parts
  Chrono timer;

  // Construct the ParMooN Databases.
  TDatabase Database;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(ParameterDatabase::default_tetgen_database(), true);

  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
#ifdef _MPI
  TDatabase::ParamDB->Comm = comm;
#endif
  TFEDatabase3D feDatabase;

  //open OUTFILE, this is where all output is written to (additionally to console)
  if(my_rank==0)
  {
    Output::set_outfile(parmoon_db["outfile"]);
  }
  Output::setVerbosity(parmoon_db["verbosity"]);

  // Choose and construct example.
  Example_TimeNSE3D example(parmoon_db);

  // Do the parameter check of the Database.
  check_parameters_consistency_NSE(parmoon_db);
  // =====================================================================
  // set the database values and generate mesh
  // =====================================================================

  //CB EXAMPLE
  //This code is specific to the WiedmeyerBatchCrystallizer Example.
  if(parmoon_db["sandwich_grid"])
  {//prepare parameters for an example specific sandwich grid
    ParameterDatabase& sw_db = parmoon_db.get_nested_database("Sandwich Grid Database");
    sw_db.add("lambda", {0.0,1.0}, "Check default_sandwich_grid_parameters for description.");

    int n_layers_inflow = sw_db["n_layers_inflow"];
    int n_layers_cone = sw_db["n_layers_cone"];
    int n_layers_outflow = sw_db["n_layers_outflow"];
    std::vector<double> lambda(n_layers_inflow + n_layers_cone + n_layers_outflow + 1);
    //TODO fill lambda somehow
    for(int i=0; i<(int)lambda.size(); ++i)
    {//fill lambda
      if(i < n_layers_inflow)
        lambda[i] = (1.0/10) * i * (1.0/n_layers_inflow);
      else if(i < n_layers_inflow + n_layers_cone)
        lambda[i] = 1.0/10 +(6.0/10) * (i - n_layers_inflow) * (1.0/n_layers_cone);
      else
        lambda[i] = 7.0/10 +(3.0/10) * (i - n_layers_inflow - n_layers_cone) * (1.0/n_layers_outflow);
    }
    //put the new lambda into the database
    sw_db["lambda"] = lambda;
  }
  //END EXAMPLE

  // Construct domain, thereby read in controls from the input file.
  Output::decreaseVerbosity(1);
  TDomain domain(parmoon_db, argv[1]);
  Output::increaseVerbosity(1);


  if(my_rank==0) //Only one process should do that.
  {
    parmoon_db.write(Output::get_outfile());
    Database.WriteParamDB(argv[0]);
    Database.WriteTimeDB();
  }
  // Initial refinement and grid collection
#ifdef _MPI
  int maxSubDomainPerDof = 0;
#endif
  std::list<TCollection* > gridCollections
     = domain.refine_and_get_hierarchy_of_collections(
       parmoon_db
#ifdef _MPI
       , maxSubDomainPerDof
#endif       
    );
  //CB EXAMPLE
  // This code is specific to the WiedmeyerBatchCrystallizer Example.
  // This will be done for only the finest grid -
  // the other grids share its vertices!
  if(parmoon_db["sandwich_grid"])
  {
    double outflow_stretch = 7.5;
    if(parmoon_db.contains("outflow_stretch"))
      outflow_stretch = parmoon_db["outflow_stretch"];
    transform_to_crystallizer_geometry(gridCollections.front(), outflow_stretch);
  }
  //END EXAMPLE

  TCollection* coll = gridCollections.front();
  TOutput3D output(0,0,0,0,std::addressof(domain),coll);
  
  //print information on the mesh partition on the finest grid
  domain.print_info("TNSE3D domain");

//  //CB DEBUG
//  // In case you want to look at the result.
//  int i=0;
//  for(auto grid : gridCollections)
//  {
//    TOutput3D output(0,0,0,0,std::addressof(domain),grid);
//    std::string my_str = "sandwich_mesh.";
//    my_str += std::to_string(i++);
//    my_str += ".vtk";
//    output.WriteVtk(my_str.c_str());
//  }
//  //END DEBUG

  // set some parameters for time stepping
  SetTimeDiscParameters(0);
  // Construct an object of the Time_NSE3D-problem type.
#ifdef _MPI
  Time_NSE3D tnse3d(gridCollections, parmoon_db, example, maxSubDomainPerDof);
#else
  Time_NSE3D tnse3d(gridCollections, parmoon_db, example);
#endif
  
  tnse3d.assemble_initial_time();
  
  double end_time = TDatabase::TimeDB->ENDTIME;
  tnse3d.current_step_ = 0;

  int n_substeps = GetN_SubSteps();

  int image = 0;
  
  LoopInfo loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1; // full verbosity
  int linear_iterations = 0; 

  timer.restart_and_print("setting up spaces, matrices and initial assembling");
  TDatabase::TimeDB->CURRENTTIME = 0.0;  
  //======================================================================
  // time iteration
  //======================================================================
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    // time measuring during every time iteration
    Chrono timer_timeit;

    tnse3d.current_step_++;

    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    for(int j = 0; j < n_substeps; ++j) // loop over substeps in one time iteration
    {
      // setting the time discretization parameters
      SetTimeDiscParameters(1);
     if( tnse3d.current_step_ == 1 && my_rank==0) // a few output, not very necessary
     {
       Output::print<1>("Theta1: ", TDatabase::TimeDB->THETA1);
       Output::print<1>("Theta2: ", TDatabase::TimeDB->THETA2);
       Output::print<1>("Theta3: ", TDatabase::TimeDB->THETA3);
       Output::print<1>("Theta4: ", TDatabase::TimeDB->THETA4);
     }
      // tau may change depending on the time discretization (adaptive time)
      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;

      if (my_rank==0)
        Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      // prepare the right hand side vector - needed only once per time step
      tnse3d.assemble_rhs();

      // assemble the nonlinear matrices
      tnse3d.assemble_nonlinear_term();

      // prepare the matrices for defect computations and solvers
      tnse3d.assemble_system();
      
      timer_timeit.restart_and_print("preparation of nonlinear iteration");

      // nonlinear iteration
     LoopInfo loop_info("nonlinear");
     loop_info.print_time_every_step = true;
     loop_info.verbosity_threshold = 1; // full verbosity
     for(unsigned int k=0; ; k++)
      {
        tnse3d.compute_residuals();

        if (my_rank==0) // some outputs
        {
          Output::print<1>("\nNONLINEAR ITERATION :", setw(3), k);
          Output::print<1>("Residuals :", tnse3d.get_residuals());
        }

        // checking residuals and stop conditions
        if(tnse3d.stop_it(k))
        {
           loop_info.finish(k, tnse3d.get_full_residual());
	   linear_iterations+=k;
	   /// @todo provide all parts of the residual 
	   /// @todo loop_info restricted to the solver only
           loop_info_time.print(linear_iterations, tnse3d.get_full_residual());
	   break;
         }
         else
           loop_info.print(k, tnse3d.get_full_residual());
 
        tnse3d.solve();

        if(tnse3d.imex_scheme(1))
          continue;

        tnse3d.assemble_nonlinear_term();

        tnse3d.assemble_system();

        timer_timeit.restart_and_print("solving and reassembling in the "
                                        "nonlinear iteration " + 
                                        std::to_string(k));
      }  // end of nonlinear loop

      timer_timeit.restart_and_print(
        "solving the time iteration " +
        std::to_string(TDatabase::TimeDB->CURRENTTIME));

      tnse3d.output(tnse3d.current_step_,image);
      timer_timeit.print_total_time(
        "time step " + std::to_string(TDatabase::TimeDB->CURRENTTIME));
    } // end of subtime loop
  } // end of time loop
  loop_info_time.finish(linear_iterations, tnse3d.get_full_residual());

  timer.print_total_time("whole solving procedure ");

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


//CB EXAMPLE
//This code is specific to the WiedmeyerBatchCrystallizer Example.
void compute_position_in_crystallizer_geometry(
    double x, double y, double z,
    double& x_trans, double& y_trans, double& z_trans,
    double outflow_stretch
    )
{// We assume that the input geometry is a cylinder with radius 1 (cm)
 // and height 50 (cm), z being the height direction.
 // The cylinder is further assumed to 'fit' the inflow of
 // the crystallizer geometry, i.e., the conical
 // part and the outflow part are gained by a stretching.
 double tol = 1e-6;
 double inflow_end = 5;
 double cone_end = 35;
 double outflow_end = 50;

 if(z < inflow_end + tol)//inflow piece,keep it unchanged.
 {
   x_trans = x;
   y_trans = y;
   z_trans = z;
 }
 else if (z < cone_end + tol)//conical piece, stretch it linearly
 {
   double lincomb = (1 - (z - inflow_end)/(cone_end - inflow_end)) * 1.0 +
                    (z - inflow_end)/(cone_end - inflow_end) * outflow_stretch;
   x_trans = x * lincomb;
   y_trans = y * lincomb;
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
void transform_to_crystallizer_geometry(TCollection *coll, double outflow_stretch)
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
        compute_position_in_crystallizer_geometry(x, y, z, x_transf, y_transf, z_transf, outflow_stretch);
        vertex->SetCoords(x_transf, y_transf, z_transf);

        //mark this vertex as treated
        vertex->SetClipBoard(1);
      }
    }
  }
}
//END EXAMPLE
