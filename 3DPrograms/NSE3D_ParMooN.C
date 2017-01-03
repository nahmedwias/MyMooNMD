/**
 * @brief Main program for solving a 3D stationary Navier Stokes equation using ParMooN.
 * @author Sashikumaar Ganesan, Clemens Bartsch
 *
 * Implementation started on 2015/01/23. Rework since 2015/12/7.
 *
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <NSE3D.h>
#include <MeshPartition.h>
#include <Chrono.h>
#include <LoopInfo.h>

//project specific
#include <CoiledPipe.h>

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

//project specific declaration
struct derived_properties;
int main(int argc, char* argv[])
{
#ifdef _MPI
  //Construct and initialise the default MPI communicator.
  MPI_Init(&argc, &argv);
#endif
  {
#ifdef _MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(my_rank==0)
    {
      Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
      Output::info("NSE3D", "MPI, using ", size, " processes");
    }
#else
    int my_rank = 0;
    Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
    Output::info("NSE3D", "SEQUENTIAL (or OMP...)");
#endif

    //start a stopwatch which measures time spent in program parts
    Chrono timer;

    // Construct the ParMooN Databases.
    TDatabase Database;
    ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
    std::ifstream fs(argv[1]);
    parmoon_db.read(fs);
    fs.close();

#ifdef _MPI
    TDatabase::ParamDB->Comm = comm;
#endif

    TFEDatabase3D feDatabase;

  //open OUTFILE, this is where all output is written to (additionally to console)
    //open OUTFILE, this is where all output is written to (addionally to console)
    if(my_rank==0)
    {
      Output::set_outfile(parmoon_db["outfile"]);
    }
    Output::setVerbosity(parmoon_db["verbosity"]);

    if(my_rank==0) //Only one process should do that.
    {
      parmoon_db.write(Output::get_outfile());
      Database.WriteParamDB(argv[0]);
    }

    // Do the old parameter check of the Database.
    Database.CheckParameterConsistencyNSE();
    
      // project specific: prepare the coiled geometry
  size_t n_twists                 = parmoon_db["twisted_pipe_n_twists"];
  size_t n_segments_per_twist     = parmoon_db["twisted_pipe_n_segments_per_twist"];
  double l_inflow                 = parmoon_db["twisted_pipe_l_inflow"];
  size_t n_segments_inflow        = parmoon_db["twisted_pipe_n_segments_inflow"];
  double l_outflow                = parmoon_db["twisted_pipe_l_outflow"];
  size_t n_segments_outflow       = parmoon_db["twisted_pipe_n_segments_outflow"];
  double tube_radius              = parmoon_db["twisted_pipe_tube_radius"];
  double twist_radius             = parmoon_db["twisted_pipe_twist_radius"];
  double space_between_twists     = parmoon_db["twisted_pipe_space_between_twists"];

  CoiledPipe::set_up_geoconsts(
      n_twists,
      n_segments_per_twist,
      l_inflow,
      n_segments_inflow,
      l_outflow,
      n_segments_outflow,
      tube_radius,
      twist_radius,
      space_between_twists
  );
  double drift_x = 0;
  double drift_y = 0;
  double drift_z =CoiledPipe::GeoConsts::l_tube;
  
  // Choose and construct example - in this project, this has to be done before
  // initiing the Domain, as the example itself inluences the geometry by setting
  // a global parameter (which is awful).
  Example_NSE3D example(parmoon_db["example"], parmoon_db);

  // Construct domain, thereby read in controls from the input file.
  TDomain domain(argv[1], parmoon_db,
                 drift_x, drift_y, drift_z, CoiledPipe::GeoConsts::segment_marks);




    // Intial refinement and grabbing of grids for multigrid.
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

    //print information on the mesh partition on the finest grid
    domain.print_info(std::string("coiled tube"));

//Now this is Felix' code - get the correct plane for grabbing values.
    size_t planeSegment = parmoon_db["twisted_pipe_planeSegment"];

  // variable zPlane
  double zPlane = 0.;

  if (planeSegment < n_segments_inflow)
  {//inflow part
	  zPlane = (double) planeSegment / n_segments_inflow * l_inflow;
  }
  else if (planeSegment < n_segments_inflow + n_twists * n_segments_per_twist )
  {//twisted part
	  size_t i_twist = planeSegment - n_segments_inflow;
	  zPlane = l_inflow + (double) i_twist / (n_twists * n_segments_per_twist)
			   * (drift_z - l_inflow - l_outflow);
  }
  else
  {//outflow part
	  size_t i_out = planeSegment - n_segments_inflow - n_twists * n_segments_per_twist;
	  zPlane = l_inflow + (drift_z - l_inflow - l_outflow)
			  + (double) i_out / n_segments_outflow * l_outflow;
  }

  // plane with position vector and normal vector
  TVertex*  posPlane = new TVertex(0., 0., zPlane);
  TVertex* normPlane = new TVertex(0., 0., 1.);

  // copy the straight pipe collection
  TCollection* straightColl = domain.GetCollection(It_Finest, 0);
  std::vector<TBaseCell*> cells = CoiledPipe::getCellsAtPlane(straightColl, posPlane, normPlane);

  // discretization of the plane
  // it will be a square of height = width = tube_radius
  int discN = 100;
  double h = 2*tube_radius/(discN-1);

  std::vector<CoiledPipe::VertexValues> vertexList;

  // some points in the plane
  for (int i = 0; i < discN; ++i)
  {
	double xCoord = -tube_radius + i*h;

	for (int j = 0; j < discN; ++j)
	{
	  double yCoord = -tube_radius + j*h;

	  // if the point is in the tube
	  if ( xCoord*xCoord + yCoord*yCoord < tube_radius*tube_radius)
	  {
		// already swapped x and z, because else the coiling of these
		// vertices is going wrong
	    CoiledPipe::VertexValues p(j+i*discN, xCoord,
	    						   yCoord, zPlane);
#ifdef _MPI
	    double x, y, z;

	    p.GetCoords(x, y, z);

	    for (TBaseCell* cell : cells)
	    {
	    	if ( cell->PointInCell(x, y, z) )
	    	{
#endif
	    		p.origZ = p.origX;
	    		p.origX = zPlane - CoiledPipe::GeoConsts::l_inflow;

	    		vertexList.push_back(p);
#ifdef _MPI
	    		break;
	    	} // end if
	    }// end for cells
#endif
	  }
	}
  }


    timer.restart_and_print("setup(domain, example, database)");
    // Construct an object of the NSE3D-problem type.
#ifdef _MPI
    NSE3D nse3d(gridCollections, parmoon_db, example, maxSubDomainPerDof);
#else
    NSE3D nse3d(gridCollections, parmoon_db, example);
#endif
    timer.restart_and_print("constructing NSE3D object");
    
    // assemble all matrices and right hand side
    nse3d.assemble_linear_terms();
    nse3d.stop_it(0);  // check initial residuals

    LoopInfo loop_info("nonlinear");
    loop_info.print_time_every_step = true;
    loop_info.verbosity_threshold = 1; // full verbosity
    if(my_rank==0)
      loop_info.print(0, nse3d.get_full_residual());

    timer.restart_and_print("assembling linear terms");

    //Timer which measures time spent in solve() method solely.
    Chrono timer_sol;
    timer_sol.stop();

    //======================================================================
    for(unsigned int k=1;; k++)
    {
      if(my_rank == 0)
        Output::print(); // new line for a new nonlinear iteration
      // solve the system
      timer_sol.start();
      nse3d.solve();
      timer_sol.stop();
    if (cells.empty())
    {
   	 Output::print<2>("No cells found at zPlane = ", zPlane);
    }
    else
    {
      Output::print<2>(cells.size(), " cells found");

      std::string filename = std::string("testVelocity.txt");

      CoiledPipe::writeVelocityOfCells(cells, vertexList,
   	  	 	 	 	 	 	 	    nse3d.get_velocity(),
									nse3d.get_pressure(),
									filename);
    }

      //no nonlinear iteration for Stokes problem
      if(parmoon_db["problem_type"].is(3))
        break;

      nse3d.assemble_non_linear_term();

      // checking residuals
      if(nse3d.stop_it(k))
      {
        loop_info.finish(k, nse3d.get_full_residual());
        break;
      }
      else
        loop_info.print(k, nse3d.get_full_residual());

    } // end for k
    timer.restart_and_print("nonlinear loop");
    
    timer_sol.print_total_time("solver only");

    nse3d.output();
    timer.restart_and_print("output");
    timer.print_total_time("NSE3D_ParMooN program");

    if(my_rank==0)
      Output::print("<<<<< ParMooN Finished: NSE3D Main Program >>>>>");

    if(my_rank == 0)
      Output::close_file();
  }

#ifdef _MPI
  MPI_Finalize();
#endif
  
  return 0;
}

