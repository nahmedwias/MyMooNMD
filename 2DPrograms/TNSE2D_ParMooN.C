#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Example_TimeNSE2D.h>
#include <Time_NSE2D.h>
#include <TimeDiscretizations.h>
#include <TimeDiscRout.h>
#include <LoopInfo.h>
#include <MeanVelocity.h>

using namespace std;


void get_output_points_inside_cells(TCollection* coll, 
				    std::vector<double> &xy_coords,  std::vector<TBaseCell *> &cells)
{
  
  xy_coords.resize(0);
  cells.resize(0);

  std::vector<double> xLines;
  xLines.resize(0);

  // x-coordinates where we look at the solution
  xLines.push_back(-49.99);
  xLines.push_back(-10.);
  xLines.push_back(0.3);
  xLines.push_back(5.);
  xLines.push_back(10.);
  xLines.push_back(16.5);
  xLines.push_back(24.2);
  xLines.push_back(29.2);
  xLines.push_back(34.2);
  xLines.push_back(39.2);
  xLines.push_back(250.);

  // define the y-values where velocity is written
  std::vector< double > yValues;
  yValues.resize(0);
  for (unsigned int j=0; j<31; j++) {
    yValues.push_back(j*0.5);
  }
  for (unsigned int j=1; j<=15; j++){
    yValues.push_back(15+j);
  }
  for (unsigned int j=1; j<=35; j++) {
    yValues.push_back(30+2*j);
  }

  //double ySpacing = 2.; // space between two sample points
  //double yMax = 100.;
  unsigned int nyPoints = yValues.size(); //yMax/ySpacing;
  
  for (unsigned int k=0; k < xLines.size(); k++)
  {
    double x = xLines.at(k);
    for (unsigned int j=0; j < nyPoints; j++)
    {
      double y = yValues.at(j); //j*ySpacing;
      
       int n_cells = coll->GetN_Cells();
      int cell_found = -1;
      
      for(int i=0; i<n_cells; i++)
      {
	TBaseCell *c = coll->GetCell(i);
	
	if(c->PointInCell(x,y))
	{
	  cells.push_back(c);
	  cell_found = 1;
	  xy_coords.push_back(x);
	  xy_coords.push_back(y);
	}
      }
	    
      if (cell_found==-1) 
	Output::print(" *** point: ",x," ",y, " no cell found ***");
      
    }
  }
}

int main(int argc, char* argv[])
{
  TDatabase Database;
  TFEDatabase2D FEDatabase2D;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  std::ifstream fs(argv[1]);
  parmoon_db.read(fs);
  fs.close();
  //std::string time_disc = parmoon_db["time_discretization"];
  //parmoon_db.merge(ParameterDatabase::default_time_database(time_disc), true);

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  TDomain Domain(parmoon_db, argv[1]);

  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(parmoon_db);

  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

  // refine grid up to the coarsest level
  size_t n_ref = Domain.get_n_initial_refinement_steps();
  for(size_t i = 0; i < n_ref; i++)
    Domain.RegRefineAll();

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);

  // set some parameters for time stepping
  //TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters(0);

  Example_TimeNSE2D example( parmoon_db );
  // create an object of Time_NSE2D class
  Time_NSE2D tnse2d(Domain, parmoon_db, example);
  
  TimeDiscretization& tss = tnse2d.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  
  // fill the array with right x and y coordinates from the file
  std::vector<double> xy_coords;
  // cells where the corrdinates belongs to
  std::vector<TBaseCell *> cells;
  if(parmoon_db["example"].is(7))
  {
    get_output_points_inside_cells(tnse2d.get_velocity_space().GetCollection(), xy_coords, cells);
    // define the list of points where we want to compute velocity
    /*
      xy_coords.resize(0);
    cells.resize(0);


    std::vector<double> xLines;
    xLines.resize(0);

    // x-coordinates where we look at the solution
    xLines.push_back(-10.);
    xLines.push_back(0.3);
    xLines.push_back(5.);
    xLines.push_back(10.);
    xLines.push_back(16.5);
    xLines.push_back(24.2);
    xLines.push_back(29.2);
    xLines.push_back(215.);

    double ySpacing = 1.; // space between two sample points
    double yMax = 60.;
    unsigned int nyPoints = yMax/ySpacing;

    for (unsigned int k=0; k < xLines.size(); k++)
    {
      double x = xLines.at(k);
      for (unsigned int j=0; j <= nyPoints; j++)
      {
	double y = j*ySpacing;
	    
	TCollection *coll = tnse2d.get_velocity_space().GetCollection();

	int n_cells = coll->GetN_Cells();
	int cell_found = -1;
	
	for(int i=0; i<n_cells; i++)
	{
	  TBaseCell *c = coll->GetCell(i);
		
	  if(c->PointInCell(x,y))
	  {
	    cells.push_back(c);
	    cell_found = 1;
	    xy_coords.push_back(x);
	    xy_coords.push_back(y);
	  }
	}
	    
	if (cell_found==-1) 
	  Output::print(" *** point: ",x," ",y, " no cell found ***");
	    
      }
    }
    */
    cout <<"*** Number of cells containing velocity-output points "<< cells.size()<< " *** " << endl;
  }


  // assemble everything at the start time
  // this includes assembling of all A's, B's
  // and M's blocks that are necessary
  tnse2d.assemble_initial_time();

  tnse2d.output(tss.current_step_);

  double end_time = TDatabase::TimeDB->ENDTIME;
  
  LoopInfo loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1;
  int linear_iteration=0;
  
  while(TDatabase::TimeDB->CURRENTTIME < end_time - 1e-10)
  {
    tss.current_step_++;

    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    // set the time parameters
    tss.set_time_disc_parameters();
    double tau = parmoon_db["time_step_length"];
    TDatabase::TimeDB->CURRENTTIME += tau;
    Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
    
    tnse2d.assemble_matrices_rhs(0);

    LoopInfo loop_info("nonlinear");
    loop_info.print_time_every_step = true;
    loop_info.verbosity_threshold = 1;
    for(unsigned int i=0;; i++)
    {
      if(tnse2d.stopIte(i))
      {
        loop_info.finish(i,tnse2d.getFullResidual());
        linear_iteration +=i;
        loop_info_time.print(linear_iteration, tnse2d.getFullResidual());
        break;
      }
      else
        loop_info.print(i, tnse2d.getFullResidual());

      tnse2d.solve();

      if(tnse2d.imex_scheme(1))
        continue;

      tnse2d.assemble_matrices_rhs(i+1);
    }
    tnse2d.output(tss.current_step_);
    
    if(parmoon_db["example"].is(7) ) 
    {
      std::string velFileName;
      velFileName = parmoon_db["output_vtk_directory"].get<std::string>() +
        "/velocityOut/" + parmoon_db["output_basename"].get<std::string>() + "_velocities.";
      MeanVelocity::compute_velocity_on_points(tnse2d, xy_coords, cells,velFileName);
    }
  }
  loop_info_time.finish(linear_iteration, tnse2d.getFullResidual());
  Output::close_file();
  return 0;
}
