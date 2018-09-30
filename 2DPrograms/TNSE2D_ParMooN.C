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


void get_corrds_cell(TDomain domain, std::vector<double> &xy_coords,  std::vector<TBaseCell *> &cells)
{
  // define the list of points where we want to compute velocity
  std::vector<double> xLines;
  xLines.resize(8);
  // windtunnel data
  xLines[0] = -10.;
  xLines[1] = 0.3;
  xLines[2] = 5.0;
  xLines[3] = 10.;
  xLines[4] = 16.5;
  xLines[5] = 24.2;
  xLines[6] = 29.2;
  // to check open boundary profile
  xLines[7] = 215.;
  double ySpacing = 0.5;
  double yMax = 60.;
  unsigned int nyPoints = yMax/ySpacing;
  double x,y;
  for (unsigned int k=0; k < xLines.size(); k++)
  {
    x = xLines[k];
    for (unsigned int j=0; j <= nyPoints; j++)
    {
      y = j*ySpacing;

      TCollection *coll = domain.GetCollection(It_Finest, 0, -4711);
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

  /*std::ifstream f("all_data.txt");
  std::string line;
  double x, y, u, v, rmsu, rmsv;

  while ((f>> x>> y >> u >>  v >> rmsu >> rmsv) )
  {
    
    TCollection *coll = Domain.GetCollection(It_Finest, 0, -4711);
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

        //Output::print(" point: ",x," ",y);
      }
      
    }
    if (cell_found==-1) {
        Output::print(" *** point: ",x," ",y, " no cell found ***");
      }
  }
  f.close();
  if (cells.size() == 0)
    {
    Output::print(" *** WARNING: no cell found. Check the input file... ***");
    Output::print(" *** Note: the file all_data.txt shall be in the     ***");
    Output::print(" *** directory where the program is started.         ***");
  }
  */
}

int main(int argc, char* argv[])
{
  TDatabase Database;
  TFEDatabase2D FEDatabase2D;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);
  
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  TDomain Domain(parmoon_db, argv[1]);

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
    get_corrds_cell(Domain, xy_coords, cells);
    cout <<"# of cells containing output points"<< cells.size()<<endl;
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
