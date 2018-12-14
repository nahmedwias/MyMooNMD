#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <Example_TimeNSE2D.h>
#include "TimeNavierStokes.h"
#include <TimeDiscretizations.h>
#include <TimeDiscRout.h>
#include <LoopInfo.h>
#include <MeanVelocity.h>

using namespace std;


void get_corrds_cell(const TDomain& domain, std::vector<double> &xy_coords,  std::vector<TBaseCell *> &cells)
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
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase2D;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);
  
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  TDomain Domain(parmoon_db);

  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(parmoon_db);

  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();

  // refine grid
  Domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);

  // set some parameters for time stepping
  SetTimeDiscParameters(0);
// fill the array with right x and y coordinates from the file
  std::vector<double> xy_coords;
  // cells where the corrdinates belongs to
  std::vector<TBaseCell *> cells;
  if(parmoon_db["example"].is(7))
  {
    get_corrds_cell(Domain, xy_coords, cells);
    cout <<"# of cells containing output points"<< cells.size()<<endl;
  }
  // create an object of TimeNavierStokes<2> class
  TimeNavierStokes<2> tnse2d(Domain, parmoon_db);
  
  TimeDiscretization& tss = tnse2d.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  
//   // fill the array with right x and y coordinates from the file
//   std::vector<double> xy_coords;
//   // cells where the corrdinates belongs to
//   std::vector<TBaseCell *> cells;
//   if(parmoon_db["example"].is(7))
//   {
//     get_corrds_cell(Domain, xy_coords, cells);
//     cout <<"# of cells containing output points"<< cells.size()<<endl;
//   }
  // assemble everything at the start time
  // this includes assembling of all A's, B's
  // and M's blocks that are necessary
  tnse2d.assemble_initial_time();

  tnse2d.output();

  
  LoopInfo loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1;
  int linear_iteration=0;
  
  double end_time = tss.get_end_time();
  TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();
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
      if(tnse2d.stop_it(i))
      {
        loop_info.finish(i,tnse2d.get_full_residual());
        linear_iteration +=i;
        loop_info_time.print(linear_iteration, tnse2d.get_full_residual());
        break;
      }
      else
        loop_info.print(i, tnse2d.get_full_residual());

      tnse2d.solve();

      if(tnse2d.imex_scheme())
        continue;

      tnse2d.assemble_matrices_rhs(i+1);
    }
    tnse2d.output();
    
    if(parmoon_db["example"].is(7) ) 
    {
      std::string velFileName;
      velFileName = parmoon_db["output_vtk_directory"].get<std::string>() +
        "/velocityOut/" + parmoon_db["output_basename"].get<std::string>() + "_velocities.";
      MeanVelocity::compute_velocity_on_points(tnse2d, xy_coords, cells,velFileName);
    }
  }
  loop_info_time.finish(linear_iteration, tnse2d.get_full_residual());
  Output::close_file();
  return 0;
}
