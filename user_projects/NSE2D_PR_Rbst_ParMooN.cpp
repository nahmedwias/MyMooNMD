#include "PresRobustNavierStokes.h"
#include "Example_NSPR_NSE2D.h"
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <TimeDiscretizations.h>
#include <TimeDiscRout.h>
#include <LoopInfo.h>


void get_velocity_pressure_orders(std::tuple<int,int,int> &velocity_pressure_orders, int d)
{
  int velocity_order = std::get<0>(velocity_pressure_orders);
  int pressure_order = std::get<1>(velocity_pressure_orders);
  int projection_order = std::get<2>(velocity_pressure_orders);
  int order;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5:
    case 12: case 13: case 14: case 15:
      if(velocity_order > 10)
        order = velocity_order-10;
      else
        order = velocity_order;
      break;
    case -1: case -2: case -3: case -4: case -5: case -101:
      order = velocity_order;
      break;
    case 100: case 201: case 302: case 403: case 504:
      if(d == 3)
        ErrThrow("velocity_order ", velocity_order, " not supported in 3D");
      order = velocity_order;
      break;
    // conforming fe spaces with bubbles on triangles
    case 22: case 23: case 24:
      order = velocity_order;
      break;
      // discontinuous spaces
    case -11: case -12: case -13:
      order = velocity_order*10;
      break;
  }
  
  switch(pressure_order)
  {
    case -4711:
    {
      switch(velocity_order)
      {
        case -1: case -2: case -3: case -4:
          // nonconforming pw (bi)linear velo/ pw constant pressure
          // conforming pw (bi)linear velo/ pw constant pressure (not stable !!!)
          pressure_order = -velocity_order-1;
          break;
        case 1: // discontinuous space
          pressure_order = 0;
          break;
        case 2: case 3: case 4: case 5:
        // standard conforming velo and continuous pressure
          pressure_order = velocity_order-1;          
          break;
          // discontinuous pressure spaces
          // standard conforming velo and discontinuous pressure
          // this is not stable on triangles !!!
        case 12: case 13: case 14: case 15:
        case -11: case -12: case -13: case -14:
          pressure_order = -(velocity_order-1)*10;
          break;
        case 22: case 23: case 24:
          pressure_order = -(velocity_order-11)*10;
          break;
        case 100: case 201: case 302: case 403: case 504:
          pressure_order = -(velocity_order%100 + 10)*10;
          break; 
      }
      break;
    }
    case 1: case 2: case 3: case 4: case 5:
      // pressure order is chosen correctly
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
    case 100: case 201: case 302: case 403: case 504:
      if(d == 3)
        ErrThrow("pressure_order ", pressure_order, " not supported in 3D");
      // pressure order is chosen correctly
      break;
    default:
      ErrThrow("pressure space is not chosen properly ", pressure_order);
  }
  //TDatabase::ParamDB->VELOCITY_SPACE = order;
  //TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  std::get<0>(velocity_pressure_orders) = velocity_order;
  std::get<1>(velocity_pressure_orders) = pressure_order;
  std::get<2>(velocity_pressure_orders) = projection_order;

  Output::print("velocity   space ", setw(6), std::get<0>(velocity_pressure_orders));
  Output::print("pressure space ", setw(6), std::get<1>(velocity_pressure_orders));
  Output::print("projection space ", setw(6), std::get<2>(velocity_pressure_orders));
}

using namespace std;

int main(int argc, char* argv[])
{
    // start a stopwatch which measures time spent in program parts
  Chrono timer;
  
  //  declaration of database, you need this in every program
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();  
  parmoon_db.read(argv[1]);
  
  //open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  bool linear_problem = (parmoon_db["problem_type"].is(3)
                         || parmoon_db["problem_type"].is(7));
  TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = linear_problem;
  
  TDomain domain(parmoon_db);
  
  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(parmoon_db);
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    domain.PS("Domain.ps", It_Finest, 0);
  
  auto collections = domain.get_grid_collections();
  TCollection *coll = collections.front();
  Example_NSPR_NSE2D ex(parmoon_db);
  int proj = parmoon_db["projection_order"];
  int velo = parmoon_db["velocity_order"];
  int pres = parmoon_db["pressure_order"];
  auto _spaces_orders_ = std::make_tuple(velo, pres, proj);
  get_velocity_pressure_orders(_spaces_orders_, 2);
  
  // create an object of the Navier-Stokes class
  PresRobustNavierStokes<2> pr_ns(domain, parmoon_db, ex, _spaces_orders_);
  
  // assemble linear matrices
  pr_ns.assemble_matrices();
  return 0;
}
