// =======================================================================
//
// Purpose:     main program for computing POD basis in TCD problems 
//
// Author:      Swetlana Giere & Alfonso Caiazzo
//
// History:     Started 15.1.19

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>

#include <TimeConvectionDiffusionPOD.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <Example_TimeCD2D.h>
#include <TimeDiscRout.h>


using namespace std;

int main(int argc, char* argv[])
{
  double t_start = GetTime();
  TDatabase Database(argv[1]);
  TFEDatabase2D FEDatabase;
  
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);
  
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  TDomain Domain(parmoon_db);
  
  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();
  
  // refine grid up to the coarsest level
  Domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);

  auto collections = Domain.get_grid_collections();
  TCollection& cellCollection = *collections.front();
  TimeConvectionDiffusionPOD<2> tcd_pod(cellCollection, parmoon_db);
  tcd_pod.compute_pod_basis();
  //pod_2d.read_basis();
  
  tcd_pod.output();
  
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();
  return 0;
} // end main
