// =======================================================================
//
// Purpose:     main program for scalar equations with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 08.08.2014

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>

#include <Time_CD2D.h>
#include <POD_TCDR2D.h>



#include <sys/stat.h>
#include <sys/types.h>

#include <Example_TimeCD2D.h>
#include <TimeDiscRout.h>


using namespace std;

int main(int argc, char* argv[])
{
  double t_start = GetTime();
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();

  for (int i=1; i<argc; ++i)
  {
	std::ifstream fs(argv[i]);
	parmoon_db.read(fs);
	fs.close();
  }

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  TDomain Domain(parmoon_db, argv[1] );
  
  Output::set_outfile(parmoon_db["outfile"]);
  Output::setVerbosity(parmoon_db["verbosity"]);

  parmoon_db.write(Output::get_outfile());
  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();
  
  // refine grid up to the coarsest level
  size_t n_ref = Domain.get_n_initial_refinement_steps();
  for(unsigned int i=0; i<n_ref; i++){
    Domain.RegRefineAll();  
  }
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    Domain.PS("Domain.ps", It_Finest, 0);

  TCollection *coll = Domain.GetCollection(It_Finest, 0);
  Example_TimeCD2D example( parmoon_db );

  POD_TCDR2D pod_2d(*coll, parmoon_db, example);
  
  pod_2d.compute_pod_basis();
  //pod_2d.read_basis();
  
  pod_2d.output();
  
  // ======================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  Output::print("used time: ", GetTime() - t_start, "s");
  // ======================================================================
  Output::close_file();
  return 0;
} // end main
