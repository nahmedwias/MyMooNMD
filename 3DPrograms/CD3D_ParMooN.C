/**
 * @brief Main program for solving a 3D stationary scalar equation using ParMooN.
 * @author Sashikumaar Ganesan, Clemens Bartsch
 *
 * Implementation started on 2015/01/23. Rework since 2015/10/19.
 *
 */
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <CD3D.h>
#include <Example_CD3D.h>
#include <MeshPartition.h>

#include <sys/stat.h>

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif


// main program
// =======================================================================
int main(int argc, char* argv[])
{
  {
#ifdef _MPI
  //Construct and initialise the default MPI communicator.
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(my_rank==0)
  {
    Output::print("<<<<< Running ParMooN: CD3D Main Program >>>>>");
    Output::info("CD3D", "MPI, using ", size, " processes");
  }
#else
  int my_rank = 0;
  Output::print("<<<<< Running ParMooN: CD3D Main Program >>>>>");
  Output::info("CD3D", "SEQUENTIAL (or OMP...)");
#endif

//  //start a stopwatch which measures time spent in program parts
//  Chrono chrono_parts;

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

  // Construct domain, thereby read in controls from the input file.
  TDomain domain(argv[1], parmoon_db);

  //open OUTFILE, this is where all output is written to (addionally to console)
  if(my_rank==0)
  {
    Output::set_outfile(parmoon_db["outfile"]);
  }
  Output::setVerbosity(parmoon_db["verbosity"]);

  if(my_rank==0) //Only one process should do that.
    Database.WriteParamDB(argv[0]);

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"] && my_rank==0)
    domain.PS("Domain.ps", It_Finest, 0);

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
  domain.print_info("cd3d domain");

  // Choose and construct example.
  Example_CD3D example(parmoon_db);

  // Construct the cd3d problem object.
#ifdef _MPI
  CD3D cd3d(gridCollections, parmoon_db, example, maxSubDomainPerDof);
#else
  CD3D cd3d(gridCollections, parmoon_db, example);
#endif

  //=========================================================================
  //Start the actual computations.
  //=========================================================================

  cd3d.assemble(); // assemble matrix and rhs

  cd3d.solve();    // solve the system

  cd3d.output();   // produce nice output

  //=========================================================================

  if(my_rank==0)
    Output::print("<<<<< ParMooN Finished: CD3D Main Program >>>>>");

  Output::close_file();

}

#ifdef _MPI
  MPI_Finalize();
#endif

  return 0;
} // end main


