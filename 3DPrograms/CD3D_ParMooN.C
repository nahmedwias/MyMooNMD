// =======================================================================
//
// Purpose:     main program for solving a 3D stationary scalar equation using ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 23.01.2015
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <CD3D.h>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef _MPI
// we need this here because for some reason (??) these are declared extern in
// e.g. TParFECommunicator3D
double bound = 0;
double timeC = 0;
#endif

// a short method to get the current time, in case the parameter is set
double get_time()
{
  if(TDatabase::ParamDB->timeprofiling)
  {
  #ifdef _MPI
    return MPI_Wtime();
  #else
    return GetTime();
  #endif
  }
  return 0.;
}

// main program
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase Database;
#ifdef _MPI
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  TDatabase::ParamDB->Comm = Comm;
#endif
  TFEDatabase3D FEDatabase;
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(argv[1]);
  
  double start_time = get_time();
  
  //set PROBLEM_TYPE to CD if not yet set
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 1;
  //open OUTFILE, this is where all output is written to (addionally to console)
  OpenFiles();
  
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
  
  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */
  domain.Init(NULL, TDatabase::ParamDB->GEOFILE); // call mesh generator
              
  // refine grid up to the coarsest level
  for(int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
  {
    domain.RegRefineAll();
    #ifdef _MPI
    domain.GenerateEdgeInfo();
    Domain_Crop(Comm, &domain);
    #endif
  }
  
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    domain.PS("Domain.ps", It_Finest, 0);
  
#ifdef _MPI
  domain.GenerateEdgeInfo();
  int MaxCpV; //MaxCpV=maximum cell per vertex
  Partition_Mesh3D(Comm, &domain, MaxCpV);
  domain.GenerateEdgeInfo();
  int size; 
  MPI_Comm_size(Comm, &size);
  int MaxSubDomainPerDof = MIN(MaxCpV, size);
  OutPut("MaxSubDomainPerDof " << MaxSubDomainPerDof << endl);
#endif
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  // choose example according to the value of TDatabase::ParamDB->EXAMPLE
  Example_CD3D example;
  
  //=========================================================================
  CD3D cd3d(&domain, &example);
  cd3d.assemble();
  cd3d.solve();
  cd3d.output();
  //=========================================================================
  
  double
  
  CloseFiles();
  return 0;
  
} // end main


