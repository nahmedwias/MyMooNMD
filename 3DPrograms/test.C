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
#include <LinAlg.h>
#include <FESpace3D.h>
#include <SystemMatScalar3D.h>
#include <Output3D.h>
#include <MainUtilities.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _MPI
#include "mpi.h"
#include <MeshPartition.h>
//#include <MeshPartition2D.h>
#include <ParFECommunicator_3D.h>
// #include <MumpsSolver.h>
// #include <ParVector3D.h>
// #include <ParVectorNSE3D.h>
// #include <Scalar_ParSolver.h>
#endif

#include <omp.h>

#define AMG 0
#define GMG 1
#define DIRECT 2

double bound = 0;
double timeC = 0;
// =======================================================================
// include current example
// =======================================================================
#include "../Examples/CD_3D/Laplace.h"
// =======================================================================

// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, N_Cells, ORDER, N_DOF,  img=1;
  int mg_type, mg_level, LEVELS;
  
  double *sol, *rhs, **Sol_array, **Rhs_array, t1, t2, errors[4];
  double start_time, end_time;
  double construction, assembling, solving, vtk, total;                         //time calc
  
  TDatabase *Database = new TDatabase();
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO, *PRM;
  char Name[] = "name";
  char Description[] = "description";
  char CString[] = "C";
  char SubID[] = "";
  
  int profiling;
  #ifdef _MPI
  int rank, size, out_rank;
  int MaxCpV, MaxSubDomainPerDof;
  
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  double time1, time2;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;
  #endif 
  
  TDomain *Domain;
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 


  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  
  //############################# TEST ############################################//
  double *local = new double[rank+1];
  for(i=0;i<rank+1;i++)
    local[i]=rank+1;
  
  double *global;
  TParFECommunicator_3D *C;
  C = new TParFECommunicator_3D;
  
  int globalSize;
  int root=0;
  
  C->GatherToRoot(global, globalSize, local, rank+1, root);
  
  if(rank == root)
  {
    printf("\nglobalSize=%d\n",globalSize);
    for(i=0;i<globalSize;i++)
    {
      printf("global[%d]=%lf\n",i,global[i]);
      global[i]++;
    }
  }
  
  C->ScatterFromRoot(global, globalSize, local, rank+1, root);
  
  int j;
  for(j=0;j<size;j++)
  {
    if(rank==j){
      printf("\nrank=%d\n",rank);
      for(i=0;i<rank+1;i++)
      {
	printf("local[%d]=%lf\n",i,local[i]);
      }
    }
    MPI_Barrier(Comm);
  }
  //############################# TEST ############################################//

      
  CloseFiles();
#ifdef _MPI
  MPI_Finalize(); 
#endif
  
  return 0;
} // end main















