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
// #include <ParFECommunicator3D.h>
// #include <MumpsSolver.h>
// #include <ParVector3D.h>
// #include <ParVectorNSE3D.h>
// #include <Scalar_ParSolver.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

#define profiling 1

double bound = 0;
// =======================================================================
// include current example
// =======================================================================
#include "../Examples_All/CD_3D/Laplace.h"
// =======================================================================

 double timeC;


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
  
  TDatabase *Database = new TDatabase();
  char SubID[] = "";
  
  if(profiling)	start_time = GetTime();
  
  #ifdef _MPI
  int rank, size, out_rank;
  int MaxCpV, MaxSubDomainPerDof;
  
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  double time1, time2;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;
    
  if(profiling){
    start_time = MPI_Wtime();
  }

  TFESpace3D **OwnScalar_Spaces;
  TCollection *own_coll;
  #endif 

  TDomain *Domain;
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TCollection *coll;
  TFESpace3D **Scalar_FeSpaces, *fesp[1];
  TFEFunction3D *Scalar_FeFunction, **Scalar_FeFunctions;
  TOutput3D *Output;
  TSystemMatScalar3D *SystemMatrix;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };

  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
  char Name[] = "name";
  char Description[] = "description";
  char CString[] = "C";
   
  std::ostringstream os;
  os << " ";   
      
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  
  OpenFiles();
  OutFile.setf(std::ios::scientific);

  #ifdef _MPI
  out_rank=TDatabase::ParamDB->Par_P0;
  //out_rank = 0;
  if(rank == out_rank)
  #endif
   {
  Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();
   }

  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  
  GEO = TDatabase::ParamDB->GEOFILE;
  Domain->Init(NULL, GEO);
   
  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  

#ifdef _MPI
  Domain->GenerateEdgeInfo();
  
  if(profiling)  t1 = MPI_Wtime();
  Partition_Mesh3D(Comm, Domain, MaxCpV);	//MaxCpV=maximum cell per vertex
  if(profiling)  t2 = MPI_Wtime(); 
  
  if(profiling){
    time2 = t2-t1;
    MPI_Reduce(&time2, &time1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
    if(rank == out_rank)
      printf("Time taken for Domain Decomposition is %e\n", time1);
  }

  Domain->GenerateEdgeInfo();
  MaxSubDomainPerDof = MIN(MaxCpV, size);
  TDatabase::ParamDB->WRITE_PS = 0;
#endif 
 
  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }
  
  if(TDatabase::ParamDB->WRITE_VTK)
   { mkdir(vtkdir, 0777); }
   
   
//=========================================================================
// set data for multigrid
//=========================================================================  
  LEVELS = TDatabase::ParamDB->LEVELS;

  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
     mg_type=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SCALAR =  mg_type;
  }
  
  if(mg_type)
   {
    mg_level =  LEVELS + 1;
    ORDER = -1;
   }
  else
   {
    mg_level = LEVELS;
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
   }
   
  if(TDatabase::ParamDB->SOLVER_TYPE==GMG)
   {
#ifdef _MPI  
    if(rank == out_rank)
    {
#endif 
    OutPut("=======================================================" << endl);
    OutPut("======           GEOMETRY  LEVEL ");
    OutPut(LEVELS-1 << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level-1 << "              ======" << endl);
    OutPut("=======================================================" << endl);   
#ifdef _MPI 
    }
#endif 
   }
    
  Scalar_FeSpaces = new TFESpace3D*[LEVELS+1];  
  Scalar_FeFunctions = new TFEFunction3D*[LEVELS+1]; 
  Sol_array = new double*[LEVELS+1];
  Rhs_array = new double*[LEVELS+1];
#ifdef _MPI    
  OwnScalar_Spaces = new TFESpace3D*[LEVELS+1];   
#endif 

//=========================================================================
// construct all finite element spaces
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
  for(i=0;i<LEVELS;i++)
   {  
    if(i)
     { Domain->RegRefineAll(); }
     
     #ifdef _MPI
     if(rank == out_rank)
       printf("Level :: %d\n\n",i);
     if(i)
     {
       Domain->GenerateEdgeInfo();
       Domain_Crop(Comm, Domain);       // removing unwanted cells in the hallo after refinement 
    }
     #endif
     
     coll = Domain->GetCollection(It_Finest, 0);
  
     // fespaces for scalar equation 
     Scalar_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);     
      
     #ifdef _MPI
     Scalar_FeSpaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
     
     own_coll = Domain->GetOwnCollection(It_Finest, 0, rank);
     OwnScalar_Spaces[i] = new TFESpace3D(own_coll,  Name, Description, BoundCondition, ORDER); 
     #endif
     
     //multilevel multigrid disc
     if(i==LEVELS-1 && i!=mg_level-1) 
      {
       ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
       Scalar_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);
       #ifdef _MPI
       Scalar_FeSpaces[mg_level-1]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
       
       OwnScalar_Spaces[mg_level-1] = new TFESpace3D(own_coll,  Name, Description, BoundCondition, ORDER); 
       #endif
      } //  if(i==LEVELS-1 && i!=mg_level-1) 
     
//======================================================================
// construct all finite element functions
//======================================================================
    N_DOF = Scalar_FeSpaces[i]->GetN_DegreesOfFreedom();
    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    Sol_array[i] = sol;
    Rhs_array[i] = rhs;   
   
    Scalar_FeFunction  = new TFEFunction3D(Scalar_FeSpaces[i], CString, CString, sol, N_DOF);  
    Scalar_FeFunctions[i] = Scalar_FeFunction;
     
    if(i==LEVELS-1 && i!=mg_level-1) 
     {  
      N_DOF = Scalar_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
      sol = new double[N_DOF];
      rhs = new double[N_DOF];
      Sol_array[mg_level-1] = sol;
      Rhs_array[mg_level-1] = rhs;

      Scalar_FeFunction = new TFEFunction3D(Scalar_FeSpaces[mg_level-1], CString, CString, sol, N_DOF);   
      Scalar_FeFunctions[mg_level-1] = Scalar_FeFunction;
     }//  if(i==LEVELS-1 && i!=mg_level-1) 
    

   }// for(i=0;i<LEVELS;i++)

   N_Cells = coll->GetN_Cells();
   OutPut("N_Cells   : " << N_Cells <<endl);
   OutPut("Dof all   : " << N_DOF  << endl);  
//    OutPut(endl);  

//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) GLS (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
    if(profiling)	t1 = GetTime();
    SystemMatrix = new TSystemMatScalar3D(mg_level, Scalar_FeSpaces, TDatabase::ParamDB->DISCTYPE, TDatabase::ParamDB->SOLVER_TYPE
    #ifdef _MPI
					  , OwnScalar_Spaces, Scalar_FeFunctions, Comm
    #endif
											);
    if(profiling){
      t2 = GetTime();
      t2 = t2-t1;
#ifdef _MPI
      MPI_Reduce(&t2, &t1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
      if(rank == out_rank)
#endif
      {
	OutPut(endl);
	OutPut( "time for constructing and initializing System Mat Components: " << t1 << endl);
      }
    } 
   
    // initilize the system matrix with the functions defined in Example file
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);

    // assemble the system matrix with given aux, sol and rhs 
    // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass it with NULL 
    if(profiling)	t1 = GetTime();
    SystemMatrix->Assemble(NULL, Sol_array, Rhs_array);
    if(profiling){
      t2 = GetTime();
      t2 = t2-t1;
#ifdef _MPI
      MPI_Reduce(&t2, &t1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
      if(rank == out_rank)
#endif
      {
	OutPut( "time for assembling: " << t1 << endl);
      }
    }
    
    //Solve the system
    if(profiling)	t1 = GetTime();
    SystemMatrix->Solve(sol, rhs);
    if(profiling){
      t2 = GetTime();
      t2 = t2-t1;
#ifdef _MPI
      MPI_Reduce(&t2, &t1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
      if(rank == out_rank)
#endif
      {
	OutPut( "time for solving: " << t1 << endl);
      }
    }
   
//======================================================================
// produce outout
//======================================================================
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput3D(2, 2, 1, 1, Domain);

    Output->AddFEFunction(Scalar_FeFunction);

     //Scalar_FeFunction->Interpolate(Exact);   
#ifdef _MPI
     double t_par1, t_par2;
     if(profiling)		t_par1 = MPI_Wtime();

     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;       
     if(profiling){
       t_par2 = MPI_Wtime();
       t_par2 = t_par2-t_par1;
       MPI_Reduce(&t_par2, &t_par1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
       if(rank == out_rank)
	printf("Time taken for writing the parvtk file %e\n", (t_par1));
     }
#else
    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      
     }
       img++;    
#endif 
     
//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] = Scalar_FeSpaces[mg_level-1];
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
     
#ifdef _MPI
       Scalar_FeFunction->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors,Comm);
#else  
      Scalar_FeFunction->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);
#endif
      
      delete aux;
      
#ifdef _MPI
      double reduced_errors[4];
      MPI_Reduce(errors, reduced_errors, 4, MPI_DOUBLE, MPI_SUM, out_rank, Comm);
      if(rank == out_rank){
	OutPut(endl);
	OutPut( "L2: " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi: " << sqrt(reduced_errors[1]) << endl);
	OutPut( "SD: " << sqrt(reduced_errors[2]) << endl);
      }
#else
      OutPut(endl);
      OutPut( "L2: " << errors[0] << endl);
      OutPut( "H1-semi: " << errors[1] << endl);
      OutPut( "SD: " << errors[2] << endl);
#endif
     } // if(TDatabase::ParamDB->MEASURE_ERRORS)

  
  
  if(profiling)	end_time = GetTime();
  #ifdef _MPI
  if(profiling){
    //printf("rank = %d out_rank = %d\n",rank,out_rank);
    end_time = MPI_Wtime();
    MPI_Reduce(&start_time, &t1, 1, MPI_DOUBLE, MPI_MIN, out_rank, Comm);
    MPI_Reduce(&end_time,   &t2, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
    if(rank == out_rank)
      OutPut( "Total time taken: " << (t2-t1) << endl);
  }
  MPI_Finalize();  
  #else 
  if(profiling)	OutPut( "Total time taken: " << (end_time-start_time) << endl);
  #endif
  CloseFiles();
  return 0;
} // end main
