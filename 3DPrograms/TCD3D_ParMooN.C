// =======================================================================
//
// Purpose:     main program for time-dependent scalar equation with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 24.01.15

// =======================================================================
 
#include <Domain.h>
#include <Database.h>
#include <SystemMatTimeScalar3D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <TimeUtilities.h>
#include <Solver.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

double bound = 0;
// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TCD_3D/Sin4.h"
// #include "../Examples/TCD_3D/ConstTSmooth.h"
//  #include "../Examples/TCD_3D/ConstT.h"
#include "../Examples/TCD_3D/amc.h"

int main(int argc, char* argv[])
{
// ======================================================================
// variable declaration
// ======================================================================
  int i, j, l, m, N_SubSteps, ORDER, LEVELS, mg_level, N_Cells, N_DOF, img=1;
  int N_Active, mg_type;
  
  double *sol, *oldsol, *rhs, *oldrhs, t1, t2, errors[5], **Sol_array, **Rhs_array;
  double tau, end_time, *defect, olderror, olderror1, hmin, hmax;
  
  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TCollection *coll;
  TFESpace3D **Scalar_FeSpaces, *fesp[1];
  TFEFunction3D *Scalar_FeFunction, **Scalar_FeFunctions;
  TOutput3D *Output;
  TSystemMatTimeScalar3D *SystemMatrix;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = {D000, D100, D010, D001};
  
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
  char Name[] = "name";
  char Description[] = "description";
  char CString[] = "C";
  double Linfty=0;

  std::ostringstream os;
  os << " ";   
    
 
  
// ======================================================================
// set the database values and generate mesh
// ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  
  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh  
  GEO = TDatabase::ParamDB->GEOFILE;
  Domain->Init(NULL, GEO);

  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
  
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
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
     mg_type=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SCALAR = mg_type;
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
    OutPut("=======================================================" << endl);
    OutPut("======           GEOMETRY  LEVEL ");
    OutPut(LEVELS << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level << "              ======" << endl);
    OutPut("=======================================================" << endl);   
   }
    
  Scalar_FeSpaces = new TFESpace3D*[mg_level];  
  Scalar_FeFunctions = new TFEFunction3D*[mg_level]; 
  Sol_array = new double*[mg_level];
  Rhs_array = new double*[mg_level];
  
//=========================================================================
// construct all finite element spaces
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
  for(i=0;i<LEVELS;i++)
   {   
    if(i)
     { Domain->RegRefineAll(); }

     coll=Domain->GetCollection(It_Finest, 0);
  
     // fespaces for scalar equation 
     Scalar_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);     

     //multilevel multigrid disc
     if(i==LEVELS-1 && mg_type==1) 
      {
       ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
       Scalar_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);
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
     
    if(i==LEVELS-1 && mg_type==1) 
     {  
      N_DOF = Scalar_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
      sol = new double[N_DOF];
      rhs = new double[N_DOF];
      Sol_array[mg_level-1] = sol;
      Rhs_array[mg_level-1] = rhs;

      Scalar_FeFunction = new TFEFunction3D(Scalar_FeSpaces[mg_level-1], CString, CString, sol, N_DOF);   
      Scalar_FeFunctions[mg_level-1] = Scalar_FeFunction;
     }//   if(i==LEVELS-1 && mg_type==1) 
    
   }// for(i=0;i<LEVELS;i++)
   
   oldrhs = new double[N_DOF];
   oldsol = new double[N_DOF];
   
   N_Cells = coll->GetN_Cells();
   OutPut("N_Cells   : " << N_Cells <<endl);
   OutPut("Dof all   : " << N_DOF  << endl);  
   OutPut(endl);   
//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    /** interpolate the initial value */
    Scalar_FeFunction->Interpolate(InitialCondition);   
   
    /** Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
     *  Solver: AMG_SOLVE (or) GMG  (or) DIRECT */
    SystemMatrix = new TSystemMatTimeScalar3D(mg_level, Scalar_FeSpaces, TDatabase::ParamDB->DISCTYPE, TDatabase::ParamDB->SOLVER_TYPE);
    
    /** initilize the system matrix with the functions defined in Example file */
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);
       
    /** assemble the system matrix with given aux, sol and rhs 
     *  aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
     *  otherwise, just pass it with NULL  */
    SystemMatrix->AssembleMRhs(NULL, Sol_array, Rhs_array); 
    
   /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
    memcpy(oldrhs, rhs, N_DOF*SizeOfDouble);     
//======================================================================
// produce outout at t=0
//======================================================================
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput3D(2, 2, 1, 1, Domain);

    Output->AddFEFunction(Scalar_FeFunction);

//     Scalar_FeFunction->Interpolate(Exact);   
    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }   

    /** measure errors to known solution */
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] = Scalar_FeSpaces[mg_level-1];       
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
     for(j=0;j<5;j++)
       errors[j] = 0;
     
      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);
     
      olderror = errors[0];
      olderror1 = errors[1]; 
     
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);     
//      Linfty=errors[0];
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  
       
       
       cout << "time " << TDatabase::TimeDB->CURRENTTIME << endl;
       
  coll->GetHminHmax(&hmin,&hmax);
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
  
  TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;
// exit(0);
//======================================================================
// time disc loop
//======================================================================    
   // parameters for time stepping scheme
   m = 0;
   N_SubSteps = GetN_SubSteps();
   end_time = TDatabase::TimeDB->ENDTIME; 

   UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
   UpdateRhs = TRUE; //check BilinearCoeffs in example file
   ConvectionFirstTime=TRUE;
   
//    cout << "init " << Ddot(N_DOF, sol, sol)<< endl;
   /** time loop starts */
   while(TDatabase::TimeDB->CURRENTTIME< end_time)
    {
     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

     for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
      {
       SetTimeDiscParameters();

      if(m==1)
       {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
       }

      /** copy the sol to old sol */  
      memcpy(oldsol, sol, N_DOF*SizeOfDouble);     
       
      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
   
      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);   

      /** unless the stiffness matrix or rhs change in time, it is enough to assemble only once at the begning   */
      if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
       {  
        if(UpdateRhs)
         { SystemMatrix->AssembleARhs(NULL, Sol_array, Rhs_array); }
        else
         { SystemMatrix->AssembleARhs(NULL, Sol_array, Rhs_array); }
        
        /**  M:= M + (tau*TDatabase::TimeDB->THETA1)*A
         *   rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs + [ M - (tau*THETA2)A]*oldsol **/
        SystemMatrix->AssembleSystMat(oldrhs, oldsol, rhs, sol);

        /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
        memcpy(oldrhs, rhs, N_DOF*SizeOfDouble); 
        ConvectionFirstTime = FALSE;
       }
     
      // solve the system matrix 
      SystemMatrix->Solve(sol);

      /** restore the mass matrix for the next time step unless the stiffness matrix 
       * or rhs change in time, it is not necessary to assemble the system matrix in every time step */
      if(UpdateStiffnessMat || UpdateRhs)
       {         
        SystemMatrix->RestoreMassMat();
       }
      
     } // for(l=0;l<N_SubSteps;l++) 
//======================================================================
// produce outout
//======================================================================
    if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      }
      
//======================================================================
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {      
      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);

      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);

      
      if(m>1)
       {      
        errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

        errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
       }
      
      olderror = errors[0];
      olderror1 = errors[1]; 
      
      if(Linfty<errors[0])
       Linfty=errors[0];

      OutPut( "Linfty " << Linfty << endl);      
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS)
//      exit(0); 
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

//======================================================================
// produce final outout
//======================================================================
  
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      }
      
  CloseFiles();
  
  return 0;
} // end main









