// =======================================================================
//
// Purpose:     main program for solving a stationary NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 27.01.2015

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <FESpace3D.h>
#include <SystemMatNSE3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <MainUtilities.h>

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
//  #include "../Examples/NSE_3D/BSExample.h" // smooth sol in unit square
// #include "../Examples/NSE_3D/AnsatzLinConst.h"
// #include "../Examples/NSE_3D/DrivenCavity3D.h"
#include "../Main_Users/Sashi/NSE_3D/DrivenCavity3D.h"
// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, j, N_Cells, N_U, N_P, N_TotalDOF, img=1, pressure_space_code;
  int Max_It, NSEType, velocity_space_code;
  int LEVELS, mg_level, mg_type;
  
  double *sol, *rhs, *defect, t1, t2, errors[4], residual, impuls_residual;
  double **Sol_array, **Rhs_array;
  double limit, u_error[6], p_error[3];
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TCollection *coll, *mortarcoll = NULL;
  TFESpace3D *velocity_space, *pressure_space, **Velocity_FeSpace, **Pressure_FeSpace, *fesp[1];
  TFEVectFunct3D **Velocity, *u;
  TFEFunction3D *p, *u1, *u2, *u3, **Pressure, *fefct[2];
  TOutput3D *Output;
  TSystemMatNSE3D *SystemMatrix;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
   
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
 
  char Name[] = "name";
  char UString[] = "u";
  char PString[] = "p";
 
  
  std::ostringstream os;
  os << " ";   
  
  mkdir(vtkdir, 0777);
    
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  
  OpenFiles();
  OutFile.setf(std::ios::scientific);
   
  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB(argv[0]);
  ExampleFile();
  
  /** needed in the new implementation */
  if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
   {  TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 1; }
 
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

//=========================================================================
// set data for multigrid
//=========================================================================  
  LEVELS = TDatabase::ParamDB->LEVELS;

  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
     mg_type=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SADDLE = mg_type;
  }
  
  if(mg_type==1)
   { mg_level =  LEVELS + 1; }
  else
   { mg_level = LEVELS; }
   
  if(TDatabase::ParamDB->SOLVER_TYPE==GMG)
   {
    OutPut("=======================================================" << endl);
    OutPut("======           GEOMETRY  LEVEL ");
    OutPut(LEVELS << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level << "              ======" << endl);
    OutPut("=======================================================" << endl);   
   }   

  Velocity_FeSpace = new TFESpace3D*[mg_level];  
  Pressure_FeSpace = new TFESpace3D*[mg_level];  
  
  Velocity = new TFEVectFunct3D*[mg_level] ;
  Pressure = new TFEFunction3D*[mg_level];
  
  Sol_array = new double*[mg_level];
  Rhs_array = new double*[mg_level];

//=========================================================================
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
  for(i=0;i<LEVELS;i++)
   {   
    if(i)
     { Domain->RegRefineAll(); }

     coll=Domain->GetCollection(It_Finest, 0);

//=========================================================================
// construct all finite element spaces
//=========================================================================  
     if(mg_type==1) // lower order FE on coarse grids 
      {
       Velocity_FeSpace[i] = new TFESpace3D(coll, Name, UString, BoundCondition, Non_USpace,1);
       Pressure_FeSpace[i] = new TFESpace3D(coll, Name, PString, BoundCondition, DiscP_PSpace,0);    
       
       if(i==LEVELS-1) // higher order on fine level
        {
         GetVelocityAndPressureSpace3D(coll, BoundCondition, velocity_space,
                                       pressure_space, &pressure_space_code,
                                       TDatabase::ParamDB->VELOCITY_SPACE,
                                       TDatabase::ParamDB->PRESSURE_SPACE); 
         Velocity_FeSpace[i+1] = velocity_space;
         Pressure_FeSpace[i+1] = pressure_space;

         // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
         TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
         velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE; 
        }
       
      }
     else
      {
       GetVelocityAndPressureSpace3D(coll, BoundCondition, velocity_space,
                                     pressure_space, &pressure_space_code,
                                     TDatabase::ParamDB->VELOCITY_SPACE,
                                     TDatabase::ParamDB->PRESSURE_SPACE); 

       Velocity_FeSpace[i] = velocity_space;
       Pressure_FeSpace[i] = pressure_space;       
       
       TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
       velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
      }
   
//======================================================================
// construct all finite element functions
//======================================================================   
     N_U = Velocity_FeSpace[i]->GetN_DegreesOfFreedom();
     N_P = Pressure_FeSpace[i]->GetN_DegreesOfFreedom();    
     N_TotalDOF = 3*N_U + N_P;    

     sol = new double[N_TotalDOF];
     memset(sol, 0, N_TotalDOF*SizeOfDouble);     
     Sol_array[i] = sol;
     
     rhs = new double[N_TotalDOF];
     memset(rhs, 0, N_TotalDOF*SizeOfDouble);     
     Rhs_array[i] = rhs;

     u = new TFEVectFunct3D(Velocity_FeSpace[i], UString,  UString,  sol, N_U, 3);
     Velocity[i] = u;
     p = new TFEFunction3D(Pressure_FeSpace[i], PString,  PString,  sol+3*N_U, N_P);    
     Pressure[i] = p;
 
     if(i==LEVELS-1 && mg_type==1)  
      {
       N_U = Velocity_FeSpace[i+1]->GetN_DegreesOfFreedom();
       N_P = Pressure_FeSpace[i+1]->GetN_DegreesOfFreedom();    
       N_TotalDOF = 3*N_U + N_P;    
    
       sol = new double[N_TotalDOF];
       memset(sol, 0, N_TotalDOF*SizeOfDouble);     
       Sol_array[i+1] = sol;
     
       rhs = new double[N_TotalDOF];
       memset(rhs, 0, N_TotalDOF*SizeOfDouble);     
       Rhs_array[i+1] = rhs;

       u = new TFEVectFunct3D(Velocity_FeSpace[i+1], UString,  UString,  sol, N_U, 3);
       Velocity[i+1] = u;
       p = new TFEFunction3D(Pressure_FeSpace[i+1], PString,  PString,  sol+3*N_U, N_P);    
       Pressure[i+1] = p;    
      }// if(i==LEVELS-1 && mg_type==1)
      
    } //  for(i=0;i<LEVELS;i++)
   
    u1 = Velocity[mg_level-1]->GetComponent(0);
    u2 = Velocity[mg_level-1]->GetComponent(1);
    u3 = Velocity[mg_level-1]->GetComponent(2);  
   
    N_Cells = coll->GetN_Cells();
    OutPut("N_Cells      : "<< setw(10) << N_Cells <<endl);
    OutPut("Dof velocity : "<< setw(10) << 3*N_U << endl);
    OutPut("Dof pressure : "<< setw(10) << N_P << endl);
    OutPut("Dof all      : "<< setw(10) << N_TotalDOF  << endl);  
//======================================================================
// produce outout
//======================================================================
   VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
   Output = new TOutput3D(2, 2, 1, 1, Domain);

   Output->AddFEVectFunct(u);
   Output->AddFEFunction(p);
   
//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    NSEType = TDatabase::ParamDB->NSTYPE;
    SystemMatrix = new TSystemMatNSE3D(mg_level, Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, 
                                       TDatabase::ParamDB->DISCTYPE, NSEType, TDatabase::ParamDB->SOLVER_TYPE);
 
    // initilize the system matrix with the functions defined in Example file
    SystemMatrix->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, U3BoundValue);
       
    
    // assemble the system matrix with given sol and rhs 
    SystemMatrix->Assemble(Sol_array, Rhs_array);
    
    // calculate the residual
    defect = new double[N_TotalDOF];
    memset(defect,0,N_TotalDOF*SizeOfDouble);
    
    SystemMatrix->GetResidual(sol, rhs, defect);
    
    //correction due to L^2_O Pressure space 
     if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20Vector3D(defect+3*N_U, N_P, pressure_space_code);
    
      residual =  Ddot(N_TotalDOF, defect, defect);
      impuls_residual = Ddot(3*N_U, defect, defect);  

      OutPut("Nonlinear iteration step   0");
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << residual-impuls_residual);
      OutPut(setw(14) << sqrt(residual) << endl);     

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
//      exit(0);
//====================================================================== 
//Solve the system
// the nonlinear iteration
//======================================================================
    limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
    Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;

    for(j=1;j<=Max_It;j++)
     {   
      // Solve the NSE system
      SystemMatrix->Solve(sol, rhs);

      //no nonlinear iteration for Stokes problem  
      if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES) 
       break;
      
      // assemble the system matrix with given aux, sol and rhs 
      SystemMatrix->AssembleNonLinear(Sol_array, Rhs_array);  
      
      // get the residual
      memset(defect,0,N_TotalDOF*SizeOfDouble);
      SystemMatrix->GetResidual(sol, rhs, defect);
      
    //correction due to L^2_O Pressure space 
     if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20Vector3D(defect+3*N_U, N_P, pressure_space_code);
    
      residual =  Ddot(N_TotalDOF, defect, defect);
      impuls_residual = Ddot(3*N_U, defect, defect);  

      OutPut("nonlinear iteration step " << setw(3) << j);
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << residual-impuls_residual);
      OutPut(setw(14) << sqrt(residual) << endl);
 
      if ((sqrt(residual)<=limit)||(j==Max_It))
       {
        break;
       }//if ((sqrt(residual)<=limit)||(j==Max_It))
       
     } //for(j=1;j<=Max_It;j
    
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        IntoL20FEFunction3D(sol+3*N_U, N_P, Pressure_FeSpace[mg_level-1]);

//  ==============================================================  
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
 
   
   
   
//    u1->Interpolate(ExactU1);
//    u2->Interpolate(ExactU2);
//    u3->Interpolate(ExactU3);
//    p->Interpolate(ExactP);
//    
//     if(TDatabase::ParamDB->WRITE_VTK)
//      {
//       os.seekp(std::ios::beg);
//        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
//          else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
//           else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
//            else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
//             else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
//       Output->WriteVtk(os.str().c_str());
//       img++;
//      }   
//  



//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {   
      SystemMatrix->MeasureErrors(ExactU1, ExactU2,ExactU3,  ExactP, u_error, p_error);

       OutPut("L2(u): " <<  sqrt(u_error[0]*u_error[0]+u_error[2]*u_error[2]+u_error[4]*u_error[4]) << endl);
       OutPut("H1-semi(u): " <<  sqrt(u_error[1]*u_error[1]+u_error[3]*u_error[3]+u_error[5]*u_error[5]) << endl);
       OutPut("L2(p): " <<  p_error[0] << endl);
       OutPut("H1-semi(p): " <<  p_error[1] << endl); 
     } // if(TDatabase::ParamDB->MEASURE_ERRORS)


  CloseFiles();
  
  return 0;
} // end main
