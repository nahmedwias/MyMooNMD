// =======================================================================
//
// Purpose:     main program for solving a stationary NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 23.08.2014

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <SystemNSE2D.h>
#include <Output2D.h>
#include <MainUtilities.h>
#include <LocalAssembling2D.h>
#include <Example_NSE2D.h>

#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase Database;
  TFEDatabase2D FEDatabase; 
  
  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain Domain(argv[1]);  
  
  //set PROBLEM_TYPE to NSE if not yet set (3 means Stokes, 5 Naver-Stokes)
  if(TDatabase::ParamDB->PROBLEM_TYPE!=3 && TDatabase::ParamDB->PROBLEM_TYPE!=5)
    TDatabase::ParamDB->PROBLEM_TYPE = 5;
  //open OUTFILE, this is where all output is written to (addionally to console)
  OpenFiles();
  
  // possibly change parameters in the database, if they are not meaningful now
  Database.CheckParameterConsistencyNSE();
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database.WriteParamDB(argv[0]);
  
  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */
  Domain.Init(NULL, TDatabase::ParamDB->GEOFILE); // call mesh generator
  
  // refine grid up to the coarsest level
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    Domain.RegRefineAll();  
  
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain.PS("Domain.ps",It_Finest,0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  Example_NSE2D example;
  
  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  // a collection is basically only an array of cells, which is needed to create
  // a finite element space
  TCollection *coll = Domain.GetCollection(It_Finest, 0);
  TCollection *mortarcoll = NULL;
  // print out some information about the mesh
  int N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
  // create finite element spaces for velocity and pressure
  TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace;
  int pressure_space_code;
  // the following function will create an inf-sup stable pair of finite element
  // spaces, if PRESSURE_SPACE is not set (-4711)
  GetVelocityAndPressureSpace(coll, example.get_bc(0), mortarcoll, 
                              Velocity_FeSpace, Pressure_FeSpace,
                              &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  
  // defaulty inf-sup pressure space will be selected based on the velocity 
  // space, so update it in database
  TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
  int velocity_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
  // print out some information on the finite element spaces
  int N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
  int N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();    
  int N_TotalDOF = 2*N_U + N_P;
  OutPut("Dof Velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("Dof Pressure : "<< setw(10) << N_P << endl);
  OutPut("Total Dof all: "<< setw(10) << N_TotalDOF  << endl);
  
  //======================================================================
  // construct all finite element functions
  //======================================================================
  double *sol = new double[N_TotalDOF];
  double *rhs = new double[N_TotalDOF];
  // set solution and right hand side vectors to zero
  memset(sol, 0, N_TotalDOF*SizeOfDouble);
  memset(rhs, 0, N_TotalDOF*SizeOfDouble);

  // create the finite element functions
  TFEVectFunct2D Velocity(Velocity_FeSpace, (char*)"u", (char*)"u",  sol, N_U,
                          2);
  TFEFunction2D Pressure(Pressure_FeSpace, (char*)"p", (char*)"p", sol+2*N_U,
                         N_P);
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  TFEFunction2D *fe_functions[3] = { Velocity.GetComponent(0),
                                     Velocity.GetComponent(1), &Pressure };
  
  //======================================================================
  // SystemMatrix construction and solution
  //======================================================================  
  TSystemMatNSE2D SystemMatrix(&Velocity, &Pressure);
 
  // initilize the system matrix with the functions defined in Example file
  SystemMatrix.Init(example.get_bd(0), example.get_bd(1));
  
  // create a local assembling objects which are needed to assemble the matrices
  LocalAssembling2D la(NSE2D_Galerkin, fe_functions, example.get_coeffs());
  LocalAssembling2D la_nonlinear(NSE2D_Galerkin_Nonlinear, fe_functions,
                                 example.get_coeffs());
  
  // assemble the system matrix with given sol and rhs 
  SystemMatrix.Assemble(la, sol, rhs);
  // if solution was not zero up to here, you should call 
  //SystemMatrix.AssembleNonLinear(la_nonlinear, sol, rhs);
  
  // calculate the residual
  double *defect = new double[N_TotalDOF];
  memset(defect,0,N_TotalDOF*SizeOfDouble);
  
  SystemMatrix.GetResidual(sol, rhs, defect);
  
  //correction due to L^2_O Pressure space 
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code);
  
  double residual =  Ddot(N_TotalDOF, defect, defect);
  double impuls_residual = Ddot(2*N_U, defect, defect);  

  OutPut("Nonlinear iteration step   0");
  OutPut(setw(14) << impuls_residual);
  OutPut(setw(14) << residual-impuls_residual);
  OutPut(setw(14) << sqrt(residual) << endl);

  //====================================================================== 
  //Solve the system
  // the nonlinear iteration
  //======================================================================
  double limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  
  for(int j=1; j<=Max_It; j++)
  {
    // Solve the NSE system
    SystemMatrix.Solve(sol, rhs);
    
    //no nonlinear iteration for Stokes problem
    if(TDatabase::ParamDB->PROBLEM_TYPE == 3)
      break;
    
    // assemble the system matrix with given aux, sol and rhs 
    SystemMatrix.AssembleNonLinear(la_nonlinear, sol, rhs);
    
    // get the residual
    memset(defect,0,N_TotalDOF*SizeOfDouble);
    SystemMatrix.GetResidual(sol, rhs, defect);
    
    //correction due to L^2_O Pressure space 
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code);
    
    residual = Ddot(N_TotalDOF, defect, defect);
    impuls_residual = Ddot(2*N_U, defect, defect);
    
    OutPut("nonlinear iteration step " << setw(3) << j);
    OutPut(setw(14) << impuls_residual);
    OutPut(setw(14) << residual-impuls_residual);
    OutPut(setw(14) << sqrt(residual) << endl);
    
    if((sqrt(residual)<=limit) || (j==Max_It))
    { // stop the nonlinear iteration
      break;
    }
  } //for(j=1;j<=Max_It;j
  
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20FEFunction(sol+2*N_U, N_P, Pressure_FeSpace, velocity_space_code, 
                      pressure_space_code);
  
  //======================================================================
  // produce outout
  //======================================================================
  TOutput2D Output(2, 2, 1, 1, &Domain);
  
  Output.AddFEVectFunct(&Velocity);
  Output.AddFEFunction(&Pressure);
  
  //fe_functions[0]->Interpolate(ExactU1);
  //fe_functions[1]->Interpolate(ExactU2);
  //fe_functions[2]->Interpolate(ExactP);
   
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
    filename += "/" + std::string(TDatabase::ParamDB->BASENAME) + ".vtk";
    Output.WriteVtk(filename.c_str());
  }   
   
  //====================================================================== 
  // measure errors to known solution
  //======================================================================    
  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    double err[4];
    TAuxParam2D NSEaux_error;
    MultiIndex2D NSAllDerivatives[3] = {D00, D10, D01};
    
    // errors in first velocity component
    fe_functions[0]->GetErrors(example.get_exact(0), 3, NSAllDerivatives, 2, 
                               L2H1Errors, NULL, &NSEaux_error, 1,
                               &Velocity_FeSpace, err);
    
    // errors in second velocity component
    fe_functions[1]->GetErrors(example.get_exact(1), 3, NSAllDerivatives, 2,
                               L2H1Errors, NULL, &NSEaux_error, 1,
                               &Velocity_FeSpace, err+2);
    OutPut("L2(u)     : " << sqrt(err[0]*err[0] + err[2]*err[2]) << endl);
    OutPut("H1-semi(u): " << sqrt(err[1]*err[1] + err[3]*err[3]) << endl);
    
    // errors in pressure
    fe_functions[2]->GetErrors(example.get_exact(2), 3, NSAllDerivatives, 2,
                               L2H1Errors, NULL, &NSEaux_error, 1,
                               &Pressure_FeSpace, err);
    OutPut("L2(p)     : " <<  err[0] << endl);
    OutPut("H1-semi(p): " <<  err[1] << endl); 
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
  
  CloseFiles();
  return 0;
} // end main
