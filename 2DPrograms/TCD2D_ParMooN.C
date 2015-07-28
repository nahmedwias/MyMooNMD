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
#include <SystemMatTimeScalar2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
#include <LocalAssembling2D.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// include current example
// =======================================================================
#include "../Examples/TCD_2D/exp.h"
// #include "../Examples/TCD_2D/SinCos1.h"
// #include "../Examples_All/TCD_2D/Time3.h"
// #include "../Examples/TCD_2D/exp_0.h"
//    #include "../Examples/TCD_2D/exp_2.h"
// #include "../Examples_All/TCD_2D/exp_1.h"

int main(int argc, char* argv[])
{
  TDatabase Database;
  TFEDatabase2D FEDatabase;
  
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  TDomain Domain(argv[1]);  
  
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 2;
  OpenFiles();

  Database.WriteParamDB(argv[0]);
  Database.WriteTimeDB();
  ExampleFile();
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh  
  Domain.ReadGeo(TDatabase::ParamDB->GEOFILE);

  // refine grid up to the coarsest level
  for(int i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    Domain.RegRefineAll();  
  
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain.PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
   
  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  TCollection *coll = Domain.GetCollection(It_Finest, 0);
  int N_Cells = coll->GetN_Cells();
  OutPut("N_Cells (space) : " << N_Cells <<endl);
  
  // fespaces for scalar equation 
  TFESpace2D *Scalar_FeSpace = new TFESpace2D(coll, (char*)"fe space", 
                                              (char*)"solution space", 
                                              BoundCondition, ORDER, NULL);
   
  int N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
  int N_Active =  Scalar_FeSpace->GetActiveBound();
  OutPut("dof all      : "<< setw(10) << N_DOF  << endl);
  OutPut("dof active   : "<< setw(10) << N_Active << endl);
   
//======================================================================
// construct all finite element functions
//======================================================================
  double *sol = new double[N_DOF];
  double *rhs = new double[N_DOF];
  double *oldrhs = new double[N_DOF];
    
  memset(sol, 0, N_DOF*SizeOfDouble);
  memset(rhs, 0, N_DOF*SizeOfDouble);

  TFEFunction2D *Scalar_FeFunction = new TFEFunction2D(
    Scalar_FeSpace, (char*)"sol", (char*)"sol", sol, N_DOF); 
  
  //interpolate the initial value
  Scalar_FeFunction->Interpolate(InitialCondition);

  //======================================================================
  // SystemMatrix construction and solution
  //======================================================================    
  // as LocalAssembling2D_type choose either
  // TCD2D_Mass_Rhs_Galerkin+TCD2D_Stiff_Rhs_Galerkin  or
  // TCD2D_Mass_Rhs_SUPG+TCD2D_Stiff_Rhs_SUPG
  LocalAssembling2D_type mass_type = TCD2D_Mass_Rhs_Galerkin;
  LocalAssembling2D_type stiff_type = TCD2D_Stiff_Rhs_Galerkin;
  if(0) // choose SUPG
  {
    mass_type = TCD2D_Mass_Rhs_SUPG;
    stiff_type = TCD2D_Stiff_Rhs_SUPG;
  }
  LocalAssembling2D la_mass_rhs(mass_type, &Scalar_FeFunction,
                                BilinearCoeffs);
  LocalAssembling2D la_stiff_rhs(stiff_type, &Scalar_FeFunction,
                                 BilinearCoeffs);
  
  // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
  // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
  TSystemMatTimeScalar2D SystemMatrix(Scalar_FeSpace);
  
  // initilize the system matrix with the functions defined in Example file
  SystemMatrix.Init(BilinearCoeffs, BoundCondition, BoundValue);
     
  // assemble the system matrix with given aux, sol and rhs 
  // aux is used to pass  addition fe functions (eg. mesh velocity) that is 
  // nedded for assembling, otherwise, just pass with NULL 
  SystemMatrix.AssembleMRhs(la_mass_rhs, sol, rhs);
  
  //======================================================================
  // produce outout at t=0
  //======================================================================
  TOutput2D Output(2, 2, 1, 1, &Domain);
  Output.AddFEFunction(Scalar_FeFunction);

  //Scalar_FeFunction->Interpolate(Exact);
  int img = 1;
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
    filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
    if(img<10) filename += ".0000";
    else if(img<100) filename += ".000";
    else if(img<1000) filename += ".00";
    else if(img<10000) filename += ".0";
    else filename += ".";
    filename += std::to_string(img) + ".vtk";
    Output.WriteVtk(filename.c_str());
    img++;
  }

  double errors[5] = { 0., 0., 0., 0., 0. };
  double Linfty; // L^infty in time of L2 in space
  double olderror, olderror1;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  // measure errors to known solution
  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    TFESpace2D *fesp[1] = { Scalar_FeSpace }; 
    TAuxParam2D aux;
    
    Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, 
                                 BilinearCoeffs, &aux, 1, fesp, errors);
   
    olderror = errors[0];
    olderror1 = errors[1]; 
   
    OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
    OutPut(" L2: " << errors[0]);
    OutPut(" H1-semi: " << errors[1] << endl);     
    Linfty = errors[0];
  } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  
  
  {
    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
  }

  // TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;
  //======================================================================
  // time disc loop
  //======================================================================    
  // parameters for time stepping scheme
  int m = 0;
  int N_SubSteps = GetN_SubSteps();
  double end_time = TDatabase::TimeDB->ENDTIME; 

  bool UpdateStiffnessMat = false; // check BilinearCoeffs in example file
  bool UpdateRhs = true; // check BilinearCoeffs in example file
  bool ConvectionFirstTime = true;

  // time loop starts
  while(TDatabase::TimeDB->CURRENTTIME < end_time)
  {
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(int l=0; l<N_SubSteps; l++) // sub steps of fractional step theta
    {
      SetTimeDiscParameters(1);

      if(m==1)
      {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
      }

      double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
       
      OutPut("\nCURRENT TIME: " << TDatabase::TimeDB->CURRENTTIME << endl);
      
      //copy rhs to oldrhs
      memcpy(oldrhs, rhs, N_DOF*SizeOfDouble); 
      
      // unless the stiffness matrix or rhs change in time, it is enough to 
      // assemble only once at the begning
      if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
      {
        SystemMatrix.AssembleARhs(la_stiff_rhs, sol, rhs);
        
        // M:= M + (tau*THETA1)*A
        // rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs +[M-(tau*THETA2)A]*oldsol
        // note! sol contains only the previous time step value, so just pass 
        // sol for oldsol
        SystemMatrix.AssembleSystMat(oldrhs, sol, rhs, sol);
        ConvectionFirstTime = false;
      }
     
      // solve the system matrix 
      SystemMatrix.Solve(sol, rhs);
    
      // restore the mass matrix for the next time step    
      // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
      if(UpdateStiffnessMat || UpdateRhs)
      {
        SystemMatrix.RestoreMassMat();
      }   
    } // for(int l=0;l<N_SubSteps;l++) 
    //======================================================================
    // produce outout
    //======================================================================
    if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
    {
      if(TDatabase::ParamDB->WRITE_VTK)
      {
        std::string filename(TDatabase::ParamDB->OUTPUTDIR);
        filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
        if(img<10) filename += ".0000";
        else if(img<100) filename += ".000";
        else if(img<1000) filename += ".00";
        else if(img<10000) filename += ".0";
        else filename += ".";
        filename += std::to_string(img) + ".vtk";
        Output.WriteVtk(filename.c_str());
        img++;
      }
    }
    
    //======================================================================
    // measure errors to known solution
    //======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      TAuxParam2D aux;
      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, &aux, 1, &Scalar_FeSpace,
                                   errors);

      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);

      errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror = errors[0];
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

      errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
      olderror1 = errors[1];   
      
      
      if(Linfty<errors[0])
      Linfty=errors[0];

      OutPut( "Linfty " << Linfty << endl);      
    } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

  CloseFiles();
  return 0;
} // end main









