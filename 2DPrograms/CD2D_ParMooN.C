// =======================================================================
//
// Purpose:     main program for solving a stationary scalar equation using ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemMatScalar2D.h>

#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>

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
// #include "../Examples/CD_2D/Hemker1996.h" // circle in a channel
#include "../Examples/CD_2D/SineLaplace.h" // smooth sol in unitsquares
// #include "../Examples/CD_2D/TwoInteriorLayers.h" // smooth sol in unitsquares
// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, N_Cells, ORDER, N_DOF,  img=1;
  
  double *sol, *rhs, t1, t2, errors[4];
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll;
  TFESpace2D *Scalar_FeSpace, *fesp[1];
  TFEFunction2D *Scalar_FeFunction;
  TOutput2D *Output;
  TSystemMatScalar2D *SystemMatrix;
  TAuxParam2D *aux;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
   
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
 
  Database->WriteParamDB(argv[0]);
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
// construct all finite element spaces
//=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
  // fespaces for scalar equation 
   Scalar_FeSpace  =  new TFESpace2D(coll, Name, Description, BoundCondition, ORDER, NULL);
   
   N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
   OutPut("dof all      : "<< setw(10) << N_DOF  << endl);
 
//======================================================================
// construct all finite element functions
//======================================================================
    sol = new double[N_DOF];
    rhs = new double[N_DOF];

    memset(sol, 0, N_DOF*SizeOfDouble);
    memset(rhs, 0, N_DOF*SizeOfDouble);

    Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, CString, CString, sol, N_DOF);   
    
//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) GLS (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
    SystemMatrix = new TSystemMatScalar2D(Scalar_FeSpace, LOCAL_PROJECTION, DIRECT);
    
    // initilize the system matrix with the functions defined in Example file
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);
       
    // assemble the system matrix with given aux, sol and rhs 
    // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass with NULL 
    SystemMatrix->Assemble(NULL, sol, rhs);
    
    //Solve the system
    t1 = GetTime();
    SystemMatrix->Solve(sol, rhs);
    t2 = GetTime();
    OutPut( "time for solving: " << t2-t1 << endl);
    
//======================================================================
// produce outout
//======================================================================
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput2D(2, 2, 1, 1, Domain);

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
     
//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] = Scalar_FeSpace;
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);

      delete aux;

      OutPut( "L2: " << errors[0] << endl);
      OutPut( "H1-semi: " << errors[1] << endl);
      OutPut( "SD: " << errors[2] << endl);

     } // if(TDatabase::ParamDB->MEASURE_ERRORS)


  CloseFiles();
  
  return 0;
} // end main
