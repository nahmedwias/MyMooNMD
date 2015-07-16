// =======================================================================
//
// Purpose:  main program for solving a stationary scalar equation using ParMooN
//
// Author:   Sashikumaar Ganesan
//
// History:  Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <SystemCD2D.h>
#include <Output2D.h>
#include <MainUtilities.h>

#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_CD2D.h>

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

  //set PROBLEM_TYPE to CD if not yet set
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 1;
  //open OUTFILE, this is where all output is written to (addionally to console)
  OpenFiles();
 
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
    Domain.PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);
  
  // choose example according to the value of TDatabase::ParamDB->EXAMPLE
  Example_CD2D example;
   
  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  // a collection is basically only an array of cells, which is needed to create
  // a finite element space
  TCollection *coll = Domain.GetCollection(It_Finest, 0);
  // print out some information about the mesh
  int N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
  // create fespace for scalar equation
  TFESpace2D Scalar_FeSpace (coll, (char*)"name", (char*)"description", 
                             example.get_bc(0), ORDER, NULL);
  // print out some information on the finite element space
  int N_DOF = Scalar_FeSpace.GetN_DegreesOfFreedom();
  OutPut("dof all      : "<< setw(10) << N_DOF  << endl);
 
  //======================================================================
  // construct all finite element functions
  //======================================================================
  double *sol = new double[N_DOF];
  double *rhs = new double[N_DOF];
  // set solution and right hand side vectors to zero
  memset(sol, 0, N_DOF*SizeOfDouble);
  memset(rhs, 0, N_DOF*SizeOfDouble);

  // create a finite element function
  TFEFunction2D Scalar_FeFunction(&Scalar_FeSpace, (char*)"C", (char*)"C", sol,
                                  N_DOF);
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  TFEFunction2D *fe_function[1] = { &Scalar_FeFunction };
  
  //======================================================================
  // SystemMatrix construction and solution
  //====================================================================== 
  TSystemMatScalar2D SystemMatrix(&Scalar_FeSpace);
  
  // initilize the system matrix with the functions defined in the example
  SystemMatrix.Init(example.get_bc(0), example.get_bd(0));
  
  // create a local assembling object which is needed to assemble the matrix
  LocalAssembling2D la(CD2D_SUPG, fe_function, example.get_coeffs());
       
  // assemble the system matrix with given local assembling, sol and rhs 
  SystemMatrix.Assemble(la, sol, rhs);
  
  //Solve the system
  {
    double t1 = GetTime();
    SystemMatrix.Solve(sol, rhs);
    double t2 = GetTime();
    OutPut( "time for solving: " << t2-t1 << endl);
  }
    
  //======================================================================
  // produce outout
  //======================================================================
  TOutput2D Output(1, 1, 0, 0, &Domain);
  Output.AddFEFunction(&Scalar_FeFunction);

  //Scalar_FeFunction.Interpolate(example.get_exact(0));
  // write solution to a vtk file
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
    filename += "/" + std::string(TDatabase::ParamDB->BASENAME) + ".vtk";
    Output.WriteVtk(filename.c_str());
  }
   
  //====================================================================== 
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  //======================================================================    
  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    double errors[4];
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
    TFESpace2D *fespace[1] = { &Scalar_FeSpace };

    Scalar_FeFunction.GetErrors(example.get_exact(0), 3, AllDerivatives, 2, 
                                 L2H1Errors, example.get_coeffs(), &aux, 1,
                                 fespace, errors);

    OutPut( "L2: " << errors[0] << endl);
    OutPut( "H1-semi: " << errors[1] << endl);
    OutPut( "SD: " << errors[2] << endl);
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)

  delete [] sol;
  delete [] rhs;
  delete coll;
  CloseFiles();
  return 0;
} // end main
