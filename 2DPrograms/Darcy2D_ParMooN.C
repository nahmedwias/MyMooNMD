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
#include <SystemDarcy2D.h>
#include <Output2D.h>
#include <MainUtilities.h>
#include <LinAlg.h>

#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <LocalAssembling2D.h>
#include <Example_Darcy2D.h>

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
  Example_Darcy2D example;
   
  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  // you have to set the two variables apropriatly to get an inf-sup stable 
  // pair. There is currently no function which checks this.
  int v_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
  int p_space_code = TDatabase::ParamDB->PRESSURE_SPACE;
  
  // a collection is basically only an array of cells, which is needed to create
  // a finite element space
  TCollection *coll = Domain.GetCollection(It_Finest, 0);
  // print out some information about the mesh
  int N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
  // create the velocity and pressure space
  TFESpace2D v_space(coll, (char*)"u", (char*)"Darcy velocity", 
                     example.get_bc(0), v_space_code, NULL);
  TFESpace2D p_space(coll, (char*)"p", (char*)"Darcy pressure", 
                     example.get_bc(1), p_space_code, NULL);
  // both spaces are discontinuous
  v_space.SetAsDGSpace(); // this is used for output
  p_space.SetAsDGSpace(); // this is used for output
  TFESpace2D *fespaces[2] = { &v_space, &p_space };
  
  // print out some information on the finite element space
  int N_U = v_space.GetN_DegreesOfFreedom();
  int N_U_Active = v_space.GetActiveBound();
  int N_P = p_space.GetN_DegreesOfFreedom();
  int N_DOF = N_U + N_P;
  
  OutPut(" dof velocity (vector-valued) : "<< setw(5) << N_U << endl);
  OutPut(" active dof velocity          : "<< setw(5) << N_U_Active << endl);
  OutPut(" dof pressure                 : "<< setw(5) << N_P << endl);
  OutPut(" dof all                      : "<< setw(5) << N_DOF << endl);
  
  //======================================================================
  // construct all finite element functions
  //======================================================================
  double *sol = new double[N_DOF];
  double *rhs = new double[N_DOF];
  // set solution and right hand side vectors to zero
  memset(sol, 0, N_DOF*SizeOfDouble);
  memset(rhs, 0, N_DOF*SizeOfDouble);

  // create finite element functions
  TFEFunction2D u(&v_space, (char*)"u", (char*)"Darcy velo", sol,N_U);
  TFEFunction2D p(&p_space, (char*)"p", (char*)"Darcy press", sol+N_U, N_P);
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  TFEFunction2D *fe_functions[2] = { &u, &p };
  
  //======================================================================
  // SystemMatrix construction and solution
  //====================================================================== 
  // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) GLS (or) SUPG (or) LOCAL_PROJECTION
  // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
  TSystemDarcy2D SystemMatrix(fespaces);
  
  // initilize the system matrix with the functions defined in the example
  SystemMatrix.Init(example.get_bc(), example.get_bd());
  
  // create a local assembling object which is needed to assemble the matrix
  LocalAssembling2D la(Darcy2D_Galerkin, fe_functions, example.get_coeffs());
       
  // assemble the system matrix with given local assembling, sol and rhs 
  SystemMatrix.Assemble(la, sol, rhs);
  
  //Solve the system
  {
    double t1 = GetTime();
    SystemMatrix.Solve(sol, rhs);
     if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       // this sould be done in TSystemDarcy2D::Solve, but we don't have 
       // access to p there
       //p.project_into_L20();
       IntoL20Vector2D(sol+N_U, N_P, p_space_code);
    double t2 = GetTime();
    OutPut( "time for solving: " << t2-t1 << endl);
  }
  /*
  for(int i = 0; i < N_DOF; ++i)
    OutPut("sol(" << i+1 << ")= " << sol[i] << ";\n");
  */
  
  //======================================================================
  // produce outout
  //======================================================================
  TOutput2D Output(1, 1, 0, 0, &Domain);
  Output.AddFEFunction(&u);
  Output.AddFEFunction(&p);

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
    DoubleFunct2D *const *Exact = &(example.get_exact())[0];
    ErrorMethod2D *L2DivH1 = L2DivH1Errors;
    double errors[6];
    u.GetErrorsForVectorValuedFunction(Exact, L2DivH1, errors);
    
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
    p.GetErrors(example.get_exact(2), 3, AllDerivatives, 2, L2H1Errors,
                example.get_coeffs(), &aux, 1, fespaces+1, errors+3);

    OutPut(" L2(u):      " << errors[0] << endl);
    OutPut(" L2(div(u)): " << errors[1] << endl);
    OutPut(" H1-semi(u): " << errors[2] << endl);
    OutPut(" L2(p):      " << errors[3] << endl);
    OutPut(" H1-semi(p): " << errors[4] << endl);
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)

  delete [] sol;
  delete [] rhs;
  delete coll;
  CloseFiles();
  return 0;
} // end main
