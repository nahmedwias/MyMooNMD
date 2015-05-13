#include <CD2D.h>
#include <LocalAssembling2D.h>
#include <Domain.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <Output2D.h>
#include <MultiGrid2D.h>
#include <MainUtilities.h> // L2H1Errors


CD2D::CD2D(TDomain *domain, const Example_CD2D* e)
 : matrix(1, NULL), rhs(1, NULL), function(1, NULL), 
   example(e != NULL ? e : new Example_CD2D()), multigrid(NULL)
{
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain->GetCollection(It_Finest, 0);
  
  // create finite elememt space, access through 'function'
  int ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  TFESpace2D* space = new TFESpace2D(coll, (char*)"scalar_space", 
                                     (char*)"description",  
                                     this->example->get_bc(0), ORDER, NULL);
  int n_dof = space->GetN_DegreesOfFreedom();
  
  // create right hand side and solution
  this->rhs[0] = new double[n_dof];
  double* sol = new double[n_dof]; // access to solution through 'function'
  // set solution and right hand side vectors to zero
  memset(sol, 0, n_dof*SizeOfDouble);
  memset(rhs[0], 0, n_dof*SizeOfDouble);
  
  this->function[0] = new TFEFunction2D(space, (char*)"solution",
                                        (char*)"solution", sol, n_dof);
  
  this->matrix[0] = new TSystemMatScalar2D(space);
  this->matrix[0]->Init(this->example->get_bc(0), this->example->get_bd(0));
  
  // print out some information
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  OutPut("N_Cells    : " << setw(12) << coll->GetN_Cells() << endl);
  OutPut("h (min,max): " << setw(12) << h_min << " " << setw(12) <<h_max<<endl);
  OutPut("dof all    : " << setw(12) << n_dof << endl);
  OutPut("dof active : " << setw(12) << space->GetN_ActiveDegrees() << endl);
}


CD2D::~CD2D()
{
  // delete matrix
  for(auto mat : this->matrix)
    delete mat;
  for(auto r : this->rhs)
    delete [] r;
  for(auto f : this->function)
    delete f;
  delete multigrid;
}


void CD2D::assemble()
{
  LocalAssembling2D_type t = CD2D_Galerkin;
  switch(TDatabase::ParamDB->DISCTYPE)
  {
    case 1: t = CD2D_Galerkin; break;
    case 2: t = CD2D_SUPG; break;
    case 6: t = CD2D_GLS; break;
    default:
      ErrMsg("currently DISCTYPE " << TDatabase::ParamDB->DISCTYPE << " is not "
             << "supported by the class CD2D");
      throw("unsupported DISCTYPE");
  }
  // create a local assembling object which is needed to assemble the matrix
  LocalAssembling2D la(t, &(this->function[0]), this->example->get_coeffs());
       
  // assemble the system matrix with given local assembling, solution and rhs 
  this->matrix[0]->Assemble(la, this->function[0]->GetValues(), rhs[0]);
}


void CD2D::solve()
{
  this->matrix[0]->Solve(this->function[0]->GetValues(), rhs[0]);
}


void CD2D::output(int i)
{
  if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;
  
  // print the value of the largest and smallest entry in the finite element 
  // vector
  this->function[0]->PrintMinMax();
  
  // write solution to a vtk file
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(1, 1, 0, 0, NULL);
    Output.AddFEFunction(this->function[0]);
    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
    filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
    if(i >= 0) filename += "_" + std::to_string(i);
    filename += ".vtk";
    Output.WriteVtk(filename.c_str());
  }
  
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    double errors[4];
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
    TFESpace2D* space = this->function[0]->GetFESpace2D();
    
    this->function[0]->GetErrors(this->example->get_exact(0), 3, AllDerivatives,
                                 2, L2H1Errors, this->example->get_coeffs(),
                                 &aux, 1, &space, errors);

    OutPut( "L2: " << errors[0] << endl);
    OutPut( "H1-semi: " << errors[1] << endl);
    OutPut( "SD: " << errors[2] << endl);
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
}
