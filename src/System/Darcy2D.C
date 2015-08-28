/** ************************************************************************ 
* @brief     source file for Darcy2D
* @author    Ulrich Wilbrandt,
* @date      15.03.15
 ************************************************************************  */
#include <Database.h>
#include <Darcy2DMixed.h>
#include <Darcy2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <AuxParam2D.h>
#include <MultiGridIte.h>
#include <MultiGrid2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <Output2D.h>
#include <MainUtilities.h>
#include <stdlib.h>
#include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>

Darcy2D::Darcy2D(TDomain *domain, Example_Darcy2D* e)
 : matrix(1, NULL), rhs(1, NULL), u(1, NULL), p(1, NULL),
   example(e != NULL ? e : new Example_Darcy2D()), multigrid(nullptr)
{
  this->set_parameters();
  // you have to set the two variables apropriatly to get an inf-sup stable 
  // pair. There is currently no function which checks this.
  int v_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
  int p_space_code = TDatabase::ParamDB->PRESSURE_SPACE;
  
  // a collection is basically only an array of cells, which is needed to create
  // a finite element space
  TCollection *coll = domain->GetCollection(It_Finest, 0);
  // print out some information about the mesh
  int N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
  // create the velocity and pressure space
  TFESpace2D * v_space = new TFESpace2D(coll, (char*)"u",
                                        (char*)"Darcy velocity",
                                        this->example->get_bc(0), v_space_code,
                                        NULL);
  TFESpace2D * p_space = new TFESpace2D(coll, (char*)"p",
                                        (char*)"Darcy pressure",
                                        this->example->get_bc(1), p_space_code,
                                        NULL);
  // both spaces are discontinuous
  v_space->SetAsDGSpace(); // this is used for output
  p_space->SetAsDGSpace(); // this is used for output
  
  // print out some information on the finite element space
  int n_u = v_space->GetN_DegreesOfFreedom();
  int n_u_active = v_space->GetActiveBound();
  int n_p = p_space->GetN_DegreesOfFreedom();
  int n_dof = n_u + n_p;
  
  OutPut(" dof velocity (vector-valued) : "<< setw(5) << n_u << endl);
  OutPut(" active dof velocity          : "<< setw(5) << n_u_active << endl);
  OutPut(" dof pressure                 : "<< setw(5) << n_p << endl);
  OutPut(" dof all                      : "<< setw(5) << n_dof << endl);
  
  this->rhs[0] = new double[n_dof];
  double* sol = new double[n_dof]; // access to solution through 'u'
  
  // set solution and right hand side vectors to zero
  memset(sol, 0, n_dof * SizeOfDouble);
  memset(this->rhs[0], 0, n_dof * SizeOfDouble);
  
  // create finite element functions
  this->u[0] = new TFEFunction2D(v_space, (char*)"u", (char*)"u", sol, n_u);
  this->p[0] = new TFEFunction2D(p_space, (char*)"u", (char*)"u", sol+n_u, n_p);
  
  this->matrix[0] = new BlockMatrixDarcy2D(v_space, p_space, example->get_bd());
}

Darcy2D::~Darcy2D()
{
  for(auto mat : this->matrix)
    delete mat;
  for(auto r : this->rhs)
    delete[] r;
  for(auto f : u)
  {
    delete [] f->GetValues();
    delete f;
  }
  for(auto f : p)
    delete f;
  delete multigrid;
}

void Darcy2D::set_parameters()
{
  if(TDatabase::ParamDB->PRESSURE_SPACE == -4711)
  {
    switch(TDatabase::ParamDB->VELOCITY_SPACE)
    {
      case 1000: // Raviart-Thomas, order 0
        TDatabase::ParamDB->PRESSURE_SPACE = 0;
        break;
      case 1001: // Raviart-Thomas, order 1
        TDatabase::ParamDB->PRESSURE_SPACE = -11;
        break;
      case 1002: // Raviart-Thomas, order 2
        TDatabase::ParamDB->PRESSURE_SPACE = -12;
        break;
      case 1003: // Raviart-Thomas, order 3
        TDatabase::ParamDB->PRESSURE_SPACE = -13;
        break;
      case 1011: // Brezzi-Douglas-Marini, order 1
        TDatabase::ParamDB->PRESSURE_SPACE = 0;
        break;
      case 1012: // Brezzi-Douglas-Marini, order 2
        TDatabase::ParamDB->PRESSURE_SPACE = -110;
        break;
      case 1013: // Brezzi-Douglas-Marini, order 3
        TDatabase::ParamDB->PRESSURE_SPACE = -120;
        break;
      default:
        ErrMsg("unknown velocity space for Darcy2D");
        throw(std::runtime_error("unknown velocity space for Darcy2D"));
        break;
    }
  }
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5 
    && TDatabase::ParamDB->SOLVER_TYPE == 1)
  {
     ErrMsg("multigrid not yet implemented for Darcy2D");
     throw(std::runtime_error("multigrid not yet implemented for Darcy2D"));
  }
}


void Darcy2D::assemble()
{
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  for(unsigned int grid = 0, n_grids = matrix.size(); grid < n_grids; ++grid)
  {
    TFEFunction2D *fe_functions[2] = { this->u[grid], this->p[grid] };
    // create a local assembling object which is needed to assemble the matrices
    LocalAssembling2D la(Darcy2D_Galerkin, fe_functions,
                         this->example->get_coeffs());
    
    this->matrix[grid]->Assemble(la, this->u[grid]->GetValues(), 
                                 this->rhs[grid]);
  }
  
  // copy Dirichlet values from rhs to solution vector (this is not really 
  // necessary in case of a direct solver)
  int n_u = this->u[0]->GetFESpace2D()->GetN_DegreesOfFreedom();
  int n_u_active = this->u[0]->GetFESpace2D()->GetActiveBound();
  double * sol = this->get_solution();
  memcpy(sol+n_u_active, rhs[0]+n_u_active, (n_u-n_u_active)*SizeOfDouble);
} // void Darcy2D::Assemble


void Darcy2D::solve()
{
   this->matrix[0]->Solve(this->get_solution(), this->rhs[0]);
   if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
   {
     // this sould be done in SystemMatDarcy2D::Solve, but we don't have 
     // access to p there
     //this->p[0]->project_into_L20();
     int n_u = this->u[0]->GetFESpace2D()->GetN_DegreesOfFreedom();
     int n_p = this->p[0]->GetFESpace2D()->GetN_DegreesOfFreedom();
     IntoL20Vector2D(this->get_solution()+n_u, n_p, 
                     TDatabase::ParamDB->PRESSURE_SPACE);
   }
}


void Darcy2D::output(int i)
{
  if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;
  
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    this->u[0]->PrintMinMax();
    this->p[0]->PrintMinMax();
  }
  
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(2, 2, 0, 0, NULL);
    Output.AddFEFunction(this->u[0]);
    Output.AddFEFunction(this->p[0]);
    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
    filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
    if(i >= 0)
      filename += "_" + std::to_string(i);
    filename += ".vtk";
    Output.WriteVtk(filename.c_str());
  }
  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    DoubleFunct2D *const *Exact = &(example->get_exact())[0];
    ErrorMethod2D *L2DivH1 = L2DivH1Errors;
    double errors[6];
    u[0]->GetErrorsForVectorValuedFunction(Exact, L2DivH1, errors);
    
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
    TFESpace2D * pointer_to_p_space = this->p[0]->GetFESpace2D();
    p[0]->GetErrors(example->get_exact(2), 3, AllDerivatives, 2, L2H1Errors,
                    example->get_coeffs(), &aux, 1, &pointer_to_p_space,
                    errors+3);

    OutPut(" L2(u):      " << errors[0] << endl);
    OutPut(" L2(div(u)): " << errors[1] << endl);
    OutPut(" H1-semi(u): " << errors[2] << endl);
    OutPut(" L2(p):      " << errors[3] << endl);
    OutPut(" H1-semi(p): " << errors[4] << endl);
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
}






