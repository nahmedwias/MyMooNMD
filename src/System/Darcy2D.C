/** ************************************************************************ 
* @brief     source file for Darcy2D
* @author    Ulrich Wilbrandt,
* @date      15.03.15
 ************************************************************************  */
#include <Database.h>
#include <Darcy2D.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <MainUtilities.h>

/** ************************************************************************ */
Darcy2D::System_per_grid::System_per_grid(const Example_Darcy2D& example,
                                          TCollection& coll)
 : velocity_space(&coll, (char*)"u", (char*)"Darcy velocity", example.get_bc(0),
                  TDatabase::ParamDB->VELOCITY_SPACE, nullptr),
   pressure_space(&coll, (char*)"p", (char*)"Darcy pressure", example.get_bc(1),
                  TDatabase::ParamDB->PRESSURE_SPACE, nullptr),
   matrix(this->velocity_space, this->pressure_space, example.get_bd()),
   rhs(this->matrix, true),
   solution(this->matrix, false),
   u(&this->velocity_space, (char*)"u", (char*)"u", this->solution.block(0),
     this->solution.length(0)),
   p(&this->pressure_space, (char*)"p", (char*)"p", this->solution.block(1),
     this->solution.length(1))
{
  velocity_space.SetAsDGSpace();
  pressure_space.SetAsDGSpace();
}

/** ************************************************************************ */
Darcy2D::Darcy2D(const TDomain& domain, int reference_id)
 : Darcy2D(domain, *(new Example_Darcy2D()), reference_id)
{
  // note that the way we construct the example above will produce a memory 
  // leak, but that class is small.
}

/** ************************************************************************ */
Darcy2D::Darcy2D(const TDomain& domain, const Example_Darcy2D& ex, 
                 int reference_id)
 : systems(), example(ex), multigrid(nullptr)
{
  // make sure all parameters in the database are set consistently
  this->set_parameters();
  
  // a collection is basically only an array of cells, which is needed to create
  // a finite element space
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  
  // create finite element spaces and functions, a matrix, rhs, and solution
  this->systems.emplace_back(example, *coll);
  
  // print out some information on the finite element space
  const TFESpace2D& v_space = this->systems.front().velocity_space;
  const TFESpace2D& p_space = this->systems.front().pressure_space;
  int n_u = v_space.GetN_DegreesOfFreedom();
  int n_u_active = v_space.GetActiveBound();
  int n_p = p_space.GetN_DegreesOfFreedom();
  int n_dof = n_u + n_p;
  
  OutPut(" dof velocity (vector-valued) : "<< setw(5) << n_u << endl);
  OutPut(" active dof velocity          : "<< setw(5) << n_u_active << endl);
  OutPut(" dof pressure                 : "<< setw(5) << n_p << endl);
  OutPut(" dof all                      : "<< setw(5) << n_dof << endl);
}

/** ************************************************************************ */
Darcy2D::~Darcy2D()
{
  // delete the collections created during the contructor
  for(auto & s : this->systems)
    delete s.velocity_space.GetCollection();
}

/** ************************************************************************ */
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

/** ************************************************************************ */
void Darcy2D::assemble()
{
  for(System_per_grid& s : this->systems)
  {
    // the class LocalAssembling2D which we will need next, requires an array of
    // pointers to finite element functions, i.e. TFEFunction2D **.
    TFEFunction2D *fe_functions[2] = { &s.u, &s.p };
    // create a local assembling object which is needed to assemble the matrices
    LocalAssembling2D la(Darcy2D_Galerkin, fe_functions,
                         this->example.get_coeffs());
    s.matrix.Assemble(la, s.rhs);
  }
  // copy Dirichlet values from rhs to solution vector (this is not really 
  // necessary in case of a direct solver)
  this->systems.front().solution.copy_nonactive(this->systems.front().rhs);
} // void Darcy2D::Assemble

/** ************************************************************************ */
void Darcy2D::solve()
{
  System_per_grid & s = this->systems.front();
  s.matrix.Solve(s.solution, s.rhs);
  
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    s.p.project_into_L20();
}

/** ************************************************************************ */
void Darcy2D::output(int i)
{
  if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;
  
  System_per_grid & s = this->systems.front();
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    s.u.PrintMinMax();
    s.p.PrintMinMax();
  }
  
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(2, 2, 0, 0, nullptr);
    Output.AddFEFunction(&s.u);
    Output.AddFEFunction(&s.p);
    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
    filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
    if(i >= 0)
      filename += "_" + std::to_string(i);
    filename += ".vtk";
    Output.WriteVtk(filename.c_str());
  }
  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    DoubleFunct2D *const *Exact = &(example.get_exact())[0];
    ErrorMethod2D *L2DivH1 = L2DivH1Errors;
    double errors[6];
    s.u.GetErrorsForVectorValuedFunction(Exact, L2DivH1, errors);
    
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
    const TFESpace2D * pointer_to_p_space = &s.pressure_space;
    s.p.GetErrors(example.get_exact(2), 3, AllDerivatives, 2, L2H1Errors,
                  example.get_coeffs(), &aux, 1, &pointer_to_p_space,
                  errors + 3);

    OutPut(" L2(u):      " << errors[0] << endl);
    OutPut(" L2(div(u)): " << errors[1] << endl);
    OutPut(" H1-semi(u): " << errors[2] << endl);
    OutPut(" L2(p):      " << errors[3] << endl);
    OutPut(" H1-semi(p): " << errors[4] << endl);
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
}

/** ************************************************************************ */




