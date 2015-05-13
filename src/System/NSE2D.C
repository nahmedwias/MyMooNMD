#include <NSE2D.h>
#include <MainUtilities.h> // GetVelocityAndPressureSpace
#include <Domain.h>
#include <Database.h>
#include <Output2D.h>
#include <LinAlg.h> // DDot

NSE2D::NSE2D(TDomain *domain, const Example_NSE2D* e)
    : matrix(1, NULL), rhs(1, NULL), u(1, NULL), p(1, NULL), u1(1, NULL),
      u2(1, NULL), example(e != NULL ? e : new Example_NSE2D()), 
      multigrid(NULL), defect(), norms_of_residuals(10, 1e10),
      initial_residual(1e10)
{
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain->GetCollection(It_Finest, 0);
  
  // create finite element spaces for velocity and pressure
  TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace;
  int pressure_space_code;
  // the following function will create an inf-sup stable pair of finite element
  // spaces, if PRESSURE_SPACE is not set (-4711)
  GetVelocityAndPressureSpace(coll, example->get_bc(0), NULL, Velocity_FeSpace,
                              Pressure_FeSpace, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  
  // defaulty inf-sup pressure space will be selected based on the velocity 
  // space, so update it in database
  TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
  
  int n_u = Velocity_FeSpace->GetN_DegreesOfFreedom();
  int n_p = Pressure_FeSpace->GetN_DegreesOfFreedom();
  int n_dof = 2 * n_u + n_p; // total number of degrees of freedom
  
  // create right hand side and solution
  this->rhs[0] = new double[n_dof];
  double* sol = new double[n_dof]; // access to solution through 'function'
  defect.resize(n_dof);
  
  // set solution and right hand side vectors to zero
  memset(sol, 0, n_dof * SizeOfDouble);
  memset(this->rhs[0], 0, n_dof * SizeOfDouble);
  
  // create the finite element functions
  this->u[0] = new TFEVectFunct2D(Velocity_FeSpace, (char*) "u", (char*) "u",
                                  sol, n_u, 2);
  this->p[0] = new TFEFunction2D(Pressure_FeSpace, (char*) "p", (char*) "p",
                                 sol + 2 * n_u, n_p);
  this->u1[0] = this->u[0]->GetComponent(0);
  this->u2[0] = this->u[0]->GetComponent(1);
  
  this->matrix[0] = new TSystemMatNSE2D(this->u[0], this->p[0]);
  this->matrix[0]->Init(example->get_bd(0), example->get_bd(1));
  
  // print out some information
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  OutPut("N_Cells     : " << setw(10) << coll->GetN_Cells() << endl);
  OutPut("h (min,max) : " << setw(10) << h_min << " " << setw(12) << h_max
         <<endl);
  OutPut("dof Velocity: " << setw(10) << 2* n_u << endl);
  OutPut("dof Pressure: " << setw(10) << n_p << endl);
  OutPut("dof all     : " << setw(10) << n_dof << endl);
}


NSE2D::~NSE2D()
{
  // delete matrix
  for(auto mat : this->matrix)
    delete mat;
  for(auto r : this->rhs)
    delete[] r;
  for(auto f : u)
  {
    delete [] f->GetValues();
    delete f;
  }
  for(auto f : u1)
    delete f;
  for(auto f : u2)
    delete f;
  for(auto f : p)
    delete f;
  delete multigrid;
}


void NSE2D::assemble()
{
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  TFEFunction2D *fe_functions[3] = { u1[0], u2[0], p[0] };
  // create a local assembling objects which are needed to assemble the matrices
  LocalAssembling2D la(NSE2D_Galerkin, fe_functions,
                       this->example->get_coeffs());
  
  this->matrix[0]->Assemble(la, this->get_solution(), this->rhs[0]);
}


void NSE2D::assemble_nonlinear_term()
{
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  TFEFunction2D *fe_functions[3] = { u1[0], u2[0], p[0] };
  LocalAssembling2D la_nonlinear(NSE2D_Galerkin_Nonlinear, fe_functions,
                                 this->example->get_coeffs());
  
  this->matrix[0]->AssembleNonLinear(la_nonlinear, this->get_solution(),
                                     this->rhs[0]);
}


bool NSE2D::stopIt(unsigned int iteration_counter)
{
  // the array norms_of_residuals has length 10, At position 'last_digit_ite' 
  // of this array is the norm of the residual from 10 iterations ago. It is
  // compared to the current norm of the residual and then replaced by it.
  const unsigned int last_digit_ite = iteration_counter%10;
  const double normOfResidual = this->normOfResidual();
  if(iteration_counter == 0)
    initial_residual = normOfResidual;
  const double oldNormOfResidual = this->norms_of_residuals[last_digit_ite];
  
  const unsigned int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  const double convergence_speed = TDatabase::ParamDB->SC_NONLIN_DIV_FACTOR;
  bool slow_conv = false;
  
  if(normOfResidual >= convergence_speed*oldNormOfResidual)
    slow_conv = true;
  this->norms_of_residuals[last_digit_ite] = normOfResidual;
  
  double limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  if (TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE)
  {
    limit *= sqrt(this->get_size());
    OutPut("stopping tolerance for nonlinear iteration " << limit << endl);
  }
  // check if the iteration has converged, or reached the maximum number of
  // iterations or if convergence is too slow. Then return true otherwise false
  if( (normOfResidual<=limit) || (iteration_counter==Max_It) || (slow_conv) )
  {
    if(slow_conv)
      OutPut(" SLOW !!! " << normOfResidual/oldNormOfResidual << endl);
    // stop iteration
    OutPut(" ITE : " << setw(4) << iteration_counter <<
           setprecision(8) << " RES : " << normOfResidual << 
           " Reduction : " << normOfResidual/initial_residual << endl);
    return true;
  }
  else
    return false;
}


double NSE2D::normOfResidual()
{
  int n_u_dof = this->u1[0]->GetLength();
  int n_p_dof = this->p[0]->GetLength();
  
  memset(&defect[0],0,(2*n_u_dof + n_p_dof)*SizeOfDouble);
  
  this->matrix[0]->GetResidual(this->get_solution(), this->rhs[0], 
                               &this->defect[0]);
  
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20Vector2D(&defect[2*n_u_dof], n_p_dof, 
                    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
  
  // square of the norms of the residual components
  double impuls_Residual = Ddot(2*n_u_dof, &this->defect[0],&this->defect[0]);
  double mass_residual = Ddot(n_p_dof, &this->defect[2*n_u_dof],
                              &this->defect[2*n_u_dof]);
  // the full residual
  double full_residual = sqrt(impuls_Residual + mass_residual);
  
  OutPut(setw(14) << sqrt(impuls_Residual) << "\t" << setw(14) 
         << sqrt(mass_residual) << "\t" << setw(14) << full_residual << endl);
  return full_residual;
}


void NSE2D::solve()
{
  this->matrix[0]->Solve(this->get_solution(), this->rhs[0]);
}


void NSE2D::output(int i)
{
  if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;
  
  // print the value of the largest and smallest entry in the finite element 
  // vector
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    this->u1[0]->PrintMinMax();
    this->u2[0]->PrintMinMax();
    this->p[0]->PrintMinMax();
  }
  
  // write solution to a vtk file
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(2, 3, 1, 0, NULL);
    Output.AddFEFunction(this->u1[0]);
    Output.AddFEFunction(this->u2[0]);
    Output.AddFEFunction(this->p[0]);
    Output.AddFEVectFunct(this->u[0]);
    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
    filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
    if(i >= 0)
      filename += "_" + std::to_string(i);
    filename += ".vtk";
    Output.WriteVtk(filename.c_str());
  }
  
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    double err[4];
    TAuxParam2D NSEaux_error;
    MultiIndex2D NSAllDerivatives[3] = {D00, D10, D01};
    TFESpace2D *velocity_space = this->u[0]->GetFESpace2D();
    TFESpace2D *pressure_space = this->p[0]->GetFESpace2D();
    
    // errors in first velocity component
    this->u1[0]->GetErrors(example->get_exact(0), 3, NSAllDerivatives, 2,
                           L2H1Errors, NULL, &NSEaux_error, 1, &velocity_space,
                           err);
    // errors in second velocity component
    this->u2[0]->GetErrors(example->get_exact(1), 3, NSAllDerivatives, 2,
                           L2H1Errors, NULL, &NSEaux_error, 1, &velocity_space,
                           err + 2);
    OutPut("L2(u)     : " << sqrt(err[0]*err[0] + err[2]*err[2]) << endl);
    OutPut("H1-semi(u): " << sqrt(err[1]*err[1] + err[3]*err[3]) << endl);
    
    // errors in pressure
    this->p[0]->GetErrors(example->get_exact(2), 3, NSAllDerivatives, 2,
                          L2H1Errors, NULL, &NSEaux_error, 1, &pressure_space,
                          err);
    OutPut("L2(p)     : " << err[0] << endl);
    OutPut("H1-semi(p): " << err[1] << endl); 
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
}

