#include <NSE2D.h>
#include <MainUtilities.h> // GetVelocityAndPressureSpace
#include <Domain.h>
#include <Database.h>
#include <Output2D.h>
#include <LinAlg.h> // DDot
#include <ItMethod.h>
#include <MultiGridIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>

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
  
  // done with the conrtuctor in case we're not using multigrid
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE != 5 
    || TDatabase::ParamDB->SOLVER_TYPE != 1)
    return;
  // else multigrid
  
  // create spaces, functions, matrices on coarser levels
  double *param = new double[2];
  param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  this->multigrid=new TNSE_MultiGrid(1, 2, param);
  // number of refinement levels for the multigrid
  int LEVELS = TDatabase::ParamDB->LEVELS;
  if(LEVELS > domain->get_ref_level() + 1)
    LEVELS = domain->get_ref_level() + 1;
  
  this->u.resize(LEVELS,nullptr);
  this->u1.resize(LEVELS,nullptr);
  this->u2.resize(LEVELS,nullptr);
  this->p.resize(LEVELS,nullptr);
  this->rhs.resize(LEVELS,nullptr);
  this->matrix.resize(LEVELS,nullptr);
  
  // the matrix and rhs side on the finest grid are already constructed 
  // now construct all matrices, rhs, and solutions on coarser grids
  for(int i=0; i<LEVELS-1; i++)
  {
    unsigned int grid = i + domain->get_ref_level() + 1 - LEVELS;
    TCollection *coll = domain->GetCollection(It_EQ, grid);
    // index of the corresponding matrix, rhs, and solution in their respective
    // vectors
    unsigned int index = LEVELS - 1 - i;
    
    GetVelocityAndPressureSpace(coll, example->get_bc(0), NULL, Velocity_FeSpace,
                              Pressure_FeSpace, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
    n_u = Velocity_FeSpace->GetN_DegreesOfFreedom();
    n_p = Pressure_FeSpace->GetN_DegreesOfFreedom();
    n_dof = 2 * n_u + n_p; 
    this->rhs[index] = new double[n_dof];
    memset(this->rhs[index],0,n_dof*SizeOfDouble);
    sol = new double[n_dof];
    memset(sol, 0, n_dof*SizeOfDouble);
    
    // create the finite element functions
    this->u[index] = new TFEVectFunct2D(Velocity_FeSpace, (char*) "u", (char*) "u",
                                    sol, n_u, 2);
    this->p[index] = new TFEFunction2D(Pressure_FeSpace, (char*) "p", (char*) "p",
                                   sol + 2 * n_u, n_p);
    this->u1.at(index) = this->u[index]->GetComponent(0);
    this->u2.at(index) = this->u[index]->GetComponent(1);
    
    this->matrix.at(index) = new TSystemMatNSE2D(this->u[index], this->p[index]);
    this->matrix.at(index)->Init(example->get_bd(0), example->get_bd(1));
    
    this->multigrid->AddLevel(this->mg_levels(i, index));
  }
  this->multigrid->AddLevel(this->mg_levels(LEVELS-1, 0));
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
  for(unsigned int grid=0, n_grids=matrix.size(); grid<n_grids; ++grid)
  {
    TFEFunction2D *fe_functions[3] = { this->u1[grid], this->u2[grid], this->p[grid] };
    // create a local assembling objects which are needed to assemble the matrices
    LocalAssembling2D la(NSE2D_Galerkin, fe_functions,
                         this->example->get_coeffs());
    
    this->matrix[grid]->Assemble(la, this->u[grid]->GetValues(), this->rhs[grid]);
  }
}


void NSE2D::assemble_nonlinear_term()
{
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  for(unsigned int grid=0, n_grids=matrix.size(); grid<n_grids; ++grid)
  {
    TFEFunction2D *fe_functions[3] = { this->u1[grid], this->u2[grid], this->p[grid] };
    LocalAssembling2D la_nonlinear(NSE2D_Galerkin_Nonlinear, fe_functions,
                                   this->example->get_coeffs());
    
    this->matrix[grid]->AssembleNonLinear(la_nonlinear, this->u[grid]->GetValues(),
                                       this->rhs[grid]);
  }
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
  if((TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE !=5)
    || (TDatabase::ParamDB->SOLVER_TYPE != 1))
  {
    this->matrix[0]->Solve(this->get_solution(), this->rhs[0]);
  }
  else
  {
    this->mg_solver();
  }
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

TNSE_MGLevel* NSE2D:: mg_levels(int i, int index)
{
  TNSE_MGLevel *mg_l;
  int n_aux;
  double alpha[2];

  int v_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
  int p_space_code = TDatabase::ParamDB->PRESSURE_SPACE;  
  
  if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
        || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE))
     n_aux=4;
  else
     n_aux=2;
  
  if (i==0)
  {
    alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
    alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
  }
  else
  {
    alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
    alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
  }
  
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      OutPut("NSTYPE: " << TDatabase::ParamDB->NSTYPE << " does not supported" << endl);
      exit(-4771);
    break;
    case 2:
      mg_l = new TNSE_MGLevel2(i, matrix[index]->get_square_matrix(0), 
                     matrix[index]->get_rectangular_matrix(0), matrix[index]->get_rectangular_matrix(1), 
                     matrix[index]->get_rectangular_matrix(2), matrix[index]->get_rectangular_matrix(3),
                     this->rhs[index], this->u[index]->GetValues(), n_aux, alpha, 
                     v_space_code, p_space_code, NULL, NULL);
    break;
    case 3:
      OutPut("NSTYPE: " << TDatabase::ParamDB->NSTYPE << " does not supported" << endl);
      exit(-4771);
    break;
    case 4:
       mg_l = new TNSE_MGLevel4(i, matrix[index]->get_square_matrix(0), 
                     matrix[index]->get_square_matrix(1), matrix[index]->get_square_matrix(2), 
                     matrix[index]->get_square_matrix(3),
                     matrix[index]->get_rectangular_matrix(0), matrix[index]->get_rectangular_matrix(1), 
                     matrix[index]->get_rectangular_matrix(2), matrix[index]->get_rectangular_matrix(3),
                     this->rhs[index], this->u[index]->GetValues(), n_aux, alpha, 
                     v_space_code, p_space_code, NULL, NULL);
    break;
    case 14:
      OutPut("WARNING: NSTYPE 14 is not supported \n");
    break;
  }
  return mg_l;
}


void NSE2D :: mg_solver()
{
  double *itmethod_rhs, *itmethod_sol;
  TItMethod *itmethod, *prec;
  int zero_start;  
  TSquareMatrix2D *sqMat[5];
  TSquareMatrix **sqmatrices = (TSquareMatrix **)sqMat;
  TMatrix2D *recMat[4];
  TMatrix **matrices = (TMatrix **)recMat;
  MatVecProc *MatVect;
  DefectProc *Defect;

  int n_dof = this->get_size();
  int n_u_dof=this->u1[0]->GetLength();
  
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      sqMat[0]=this->matrix[0]->get_square_matrix(0);
      recMat[0]=this->matrix[0]->get_rectangular_matrix(0);
      recMat[1]=this->matrix[0]->get_rectangular_matrix(1);
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
      break;
    case 2:
      sqMat[0]=this->matrix[0]->get_square_matrix(0);
      
      recMat[0]=this->matrix[0]->get_rectangular_matrix(0);
      recMat[1]=this->matrix[0]->get_rectangular_matrix(1);
      recMat[2]=this->matrix[0]->get_rectangular_matrix(2);
      recMat[3]=this->matrix[0]->get_rectangular_matrix(3);
      
      MatVect = MatVect_NSE2;
      Defect = Defect_NSE2;
      break;
    case 3:
      sqMat[0]=this->matrix[0]->get_square_matrix(0);
      sqMat[1]=this->matrix[0]->get_square_matrix(1);
      sqMat[2]=this->matrix[0]->get_square_matrix(2);
      sqMat[3]=this->matrix[0]->get_square_matrix(3);
      recMat[0]=this->matrix[0]->get_rectangular_matrix(0);
      recMat[1]=this->matrix[0]->get_rectangular_matrix(1);
      MatVect = MatVect_NSE3;
      Defect = Defect_NSE3;
      break;
    case 4:
      sqMat[0]=this->matrix[0]->get_square_matrix(0);
      sqMat[1]=this->matrix[0]->get_square_matrix(1);
      sqMat[2]=this->matrix[0]->get_square_matrix(2);
      sqMat[3]=this->matrix[0]->get_square_matrix(3);
      recMat[0]=this->matrix[0]->get_rectangular_matrix(0);
      recMat[1]=this->matrix[0]->get_rectangular_matrix(1);
      recMat[2]=this->matrix[0]->get_rectangular_matrix(2);
      recMat[3]=this->matrix[0]->get_rectangular_matrix(3);
      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
      break;
    case 14:
      OutPut("WARNING: NSTYPE 14 is not fully supported, take NSTYPE 4\n");
      exit(-4711);
      break;
  }


  if(TDatabase::ParamDB->SOLVER_TYPE ==1)
  {
    switch(TDatabase::ParamDB->SC_SOLVER_SADDLE)
    {
      case 11:
        zero_start=1;
        break;
      case 16:
        zero_start=0;
        break;
    }
    switch(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE)
    {
      case 5:
        prec = new TMultiGridIte(MatVect, Defect, NULL, 0, 
                                 n_dof, this->multigrid, zero_start);
        break;
      default:
        OutPut("Unknown preconditioner !!!" << endl);
        exit(4711);
    }
    
    if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
    {
      itmethod_sol = new double[n_dof];
      itmethod_rhs = new double[n_dof];
      
      memcpy(itmethod_sol, this->u[0]->GetValues(), n_dof*SizeOfDouble);
      
      memcpy(itmethod_rhs, this->getRhs(), n_dof*SizeOfDouble);
    }
    else
    {
      itmethod_sol = this->u[0]->GetValues();
      itmethod_rhs = this->getRhs();
    }

    switch(TDatabase::ParamDB->SC_SOLVER_SADDLE)
    {
      case 11:
        itmethod = new TFixedPointIte(MatVect, Defect, prec, 0, n_dof, 0);
        break;
      case 16:
        itmethod = new TFgmresIte(MatVect, Defect, prec, 0, n_dof, 0);
        break;
      default:
        OutPut("Unknown solver !!!!! " << endl);
        exit(4711);
    }
  }
    
  switch(TDatabase::ParamDB->SOLVER_TYPE)
  {
    case 1:
      itmethod->Iterate(sqmatrices,matrices,itmethod_sol,itmethod_rhs);
      break;
    case 2:
      break;
  }
  
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
  {
    memcpy(this->get_solution(), itmethod_sol, n_dof*SizeOfDouble);
    
    memcpy(this->getRhs(), itmethod_rhs, n_dof*SizeOfDouble);
    
    delete itmethod; delete prec;
    delete [] itmethod_rhs;
    delete [] itmethod_sol;
  }
}
