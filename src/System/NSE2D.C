#include <NSE2D.h>
#include <MainUtilities.h> // GetVelocityAndPressureSpace
#include <Database.h>
#include <Output2D.h>
#include <LinAlg.h> // DDot
#include <MultiGridIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <DirectSolver.h>


/** ************************************************************************ */
NSE2D::System_per_grid::System_per_grid (const Example_NSE2D& example,
               TCollection& coll, std::pair<int,int> velocity_pressure_orders)
 : velocity_space(&coll, (char*)"u", (char*)"Darcy velocity", example.get_bc(0),
                  velocity_pressure_orders.first, nullptr),
   pressure_space(&coll, (char*)"p", (char*)"Darcy pressure", example.get_bc(2),
                  velocity_pressure_orders.second, nullptr),
   matrix(this->velocity_space, this->pressure_space, example.get_bd()),
   rhs(this->matrix, true),
   solution(this->matrix, false),
   u(&this->velocity_space, (char*)"u", (char*)"u", this->solution.block(0),
     this->solution.length(0), 2),
   p(&this->pressure_space, (char*)"p", (char*)"p", this->solution.block(2),
     this->solution.length(2))
{

}

/** ************************************************************************ */
NSE2D::NSE2D(const TDomain& domain, int reference_id)
 : NSE2D(domain, *(new Example_NSE2D()), reference_id)
{
  // note that the way we construct the example above will produce a memory 
  // leak, but that class is small.
}

/** ************************************************************************ */
NSE2D::NSE2D(const TDomain & domain, const Example_NSE2D & e,
             unsigned int reference_id)
    : systems(), example(e), multigrid(), defect(), oldResiduals(),
      initial_residual(1e10)
{
  std::pair <int,int> 
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and preesure spaces
  // this function returns a pair which consists of 
  // velocity and pressure order
  this->get_velocity_pressure_orders(velocity_pressure_orders);
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  
  this->systems.emplace_back(example, *coll, velocity_pressure_orders);
  
  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);
  
  // print out some information  
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_u_active = this->get_velocity_space().GetN_ActiveDegrees();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 2 * n_u + n_p; // total number of degrees of freedom
  
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  Output::print<1>("N_Cells            : ", setw(10), coll->GetN_Cells());
  Output::print<1>("h (min,max)        : ", setw(10), h_min, " ", setw(12),
                   h_max);
  Output::print<1>("dof velocity       : ", setw(10), 2* n_u);
  Output::print<1>("dof velocity active: ", setw(10), 2* n_u_active);
  Output::print<1>("dof pressure       : ", setw(10), n_p);
  Output::print<1>("dof all            : ", setw(10), n_dof);
  
  // done with the conrtuctor in case we're not using multigrid
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE != 5 
    || TDatabase::ParamDB->SOLVER_TYPE != 1)
    return;
  // else multigrid
  
  // create spaces, functions, matrices on coarser levels
  double *param = new double[2];
  param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  this->multigrid.reset(new TNSE_MultiGrid(1, 2, param));
  // number of refinement levels for the multigrid
  int LEVELS = TDatabase::ParamDB->LEVELS;
  if(LEVELS > domain.get_ref_level() + 1)
    LEVELS = domain.get_ref_level() + 1;
  
  // the matrix and rhs side on the finest grid are already constructed 
  // now construct all matrices, rhs, and solutions on coarser grids
  for(int i = LEVELS - 2; i >= 0; i--)
  {
    unsigned int grid = i + domain.get_ref_level() + 1 - LEVELS;
    TCollection *coll = domain.GetCollection(It_EQ, grid, reference_id);
    this->systems.emplace_back(example, *coll, velocity_pressure_orders);
  }
  
  // create multigrid-level-objects, must be coarsest first
  unsigned int i = 0;
  for(auto it = this->systems.rbegin(); it != this->systems.rend(); ++it)
  {
    this->multigrid->AddLevel(this->mg_levels(i, *it));
    i++;
  }
}

/** ************************************************************************ */
NSE2D::~NSE2D()
{
}

/** ************************************************************************ */
void NSE2D::get_velocity_pressure_orders(std::pair <int,int> 
                 &velocity_pressure_orders)
{
  int velocity_order = velocity_pressure_orders.first;
  int pressure_order = velocity_pressure_orders.second;
  int order = 0;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5:
    case 12: case 13: case 14: case 15:
      if(velocity_order > 10)
        order = velocity_order-10;
      else
        order = velocity_order;
      break;
    case -1: case -2: case -3: case -4: case -5:
    case -101:
      order = velocity_order;
      break;
    // conforming fe spaces with bubbles on triangles 
    case 22: case 23: case 24:
      order = velocity_order;
      break;
      // discontinuous spaces 
    case -11: case -12: case -13:
      order = velocity_order*10;
      break;
  }
  TDatabase::ParamDB->VELOCITY_SPACE = order;
  velocity_pressure_orders.first = order;
  switch(pressure_order)
  {
    case -4711:
      switch(velocity_order)
      {
        case -1:
        case -2:
        case -3:
        case -4:
          // nonconforming pw (bi)linear velo/ pw constant pressure
          // conforming pw (bi)linear velo/ pw constant pressure (not stable !!!)
          pressure_order = -velocity_order-1;
          break; 
        case 1: // discontinuous space 
          pressure_order = 0;
          break;
        case 2: case 3: case 4: case 5:
        // standard conforming velo and continuous pressure
          pressure_order = velocity_order-1;
          break;
          // discontinuous pressure spaces 
          // standard conforming velo and discontinuous pressure
          // this is not stable on triangles !!!
        case 12: case 13: case 14: case 15:
          pressure_order = -(velocity_order-1)*10;
          break;
        case 22: case 23: case 24:
          pressure_order = -(velocity_order-11)*10;
          break;
      }
      break;
    // continuous pressure spaces
    case 1: case 2: case 3: case 4: case 5:
      pressure_order = 1;
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
  }
  TDatabase::ParamDB->PRESSURE_SPACE  = pressure_order;
  velocity_pressure_orders.second = pressure_order;
  
  Output::print("velocity space", setw(10), TDatabase::ParamDB->VELOCITY_SPACE);
  Output::print("pressure space", setw(10), TDatabase::ParamDB->PRESSURE_SPACE);
  
  // projection spaces for reconstructions
  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case 22:
      TDatabase::ParamDB->PROJECTION_SPACE = 1012;
      break;
    case 3:
      TDatabase::ParamDB->PROJECTION_SPACE = 1013;
      break;
    case 4:
      TDatabase::ParamDB->PROJECTION_SPACE = 1014;
      break;
  }
}

/** ************************************************************************ */
void NSE2D::assemble()
{
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  for(System_per_grid& s : this->systems)
  {
    TFEFunction2D *fe_functions[3] = 
      { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    LocalAssembling2D la(NSE2D_Galerkin, fe_functions,
                         this->example.get_coeffs());
    s.matrix.Assemble(la, s.rhs);
    // set rhs for Dirichlet nodes
    s.solution.copy_nonactive(s.rhs);
    delete fe_functions[0];
    delete fe_functions[1];
  }
}

/** ************************************************************************ */
void NSE2D::assemble_nonlinear_term()
{
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  for(System_per_grid& s : this->systems)
  {
    TFEFunction2D *fe_functions[3] = 
      { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    LocalAssembling2D la_nonlinear(NSE2D_Galerkin_Nonlinear, fe_functions,
                                   this->example.get_coeffs());
    s.matrix.AssembleNonLinear(la_nonlinear);
    // set rhs for Dirichlet nodes (is this necessary here?)
    s.solution.copy_nonactive(s.rhs);
  }
}

/** ************************************************************************ */
bool NSE2D::stopIt(unsigned int iteration_counter)
{
  // compute the residuals with the current matrix and solution
  this->normOfResidual();
  // the current norm of the residual
  const double normOfResidual = this->getFullResidual();
  // store initial residual, so later we can print the overall reduction
  if(iteration_counter == 0)
    initial_residual = normOfResidual;
  // the residual from 10 iterations ago
  const double oldNormOfResidual = this->oldResiduals.front().fullResidual;
  
  const unsigned int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  const double convergence_speed = TDatabase::ParamDB->SC_NONLIN_DIV_FACTOR;
  bool slow_conv = false;
  
  
  if(normOfResidual >= convergence_speed*oldNormOfResidual)
    slow_conv = true;
  
  double limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  if (TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE)
  {
    limit *= sqrt(this->get_size());
    Output::print<1>("stopping tolerance for nonlinear iteration ", limit);
  }

  // check if the iteration has converged, or reached the maximum number of
  // iterations or if convergence is too slow. Then return true otherwise false
  if( (normOfResidual<=limit) || (iteration_counter==Max_It) || (slow_conv) )
  {
    if(slow_conv)
      Output::print<1>(" SLOW !!! ", normOfResidual/oldNormOfResidual);
    // stop iteration
    Output::print<1>(" ITE : ", setw(4), iteration_counter, setprecision(8),
                     " RES : ", normOfResidual, " Reduction : ",
                     normOfResidual/initial_residual);
    return true;
  }
  else
    return false;
}

/** ************************************************************************ */
void NSE2D::normOfResidual()
{
  System_per_grid& s = this->systems.front();
  unsigned int n_u_dof = s.solution.length(0);
  unsigned int n_p_dof = s.solution.length(2);
  
  // copy rhs to defect
  this->defect = s.rhs;
  s.matrix.apply_scaled_add(s.solution.get_entries(), defect.get_entries(),-1.);
  
  //this->matrix[0]->GetResidual(this->get_solution(), this->rhs[0], 
  //                             &this->defect[0]);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    IntoL20FEFunction(&defect[2*n_u_dof], n_p_dof, &this->get_pressure_space(),
                      TDatabase::ParamDB->VELOCITY_SPACE, 
                      TDatabase::ParamDB->PRESSURE_SPACE);
  }
  
  // square of the norms of the residual components
  double impuls_Residual = Ddot(2*n_u_dof, &this->defect[0],&this->defect[0]);
  double mass_residual = Ddot(n_p_dof, &this->defect[2*n_u_dof],
                              &this->defect[2*n_u_dof]);
  
  Residuals currentResiduals(impuls_Residual, mass_residual);
  oldResiduals.add(currentResiduals);
}

/** ************************************************************************ */
void NSE2D::solve()
{
  System_per_grid& s = this->systems.front();
  if((TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE !=5)
    || (TDatabase::ParamDB->SOLVER_TYPE != 1))
  {
    s.matrix.Solve(s.solution.get_entries(), s.rhs.get_entries());
    /*
    std::shared_ptr<TMatrix> m = s.matrix.get_combined_matrix();
    TStructure*s = new TStructure(m->GetN_Rows(), m->GetN_Entries(),
                                              m->GetKCol(), m->GetRowPtr());
    s->Sort();
    TSquareMatrix * sm = new TSquareMatrix(s);
    sm->setEntries(m->GetEntries());
    
    DirectSolver(sm,  s.rhs.get_entries(), s.solution.get_entries());
    */
  }
  else
  {
    this->mg_solver();
  }
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    s.p.project_into_L20();
}

/** ************************************************************************ */
void NSE2D::output(int i)
{
  if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;
  
  System_per_grid& s = this->systems.front();
  TFEFunction2D* u1 = s.u.GetComponent(0);
  TFEFunction2D* u2 = s.u.GetComponent(1);
  
  // print the value of the largest and smallest entry in the finite element 
  // vector
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    s.p.PrintMinMax();
  }
  
  // write solution to a vtk file
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(2, 3, 1, 0, NULL);
    Output.AddFEFunction(&s.p);
    Output.AddFEVectFunct(&s.u);
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
    const TFESpace2D *velocity_space = &this->get_velocity_space();
    const TFESpace2D *pressure_space = &this->get_pressure_space();
    
    // errors in first velocity component
    u1->GetErrors(example.get_exact(0), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &velocity_space, err);
    // errors in second velocity component
    u2->GetErrors(example.get_exact(1), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &velocity_space, err + 2);
    Output::print<1>("L2(u)     : ", sqrt(err[0]*err[0] + err[2]*err[2]));
    Output::print<1>("H1-semi(u): ", sqrt(err[1]*err[1] + err[3]*err[3]));
    
    // errors in pressure
    s.p.GetErrors(example.get_exact(2), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &pressure_space, err);
    Output::print<1>("L2(p)     : ", err[0]);
    Output::print<1>("H1-semi(p): ", err[1]); 
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
  delete u1;
  delete u2;
}

/** ************************************************************************ */
TNSE_MGLevel* NSE2D::mg_levels(int i, System_per_grid& s)
{
  TNSE_MGLevel *mg_l = nullptr;
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
      ErrThrow("NSE2D::mg_levels: NSTYPE 1 is not supported");
      break;
    
    case 2:
      mg_l = new TNSE_MGLevel2(i, s.matrix.get_A_block(0), 
                               s.matrix.get_B_block(0), 
                               s.matrix.get_B_block(1), 
                               s.matrix.get_BT_block(0),
                               s.matrix.get_BT_block(1),
                               s.rhs.get_entries(), 
                               s.solution.get_entries(), 
                               n_aux, alpha, v_space_code, p_space_code, 
                               nullptr, nullptr);
      break;
      
    case 3:
      ErrThrow("NSE2D::mg_levels: NSTYPE 3 is not supported");
      break;
      
    case 4:
    case 14:
       mg_l = new TNSE_MGLevel4(i, s.matrix.get_A_block(0), 
                                s.matrix.get_A_block(1),
                                s.matrix.get_A_block(2),
                                s.matrix.get_A_block(3),
                                s.matrix.get_B_block(0),
                                s.matrix.get_B_block(1),
                                s.matrix.get_BT_block(0),
                                s.matrix.get_BT_block(1),
                                s.rhs.get_entries(), 
                                s.solution.get_entries(), 
                                n_aux, alpha, v_space_code, p_space_code, 
                                nullptr, nullptr);
    if(TDatabase::ParamDB->NSTYPE == 14)
      Output::print<1>("WARNING: NSTYPE 14 is not supported. C block is ignored");
    break;
  }
  return mg_l;
}

/** ************************************************************************ */
void NSE2D :: mg_solver()
{
  System_per_grid& s = this->systems.front();
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
  
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      sqMat[0] = s.matrix.get_A_block(0);
      recMat[0] = s.matrix.get_B_block(0);
      recMat[1] = s.matrix.get_B_block(1);
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
      break;
    case 2:
      sqMat[0] = s.matrix.get_A_block(0);
      
      recMat[0] = s.matrix.get_B_block(0);
      recMat[1] = s.matrix.get_B_block(1);
      recMat[2] = s.matrix.get_BT_block(0);
      recMat[3] = s.matrix.get_BT_block(1);
      
      MatVect = MatVect_NSE2;
      Defect = Defect_NSE2;
      break;
    case 3:
      sqMat[0] = s.matrix.get_A_block(0);
      sqMat[1] = s.matrix.get_A_block(1);
      sqMat[2] = s.matrix.get_A_block(2);
      sqMat[3] = s.matrix.get_A_block(3);
      recMat[0] = s.matrix.get_B_block(0);
      recMat[1] = s.matrix.get_B_block(1);
      MatVect = MatVect_NSE3;
      Defect = Defect_NSE3;
      break;
    case 4:
      sqMat[0] = s.matrix.get_A_block(0);
      sqMat[1] = s.matrix.get_A_block(1);
      sqMat[2] = s.matrix.get_A_block(2);
      sqMat[3] = s.matrix.get_A_block(3);
      recMat[0] = s.matrix.get_B_block(0);
      recMat[1] = s.matrix.get_B_block(1);
      recMat[2] = s.matrix.get_BT_block(0);
      recMat[3] = s.matrix.get_BT_block(1);
      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
      break;
    case 14:
      Output::print<1>("WARNING: NSTYPE 14 is not fully supported, take NSTYPE 4");
      exit(-4711);
      break;
  }


  if(TDatabase::ParamDB->SOLVER_TYPE ==1)
  {
    switch(TDatabase::ParamDB->SC_SOLVER_SADDLE)
    {
      case 11:
        zero_start = 1;
        break; 
      case 16:
        zero_start = 0;
        break;
    }
    switch(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE)
    {
      case 5:
        prec = new TMultiGridIte(MatVect, Defect, nullptr, 0, n_dof, 
                                 this->multigrid.get(), zero_start);
        break;
      default:
        ErrThrow("Unknown preconditioner !!!");
    }
    
    if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
    {
      itmethod_sol = new double[n_dof];
      itmethod_rhs = new double[n_dof];
      
      memcpy(itmethod_sol, s.solution.get_entries(), n_dof*SizeOfDouble);
      memcpy(itmethod_rhs, s.rhs.get_entries(), n_dof*SizeOfDouble);
    }
    else
    {
      itmethod_sol = s.solution.get_entries();
      itmethod_rhs = s.rhs.get_entries();
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
        ErrThrow("Unknown preconditioner !!!");
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
    memcpy(s.solution.get_entries(), itmethod_sol, n_dof*SizeOfDouble);
    memcpy(s.rhs.get_entries(), itmethod_rhs, n_dof*SizeOfDouble);
    
    delete itmethod; delete prec;
    delete [] itmethod_rhs;
    delete [] itmethod_sol;
  }
}

/** ************************************************************************ */
const NSE2D::Residuals& NSE2D::getResiduals() const
{
  return this->oldResiduals.back();
}

/** ************************************************************************ */
double NSE2D::getImpulsResidual() const
{
  return this->oldResiduals.back().impulsResidual;
}

/** ************************************************************************ */
double NSE2D::getMassResidual() const
{
  return this->oldResiduals.back().massResidual;
}

/** ************************************************************************ */
double NSE2D::getFullResidual() const
{
  return this->oldResiduals.back().fullResidual;
}

/** ************************************************************************ */
NSE2D::Residuals::Residuals()
 : impulsResidual(1e10), massResidual(1e10), fullResidual(1e10)
{}

/** ************************************************************************ */
NSE2D::Residuals::Residuals(double imR, double maR)
 : impulsResidual(sqrt(imR)), massResidual(sqrt(maR)),
   fullResidual(sqrt(imR + maR))
{}

/** ************************************************************************ */
std::ostream& operator<<(std::ostream& s, const NSE2D::Residuals& n)
{
  s << setw(14) << n.impulsResidual << "\t" << setw(14)
    << n.massResidual << "\t" << setw(14) << n.fullResidual;
  return s;
}

