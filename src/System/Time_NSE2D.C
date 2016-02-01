#include <Time_NSE2D.h>
#include <Database.h>
#include <Assemble2D.h>
#include <LinAlg.h>
#include <ItMethod.h>
#include <MultiGridIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <Output2D.h>
#include <DirectSolver.h>

/**************************************************************************** */
Time_NSE2D::System_per_grid::System_per_grid(const Example_NSE2D& example, 
                              TCollection& coll, std::pair<int,int> velo_pres_order)
 : velocity_space(&coll, (char*)"u", (char*)"velocity space",  example.get_bc(0),
                  velo_pres_order.first, nullptr),
   pressure_space(&coll, (char*)"p", (char*)"Darcy pressure", example.get_bc(1),
                  velo_pres_order.second, nullptr),
   matrix(this->velocity_space, this->pressure_space, example.get_bd()),
   Mass_Matrix(this->velocity_space, this->pressure_space, example.get_bd(), true),
   rhs(this->matrix, true),
   solution(this->matrix, false),
   u(&this->velocity_space, (char*)"u", (char*)"u", this->solution.block(0),
     this->solution.length(0), 2),
   p(&this->pressure_space, (char*)"p", (char*)"p", this->solution.block(2),
     this->solution.length(2))
{
  
}

/**************************************************************************** */
Time_NSE2D::Time_NSE2D(const TDomain& domain, int reference_id)
  : Time_NSE2D(domain, *(new Example_NSE2D()), reference_id)
{
  
}

/**************************************************************************** */
Time_NSE2D::Time_NSE2D(const TDomain& domain, const Example_NSE2D& ex, 
                       int reference_id)
 : systems(), example(ex), multigrid(), defect(), 
   oldResidual(0), initial_residual(1e10), errors(6,0.), oldtau(0.0)
{
  this->set_parameters();
  
  std::pair <int,int> velo_pres_order(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  this->get_velocity_pressure_orders(velo_pres_order);
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  
  this->systems.emplace_back(example, *coll, velo_pres_order);
  
  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);
  
  // print out some information  
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 2 * n_u + n_p; // total number of degrees of freedom
  int nActive = this->get_velocity_space().GetN_ActiveDegrees();  
  
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  Output::print<1>("N_Cells     : ", setw(10), coll->GetN_Cells());
  Output::print<1>("h (min,max) : ", setw(10), h_min ," ", setw(12), h_max);
  Output::print<1>("dof Velocity: ", setw(10), 2* n_u);
  Output::print<1>("dof Pressure: ", setw(10), n_p   );
  Output::print<1>("dof all     : ", setw(10), n_dof );
  Output::print<1>("active dof  : ", setw(10), 2*nActive);
  
  TFEFunction2D * u1 = this->systems.front().u.GetComponent(0);
  TFEFunction2D * u2 = this->systems.front().u.GetComponent(1);
  
  u1->Interpolate(example.get_initial_cond(0));
  u2->Interpolate(example.get_initial_cond(1));
  
  // done with the conrtuctor in case we're not using multigrid
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE!= 5 
    || TDatabase::ParamDB->SOLVER_TYPE != 1)
    return;
  // else multigrid
  
  // create spaces, functions, matrices on coarser levels
  double *param = new double[10];
  param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  param[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
  param[9] = 0;
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
    this->systems.emplace_back(example, *coll, velo_pres_order);
  }
  
  // create multigrid-level-objects, must be coarsest first
  unsigned int i = 0;
  for(auto it = this->systems.rbegin(); it != this->systems.rend(); ++it)
  {
    ErrThrow("NSE2D-multigrid needs to be checked");
    this->multigrid->AddLevel(this->mg_levels(i, *it));
    i++;
  }
}

/**************************************************************************** */
void Time_NSE2D::set_parameters()
{
  if(TDatabase::ParamDB->EXAMPLE < 101)
  {
    ErrMsg("Example " << TDatabase::ParamDB->EXAMPLE 
    <<"does not supported for time dependent problem");
    exit(1);
  }
  
  if(TDatabase::TimeDB->TIME_DISC == 0)
  {
    ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC 
          << " does not supported");
    throw("TIME_DISC: 0 does not supported");
  }  
}

/**************************************************************************** */
void Time_NSE2D::get_velocity_pressure_orders(std::pair< int, int > &velo_pres_order)
{
  int velocity_order = velo_pres_order.first;
  int pressure_order = velo_pres_order.second;
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
  velo_pres_order.first = order;
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
    // discontinuous spaces
    case 1:case 2: case 3: case 4: case 5:
      pressure_order = -(velocity_order-1)*10;
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
  }
  TDatabase::ParamDB->PRESSURE_SPACE  = pressure_order;
  velo_pres_order.second = pressure_order;
  
  Output::print("velocity space", setw(10), velo_pres_order.first);
  Output::print("pressure space", setw(10), velo_pres_order.second);
}

/**************************************************************************** */
void Time_NSE2D::assemble_initial_time()
{
  LocalAssembling2D_type nsall_rhs = LocalAssembling2D_type::TNSE2D;
  
  for(auto &s: this->systems)
  {
    TFEFunction2D *fe_functions[3] = 
       { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    
    
    LocalAssembling2D la_nsall_rhs(nsall_rhs, fe_functions, 
                              this->example.get_coeffs());
    
    s.matrix.Assemble(la_nsall_rhs, s.rhs, &s.Mass_Matrix);
    // set rhs for Dirichlet dofs
    s.solution.copy_nonactive(s.rhs);
  }
    // copy the current right hand side vector to the old_rhs 
  this->old_rhs = this->systems.front().rhs;  
}

/**************************************************************************** */
void Time_NSE2D::assemble_rhs()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  const double theta2 = TDatabase::TimeDB->THETA2;
  const double theta3 = TDatabase::TimeDB->THETA3;
  const double theta4 = TDatabase::TimeDB->THETA4;
  
  System_per_grid& s = this->systems.front();
  s.rhs.reset();

  TFEFunction2D *fe_functions[3] = 
       { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
  
  LocalAssembling2D la_rhs(TNSE2D_Rhs, fe_functions,
                           this->example.get_coeffs());
  this->AssembleRhs(la_rhs, s.rhs);
  
  // s.rhs.print("c");
  // now it is this->systems[i].rhs = f^k
  // scale by time step length and theta4 (only active dofs)  
  s.rhs.scaleActive(tau*theta4);
  // add rhs from previous time step 
  s.rhs.addScaledActive((this->old_rhs), tau*theta3);
  
  // now it is this->systems[i].rhs = tau*theta3*f^{k-1} + tau*theta4*f^k
  // next we want to set old_rhs to f^k (to be used in the next time step)
  this->old_rhs.addScaledActive(s.rhs, -1./(tau*theta3));
  this->old_rhs.scaleActive(-theta3/theta4);

  s.Mass_Matrix.applyScaledAddActive(s.solution.get_entries(), 
					s.rhs.get_entries(), 1.0);
  
  double factor = tau*theta2;
  s.matrix.applyScaledAddActive(s.solution.get_entries(), 
				   s.rhs.get_entries(), -factor);
  
  // scale the BT blocks with time step length
  for(System_per_grid& s : this->systems)
  {
    if(tau != this->oldtau)
    {
      factor = tau;
      if(this->oldtau != 0.0)
      {
        factor /= this->oldtau;
        Output::print<1>("change in tau", this->oldtau, "->", tau);
      }
      // scale the gradient matrices with step length tau
      s.matrix.get_BT_block(0)->operator*=(factor);
      s.matrix.get_BT_block(1)->operator*=(factor);
      if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
      {
        s.matrix.get_B_block(0)->operator*=(factor);
        s.matrix.get_B_block(1)->operator*=(factor);
      }
    }
  }
  this->oldtau = tau;
  // copy non active from rhs into solution vector
  s.solution.copy_nonactive(s.rhs);
  
  if(TDatabase::ParamDB->SOLVER_TYPE == GMG
     && TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
     this->multigrid->RestrictToAllGrids();
  Output::print<3>("assembling system matrix ");
}

/**************************************************************************** */
void Time_NSE2D::AssembleRhs(LocalAssembling2D& la, BlockVector& rhs)
{
  int N_Rhs = 3;
  const TFESpace2D * v_space = &this->get_velocity_space();
  const TFESpace2D * p_space = &this->get_pressure_space();
  
  rhs.reset();
  double *RHSs[3] = {rhs.block(0), rhs.block(1), rhs.block(2)};
  
  const TFESpace2D *fespmat[2] = {v_space, p_space};
  const TFESpace2D *fesprhs[3] = {v_space, v_space, p_space};
  
  BoundCondFunct2D * boundary_conditions[3] = {
             v_space->GetBoundCondition(), v_space->GetBoundCondition(), 
              p_space->GetBoundCondition() };

   std::array<BoundValueFunct2D*, 3> non_const_bound_values;
   non_const_bound_values[0] = this->example.get_bd(0);
   non_const_bound_values[1] = this->example.get_bd(1);
   non_const_bound_values[2] = this->example.get_bd(2);

   Assemble2D(1, fespmat, 0, nullptr,
              0, nullptr, N_Rhs, RHSs, fesprhs,
              boundary_conditions, non_const_bound_values.data(), la);

  //rhs.print("hh");exit(0);
}

/**************************************************************************** */
void Time_NSE2D::assemble_system()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;
  
  for(System_per_grid& s : this->systems)
  {
    s.matrix.scaleActive(factor);
    s.matrix.addScaledActive(s.Mass_Matrix,1.0);
  }
}

/**************************************************************************** */
void Time_NSE2D::assemble_nonlinear_term()
{
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  for(System_per_grid& s : this->systems)
  {
    TFEFunction2D *fe_functions[3] = 
      { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    LocalAssembling2D la_nonlinear(TNSE2D_NL, fe_functions,
                                   this->example.get_coeffs());
    s.matrix.AssembleNonLinear(la_nonlinear);
  }

  
  Output::print<3>("Assembling nonlinear matrices");
}

/**************************************************************************** */
bool Time_NSE2D::stopIte(unsigned int it_counter)
{
  System_per_grid& s = this->systems.front();
  unsigned int nuDof = s.solution.length(0);
  unsigned int npDof = s.solution.length(2);
  
  this->defect = s.rhs; 
  s.matrix.apply_scaled_add(s.solution.get_entries(),
                            defect.get_entries(),-1.);
  // 
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20FEFunction(&defect[2*nuDof], npDof, &this->get_pressure_space(),
                      TDatabase::ParamDB->VELOCITY_SPACE, 
                      TDatabase::ParamDB->PRESSURE_SPACE);
  double residual =  Ddot(2*nuDof+npDof, &this->defect[0], &this->defect[0]);
  double impulse_residual = Ddot(2*nuDof, &this->defect[0],
         &this->defect[0]);
  double mass_residual    = Ddot(npDof,&this->defect[2*nuDof],
         &this->defect[2*nuDof]);
  
  Output::print("nonlinear step  :  " , setw(3), it_counter);
  Output::print("impulse_residual:  " , setw(3), impulse_residual);
  Output::print("mass_residual   :  " , setw(3), mass_residual);
  Output::print("residual        :  " , setw(3), sqrt(residual));
  
  if (it_counter>0)
  {
  Output::print("rate:           :  " , setw(3), sqrt(residual)/oldResidual);
  }
  
  oldResidual = sqrt(residual);
  if(it_counter == 0)
    initial_residual = sqrt(residual);
  
  int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  double limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  if (TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE)
  {
    limit *= sqrt(this->get_size());
    Output::print("stopping tolerance for nonlinear iteration ", limit);
  }
  
  if ((((sqrt(residual)<=limit)||(it_counter==Max_It)))
   && (it_counter>=TDatabase::ParamDB->SC_MINIT))
   {
     Output::print("ITE : ", setw(3), it_counter, "  RES : ", sqrt(residual), 
                   " Reduction : ",  sqrt(residual)/initial_residual);
     // descale the matrices, since only the diagonal A block will 
     // be reassembled in the next time step
     this->deScaleMatrices();     
     return true;
   }
   else
     return false;
}

/**************************************************************************** */
void Time_NSE2D::solve()
{
  System_per_grid& s = this->systems.front();

  if((TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE !=5)
    || (TDatabase::ParamDB->SOLVER_TYPE !=1 ))
  {
    if(TDatabase::ParamDB->SOLVER_TYPE != 2)
      ErrThrow("only the direct solver is supported currently");
    
    /// @todo consider storing an object of DirectSolver in this class
    // use keyword class here, until all methods with the same name are removed
    class DirectSolver direct_solver(s.matrix, 
                                     DirectSolver::DirectSolverTypes::umfpack);
    direct_solver.solve(s.rhs, s.solution);
  }
  else
    this->mg_solver();
  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11 and A22 matrices are 
  // reset and assembled again but the A12 and A21 are scaled, so
  // for the next iteration we have to descale, see assemble_system()
    this->deScaleMatrices();

  Output::print<3>("solver done");
}

/**************************************************************************** */
void Time_NSE2D::deScaleMatrices()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;  
  for(System_per_grid& s : this->systems)
  {
    s.matrix.addScaledActive(s.Mass_Matrix, -1.0);
    s.matrix.scaleActive(1./factor);
  }
}

/**************************************************************************** */
TNSE_MGLevel* Time_NSE2D::mg_levels(int i, Time_NSE2D::System_per_grid& s)
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
    break;
    case 14:
      mg_l = new TNSE_MGLevel14(i, s.matrix.get_A_block(0), 
                                s.matrix.get_A_block(1),
                                s.matrix.get_A_block(2),
                                s.matrix.get_A_block(3),
                                s.matrix.get_C_block(),
                                s.matrix.get_B_block(0),
                                s.matrix.get_B_block(1),
                                s.matrix.get_BT_block(0),
                                s.matrix.get_BT_block(1),
                                s.rhs.get_entries(), 
                                s.solution.get_entries(), 
                                n_aux, alpha, v_space_code, p_space_code, 
                                nullptr, nullptr);
    break;
  }
  return mg_l;
}

/**************************************************************************** */
void Time_NSE2D::mg_solver()
{
  System_per_grid& s = this->systems.front(); 
  TSquareMatrix2D *sqMat[5];
  TSquareMatrix **sqmatrices = (TSquareMatrix **)sqMat;
  TMatrix2D *recMat[4];
  TMatrix **matrices = (TMatrix **)recMat;
  MatVecProc *MatVect;
  DefectProc *Defect;  
  
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      sqMat[0]=s.matrix.get_A_block(0);
      recMat[0] = s.matrix.get_B_block(0);
      recMat[1] = s.matrix.get_B_block(1);
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
      break;
    case 2:
      sqMat[0]=s.matrix.get_A_block(0);
      
      recMat[0] = s.matrix.get_B_block(0);
      recMat[1] = s.matrix.get_B_block(1);
      recMat[2] = s.matrix.get_BT_block(0);
      recMat[3] = s.matrix.get_BT_block(1);
      MatVect = MatVect_NSE2;
      Defect = Defect_NSE2;
      break;
    case 3:
      sqMat[0]=s.matrix.get_A_block(0);
      sqMat[1]=s.matrix.get_A_block(1);
      sqMat[2]=s.matrix.get_A_block(2);
      sqMat[3]=s.matrix.get_A_block(3);
      
      recMat[0] = s.matrix.get_B_block(0);
      recMat[1] = s.matrix.get_B_block(1);
      MatVect = MatVect_NSE3;
      Defect = Defect_NSE3;
      break;
    case 4:
      sqMat[0]=s.matrix.get_A_block(0);
      sqMat[1]=s.matrix.get_A_block(1);
      sqMat[2]=s.matrix.get_A_block(2);
      sqMat[3]=s.matrix.get_A_block(3);
      
      recMat[0] = s.matrix.get_B_block(0);
      recMat[1] = s.matrix.get_B_block(1);
      recMat[2] = s.matrix.get_BT_block(0);
      recMat[3] = s.matrix.get_BT_block(1);
      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
      break;
  }
 
  int zero_start;  
  int nDof = this->get_size();
  double *itmethod_rhs, *itmethod_sol;
  TItMethod *itmethod, *prec;
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
        prec = new TMultiGridIte(MatVect, Defect, nullptr, 0, nDof, 
                                 this->multigrid.get(), zero_start);
        break;
      default:
        ErrThrow("Unknown preconditioner !!!");
    }
    
    if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
    {
      itmethod_sol = new double[nDof];
      itmethod_rhs = new double[nDof];
      
      memcpy(itmethod_sol, s.solution.get_entries(), nDof*SizeOfDouble);
      memcpy(itmethod_rhs, s.rhs.get_entries(), nDof*SizeOfDouble);
    }
    else
    {
      itmethod_sol = s.solution.get_entries();
      itmethod_rhs = s.rhs.get_entries();
    }

    switch(TDatabase::ParamDB->SC_SOLVER_SADDLE)
    {
      case 11:
        itmethod = new TFixedPointIte(MatVect, Defect, prec, 0, nDof, 0);
        break;
      case 16:
        itmethod = new TFgmresIte(MatVect, Defect, prec, 0, nDof, 0);
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
    memcpy(s.solution.get_entries(), itmethod_sol, nDof*SizeOfDouble);
    memcpy(s.rhs.get_entries(), itmethod_rhs, nDof*SizeOfDouble);
    
    delete itmethod; delete prec;
    delete [] itmethod_rhs;
    delete [] itmethod_sol;
  }
}

/**************************************************************************** */
void Time_NSE2D::output(int m, int& image)
{
  if(!TDatabase::ParamDB->WRITE_VTK 
    && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;

  System_per_grid& s = this->systems.front();
  TFEFunction2D * u1 = s.u.GetComponent(0);
  TFEFunction2D * u2 = s.u.GetComponent(1);
  
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       s.p.project_into_L20();

  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    s.p.PrintMinMax();
  }

  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    double locerr[8];
    MultiIndex2D allderiv[3]= {D00, D10, D01};
    const TFESpace2D *v_sp = &this->get_velocity_space();
    const TFESpace2D *p_sp = &this->get_pressure_space();
    TAuxParam2D aux;
    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    
    u1->GetErrors(example.get_exact(0), 3, allderiv, 2, L2H1Errors,nullptr,
                  &aux,1, &v_sp,locerr);
    
    u2->GetErrors(example.get_exact(1), 3, allderiv, 2, L2H1Errors,nullptr,
                  &aux,1, &v_sp,locerr+2);

    errors[0] += (locerr[0]*locerr[0]+locerr[2]*locerr[2] 
                  + this->errors[1])*tau*0.5;
    errors[1] = locerr[0]*locerr[0]+locerr[2]*locerr[2];
    errors[2] += (locerr[1]*locerr[1]+locerr[3]*locerr[3] 
                  + this->errors[3])*tau*0.5;
    errors[3] = locerr[1]*locerr[1]+locerr[3]*locerr[3];  

    Output::print<1>("L2(u) : ", sqrt(this->errors[1]));
    Output::print<1>("H1-semi(u) : ", sqrt(this->errors[3]));

    Output::print<1>("L2(0,t,L2) : ", sqrt(this->errors[0]));
    Output::print<1>("L2(0,t,H1) : ", sqrt(this->errors[2]));

    s.p.GetErrors(example.get_exact(2), 3, allderiv, 2, L2H1Errors, 
                  nullptr, &aux, 1, &p_sp, locerr);

    Output::print<1>("L2(p) : ", locerr[0]);
    Output::print<1>("H1-semi(p)) : " , locerr[1] );

    errors[4] += (locerr[0]*locerr[0] + this->errors[5])*tau*0.5;
    errors[5] = locerr[0]*locerr[0];
    Output::print<1>("L2(0,t,L2) : ", sqrt(errors[4]) );    
  }
   delete u1;
   delete u2;
  
  if((m==0) || (m%TDatabase::TimeDB->STEPS_PER_IMAGE) )
  {
    if(TDatabase::ParamDB->WRITE_VTK)
    {
      TOutput2D output(2, 3, 1, 0, NULL);
      output.AddFEFunction(&s.p);
      output.AddFEVectFunct(&s.u);
      std::string filename(TDatabase::ParamDB->OUTPUTDIR);
      filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
      if(image<10) filename += ".0000";
      else if(image<100) filename += ".000";
      else if(image<1000) filename += ".00";
      else if(image<10000) filename += ".0";
      else filename += ".";
      filename += std::to_string(image) + ".vtk";
      output.WriteVtk(filename.c_str());
      image++;
    }
  }
}
