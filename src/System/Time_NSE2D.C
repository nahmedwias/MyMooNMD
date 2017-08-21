#include <Time_NSE2D.h>
#include <Database.h>
#include <Assemble2D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <BoundaryAssembling2D.h>
#include <GridTransfer.h>
#include <Domain.h>

/* *************************************************************************** */
  //TODO  So far of this object only the nonlin it stuff is used - switch entirely!
ParameterDatabase get_default_TNSE2D_parameters()
{
  Output::print<5>("creating a default TNSE2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default TNSE2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TNSE2D parameter database");

  //Time_NSE2D requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  return db;
}
/* *************************************************************************** */

/**************************************************************************** */
Time_NSE2D::System_per_grid::System_per_grid(const Example_TimeNSE2D& example, 
                  TCollection& coll, std::pair< int, int > order, 
                  Time_NSE2D::Matrix type)
 : velocity_space(&coll, (char*)"u", (char*)"velocity space",  example.get_bc(0),
                  order.first, nullptr),
   pressure_space(&coll, (char*)"p", (char*)"pressure space", example.get_bc(2),
                  order.second, nullptr),
   matrix({&velocity_space, &velocity_space, &pressure_space}),
   Mass_Matrix({&velocity_space, &velocity_space}),
   rhs(matrix, true),
   solution(matrix, false),
   u(&velocity_space, (char*)"u", (char*)"u", solution.block(0), 
     solution.length(0), 2),
   p(&pressure_space, (char*)"p", (char*)"p", this->solution.block(2),
     solution.length(2)), 
   MatrixK({&velocity_space, &velocity_space, &pressure_space}),
   solution_m1(matrix, false),
   u_m1(&velocity_space, (char*)"u", (char*)"u", solution_m1.block(0), 
        solution_m1.length(0), 2),
   solution_m2(matrix, false), 
   u_m2(&velocity_space, (char*)"u", (char*)"u", solution_m2.block(0), 
        solution_m2.length(0), 2),
   combined_old_sols(matrix, false),
   comb_old_u(&velocity_space, (char*)"u", (char*)"u", combined_old_sols.block(0), 
        combined_old_sols.length(0), 2), 
   extrapolate_sol(matrix, false),
   extrapolate_u(&velocity_space, (char*)"u", (char*)"u", extrapolate_sol.block(0), 
        extrapolate_sol.length(0), 2)
{
  // Mass Matrix
  // Output::increaseVerbosity(5);
  //solution_m2 = solution_m1;
  Mass_Matrix = BlockFEMatrix::Mass_NSE2D(velocity_space);
      
  switch(type)
  {
    case Time_NSE2D::Matrix::Type1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
      break;
    case Time_NSE2D::Matrix::Type2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
      break;
    case Time_NSE2D::Matrix::Type3:
      matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
      break;
    case Time_NSE2D::Matrix::Type4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      if(TDatabase::ParamDB->DISCTYPE == RESIDUAL_VMS)
      {
        MatrixK = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      }
      break;
    case Time_NSE2D::Matrix::Type14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
      if(TDatabase::ParamDB->DISCTYPE == RESIDUAL_VMS)
      {
        MatrixK = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      }
      break;
  }
}

/**************************************************************************** */
Time_NSE2D::Time_NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
                       int reference_id)
  : Time_NSE2D(domain, param_db, Example_TimeNSE2D(param_db), reference_id)
{
  
}

/**************************************************************************** */
Time_NSE2D::Time_NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
                       const Example_TimeNSE2D& ex, int reference_id)
 : db(get_default_TNSE2D_parameters()), outputWriter(param_db), systems(),
   example(ex), solver(param_db), defect(), oldResiduals(), 
   initial_residual(1e10), errors(10,0.), oldtau(0.0)
{
  db.merge(param_db);
  this->set_parameters();
  
  std::pair <int,int> velo_pres_order(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  this->get_velocity_pressure_orders(velo_pres_order);
  
  Time_NSE2D::Matrix type;
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case  1: type = Matrix::Type1;  break;
    case  2: type = Matrix::Type2;  break;
    case  3: type = Matrix::Type3;  break;
    case  4: type = Matrix::Type4;  break;
    case 14: type = Matrix::Type14; break;
    default:
      ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
               " That NSE Block Matrix Type is unknown to class NSE2D.");
  }
  
  bool usingMultigrid = this->solver.is_using_multigrid();
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  // create finite element space, functions, matrices, rhs and solution
  // at the finest grid
  this->systems.emplace_back(example, *coll, velo_pres_order, type);
  
  TFEFunction2D * u1 = this->systems.front().u.GetComponent(0);
  TFEFunction2D * u2 = this->systems.front().u.GetComponent(1);
  
  u1->Interpolate(example.get_initial_cond(0));
  u2->Interpolate(example.get_initial_cond(1));

  if(usingMultigrid)
  {
    auto mg = solver.get_multigrid();
    bool mdml = mg->is_using_mdml();
    if(mdml)
    {
      // change the discretization on the coarse grids to lowest order 
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velo_pres_order = {-1, -4711};
      this->get_velocity_pressure_orders(velo_pres_order);
    }
    //determine multigrid levels
    int second_grid;
    int coarsest_grid = domain.get_ref_level() -
                         mg->get_n_geometric_levels() + 1;
    if(mdml)
      second_grid = domain.get_ref_level();
    else 
      second_grid = domain.get_ref_level() - 1;
    
    if(coarsest_grid < 0 )
    {
      ErrThrow("The domain has not been refined often enough to do multigrid "
          "on ", mg->get_n_geometric_levels(), " geometric levels (",
          mg->get_n_algebraic_levels()," algebraic levels). There are"
          " only ", domain.get_ref_level() + 1, " geometric grid levels.");
    }
    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    //prepare input argument for multigrid object
    matrices.push_back(&systems.back().matrix);
    for(int grid_no = second_grid; grid_no >= coarsest_grid; --grid_no)
    {
      TCollection *coll = domain.GetCollection(It_EQ, grid_no, reference_id);
      systems.emplace_back(example, *coll, velo_pres_order,
                            type);
      matrices.push_front(&systems.back().matrix);
    }
    mg->initialize(matrices);
  }
  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);
  
  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  
  // print out the information (cells, dofs, etc)
  this->output_problem_size_info();
}

/**************************************************************************** */
void Time_NSE2D::set_parameters()
{
  if(!db["problem_type"].is(6))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 6;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to Time_NSE."
          "It is now reset to the correct value for Time_NSE (=6).");
      db["problem_type"] = 6;
    }
  }
  if(TDatabase::TimeDB->TIME_DISC == 0)
  {
    ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC 
          << " does not supported");
    throw("TIME_DISC: 0 is not supported");
  }
  if(TDatabase::ParamDB->DISCTYPE == SUPG)
  {
    if(TDatabase::TimeDB->TIME_DISC !=5 && TDatabase::TimeDB->TIME_DISC !=1)
    {
      ErrThrow("TIME_DISC: " , TDatabase::TimeDB->TIME_DISC, 
               " does not supported for SUPG method" , 
               " Only (BDF2) TIME_DISC: 5 can be used");
    }
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
    case 1:case 2: case 3: case 4: case 5:
      // pressure order is chosen correctly
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
    default:
      ErrThrow("pressure space is not chosen properly ", pressure_order);
  }
  TDatabase::ParamDB->PRESSURE_SPACE  = pressure_order;
  velo_pres_order.second = pressure_order;
  
  Output::print("velocity space", setw(10), velo_pres_order.first);
  Output::print("pressure space", setw(10), velo_pres_order.second);
}

/**************************************************************************** */
void Time_NSE2D::assemble_initial_time()
{
  for(auto &s : this->systems)
  {
    call_assembling_routine(s, TNSE2D);
    // copy nonactives
    s.solution.copy_nonactive(s.rhs);
    //
    s.solution_m1 = s.solution;
    s.solution_m2 = s.solution;
  }
  // copy the current right hand side vector to the old_rhs 
  this->old_rhs = this->systems.front().rhs; 
  // set the solution vectors
  this->old_solution = this->systems.front().solution;
}

/**************************************************************************** */
void Time_NSE2D::assemble_rhs()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  const double theta2 = TDatabase::TimeDB->THETA2;
  const double theta3 = TDatabase::TimeDB->THETA3;
  const double theta4 = TDatabase::TimeDB->THETA4;
  System_per_grid& s = this->systems.front();
  
  // assemble right hand side 
  call_assembling_routine(s, TNSE2D_Rhs);
  // copy the non active to the solution vector
  // since the rhs vector will be passed to the solver
  // and is modified with matrix vector multiplication
  // which also uses the non-actives
  s.solution.copy_nonactive(s.rhs);
  
  // now it is this->systems[i].rhs = f^k
  // scale by time step length and theta4 (only active dofs)  
  s.rhs.scaleActive(tau*theta4);

  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
    this->modify_slip_bc(true);
   // add rhs from previous time step 
  if(theta3 != 0)
  {
    s.rhs.addScaledActive((this->old_rhs), tau*theta3);
    // now it is this->systems[i].rhs = tau*theta3*f^{k-1} + tau*theta4*f^k
    this->old_rhs.addScaledActive(s.rhs, -1./(tau*theta3));
    this->old_rhs.scaleActive(-theta3/theta4);
    this->old_rhs.copy_nonactive(s.rhs);
  }
  // FIXME FInd other solution than this submatrix method.
  // M u^{k-1}
  s.Mass_Matrix.apply_scaled_submatrix(old_solution, s.rhs, 2, 2, 1.0);
    
  // -tau*theta2 * A u^{k-1}
  double factor = -tau*theta2;
  s.matrix.apply_scaled_submatrix(old_solution, s.rhs, 2, 2, factor);
  
  // scale the BT blocks with time step length
  for(System_per_grid& s : this->systems)
  {
    if(tau != oldtau)
    {
      // TODO: change the factor to be THETA1*tau;
      factor = /*TDatabase::TimeDB->THETA1**/tau;
      if(this->oldtau != 0.0)
      {
        factor /= this->oldtau;
        Output::print<1>("change in tau", this->oldtau, "->", tau);
      }
      // scale the BT transposed blocks with the current time step
      const std::vector<std::vector<size_t>> cell_positions = {{0,2}, {1,2}};
	s.matrix.scale_blocks(factor, cell_positions);      
      if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
      {
        const std::vector<std::vector<size_t>> cell_positions_t = {{2,0}, {2,1}};
	s.matrix.scale_blocks(factor, cell_positions_t);
      }      
    }
  }
  this->oldtau = tau;
  // copy non active from solution into rhs vector
  s.rhs.copy_nonactive(s.solution);  

  Output::print<5>("assembled the system right hand side ");  
}

/**************************************************************************** */

void Time_NSE2D::assemble_system()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;
  
  for(System_per_grid& s : this->systems)
  {
    const std::vector<std::vector<size_t>>
      cell_positions = {{0,0}, {0,1}, {1, 0}, {1, 1}};
    // note: declaring the auxiliary cell_positions is needed by the compiler
    // to sort out the overriding of the function scale_blocks_actives(...,...)
    s.matrix.scale_blocks_actives(factor, cell_positions);
    if(TDatabase::ParamDB->NSTYPE == 14) {
      const std::vector<std::vector<size_t>>
	cell_positions_2_2 = {{2,2}};
      s.matrix.scale_blocks_actives(factor, cell_positions_2_2);
    }
    const FEMatrix& mass = *s.Mass_Matrix.get_blocks().at(0).get();
    s.matrix.add_matrix_actives(mass, 1.0, {{0,0}, {1,1}}, {false, false});
  }
  Output::print<5>("Assembled the system matrix which will be passed to the ", 
                   "solver");
}
/**************************************************************************** */

void Time_NSE2D::assemble_nonlinear_term()
{
  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(systems.size() > 1)
  {
    for( int block = 0; block < 2 ;++block)
    {
      std::vector<const TFESpace2D*> spaces;
      std::vector<double*> u_entries;
      std::vector<size_t> u_ns_dofs;
      for(auto &s : systems )
      {
        spaces.push_back(&s.velocity_space);
        u_entries.push_back(s.solution.block(block));
        u_ns_dofs.push_back(s.solution.length(block));
      }
      GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
    }
  }
  
  for(System_per_grid& s : this->systems)
    call_assembling_routine(s, TNSE2D_NL);
  
  if(TDatabase::ParamDB->DISCTYPE == SUPG)
  {
    // nonlinear right hand side only for the SUPG, residual based VMS methods
    this->assemble_rhs_supg();
  }

  Output::print<5>("Assembled the nonlinear matrix only ");
}
/**************************************************************************** */

bool Time_NSE2D::stopIte(unsigned int it_counter)
{
  //TODO This has no "slow convergence criterion yet!"
  System_per_grid& s = this->systems.front();
  unsigned int nuDof = s.solution.length(0);
  unsigned int npDof = s.solution.length(2);
  unsigned int sc_minit = db["nonlinloop_minit"];
  
  this->defect = s.rhs; 
  s.matrix.apply_scaled_add(s.solution, defect,-1.);
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
  Output::print<1>("nonlinear step : ", it_counter);
  if(it_counter > 0)
  {
    Output::print("Residuals  :  " , setw(15), impulse_residual, 
                setw(15), mass_residual, setw(15), sqrt(residual),
                setw(15), sqrt(residual)/oldResidual);
  }
  else
  {
    Output::print("Residuals  :  " , setw(15), impulse_residual, 
                setw(15), mass_residual, setw(15), sqrt(residual));
  }
  
  oldResidual = sqrt(residual);
  if(it_counter == 0)
    initial_residual = sqrt(residual);
  
  size_t Max_It = db["nonlinloop_maxit"];
  double limit = db["nonlinloop_epsilon"];
  if (db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= sqrt(this->get_size());
    Output::print("stopping tolerance for nonlinear iteration ", limit);
  }
  
  if ( (((sqrt(residual)<=limit)||(it_counter==Max_It))) )
   {
     for(System_per_grid& s: this->systems)
     {
       s.solution_m2 = s.solution_m1;
       s.solution_m1 = s.solution;
     }
     this->old_solution = s.solution;
     
     Output::print("ITE : ", setw(3), it_counter, "  RES : ", sqrt(residual), 
                   " Reduction : ",  sqrt(residual)/initial_residual);
     
     if(imex_scheme(0) && it_counter >0)
     {
       return true;
     }
     else
     {
       // descale the matrices, since only the diagonal A block will 
       // be reassembled in the next time step
       this->deScaleMatrices();
       return true;
     }     
   }
   else
     return false;  
}
/**************************************************************************** */
void Time_NSE2D::solve()
{
  System_per_grid& s = this->systems.front();
//  s.rhs.print("r");exit(0);
  solver.solve(s.matrix, s.rhs, s.solution);
  
  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11 and A22 matrices are 
  // reset and assembled again but the A12 and A21 are scaled, so
  // for the next iteration we have to descale, see assemble_system()
  this->deScaleMatrices();

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       s.p.project_into_L20();
}

/**************************************************************************** */
void Time_NSE2D::deScaleMatrices()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;  
  for(System_per_grid& s : this->systems)
  {
    const FEMatrix& mass_bloks = *s.Mass_Matrix.get_blocks().at(0).get();
    s.matrix.add_matrix_actives(mass_bloks, -1.0, {{0,0}, {1,1}}, {false, false});
    const std::vector<std::vector<size_t>>
    cell_positions = {{0,0}, {0,1}, {1, 0}, {1, 1}};
    // note: declaring the auxiliary cell_positions is needed by the compiler
    // to sort out the overriding of the function scale_blocks_actives(...,...)
    s.matrix.scale_blocks_actives(1./factor, cell_positions);
    
    if(TDatabase::ParamDB->NSTYPE == 14)
    {
      const std::vector<std::vector<size_t>>
	cell_positions_2_2 = {{2,2}};
      s.matrix.scale_blocks(1./factor, cell_positions_2_2);
    }
  }
}

/**************************************************************************** */
void Time_NSE2D::output(int m)
{
  bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
  if(no_output)
    return;

  System_per_grid& s = this->systems.front();
  TFEFunction2D * u1 = s.u.GetComponent(0);
  TFEFunction2D * u2 = s.u.GetComponent(1);

  if((size_t)db["verbosity"]> 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    s.p.PrintMinMax();
  }

  if(db["output_compute_errors"])
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

    Output::print<1>("L2(u) : ", setprecision(10), sqrt(this->errors[1]));
    Output::print<1>("H1-semi(u) : ", setprecision(10), sqrt(this->errors[3]));

    Output::print<1>("L2(0,t,L2(u)) : ", sqrt(this->errors[0]));
    Output::print<1>("L2(0,t,H1-semi(u)) : ", sqrt(this->errors[2]));

    s.p.GetErrors(example.get_exact(2), 3, allderiv, 2, L2H1Errors, 
                  nullptr, &aux, 1, &p_sp, locerr);

    Output::print<1>("L2(p) : ", setprecision(10), locerr[0]);
    Output::print<1>("H1-semi(p)) : " , setprecision(10), locerr[1] );

    errors[4] += (locerr[0]*locerr[0] + this->errors[5])*tau*0.5;
    errors[5] = locerr[0]*locerr[0];
    Output::print<1>("L2(0,t,L2(p)) : ", sqrt(errors[4]) );
    
    errors[6] += (locerr[1]*locerr[1] + this->errors[7])*tau*0.5;
    errors[7] = locerr[1]*locerr[1];
    Output::print<1>("L2(0,t,H1-semi(p)) : ", sqrt(errors[6]) );
  }
  delete u1;
  delete u2;
  double val=0.;
  // example.do_post_processing(*this, val);
  
  if((m==0) || (m/TDatabase::TimeDB->STEPS_PER_IMAGE) )
  {
    if(db["output_write_vtk"])
    {
      outputWriter.write(TDatabase::TimeDB->CURRENTTIME);
    }
  }
}
/**************************************************************************** */
void Time_NSE2D::output_problem_size_info() const
{
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_u_active = this->get_velocity_space().GetN_ActiveDegrees();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 2 * n_u + n_p; // total number of degrees of freedom
  
  double h_min, h_max;
  TCollection * coll = this->get_velocity_space().GetCollection();
  coll->GetHminHmax(&h_min, &h_max);
  Output::stat("NSE2D", "Mesh data and problem size");
  Output::dash("cells              :  ", setw(10), coll->GetN_Cells());
  Output::dash("h (min, max)       :  ", setw(10), h_min, setw(10), " ", h_max);
  Output::dash("dof velocity       :  ", setw(10), 2*n_u );
  Output::dash("dof velocity active:  ", setw(10), 2*n_u_active);
  Output::dash("dof pressure       :  ", setw(10), n_p);
  Output::dash("dof all            :  ", setw(10), n_dof);
}

/**************************************************************************** */
std::array< double, int(6) > Time_NSE2D::get_errors()
{
  std::array<double, int(6)> error_at_time_points;
  error_at_time_points[0] = sqrt(errors[1]); // L2 velocity error
  error_at_time_points[1] = sqrt(errors[3]); // H1 velocity error
  error_at_time_points[2] = sqrt(errors[5]); // L2 pressure error
  error_at_time_points[3] = sqrt(errors[7]); // H1 pressure error
  
  return error_at_time_points;
}

/**************************************************************************** */
void Time_NSE2D::call_assembling_routine(Time_NSE2D::System_per_grid& s, 
                               LocalAssembling2D_type type)
{
  // set arrays of spaces for matrices and rhs 
  std::vector<const TFESpace2D*> spaces_mat;
  std::vector<const TFESpace2D*> spaces_rhs;
  std::vector<TFEFunction2D*> fefunctios;
  // call to routine to set arrays 
  set_arrays(s, spaces_mat, spaces_rhs, fefunctios);
  
  // prepare matrices and rhs for assembling 
  std::vector<TSquareMatrix2D*> sqMatrices;
  std::vector<TMatrix2D*> rectMatrices;
  std::vector<double*> rhs_array;
  // call the routine to prepare the matrices
  set_matrices_rhs(s, type, sqMatrices, rectMatrices, rhs_array);
  // boundary conditions and boundary values array
  // boundary conditions:
  BoundCondFunct2D* bc[3] = {
    s.velocity_space.GetBoundCondition(),
    s.velocity_space.GetBoundCondition(),
    s.pressure_space.GetBoundCondition()};
  
  // boundary values:
  std::vector<BoundValueFunct2D*>bv(3);
  bv[0]=example.get_bd(0);
  bv[1]=example.get_bd(1);
  bv[2]=example.get_bd(2);

  // local assembling settings
  LocalAssembling2D la(type, fefunctios.data(),
                           this->example.get_coeffs());

  // assemble all the matrices and right hand side 
  Assemble2D(spaces_mat.size(), spaces_mat.data(), 
               sqMatrices.size(), sqMatrices.data(), 
               rectMatrices.size(), rectMatrices.data(), 
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(), 
               bc, bv.data(), la);

  // add boundary integrals to stabilize backflow at open boundaries
  ///@todo this should be controlled by input parameters
  bool add_backflow_stab = false;
  if (add_backflow_stab) {
    BoundaryAssembling2D boundary_integral;
    int neumann_boundary_component = 17;
    double beta_backflow_stab = 0.4;
    const TFESpace2D * v_space = &s.velocity_space;
    std::vector< TFEFunction2D* > u_conv;
    u_conv.resize(2);
    u_conv[0] = s.u.GetComponent(0);
    u_conv[1] = s.u.GetComponent(1);
    
    boundary_integral.matrix_u_v_backflow_stab(s.matrix,
					       v_space,
					       u_conv,
					       neumann_boundary_component,   
					       beta_backflow_stab);     
  }
  
}

/**************************************************************************** */
void Time_NSE2D::set_matrices_rhs(Time_NSE2D::System_per_grid& s, 
                   LocalAssembling2D_type type, std::vector< TSquareMatrix2D* >& sqMat, 
                   std::vector< TMatrix2D* >& reMat, std::vector< double* >& rhs_array)
{
  rhs_array.resize(0);
  sqMat.resize(0);
  reMat.resize(0);
  
  std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix.get_blocks_uniquely();
         
  switch(type)
  {
    case TNSE2D:
    {
      // right hand side: for NSTYPE: 1,2 and 3, size is 2
        rhs_array.resize(2);
        rhs_array[0] = s.rhs.block(0);
        rhs_array[1] = s.rhs.block(1);
      // get the blocks of the mass matrix 
      std::vector<std::shared_ptr<FEMatrix>> mass_blocks
        = s.Mass_Matrix.get_blocks_uniquely();
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          if(blocks.size() != 3)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size());
          }
          sqMat.resize(2);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
          // rectangular matrices
          reMat.resize(2);
          reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get());
          reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
          break;
        case 2:
          if(blocks.size() != 5)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size());
          }
          sqMat.resize(2);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
          // rectangular matrices
          reMat.resize(4);
          reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(3).get());//first the lying B blocks
          reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(4).get());
          reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get()); //than the standing B blocks
          reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
          break;
        case 3:
          if(blocks.size() != 6)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size());
          }
          sqMat.resize(5);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
          sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());  
          sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
          // mass matrix
          sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
          // rectangular matrices
          reMat.resize(2);
          reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //first the lying B blocks
          reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
          break;
        case 4:
        case 14:
          if(TDatabase::ParamDB->NSTYPE==14 && blocks.size() != 9)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size());
          }
          else if(TDatabase::ParamDB->NSTYPE==4 && blocks.size() != 8)
          {
            ErrThrow("Wrong blocks.size() ", blocks.size());
          }
          sqMat.resize(5);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
          sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());  
          sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
          // mass matrix
          sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
          // rectangular matrices
          reMat.resize(4);
          reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
          reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
          reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
          reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
          
          // right hand side 
          rhs_array.resize(2);
          rhs_array[0] = s.rhs.block(0);
          rhs_array[1] = s.rhs.block(1);
          if(TDatabase::ParamDB->NSTYPE == 14)
          {
            // C block
            sqMat.resize(6);
            sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
            // additional right hand sides
            rhs_array.resize(3);
            rhs_array[2] = s.rhs.block(2);
          }                    
          break;        
      }
      // right hand sides are assembled for the initial time step
      // for the remaining time steps, they are assembled in 
      // another function. so reset to zero here
      s.rhs.reset();
      break;
    }
    case TNSE2D_NL:
    {
      std::vector<std::shared_ptr<FEMatrix>> blocks 
         = s.matrix.get_blocks_uniquely();
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          sqMat.resize(1);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          reMat.resize(0);
          break;
        case 3:
        case 4:
          sqMat.resize(2);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
          reMat.resize(0);
          // right hand side
          rhs_array.resize(0);
          if(TDatabase::ParamDB->DISCTYPE==SUPG)
          {
            sqMat.resize(3);
            std::vector<std::shared_ptr<FEMatrix>> mass_blocks
               = s.Mass_Matrix.get_blocks_uniquely();
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
            
            rhs_array.resize(2);            
            rhs_array[0] = s.rhs.block(0);
            rhs_array[1] = s.rhs.block(1);
            
            reMat.resize(2);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //first the lying B blocks
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
            s.rhs.reset();
          }
          break;
        case 14:
          // we need to re-assemble all the matrices due to the solution
          // dependency of the stabilization parameters
          sqMat.resize(6);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
          sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());  
          sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
          // mass matrix
          std::vector<std::shared_ptr<FEMatrix>> mass_blocks
               = s.Mass_Matrix.get_blocks_uniquely();
          sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
          // pressure-pressure block
          sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
          // rectangular matrices
          reMat.resize(4);
          reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); 
          reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
          reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); 
          reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
          
          rhs_array.resize(3);
          rhs_array[0] = s.rhs.block(0);
          rhs_array[1] = s.rhs.block(1);
          rhs_array[2] = s.rhs.block(2);
          s.rhs.reset();
          break;
      }// endswitch NSTYPE
      break;
    }
    case TNSE2D_Rhs:
    {
      // no matrices to be assembled
      sqMat.resize(0);
      reMat.resize(0);
      // right hand side      
      rhs_array.resize(2);
      rhs_array[0]=s.rhs.block(0);
      rhs_array[1]=s.rhs.block(1);
      if(TDatabase::ParamDB->NSTYPE == 14)
      {
        rhs_array.resize(3);
        rhs_array[2] = s.rhs.block(2); // pressure block
      }
      // reset them to zero
      s.rhs.reset();
      break;
    }
  }
  // reset matrices
  for(auto sm : sqMat)
    sm->reset();
  for(auto rm : reMat)
    rm->reset();
}

/**************************************************************************** */
void Time_NSE2D::set_arrays(Time_NSE2D::System_per_grid& s, 
                            std::vector< const TFESpace2D* >& spaces, 
                            std::vector< const TFESpace2D* >& spaces_rhs,
                            std::vector< TFEFunction2D* >& functions)
{
  spaces.resize(2);
  spaces_rhs.resize(2);
  
  spaces[0] = &s.velocity_space;
  spaces[1] = &s.pressure_space;
  
  spaces_rhs[0] = &s.velocity_space;
  spaces_rhs[1] = &s.velocity_space;
  if(TDatabase::ParamDB->NSTYPE == 14)
  {
    spaces_rhs.resize(3);
    spaces_rhs[2] = &s.pressure_space;
  }
  
  bool is_imex = imex_scheme(0);
  functions.resize(3);
  if(!is_imex)
  {
    functions[0] = s.u.GetComponent(0);
    functions[1] = s.u.GetComponent(1);
    functions[2] = &s.p;
  }
  else
  {
    // extrapolate velocity and pressure if needed
    s.extrapolate_sol.reset();
    s.extrapolate_sol = s.solution_m1;
    s.extrapolate_sol.scale(2.);
    s.extrapolate_sol.add_scaled(s.solution_m2, -1.);   
    
    functions[0] = s.extrapolate_u.GetComponent(0);
    functions[1] = s.extrapolate_u.GetComponent(1);
    functions[2] = &s.p;
  }
  
  if(TDatabase::ParamDB->NSTYPE == 14)
  {
    if((TDatabase::TimeDB->TIME_DISC==1) || (current_step_ <= 1))
    {
      functions.resize(4);
      functions[2] = s.u_m1.GetComponent(0);
      functions[3] = s.u_m1.GetComponent(1);

    }
    else
    {
      // reset 
      s.combined_old_sols.reset();
      // copy and scale the solution at previous time step with 
      // factor 2
      s.combined_old_sols = s.solution_m1;
      s.combined_old_sols.scale(2.);
      // subtract with right factor the solution at pre-previous 
      // solution
      s.combined_old_sols.add_scaled(s.solution_m2, -1./2.);

      functions.resize(4);
      functions[2] = s.comb_old_u.GetComponent(0);
      functions[3] = s.comb_old_u.GetComponent(1);
    }
  }
}
/**************************************************************************** */
bool Time_NSE2D::imex_scheme(bool print_info)
{
  
  if(!TDatabase::TimeDB->EXTRAPOLATE_VELOCITY)
    return false;
  bool interruption_condition  = 
    TDatabase::TimeDB->EXTRAPOLATE_VELOCITY * (current_step_ >= 3);
    
  // change maximum number of nonlin_iterations to 1 in IMEX case
  if (interruption_condition)
  {
    db["nonlinloop_maxit"] = 1;
    if(print_info) // condition is here just to print it once
        Output::info<1>("Nonlinear Loop MaxIteration",
                        "The parameter 'nonlinloop_maxit' was changed to 1."
                        " Only one non-linear iteration is done, because the IMEX scheme was chosen.\n");
  }
  return interruption_condition;
}

/**************************************************************************** */
void Time_NSE2D::perform_bdf1_first()
{
  TDatabase::TimeDB->TIME_DISC = 1;
  db["time_discretization"] = 1;
  // assemble rhs
  this->assemble_rhs();
  // assemble nonlinear matrices
  this->assemble_nonlinear_term();
  // prepare matrices for solver
  this->assemble_system();
  //
  for(unsigned int i=0; ; i++)
  {
    if(stopIte(i))
      break;
    this->solve();

    this->assemble_nonlinear_term();
    this->assemble_system();
  }
  output(1);

  // declare the cell_position vectors
  // Mac compiler cannot distinguish between the
  // two variants of BlockFEMatrix::scale_blocks
  const std::vector<std::vector<size_t>>
    cell_positions_0_2_1_2 = {{0,2},{1,2}};
  const std::vector<std::vector<size_t>>
    cell_positions_2_0_2_1 = {{2,0},{2,1}};
  
  // descale the B and BT blocks:: these needs 
  // to be scaled with the correct scales for 
  // the BDF2 method
  
  
  for(System_per_grid& s : this->systems)
  {
    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
     
     
     s.matrix.scale_blocks(1./tau, cell_positions_0_2_1_2);
    if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
    {
      
      s.matrix.scale_blocks(1./tau, cell_positions_2_0_2_1);
    }
  }
  // scale the B and BT blocks with the factor=2./3.*tau
  for(System_per_grid& s : this->systems)
  {
    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    s.matrix.scale_blocks(2./3. * tau, cell_positions_0_2_1_2);
    if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
    {
      s.matrix.scale_blocks(2./3. * tau, cell_positions_2_0_2_1);
    }
  }
  // set the time disc back to BDF2 for the computation of the next
  // time steps
  TDatabase::TimeDB->TIME_DISC = 5;
  db["time_discretization"] = 5;
  // reset the thetas for the BDF2 methods
  TDatabase::TimeDB->THETA1 = 2./3.;
  TDatabase::TimeDB->THETA2 = 0.;
  TDatabase::TimeDB->THETA3 = 0.;
  TDatabase::TimeDB->THETA4 = 2./3.;
  
  Output::print<1>("the new parameters for the BDF's schemes");
  Output::print<1>("Theta1: ", TDatabase::TimeDB->THETA1);
  Output::print<1>("Theta2: ", TDatabase::TimeDB->THETA2);
  Output::print<1>("Theta3: ", TDatabase::TimeDB->THETA3);
  Output::print<1>("Theta4: ", TDatabase::TimeDB->THETA4);
}

/**************************************************************************** */
void Time_NSE2D::bdf_assemble_rhs()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  System_per_grid& s = this->systems.front();

  // assemble right hand side
  if(TDatabase::ParamDB->DISCTYPE != SUPG)
  {
    s.rhs.reset();
    call_assembling_routine(s, TNSE2D_Rhs);
  }
  s.solution.copy_nonactive(s.rhs);

  // scale the right-hand side with 2*tau, only actives
  double factor = TDatabase::TimeDB->THETA4*tau;
  s.rhs.scaleActive(factor);

  s.Mass_Matrix.apply_scaled_submatrix(s.solution_m1, s.rhs, 2, 2, 4./3.);
  s.Mass_Matrix.apply_scaled_submatrix(s.solution_m2, s.rhs, 2, 2, -1./3.);
  s.rhs.copy_nonactive(s.solution);
  s.solution.copy_nonactive(s.rhs);
}

/**************************************************************************** */
void Time_NSE2D::assemble_rhs_supg()
{
  // declare the cell_position vectors
  // Mac compiler cannot distinguish between the
  // two variants of BlockFEMatrix::scale_blocks
  const std::vector<std::vector<size_t>>
       cell_positions_0_2_1_2 = {{0,2},{1,2}};
  const std::vector<std::vector<size_t>>
    cell_positions_2_0_2_1 = {{2,0},{2,1}};
  
  // cout<<TDatabase::TimeDB->THETA1<<endl;
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  if((TDatabase::TimeDB->TIME_DISC==1) || (current_step_ == 1))
  {
    System_per_grid& s = this->systems.front();
    s.solution.copy_nonactive(s.rhs);
    s.rhs.scaleActive(tau);

    s.Mass_Matrix.apply_scaled_submatrix(old_solution, s.rhs, 2, 2, 1.0);
    // rescale the BT blocks, because they were 
    // reassembled during nonlinear iteration
    for(System_per_grid& s : this->systems)
    {
      s.matrix.scale_blocks(TDatabase::TimeDB->THETA1*tau, cell_positions_0_2_1_2);
      if(TDatabase::ParamDB->NSTYPE==14)
        s.matrix.scale_blocks(TDatabase::TimeDB->THETA1*tau, cell_positions_2_0_2_1);
    }      
    s.rhs.copy_nonactive(s.solution);
  }
  else
  {
    bdf_assemble_rhs();
    for(System_per_grid& s : this->systems)
    {
      s.matrix.scale_blocks(TDatabase::TimeDB->THETA1*tau, cell_positions_0_2_1_2);
      if(TDatabase::ParamDB->NSTYPE==14)
        s.matrix.scale_blocks(TDatabase::TimeDB->THETA1*tau, cell_positions_2_0_2_1);
    }
  }
}

/**************************************************************************** */
void Time_NSE2D::modify_slip_bc(bool BT_Mass)
{
  // modification of the matrices due to the 
  // slip type boundary conditions: If the mass matrices,
  // the off-diagonal A-blocks , and the BT's block,
  // are unchanged during the time iteration, then this modification
  // is done only once in the time loop. However, in the SUPG
  // and residual based VMS method these matrices are also 
  // updated during the time steps, so modification of all 
  // of them including the right-hand side is necessary.The 
  // modification of the diagonal A-blocks are necessary 
  // in any case.
  if(TDatabase::ParamDB->NSTYPE < 4)
  {
    ErrThrow("Slip with friction b.c. is only implemented for NSTYPE 4");
  }
  std::vector<const TFESpace2D*> spaces_mat(1);
  std::vector<double*> rhs_array(2);
  std::vector<const TFESpace2D*> rhs_space(2);
  
  for(System_per_grid& s: this->systems)
  {
    spaces_mat[0] = &s.velocity_space;
    rhs_space[0] = spaces_mat[0];
    rhs_space[1] = spaces_mat[0];
    
    rhs_array[0] = s.rhs.block(0);
    rhs_array[1] = s.rhs.block(1);
    
    std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix.get_blocks_uniquely();
    
    std::vector<std::shared_ptr<FEMatrix>> mass_blocks;
    
    mass_blocks = s.Mass_Matrix.get_blocks_uniquely();
    if(db["disctype"].is("residual_based_vms"))
    {
      mass_blocks = s.MatrixK.get_blocks_uniquely();
    }
    
    /*std::vector<const BoundCondFunct2D*> bc(3);
    bc[0]=s.velocity_space.GetBoundCondition();
    bc[1]=bc[0];
    bc[2]=s.pressure_space.GetBoundCondition();
    */
    BoundCondFunct2D* bc[3] = {
    s.velocity_space.GetBoundCondition(),
    s.velocity_space.GetBoundCondition(),
    s.pressure_space.GetBoundCondition()};
    
    // boundary values:
    std::vector<BoundValueFunct2D*>bv(3);
    bv[0]=example.get_bd(0);
    bv[1]=example.get_bd(1);
    bv[2]=example.get_bd(2);
    
    std::vector<TSquareMatrix2D*> sqMat;
    std::vector<TMatrix2D*> reMat;
    sqMat.resize(5);
    sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
    sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());  
    sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());

     if(db["disctype"].is("residual_based_vms"))
     {
       sqMat.resize(8);
       sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
       sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(1).get());
       sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(3).get());  
       sqMat[7] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());
     }
     if(BT_Mass)
       sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());

     reMat.resize(0);
     if(BT_Mass)
     {
       reMat.resize(2);    
       reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
       reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
     }
     // update the matrices and right hand side
     Assemble2DSlipBC(spaces_mat.size(), spaces_mat.data(), 
                   sqMat.size(), sqMat.data(), reMat.size(), reMat.data(), 
                   rhs_array.size(), rhs_array.data(), rhs_space.data(), 
                   bc, bv.data(), s.u.GetComponent(0), s.u.GetComponent(1));
  }

}
