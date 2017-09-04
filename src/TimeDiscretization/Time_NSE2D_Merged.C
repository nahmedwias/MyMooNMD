#include <Time_NSE2D_Merged.h>
#include <Database.h>
#include <Assemble2D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <BoundaryAssembling2D.h>
#include <GridTransfer.h>
#include <Domain.h>
#include <LocalProjection.h>

/* *************************************************************************** */
  //TODO  So far of this object only the nonlin it stuff is used - switch entirely!
ParameterDatabase get_default_TNSE2D_parameters()
{
  Output::print<5>("creating a default TNSE2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default TNSE2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TNSE2D parameter database");

  //Time_NSE2D_Merged requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  return db;
}
/* *************************************************************************** */

/**************************************************************************** */
Time_NSE2D_Merged::System_per_grid::System_per_grid(const Example_TimeNSE2D& example,
                  TCollection& coll, std::pair< int, int > order,
                  Time_NSE2D_Merged::Matrix type)
 : velocity_space(&coll, (char*)"u", (char*)"velocity space",  example.get_bc(0),
                  order.first, nullptr),
   pressure_space(&coll, (char*)"p", (char*)"pressure space", example.get_bc(2),
                  order.second, nullptr),
   matrix({&velocity_space, &velocity_space, &pressure_space}),
   mass_matrix({&velocity_space, &velocity_space}),
   rhs(matrix, true),
   solution(matrix, false),
   u(&velocity_space, (char*)"u", (char*)"u", solution.block(0),
     solution.length(0), 2),
   p(&pressure_space, (char*)"p", (char*)"p", this->solution.block(2),
     solution.length(2)),
   solution_m1(matrix, false),
   u_m1(&velocity_space, (char*)"u", (char*)"u", solution_m1.block(0),
        solution_m1.length(0), 2),
   p_old(&pressure_space, (char*)"p", (char*)"p", this->solution_m1.block(2),
     solution_m1.length(2)),
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
  mass_matrix = BlockFEMatrix::Mass_Matrix_NSE2D(velocity_space, pressure_space);

  switch(type)
  {
    case Time_NSE2D_Merged::Matrix::Type1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
      break;
    case Time_NSE2D_Merged::Matrix::Type2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
      break;
    case Time_NSE2D_Merged::Matrix::Type3:
      matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
      break;
    case Time_NSE2D_Merged::Matrix::Type4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      break;
    case Time_NSE2D_Merged::Matrix::Type14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
      break;
  }
}

/**************************************************************************** */
Time_NSE2D_Merged::Time_NSE2D_Merged(const TDomain& domain,
           const ParameterDatabase& param_db, int reference_id)
: Time_NSE2D_Merged(domain, param_db, Example_TimeNSE2D(param_db), reference_id)
{
}

/**************************************************************************** */
Time_NSE2D_Merged::Time_NSE2D_Merged(const TDomain& domain,
 const ParameterDatabase& param_db, const Example_TimeNSE2D& ex, int reference_id)
 : db(get_default_TNSE2D_parameters()), outputWriter(param_db), systems(),
   example(ex), solver(param_db), defect(), oldResiduals(),
   initial_residual(1e10), errors(10,0.), oldtau(0.0),
   time_stepping_scheme(param_db), is_rhs_and_mass_matrix_nonlinear(false)
{
  db.merge(param_db);
  this->set_parameters();

  std::pair <int,int> velo_pres_order(TDatabase::ParamDB->VELOCITY_SPACE,
                               TDatabase::ParamDB->PRESSURE_SPACE);
  this->get_velocity_pressure_orders(velo_pres_order);

  Time_NSE2D_Merged::Matrix type;
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

  // post-processing
  this->prepared_postprocessing(coll);
  // stabilization parameter
  if(TDatabase::ParamDB->NSTYPE==14)
  {
    stab_space = std::make_shared<TFESpace2D>(coll, (char*)"stab space", 
                       (char*)"stab space", example.get_bc(2), DiscP_PSpace, 0, nullptr);
    
    this->stab_param.resize(coll->GetN_Cells(), 0);
    
    stab_param_function=
   std::make_shared<TFEFunction2D>(stab_space.get(), (char*)"stab", (char*)"stab",
                                   stab_param.data(), coll->GetN_Cells());
  }

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

      TFEFunction2D * u1 = systems.back().u.GetComponent(0);
      TFEFunction2D * u2 = systems.back().u.GetComponent(1);
      u1->Interpolate(example.get_initial_cond(0));
      u2->Interpolate(example.get_initial_cond(1));
    }
    mg->initialize(matrices);
  }
  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);

  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  
  if(TDatabase::ParamDB->NSTYPE==14)
  {
    outputWriter.add_fe_function(stab_param_function.get());
    outputWriter.setCellValues(stab_param.data());
    std::string name = "param";
    outputWriter.setCellValuesName(name);
  }

  // print out the information (cells, dofs, etc)
  this->output_problem_size_info();
}

/**************************************************************************** */
void Time_NSE2D_Merged::set_parameters()
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

  // set the discretization parameters
  // standard Galerkin
  if(db["disctype"].is("galerkin"))
  {
    TDatabase::ParamDB->DISCTYPE = 1;
    /// set scaling factor for B, BT's block
    time_stepping_scheme.n_scale_block = 4;
    time_stepping_scheme.b_bt_linear_nl = "linear";
  }
  if(db["disctype"].is("supg"))
  {
    if(db["extrapolate_velocity"].is("constant_extrapolate") 
      || db["extrapolate_velocity"].is("linear_extrapolate") )
    {
      TDatabase::ParamDB->DISCTYPE = -2;
    }
    else
      TDatabase::ParamDB->DISCTYPE = 2;
    /// set scaling factor for B, BT's block
    // depends on how to deal the nonlinearity in the 
    // test function: fully implicit case
    time_stepping_scheme.b_bt_linear_nl = "nonlinear";
    time_stepping_scheme.n_scale_block = 2;
    if(TDatabase::ParamDB->NSTYPE==14)
      time_stepping_scheme.n_scale_block = 5;
  }
  if(db["disctype"].is("residual_based_vms"))
  {
    TDatabase::ParamDB->DISCTYPE = -101;
    // here the explicit version is used to handle 
    // the nonlinearity due to tests
    time_stepping_scheme.b_bt_linear_nl = "solution_dependent";
    time_stepping_scheme.n_scale_block = 4;
    if(TDatabase::ParamDB->NSTYPE==14)
      time_stepping_scheme.n_scale_block = 5;
  }
  
  // Smagorinsky
  if(db["disctype"].is("smagorinsky"))
  {
    TDatabase::ParamDB->DISCTYPE = 4;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE =1;
  }

  if(db["disctype"].is("local_projection"))
  {
     TDatabase::ParamDB->DISCTYPE = 1;
  }

  // the only case where one have to re-assemble the right hand side
  if(db["disctype"].is("supg"))
  {
    if(!db["extrapolate_velocity"].is("constant_extrapolate")
        && !db["extrapolate_velocity"].is("linear_extrapolate")  )
    {       
       is_rhs_and_mass_matrix_nonlinear = true;
    }
  }
  compute_param = false;
}

/**************************************************************************** */
void Time_NSE2D_Merged::get_velocity_pressure_orders(std::pair< int, int > &velo_pres_order)
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
    case -1: case -2: case -3: case -4: case -5: case -101:
    case 100: case 201: case 302:
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
        case -1: case -2: case -3: case -4:
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
        case -11: case -12: case -13: case -14:
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
void Time_NSE2D_Merged::assemble_initial_time()
{
  for(auto &s : this->systems)
  {
    call_assembling_routine(s, TNSE2D);
    //update matrices for local projection stabilization
    if(db["disctype"].is("local_projection"))
      update_matrices_lps(s);
    // copy nonactives
    s.solution.copy_nonactive(s.rhs);

    s.solution_m1 = s.solution;
    s.solution_m2 = s.solution;
  }
  // copy the current right hand side vector to the old_rhs
  this->old_rhs = this->systems.front().rhs;
  // set the solution vectors
  this->old_solution = this->systems.front().solution;  
}

/**************************************************************************** */
void Time_NSE2D_Merged::assemble_matrices_rhs(unsigned int it_counter)
{
  // case 1. standard Galerkin case:: Mass matrix, B's blocks
  // and right hand side is linear 
  // case 2. SUPG method: all matrices and right hand sides are
  // nonlinear except B-matrices
  if(it_counter == 0 /*&& time_stepping_scheme.current_step_==1*/)
  {
    // initialize the rhs from the time discretization
    rhs_from_time_disc = this->systems.front().rhs;
    rhs_from_time_disc.reset();
    System_per_grid& s = this->systems.front();
    // only assembles the right-hand side
    if(db["disctype"].is("residual_based_vms"))
    {
      for(System_per_grid& sys : this->systems)
        call_assembling_routine(sys, LocalAssembling2D_type::TNSE2D_Rhs);
    }
    else
      call_assembling_routine(s, LocalAssembling2D_type::TNSE2D_Rhs);
    //BEGIN DEBUG
    // s.rhs.print("rhs");
    //END DEBUG
    // copy the non active to the solution vector
    // since the rhs vector will be passed to the solver
    // and is modified with matrix vector multiplication
    // which also uses the non-actives
    s.solution.copy_nonactive(s.rhs);
    // copy the right hand side to the "rhs_from_time_disc"
    rhs_from_time_disc = s.rhs;
    // all matrices from the previous time step are available
    unsigned int n_sols = time_stepping_scheme.n_old_solutions();
    std::vector<BlockVector> oldsolutions(n_sols);
    oldsolutions[0] = s.solution_m1;
    if(oldsolutions.size() == 2)
      oldsolutions[1] = s.solution_m2;
    // one needs two right hand sides only for the crank-Nicolson
    // and fractional step theta schemes
    std::vector<BlockVector> rhs_(2);
    rhs_[0] = rhs_from_time_disc; // current rhs
    rhs_[1] = old_rhs; // old right hand side is needed for the Crank-Nicolson time stepping
    // modification of the matrices due to slip b.c
    if((time_stepping_scheme.current_step_ == 1 || db["disctype"].is("residual_based_vms")) 
        && TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1)
        this->modify_slip_bc(true, true);
    //NOTE: scale the B blocks only at the first iteration
    for(System_per_grid& sys : this->systems)
      time_stepping_scheme.scale_descale_all_b_blocks(sys.matrix, "scale");
    //BEGIN DEBUG
     // s.matrix.get_blocks().at(2)->Print("B1T");exit(0);
    //END DEBUG
    // prepare the right hand side for the solver
    time_stepping_scheme.prepare_rhs_from_time_disc(s.matrix, s.mass_matrix,
                     rhs_, oldsolutions);
    rhs_from_time_disc=rhs_[0];
    old_rhs=s.rhs;
    // copy the non-actives
    rhs_from_time_disc.copy_nonactive(s.solution);
    //BEGIN DEBUG
    if(TDatabase::TimeDB->CURRENTTIME==-0.0025)
    {
      s.matrix.get_blocks().at(7)->Print("B1T");exit(0);
    }
    //END DEBUG
  }
  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(this->systems.size() > 1)
    this->restrict_function();
  // assemble the nonlinear matrices
  compute_param = true;
  for(System_per_grid & s : systems)
  {
    call_assembling_routine(s, LocalAssembling2D_type::TNSE2D_NL);
    compute_param = false;
    //BEGIN DEBUG
    // cout <<db["time_discretization"]<<endl;
    // s.solution.print("s");
    // s.matrix.get_blocks().at(6)->Print("B1T");exit(0);
    //END DEBUG
    if(db["disctype"].is("local_projection"))
      update_matrices_lps(s);
  }
  // slip boundary modification of matrices
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1 )
  {
    if(is_rhs_and_mass_matrix_nonlinear)
      this->modify_slip_bc(true, true);
    else if(db["disctype"].is("residual_based_vms"))
      this->modify_slip_bc(false, true);
    else
      this->modify_slip_bc();
  }
  // prepare the right-hand side if it is nonlinear, assembling already 
  // done together with the matrices
  if(is_rhs_and_mass_matrix_nonlinear)
  {
    this->assemble_rhs_nonlinear();
  }
  //
  // also prepare the system matrices for the solver
  for(System_per_grid& s : this->systems)
  {
    // call the preparing method
    time_stepping_scheme.prepare_system_matrix(s.matrix, s.mass_matrix, it_counter);
    if(db["disctype"].is("supg") || db["disctype"].is("residual_based_vms"))
      time_stepping_scheme.scale_nl_b_blocks(s.matrix);
  }
  Output::print<5>("Assembling of matrices and right hand side is done");
}

/**************************************************************************** */
void Time_NSE2D_Merged::restrict_function()
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

/**************************************************************************** */
bool Time_NSE2D_Merged::stopIte(unsigned int it_counter)
{
  //TODO This has no "slow convergence criterion yet!"
  System_per_grid& s = this->systems.front();
  unsigned int nuDof = s.solution.length(0);
  unsigned int npDof = s.solution.length(2);
  // unsigned int sc_minit = db["nonlinloop_minit"];
  //BEGIN DEBUG
  // s.matrix.get_blocks().at(2)->Print("A1");
  //END DEBUG

  this->defect = rhs_from_time_disc;
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
  if(it_counter > 0)
  {
    Output::print("nonlinear step ", setw(3), it_counter , setw(14), impulse_residual,
                setw(14), mass_residual, setw(14), sqrt(residual),
                setw(14), sqrt(residual)/oldResidual);
  }
  else
  {
    Output::print("nonlinear step " , setw(3), it_counter, setw(14), impulse_residual,
                setw(14), mass_residual, setw(14), sqrt(residual));
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
     // project the solution
     s.p.project_into_L20();

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
       for(System_per_grid & s : this->systems)
       {
         time_stepping_scheme.reset_linear_matrices(s.matrix, s.mass_matrix);
         // descale if it's rescaled at the next time step for bdf schemes
         time_stepping_scheme.scale_descale_all_b_blocks(s.matrix, "descale");
       }
       return true;
     }
   }
   else
     return false;
}

/**************************************************************************** */
void Time_NSE2D_Merged::solve()
{
  System_per_grid& s = this->systems.front();
  
  solver.solve(s.matrix, rhs_from_time_disc, s.solution);
    if(TDatabase::TimeDB->CURRENTTIME==-0.0025)
    {
      Output::print(db["disctype"]);
      systems.front().matrix.get_combined_matrix()->Print("B1T");exit(0);
    }

  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11 and A22 matrices are
  // reset and assembled again but the A12 and A21 are scaled, so
  // for the next iteration we have to descale, see assemble_system()
  for(System_per_grid & s : this->systems)
    time_stepping_scheme.reset_linear_matrices(s.matrix, s.mass_matrix);
}

/**************************************************************************** */
void Time_NSE2D_Merged::output(int m)
{
  System_per_grid& s = this->systems.front();

  bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
  if(no_output)
    return;

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
    double t=TDatabase::TimeDB->CURRENTTIME;
    if(db["example"].is(3))
    {
      u1->GetErrors(ExactNull,3, allderiv, 2, L2H1Errors, nullptr,
                  &aux,1, &v_sp,locerr);
      double locerrorKE = locerr[0]*locerr[0];
      u2->GetErrors(ExactNull,3, allderiv, 2, L2H1Errors, nullptr,
                  &aux,1, &v_sp,locerr);
      locerrorKE += locerr[0]*locerr[0];
      Output::print( t, " kinetic energy ", locerrorKE/2 );
    }
    else
    {
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
      
      
      Output::print<1>(t, " L2(u) : ", setprecision(10), sqrt(this->errors[1]));
      Output::print<1>(t, " H1-semi(u) : ", setprecision(10), sqrt(this->errors[3]));
      
      Output::print<1>(t, " L2(0,t,L2(u)) : ", sqrt(this->errors[0]));
      Output::print<1>(t, " L2(0,t,H1-semi(u)) : ", sqrt(this->errors[2]));
      
      s.p.GetErrors(example.get_exact(2), 3, allderiv, 2, L2H1Errors,
                    nullptr, &aux, 1, &p_sp, locerr);
      
      Output::print<1>(t, " L2(p) : ", setprecision(10), locerr[0]);
      Output::print<1>(t, " H1-semi(p)) : " , setprecision(10), locerr[1] );
      
      errors[4] += (locerr[0]*locerr[0] + this->errors[5])*tau*0.5;
      errors[5] = locerr[0]*locerr[0];
      Output::print<1>(t, " L2(0,t,L2(p)) : ", sqrt(errors[4]) );
      
      errors[6] += (locerr[1]*locerr[1] + this->errors[7])*tau*0.5;
      errors[7] = locerr[1]*locerr[1];
      Output::print<1>(t, " L2(0,t,H1-semi(p)) : ", sqrt(errors[6]) );
    }
  }

  
  int n= s.solution.length(0);
  double *sol = s.solution.get_entries();
  StreamFunction(&s.velocity_space, sol,sol+n,
                     stream_function_space.get(), psi.data());
  if(db["example"].is(3))// mixing layer example
  {
    ComputeVorticityDivergence(&s.velocity_space,u1, u2, vorticity_space.get(),
                               vorticity_funct->GetValues(), divergence->GetValues());
    example.do_post_processing(*this, zero_vorticity);
  }
  else
  {
    double temp=0;
    example.do_post_processing(*this, temp);
  }

  if(time_stepping_scheme.current_step_ % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
  {
    if(db["output_write_vtk"])
    {
      outputWriter.write(TDatabase::TimeDB->CURRENTTIME);
    }
  }
  delete u1;
  delete u2;
}

/**************************************************************************** */
void Time_NSE2D_Merged::output_problem_size_info() const
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
std::array< double, int(6) > Time_NSE2D_Merged::get_errors()
{
  std::array<double, int(6)> error_at_time_points;
  error_at_time_points[0] = sqrt(errors[1]); // L2 velocity error
  error_at_time_points[1] = sqrt(errors[3]); // H1 velocity error
  error_at_time_points[2] = sqrt(errors[5]); // L2 pressure error
  error_at_time_points[3] = sqrt(errors[7]); // H1 pressure error

  return error_at_time_points;
}

/**************************************************************************** */
void Time_NSE2D_Merged::call_assembling_routine(Time_NSE2D_Merged::System_per_grid& s,
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

  // local assembling settings
  LocalAssembling2D la(type, fefunctios.data(),
                           this->example.get_coeffs());

  double *temp;
  if(compute_param)
    temp = stab_param.data();
  else
    temp = nullptr;
  // assemble all the matrices and ru_coight hand side
  Assemble2D(spaces_mat.size(), spaces_mat.data(),
               sqMatrices.size(), sqMatrices.data(),
               rectMatrices.size(), rectMatrices.data(),
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
               bc, bv.data(), la, temp);

  // add boundary integrals to stabilize backflow at open boundaries
  const TFESpace2D * v_space = &s.velocity_space;
  std::vector< TFEFunction2D* > u_conv;
  u_conv.resize(2);
  u_conv[0] = s.u.GetComponent(0);
  u_conv[1] = s.u.GetComponent(1);
  for (int k=0;k<TDatabase::ParamDB->n_stab_backflow_boundary;k++)
  {
    Output::print(" Backflow stab on boundary ",
    		  TDatabase::ParamDB->stab_backflow_boundary_id[k],
    		  " beta = ", TDatabase::ParamDB->stab_backflow_boundary_beta[k]);
    BoundaryAssembling2D boundary_integral;
    int neumann_boundary_component = TDatabase::ParamDB->stab_backflow_boundary_id[k];
    double beta_backflow_stab = TDatabase::ParamDB->stab_backflow_boundary_beta[k];
    boundary_integral.matrix_u_v_backflow_stab(s.matrix,
				 v_space,
				 u_conv,
				 neumann_boundary_component,   
				 beta_backflow_stab);     
    
  }
  
  // we are assembling only one mass matrix M11, but in general we need the
  // diagonal block M22 to be the same
  // For SUPG method, mass matrix is non-linear and will be changed during 
  // the nonlinear assembling, therefor also copy the M11 to M22.
  if((!db["disctype"].is("residual_based_vms") && (time_stepping_scheme.current_step_ == 0) )
    || db["disctype"].is("supg"))
  {
    s.mass_matrix.replace_blocks(*s.mass_matrix.get_blocks().at(0).get(), {{1,1}}, {false});
  }
}

/**************************************************************************** */
void Time_NSE2D_Merged::set_matrices_rhs(Time_NSE2D_Merged::System_per_grid& s,
                   LocalAssembling2D_type type, std::vector< TSquareMatrix2D* >& sqMat,
                   std::vector< TMatrix2D* >& reMat, std::vector< double* >& rhs_array)
{
  rhs_array.resize(0);
  sqMat.resize(0);
  reMat.resize(0);

  std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix.get_blocks_uniquely();
        // get the blocks of the mass matrix
  //std::vector<std::shared_ptr<FEMatrix>> mass_blocks
  //      = s.mass_matrix.get_blocks_uniquely();
  std::vector<std::shared_ptr<FEMatrix>> mass_blocks
              = s.mass_matrix.get_blocks_uniquely(true);
  switch(type)
  {
    case TNSE2D:
    {
      // right hand side: for NSTYPE: 1,2 and 3, size is 2
      rhs_array.resize(2);
      rhs_array[0] = s.rhs.block(0);
      rhs_array[1] = s.rhs.block(1);
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
          if(db["disctype"].is("residual_based_vms"))
          {
            sqMat.resize(8);
            sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
            sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(1).get());
            sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(3).get());
            sqMat[7] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());

            if(TDatabase::ParamDB->NSTYPE == 14)
            {// C block
              sqMat.resize(9);
              sqMat[8] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
            }
          }
          else
          {
            sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
            if(TDatabase::ParamDB->NSTYPE == 14)
            {// C block
              sqMat.resize(6);
              sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
            }
          }
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
    }// case TNSE2D
// case TNSE2D       
    case TNSE2D_NL:
    {
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
          if(db["disctype"].is("smagorinsky"))
          {
            sqMat.resize(4);
            sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
            sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
          }  
          // In the case of SUPG: together with the other contributions to 
          // the viscous and nonlinear terms, additional Mass matrix, BT-block, 
          // and right-had-side needs to be assembled during the nonlinear 
          // iteration due to the weighted test function.
          if(db["disctype"].is("supg"))
          {
            sqMat.resize(3);
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
            reMat.resize(2); 
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); 
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
            rhs_array.resize(2);
            rhs_array[0] = s.rhs.block(0);
            rhs_array[1] = s.rhs.block(1);
            s.rhs.reset(); // reset to zero
          }
          if(db["disctype"].is("residual_based_vms"))
          {
            sqMat.resize(4);
            sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
            sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
            rhs_array.resize(0);
          }
          break;
        case 14:
          if(!db["disctype"].is("supg") && !db["disctype"].is("residual_based_vms"))
          {
            ErrThrow("NSTYPE 14 only supports SUPG or RBVMS");
          }
          // we need to re-assemble all the matrices due to the solution
          // dependency of the stabilization parameters
          sqMat.resize(4);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
          sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
          sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
          
          if(db["disctype"].is("supg"))
          {
            sqMat.resize(6);
            sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
            // pressure-pressure block
            sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
            
            reMat.resize(4);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
            reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
            if(!db["extrapolate_velocity"].is("constant_extrapolate")
              && !db["extrapolate_velocity"].is("linear_extrapolate") )            
            {
              rhs_array.resize(3);
              rhs_array[0] = s.rhs.block(0);
              rhs_array[1] = s.rhs.block(1);
              rhs_array[2] = s.rhs.block(2);
              s.rhs.reset();
            }
          }
          if(db["disctype"].is("residual_based_vms"))
          {
            reMat.resize(2);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get());
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            rhs_array.resize(0);
          }
          break;
      }// endswitch NSTYPE
      break;
    }
//---------------------------    
    case TNSE2D_Rhs:
    {
      // no matrices to be assembled
      sqMat.resize(0);
      reMat.resize(0);
      // right hand side
      rhs_array.resize(2);
      rhs_array[0]=s.rhs.block(0);
      rhs_array[1]=s.rhs.block(1);
      if(db["disctype"].is("residual_based_vms"))
      {
        // we need to assemble the mass matrices and the B-blocks
        // due to the extrapolation of velocity. We treat the nonlinearity
        // of the tests using velocity-extrapolation.
        sqMat.resize(4);
        sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(1).get());
        sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(3).get());  
        sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());
        if(TDatabase::ParamDB->NSTYPE == 14)
        {
          // this is because the stabilization parameter depends on solution
          // vector which is also treated by extrapolation
          sqMat.resize(5);
          sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
        }
        
        reMat.resize(4);
        reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); 
        reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); 
        reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());  
      }
      
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
void Time_NSE2D_Merged::set_arrays(Time_NSE2D_Merged::System_per_grid& s,
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

  // bool is_imex = imex_scheme(0);
  functions.resize(3);
  functions[0] = s.u.GetComponent(0);
  functions[1] = s.u.GetComponent(1);
  functions[2] = &s.p;
  
  // finite element functions for the supg method: The nonlinear
  // version of the SUPG method.
  if(db["disctype"].is("supg") && TDatabase::ParamDB->NSTYPE == 14)
  {
    // part of the time derivative tested with pressue needs 
    // to be assembled together with the right-hand side (rhs3)
    if(db["time_discretization"].is("backward_euler") || 
       (time_stepping_scheme.pre_stage_bdf) )
    {
      functions.resize(4);
      functions[2] = s.u_m1.GetComponent(0);
      functions[3] = s.u_m1.GetComponent(1);
      
      if(db["extrapolate_velocity"].is("constant_extrapolate")
        || db["extrapolate_velocity"].is("linear_extrapolate")
      )
      {
        s.extrapolate_sol.reset();
        s.extrapolate_sol = s.solution_m1;
        
        functions.resize(6);
        functions[4] = s.extrapolate_u.GetComponent(0);
        functions[5] = s.extrapolate_u.GetComponent(1);
      }
      else
      {
        Output::print<5>("Fully nonlinear version of SUPG method");
      }
    }
    if(db["time_discretization"].is("bdf_two") && !(time_stepping_scheme.pre_stage_bdf))
    {
      //BEGIN DEBUG
      // cout << time_stepping_scheme.current_step_<<endl;
      // cout << "its not first second step yet "<< endl;exit(0);
      //END DEBUG
      s.combined_old_sols.reset();
      // copy and scale the solution at previous time step with factor 2
      s.combined_old_sols = s.solution_m1;
      s.combined_old_sols.scale(2.);
      // subtract with right factor the solution at pre-previous solution
      s.combined_old_sols.add_scaled(s.solution_m2, -1./2.);

      functions.resize(4);
      functions[2] = s.comb_old_u.GetComponent(0);
      functions[3] = s.comb_old_u.GetComponent(1);
      
      if(db["extrapolate_velocity"].is("constant_extrapolate"))
      {
        s.extrapolate_sol.reset();
        s.extrapolate_sol = s.solution_m1;
        
        functions.resize(6);
        functions[4] = s.extrapolate_u.GetComponent(0);
        functions[5] = s.extrapolate_u.GetComponent(1);
      }
      else if(db["extrapolate_velocity"].is("linear_extrapolate"))
      {
        s.extrapolate_sol.reset();
        s.extrapolate_sol = s.solution_m1;
        s.extrapolate_sol.scale(2.);
        s.extrapolate_sol.add_scaled(s.solution_m2, -1.);
       
        functions.resize(6);
        functions[4] = s.extrapolate_u.GetComponent(0);
        functions[5] = s.extrapolate_u.GetComponent(1);
      }
      else
      {
        Output::print<5>("Fully nonlinear version of SUPG method");
      }
    }
  }
  // finite element functions for the residual-based VMS method. In 
  // the first case the velocity in the test is treated by 
  // extrapolation
  if(db["disctype"].is("residual_based_vms"))
  {
    double tau = time_stepping_scheme.get_step_length();
    if(db["time_discretization"].is("backward_euler") 
       || (time_stepping_scheme.pre_stage_bdf))
    {
      //BEGIN DEBUG
      // cout << "time_stepping_scheme.current_step_ " 
      // << time_stepping_scheme.current_step_ << "  " << tau<< endl;
      //END DEBUG
      // for the residual in the nonlinear term
      // constant or linear extrapolation of the solution
      // is used
      s.extrapolate_sol.reset();
      s.extrapolate_sol = s.solution_m1;
      functions.resize(9);
      functions[2] = s.extrapolate_u.GetComponent(0);
      functions[3] = s.extrapolate_u.GetComponent(1);
      functions[4] = &s.p_old;
      
      //For the residual based vms scheme fill this
      //with the time derivative which is used in the residual 
      //computation 
      s.combined_old_sols.reset();
      s.combined_old_sols = s.solution;
      s.combined_old_sols.add_scaled(s.solution_m1, -1.);
      s.combined_old_sols.scale(1./tau);
      
      functions[5] = s.comb_old_u.GetComponent(0);
      functions[6] = s.comb_old_u.GetComponent(1);
      
      // also the solution from previous time steps
      // is used for the computation of right-hand
      // side fro pressure term
      s.combined_old_sols.reset();
      s.combined_old_sols = s.solution_m1;
      functions[7] = s.comb_old_u.GetComponent(0);
      functions[8] = s.comb_old_u.GetComponent(1);
    }
    // the first time step is treated by the backward Euler scheme
    if(db["time_discretization"].is("bdf_two") 
       && !(time_stepping_scheme.pre_stage_bdf))
    {
      // BDF2 consider the linear extrapolation of the 
      // test function "u"
      functions.resize(9);
      //BEGIN DEBUG
      if(time_stepping_scheme.current_step_ == 1)
      {
        cout<<"Hey here " << time_stepping_scheme.current_step_<<endl;
        exit(0);
      }
      //END DEBUG
      s.extrapolate_sol.reset();
      s.extrapolate_sol = s.solution_m1;
      s.extrapolate_sol.scale(2.);
      s.extrapolate_sol.add_scaled(s.solution_m2, -0.1);
      
      functions[2] = s.extrapolate_u.GetComponent(0);
      functions[3] = s.extrapolate_u.GetComponent(1);
      functions[4] = &s.p_old;
      
      //For the residual based vms scheme fill this
      //with the time derivative which is used in the residual 
      //computation 
      s.combined_old_sols.reset();
      s.combined_old_sols = s.solution_m1;
      s.combined_old_sols.scale(1./tau);
      s.combined_old_sols.add_scaled(s.solution_m2, -.5/tau);
      
      functions[5] = s.comb_old_u.GetComponent(0);
      functions[6] = s.comb_old_u.GetComponent(1);
      
      // also the solution from previous time steps
      // is used for the computation of right-hand
      // side fro pressure term
      s.combined_old_sols.reset();
      s.combined_old_sols = s.solution_m1;
      s.combined_old_sols.scale(2.);
      s.combined_old_sols.add_scaled(s.solution_m2, -0.5);
      functions[7] = s.comb_old_u.GetComponent(0);
      functions[8] = s.comb_old_u.GetComponent(1);
    }
  }
  
}

/**************************************************************************** */
bool Time_NSE2D_Merged::imex_scheme(bool print_info)
{

  if(!TDatabase::TimeDB->EXTRAPOLATE_VELOCITY)
    return false;
  bool interruption_condition  =
    TDatabase::TimeDB->EXTRAPOLATE_VELOCITY *
          (time_stepping_scheme.current_step_ >= 3);

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
void Time_NSE2D_Merged::modify_slip_bc(bool BT_Mass, bool slip_A_nl)
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

    std::vector<std::shared_ptr<FEMatrix>> blocks;
    blocks = s.matrix.get_blocks_uniquely();

    std::vector<std::shared_ptr<FEMatrix>> mass_blocks;
    mass_blocks = s.mass_matrix.get_blocks_uniquely(true);
    // s.mass_matrix.print_coloring_pattern("M", true);exit(0);
    //cout<<mass_blocks.size()<<endl;exit(0);
    /*std::vector<const BoundCondFunct2D*> bc(3);
    bc[0]=s.velocity_space.GetBoundCondition();
    bc[1]=bc[0];
    bc[2]=s.pressure_space.GetBoundCondition();*/
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
    sqMat.resize(2);
    // all 4 A blocks at the first time step
    // and only the first 2 within the nonlinear loop
    sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());//a11
    sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());//a22

    // if the off-diagonal are not changing within the non-linear loop
    // then dont need to assemble them again
    if(slip_A_nl)
    {
      sqMat.resize(4);
      sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());//a12
      sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());//a21
    }

    // either at the first time step
    // or every time step if M and B's are changing
    reMat.resize(0);
    if(BT_Mass)
    {
      sqMat.resize(8);
      sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
      sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());
      sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(1).get());
      sqMat[7] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(3).get());
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

/**************************************************************************** */
void Time_NSE2D_Merged::prepared_postprocessing(TCollection *coll)
{
  stream_function_space
     = std::make_shared<TFESpace2D>(coll, (char*)"stream function space",
                         (char*)"stream function space", example.get_bc(0), 1, nullptr);
  n_psi = stream_function_space->GetN_DegreesOfFreedom();
  psi.resize(n_psi, 0.);
  stream_function
     = std::make_shared<TFEFunction2D>(stream_function_space.get(), (char*)"streamfunction",
                          (char*)"streamfunction", psi.data(), n_psi);
  // add to the wrapper
  outputWriter.add_fe_function(stream_function.get());

  if(db["example"].is(3))
  {
    zero_vorticity = -4711;
    vorticity_space
      =std::make_shared<TFESpace2D>(coll, (char*)"vorticity space", (char*)"vorticity space",
                                example.get_bc(0), ContP_USpace, 1, nullptr);
    n_vort_dofs = vorticity_space->GetN_DegreesOfFreedom();
    vorticity.resize(2*n_vort_dofs, 0.);
    vorticity_funct 
      = std::make_shared<TFEFunction2D>(vorticity_space.get(), (char*)"voritcity",
                              (char*)"vorticity", vorticity.data()+n_vort_dofs, n_vort_dofs);
    divergence
      = std::make_shared<TFEFunction2D>(vorticity_space.get(), (char*)"dievergence",
                            (char*)"divergence", vorticity.data(), n_vort_dofs);
    outputWriter.add_fe_function(vorticity_funct.get());
    outputWriter.add_fe_function(divergence.get());
  }
}

void Time_NSE2D_Merged::update_matrices_lps(System_per_grid &s)
{
  std::vector<std::shared_ptr<FEMatrix>> blocks;
  blocks = s.matrix.get_blocks_uniquely();
  if(TDatabase::ParamDB->NSTYPE==3 || TDatabase::ParamDB->NSTYPE==4)
  {
    //update matrices for local projection stabilization
    std::vector< TSquareMatrix2D* > sqMat(2);
    sqMat[0]=reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    sqMat[1]=reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
    UltraLocalProjection(sqMat[0], FALSE);
    UltraLocalProjection(sqMat[1], FALSE);
  }
  else
  {
    std::vector< TSquareMatrix2D* > sqMat(1);
    sqMat[0]=reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    UltraLocalProjection(sqMat[0], FALSE);
  }
  
}

void Time_NSE2D_Merged::assemble_rhs_nonlinear()
{
  // initialize the rhs from the time discretization
  rhs_from_time_disc = this->systems.front().rhs;
  rhs_from_time_disc.reset();
  System_per_grid& s = this->systems.front();
  
  // copy the right hand side to the "rhs_from_time_disc"
  rhs_from_time_disc = s.rhs;
  // all matrices from the previous time step are available
  unsigned int n_sols = time_stepping_scheme.n_old_solutions();
  std::vector<BlockVector> oldsolutions(n_sols);
  oldsolutions[0] = s.solution_m1;
  if(oldsolutions.size() == 2)
    oldsolutions[1] = s.solution_m2;
  // one needs two right hand sides only for the crank-Nicolson
  // and fractional step theta schemes
  std::vector<BlockVector> rhs_(2);
  rhs_[0] = rhs_from_time_disc; // current rhs
  rhs_[1] = old_rhs; // old right hand side is needed for the Crank-Nicolson time stepping

  // prepare the right hand side for the solver
  time_stepping_scheme.prepare_rhs_from_time_disc(s.matrix, s.mass_matrix,
                   rhs_, oldsolutions);
  rhs_from_time_disc=rhs_[0];
  old_rhs=s.rhs;
  // copy the non-actives
  rhs_from_time_disc.copy_nonactive(s.rhs);
  s.solution.copy_nonactive(s.rhs);
  // CB SEBUG
  if(db["disctype"].is("residual_based_vms"))
  {
    ErrThrow("rhs for Disctype: ", db["disctype"] , " is linear ");
  }
  // END DEBUG
}
