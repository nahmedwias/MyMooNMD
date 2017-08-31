#include <Time_NSE3D_Merged.h>

#include <Database.h>
#include <Assemble3D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <GridTransfer.h>
#include <Domain.h>
#include <Variational_MultiScale3D.h>

/* *************************************************************************** */
  //TODO  So far of this object only the nonlin it stuff is used - switch entirely!
ParameterDatabase get_default_TNSE3D_parameters()
{
  Output::print<5>("creating a default TNSE2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default TNSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TNSE3D parameter database");

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
Time_NSE3D_Merged::System_per_grid::System_per_grid(const Example_TimeNSE3D& example,
                  TCollection& coll, std::pair< int, int > order,
                  Time_NSE3D_Merged::Matrix type
#ifdef _MPI
                  , int maxSubDomainPerDof
#endif
)
 : velocity_space(&coll, (char*)"u", (char*)"velocity space",  example.get_bc(0),
                  order.first),
   pressure_space(&coll, (char*)"p", (char*)"pressure space", example.get_bc(3),
                  order.second),
   matrix({&velocity_space, &velocity_space, &pressure_space}),
   mass_matrix({&velocity_space, &velocity_space}),
   rhs(matrix, true),
   solution(matrix, false),
   u(&velocity_space, (char*)"u", (char*)"u", solution.block(0),
     solution.length(0), 3),
   p(&pressure_space, (char*)"p", (char*)"p", this->solution.block(2),
     solution.length(3)),
   solution_m1(matrix, false),
   u_m1(&velocity_space, (char*)"u", (char*)"u", solution_m1.block(0),
        solution_m1.length(0), 3),
   p_old(&pressure_space, (char*)"p", (char*)"p", this->solution_m1.block(2),
     solution_m1.length(3)),
   solution_m2(matrix, false),
   u_m2(&velocity_space, (char*)"u", (char*)"u", solution_m2.block(0),
        solution_m2.length(0), 3),
   combined_old_sols(matrix, false),
   comb_old_u(&velocity_space, (char*)"u", (char*)"u", combined_old_sols.block(0),
        combined_old_sols.length(0), 3),
   extrapolate_sol(matrix, false),
   extrapolate_u(&velocity_space, (char*)"u", (char*)"u", extrapolate_sol.block(0),
        extrapolate_sol.length(0), 3)
{
  mass_matrix = BlockFEMatrix::Mass_Matrix_NSE3D(velocity_space, pressure_space);
      
  switch(type)
  {
    case Time_NSE3D_Merged::Matrix::Type1:
      matrix = BlockFEMatrix::NSE3D_Type1(velocity_space, pressure_space);
      break;
    case Time_NSE3D_Merged::Matrix::Type2:
      matrix = BlockFEMatrix::NSE3D_Type2(velocity_space, pressure_space);
      break;
    case Time_NSE3D_Merged::Matrix::Type3:
      matrix = BlockFEMatrix::NSE3D_Type3(velocity_space, pressure_space);
      break;
    case Time_NSE3D_Merged::Matrix::Type4:
      matrix = BlockFEMatrix::NSE3D_Type4(velocity_space, pressure_space);
      break;
    case Time_NSE3D_Merged::Matrix::Type14:
      matrix = BlockFEMatrix::NSE3D_Type14(velocity_space, pressure_space);
      break;
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }
#ifdef _MPI
  velocity_space.initialize_parallel(maxSubDomainPerDof);
  pressure_space.initialize_parallel(maxSubDomainPerDof);
  //print some information on the parallel infrastructure
  velocity_space.get_communicator().print_info();
  pressure_space.get_communicator().print_info();
#endif
}

Time_NSE3D_Merged::Time_NSE3D_Merged(std::list< TCollection* > collections_, const ParameterDatabase& param_db, 
                       const Example_TimeNSE3D& ex
#ifdef _MPI
, int maxSubDomainPerDof
#endif  
                   )
: db_(get_default_TNSE3D_parameters()), systems_(), example_(ex),
   solver_(param_db), defect_(), old_residual_(), 
   initial_residual_(1e10), errors_(), oldtau_(), time_stepping_scheme(param_db)
{
  db_.merge(param_db, false);
  this->check_and_setparameters();
  
  // get the velocity and pressure orders
  std::pair <int,int>
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                               TDatabase::ParamDB->PRESSURE_SPACE);
  this->get_velocity_pressure_orders(velocity_pressure_orders);
  
  // get matrix type
  Time_NSE3D_Merged::Matrix type;
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case  1: type = Matrix::Type1;  break;
    case  2: type = Matrix::Type2;  break;
    case  3: type = Matrix::Type3;  break;
    case  4: type = Matrix::Type4;  break;
    case 14: type = Matrix::Type14; break;
    default:
      ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
               " That NSE Block Matrix Type is unknown to class Time_NSE3D.");
  }
  // check for supg method
  if((db_["disctype"].is("supg") || db_["disctype"].is("residual_based_vms"))
    && (type != Matrix::Type14 && type !=Matrix::Type4) )
  {
    ErrThrow("The SUPG method is only implemented for NSTYPE 4 and 14");
  }
  // multigrid method 
  bool usingMultigrid = solver_.is_using_multigrid();
  TCollection *coll = collections_.front(); // finest grid collection
  // create finite element space and function, a matrix, rhs, and solution
#ifdef _MPI
  systems_.emplace_back(example_, *coll, velocity_pressure_orders, type,
                        maxSubDomainPerDof);
#else
  systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);
#endif
  // projection-based VMS, only on the finest level 
  if(db_["disctype"].is("vms_projection"))
  {
    this->set_matrices_vms(collections_.front());
  }
  // prepartaion of the multigrid 
  if(usingMultigrid)
  {
    auto mg = this->solver_.get_multigrid();
    size_t n_multigrid_levels = mg->get_n_geometric_levels();
    size_t n_grids=collections_.size();
    if(n_multigrid_levels != n_grids)
    {
      ErrThrow("Wrong number of grids for multigrid! expecting ",
               n_multigrid_levels, " geometric grids but got", n_grids,".");
    }
    
    if(mg->is_using_mdml())
    {
      // change the discretization on the coarse grids to lowest order 
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velocity_pressure_orders={-1, -4711};
      this->get_velocity_pressure_orders(velocity_pressure_orders);
    }
    else
    {
      // for standard multigrid, pop the finest collection - it was already
      // used to construct a space before the "if(usingMultigrid)" clause
      // and will not (as in mdml) be used a second time with a different discretization
      collections_.pop_front();
    }
    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    // matrix on finest grid is already constructed
    matrices.push_back(&systems_.back().matrix);
    
    for(auto coll : collections_)
    {
#ifdef _MPI
      systems_.emplace_back(example_, *coll, velocity_pressure_orders, type, 
                            maxSubDomainPerDof);
#else
      systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);
#endif
      // prepare input argument for multigrid 
      matrices.push_front(&systems_.back().matrix);
    }
    // initialize multigrid object with matrices on all levels
    mg->initialize(matrices);
  }
  this->output_problem_size_info();
  // initialize the defect of the system. It has the same structure as
  // the rhs (and as the solution)
  this->defect_.copy_structure(this->systems_.front().rhs);
  for(System_per_grid& s : this->systems_)
  {
    s.u.GetComponent(0)->Interpolate(example_.get_initial_cond(0));
    s.u.GetComponent(1)->Interpolate(example_.get_initial_cond(1));
    s.u.GetComponent(2)->Interpolate(example_.get_initial_cond(2));
    
    s.solution_m1 = s.solution;
    s.solution_m2= s.solution;
  }
}

///**************************************************************************** */
void Time_NSE3D_Merged::check_and_setparameters()
{
 // Check problem_type
 if(!db_["problem_type"].is(6))
 {
   if (db_["problem_type"].is(0))
   {
     db_["problem_type"] = 6;
   }
   else
   {
     Output::warn<2>("The parameter problem_type doesn't correspond to Time_NSE."
         "It is now reset to the correct value for Time_NSE (=6).");
     db_["problem_type"] = 6;
   }
 }

 // Tell the user he is using IMEX
 if(db_["extrapolate_velocity"].is("linear_extrapolate"))
 {
   if(solver_.is_using_multigrid())
   {
     ErrThrow("Multigrid with IMEX-scheme is not implemented yet");
   }
   else
   {
     Output::info<1>("check_and_setparameters",
                     "The IMEX scheme has been chosen as a time discretization scheme!\n");
   }
 }

 if(TDatabase::TimeDB->TIME_DISC == 0)
 {
   ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC
         << " does not supported");
   throw("TIME_DISC: 0 is not supported");
 }
 
 if(db_["disctype"].is("galerkin"))
 {
    TDatabase::ParamDB->DISCTYPE = 1;
    /// set scaling factor for B, BT's block
    time_stepping_scheme.n_scale_block = 6;
    time_stepping_scheme.b_bt_linear_nl = "linear";
 }
 
 if(db_["disctype"].is("vms_projection"))
 {
   TDatabase::ParamDB->DISCTYPE = 9;
 }
 
 if(db_["disctype"].is("supg"))
 {
   if(db_["extrapolate_velocity"].is("constant_extrapolate") 
     || db_["extrapolate_velocity"].is("linear_extrapolate") )
   {
     TDatabase::ParamDB->DISCTYPE = -2;
   }
   else
     TDatabase::ParamDB->DISCTYPE = 2;
   /// set scaling factor for B, BT's block
   // depends on how to deal the nonlinearity in the 
   // test function: fully implicit case
   time_stepping_scheme.b_bt_linear_nl = "nonlinear";
   time_stepping_scheme.n_scale_block = 3;
   if(TDatabase::ParamDB->NSTYPE==14)
     time_stepping_scheme.n_scale_block = 7;
 }
 
}

///**************************************************************************** */
void Time_NSE3D_Merged::set_matrices_vms(TCollection* coll_)
{
  // VMS projection order 
  int projection_order = db_["vms_projection_space_order"];
  if(projection_order < 0)
    projection_order = 0;
    
  // projection space
  projection_space_ = 
         std::make_shared<TFESpace3D>(coll_, (char*)"L",  
                           (char*)"vms projection space", example_.get_bc(2), 
                           DiscP_PSpace, projection_order);
    int ndofP = projection_space_->GetN_DegreesOfFreedom();
   
    // create vector for vms projection, used if small resolved scales are needed
    this->vms_small_resolved_scales.resize(6*ndofP);
    // finite element vector function for vms projection
    this->vms_small_resolved_scales_fefct = 
       std::make_shared<TFEVectFunct3D>(projection_space_.get(), (char*)"v", (char*)"v", 
                                        &vms_small_resolved_scales[0], ndofP, 6);
       
    // create the label space, used in the adaptive method   
    int n_cells = coll_->GetN_Cells();
    // initialize the piecewise constant vector
    if(projection_order == 0)
      this->label_for_local_projection.resize(n_cells, 0);
    else if(projection_order == 1)
      this->label_for_local_projection.resize(n_cells, 1);
    else
      ErrThrow("local projection space is not defined");
    
    // create fefunction for the labels such that they can be passed to the assembling routines
    label_for_local_projection_space_ = 
      std::make_shared<TFESpace3D>(coll_, (char*)"label_for_local_projection", 
                                   (char*)"label_for_local_projection", example_.get_bc(2), DiscP_PSpace, 0);
    // finite element function for local projection space
    this->label_for_local_projection_fefct = 
      std::make_shared<TFEFunction3D>(label_for_local_projection_space_.get(), (char*)"vms_local_projection_space_fefct", 
                                      (char*)"vms_local_projection_space_fefct", &label_for_local_projection[0], 
                                      n_cells);
    // matrices for the vms method needs to assemble only on the 
    // finest grid. On the coarsest grids, the SMAGORINSKY model 
    // is used and for that matrices are available on all grids:  
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1: case 2:
        ErrThrow("VMS projection cannot be supported for NSTYPE  ", 
                 TDatabase::ParamDB->NSTYPE);
        break;
      case 3: case 4: case 14:
        const TFESpace3D& velocity_space = this->get_velocity_space();
        // matrices G_tilde
        matrices_for_turb_mod.at(0) = std::make_shared<FEMatrix>(&velocity_space, projection_space_.get());
        matrices_for_turb_mod.at(1) = std::make_shared<FEMatrix>(&velocity_space, projection_space_.get());
        matrices_for_turb_mod.at(2) = std::make_shared<FEMatrix>(&velocity_space, projection_space_.get());
        // matrices G
        matrices_for_turb_mod.at(3) = std::make_shared<FEMatrix>(projection_space_.get(), &velocity_space);
        matrices_for_turb_mod.at(4) = std::make_shared<FEMatrix>(projection_space_.get(), &velocity_space);
        matrices_for_turb_mod.at(5) = std::make_shared<FEMatrix>(projection_space_.get(), &velocity_space);
        // mass matrix 
        matrices_for_turb_mod.at(6) = std::make_shared<FEMatrix>(projection_space_.get(), projection_space_.get());
        break;
    }
}

/**************************************************************************** */
void Time_NSE3D_Merged::get_velocity_pressure_orders(std::pair< int, int >& velocity_pressure_orders)
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
          Output::print<1>("Warning: The P1/P0 element pair (Q1/Q0 on hexa) is "
              " not stable. Make sure to use stabilization!");
          break;
        case 2: case 3: case 4: case 5:
        // standard conforming velocity and continuous pressure
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
      // pressure order is chosen correctly
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
  }
  TDatabase::ParamDB->PRESSURE_SPACE  = pressure_order;
  velocity_pressure_orders.second = pressure_order;

  Output::print("velocity space", setw(10), velocity_pressure_orders.first);
  Output::print("pressure space", setw(10), velocity_pressure_orders.second);
}

/**************************************************************************** */
void Time_NSE3D_Merged::assemble_initial_time()
{
  std::vector<TSquareMatrix3D*> sqMatrices;
  std::vector<TMatrix3D*>       rectMatrices;
  std::vector<double*> rhsArray(3);
  // from front to back (fine to coarse grid)
  // from front to back (fine to coarse grid)
  for(System_per_grid& s : this->systems_) 
  {
    call_assembling_routine(s, LocalAssembling3D_type::TNSE3D_LinGAL);
    /** manage dirichlet condition by copying non-actives DoFs
     * from rhs to solution of front grid (=finest grid)
     * Note: this operation can also be done inside the loop, so that
     * the s.solution is corrected on every grid. This is the case in
     * TNSE2D.
     * TODO: CHECK WHAT IS THE DIFFERENCE BETWEEN doing this on every grid
     * and doing it only on the finest grid!
     **/
    s.solution.copy_nonactive(s.rhs);    
    // copy the solution to the old solution for the residual computations
    // used in the RESIDUAL_VMS method
    s.solution_m1.copy_nonactive(s.rhs);
    s.solution_m2.copy_nonactive(s.rhs);
    if(db_["disctype"].is("vms_projection"))
    {
      std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix.get_blocks_uniquely();
      // update mass matrix of projection
      LumpMassMatrixToDiagonalMatrix3D(matrices_for_turb_mod.at(6));
      // update stiffness matrix
      VMS_ProjectionUpdateMatrices3D(blocks, matrices_for_turb_mod);
      // reset flag for projection-based VMS method such that Smagorinsky LES method
      // is used on coarser grids 
      db_["disctype"] = "smagorinsky_coarse";
      TDatabase::ParamDB->DISCTYPE = SMAGORINSKY_COARSE;
    }
  }// end for system per grid - the last system is the finer one (front)
  // reset   DISCTYPE to VMS_PROJECTION to be correct in the next assembling
  if(db_["disctype"].is("smagorinsky_coarse"))
    db_["disctype"] = "vms_projection";
  
  #ifdef _MPI
    double *u1  = this->systems_.front().solution_.block(0);
    double *u2  = this->systems_.front().solution_.block(1);
    double *u3  = this->systems_.front().solution_.block(2);
    double *p   = this->systems_.front().solution_.block(3);
    this->systems_.front().velocity_space.get_communicator().consistency_update(u1, 3);
    this->systems_.front().velocity_space.get_communicator().consistency_update(u2, 3);
    this->systems_.front().velocity_space.get_communicator().consistency_update(u3, 3);
    this->systems_.front().velocity_space.get_communicator().consistency_update(p, 3);
  #endif
  // copy the last right hand side and solution vectors to the old ones
  this->old_rhs_      = this->systems_.front().rhs;
  this->old_solution_ = this->systems_.front().solution;
  // this->systems_.front().solution_minus2 = this->old_solution_;
  this->old_rhs_w_old_sol = this->old_rhs_;
}

/**************************************************************************** */
void Time_NSE3D_Merged::call_assembling_routine(Time_NSE3D_Merged::System_per_grid& s, 
                                                LocalAssembling3D_type la_type)
{
  // spaces for matrices and right hand side, 
  // FE functions for nonlinearity, boundary conditions
  // and boundary values for assembling
  std::vector<const TFESpace3D*> spaces_mat;
  std::vector<const TFESpace3D*> spaces_rhs;
  std::vector<TFEFunction3D*> fefunctions;
  std::vector<BoundCondFunct3D*>  bound_cond;
  std::vector<BoundValueFunct3D*> bound_valu;
  
  set_arrays(s, spaces_mat, spaces_rhs, fefunctions, bound_cond, bound_valu);
  
  // prepare matrices and rhs for assembling 
  std::vector<TSquareMatrix3D*> sqMatrices;
  std::vector<TMatrix3D*> rectMatrices;
  std::vector<double*> rhs_array;
  
  prepare_matrices_rhs(s, la_type, sqMatrices, rectMatrices, rhs_array);
  
  // local assembling object - used in Assemble3D
  const LocalAssembling3D la(la_type, fefunctions.data(),this->example_.get_coeffs());
  
  // assemble all the matrices and right hand side
  Assemble3D(spaces_mat.size(), spaces_mat.data(),
             sqMatrices.size(), sqMatrices.data(),
             rectMatrices.size(), rectMatrices.data(),
             rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
             bound_cond.data(), bound_valu.data(), la);
}

/**************************************************************************** */
void Time_NSE3D_Merged::set_arrays(Time_NSE3D_Merged::System_per_grid& s, 
     std::vector< const TFESpace3D* >& spaces_mat,  std::vector< const TFESpace3D* >& spaces_rhs, 
     std::vector< TFEFunction3D* >& functions, std::vector<BoundCondFunct3D*> &bc, 
     std::vector<BoundValueFunct3D*> &bv)
{
  spaces_mat.resize(2);
  spaces_rhs.resize(3);
  // set arrays of spaces used to 
  // assemble matrices and right hand side 
  spaces_mat[0] = &s.velocity_space;
  spaces_mat[1] = &s.pressure_space;
  // spaces for right hand side
  spaces_rhs[0] = &s.velocity_space;
  spaces_rhs[1] = &s.velocity_space;
  spaces_rhs[2] = &s.velocity_space;
  
  if(TDatabase::ParamDB->NSTYPE == 14)
  {
    spaces_rhs.resize(4);
    spaces_rhs[3] = &s.pressure_space;
  }
  
  functions.resize(3);
  functions[0] = s.u.GetComponent(0);
  functions[1] = s.u.GetComponent(1);
  functions[2] = s.u.GetComponent(2);
  
  // assign boundary conditions
  bc.resize(3); 
  bc.at(0) = s.velocity_space.getBoundCondition();
  bc.at(1) = s.velocity_space.getBoundCondition();
  bc.at(2) = s.velocity_space.getBoundCondition();
  // bc.at(3) = s.pressureSpace_.getBoundCondition();
  // assign boundary values  
  bv.resize(3);
  bv.at(0) = example_.get_bd(0);
  bv.at(1) = example_.get_bd(1);
  bv.at(2) = example_.get_bd(2);
  
  if(db_["disctype"].is("vms_projection"))
  {
    spaces_mat.resize(4);
    // projection space
    spaces_mat[2]=projection_space_.get();
    // label space that indicates local projection space
    spaces_mat[3]=label_for_local_projection_space_.get();
    // append function for the labels of the local projection
    functions.resize(4);
    functions[3]=label_for_local_projection_fefct.get();
  }
  
  // finite element functions for the supg method: The nonlinear
  // version of the SUPG method.
  if(db_["disctype"].is("supg") && TDatabase::ParamDB->NSTYPE == 14)
  {
    // part of the time derivative tested with pressue needs 
    // to be assembled together with the right-hand side (rhs3)
    if(db_["time_discretization"].is("backward_euler") || 
       (time_stepping_scheme.pre_stage_bdf) )
    {
      functions.resize(4);
      functions[2] = s.u_m1.GetComponent(0);
      functions[3] = s.u_m1.GetComponent(1);
      
      if(db_["extrapolate_velocity"].is("constant_extrapolate")
        || db_["extrapolate_velocity"].is("linear_extrapolate")
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
    if(db_["time_discretization"].is("bdf_two") && !(time_stepping_scheme.pre_stage_bdf))
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
      
      if(db_["extrapolate_velocity"].is("constant_extrapolate"))
      {
        s.extrapolate_sol.reset();
        s.extrapolate_sol = s.solution_m1;
        
        functions.resize(6);
        functions[4] = s.extrapolate_u.GetComponent(0);
        functions[5] = s.extrapolate_u.GetComponent(1);
      }
      else if(db_["extrapolate_velocity"].is("linear_extrapolate"))
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
  }// endif supg method  
  
}

/**************************************************************************** */
void Time_NSE3D_Merged::prepare_matrices_rhs(Time_NSE3D_Merged::System_per_grid& s, 
     LocalAssembling3D_type la_type, std::vector< TSquareMatrix3D* >& sqMatrices, 
     std::vector< TMatrix3D* >& rectMatrices, std::vector< double* >& rhs_array)
{
  rectMatrices.resize(0);
  sqMatrices.resize(0);
  rhs_array.resize(0);
  
  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();
  int ns_type = TDatabase::ParamDB->NSTYPE;
  // get the blocks of the mass matrix 
  std::vector<std::shared_ptr<FEMatrix>> mass_blocks
         = s.mass_matrix.get_blocks_uniquely();
  switch(la_type)
  {
    // TODO: remove GAL because that will be for all, the difference
    // withing different methods could be done in LocalAssembling3D routines
    case LocalAssembling3D_type::TNSE3D_LinGAL: 
    {
      // right had side which is common for all nstypes
      rhs_array.resize(3);
      rhs_array[0]=s.rhs.block(0);
      rhs_array[1]=s.rhs.block(1);
      rhs_array[2]=s.rhs.block(2);
      
      switch(ns_type)
      {
        case 1:
          if(blocks.size() != 4)
          {
            ErrThrow("NSTYPE 1: Wrong blocks.size() ", blocks.size(), " instead of 4.");
          }
          sqMatrices.resize(2);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          // mass matrix
          sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
          
          rectMatrices.resize(3);
          rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(1).get());
          rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
          rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get());
          break;
        case 2:
          if(blocks.size() != 7)
          {
            ErrThrow("NSTYPE 2: Wrong blocks.size() ", blocks.size(), " instead of 7.");
          }
          sqMatrices.resize(2);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          // mass matrix
          sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
          
          rectMatrices.resize(6);
          rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(4).get()); //first the lying B blocks
          rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(5).get());
          rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(6).get());
          rectMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks.at(1).get()); //then the standing B blocks
          rectMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
          rectMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get());
          break;
        case 3:
          if(blocks.size() != 12)
          {
            ErrThrow("NSTYPE 3: Wrong blocks.size() ", blocks.size(), " instead of 12.");
          }
          sqMatrices.resize(10);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
          sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
          sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
          sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
          sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
          sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
          sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
          sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
          // mass matrix 
          sqMatrices[9]=reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
          
          // rectangular matrices 
          rectMatrices.resize(3);
          rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get());  // standing B blocks
          rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
          rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
          break;
        case 4: case 14:
          if(ns_type==4 && blocks.size() != 15)
          {
            ErrThrow("NSTYPE 4: Wrong blocks.size() ", blocks.size(), " instead of 15.");
          }
          if(ns_type==14 && blocks.size() != 16)
          {
            ErrThrow("NSTYPE 14: Wrong blocks.size() ", blocks.size(), " instead of 16.");
          }
          // square matrices
          sqMatrices.resize(10);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
          sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
          sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
          sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
          sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
          sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
          sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
          sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
          // mass matrix 
          sqMatrices[9] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
          // rectangular matrices 
          rectMatrices.resize(6);
          rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(12).get()); // standing B blocks
          rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(13).get());
          rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(14).get());
          rectMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get()); //than the standing B blocks
          rectMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
          rectMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
          if(db_["disctype"].is("vms_projection"))
          {
            sqMatrices.resize(11);
            sqMatrices[10]=reinterpret_cast<TSquareMatrix3D*>(
                                            matrices_for_turb_mod.at(6).get());
            rectMatrices.resize(12);
            // matrices  \tilde G_11, \tilde G_24, \tilde G_36
            rectMatrices[6]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(0).get());
            rectMatrices[7]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(1).get());
            rectMatrices[8]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(2).get());
            // matrices G_11, G_42, G_63
            rectMatrices[9]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(3).get());
            rectMatrices[10]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(4).get());
            rectMatrices[11]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(5).get());
          }
          if(ns_type ==14 )
          {
            rhs_array.resize(4);
            rhs_array[3] = s.rhs.block(3);
          }
          break;
      }// swith over different nstypes
    }// case LocalAssembling3D_type::TNSE3D_LinGAL:
    break;
    case LocalAssembling3D_type::TNSE3D_NLGAL: 
    {
      switch(ns_type)
      {
        case 1: 
          break;
        case 2:
          break;
        case 3:
          break;
        case 4: case 14:
          break;
      }// swith over different nstypes
    }// case LocalAssembling3D_type::TNSE3D_NLGAL:
    break;
    case LocalAssembling3D_type::TNSE3D_Rhs: 
    {
    }// case LocalAssembling3D_type::TNSE3D_Rhs:
    break;
  }
}


