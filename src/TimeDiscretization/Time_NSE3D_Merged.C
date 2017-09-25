#include <Time_NSE3D_Merged.h>

#include <Database.h>
#include <Assemble3D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <GridTransfer.h>
#include <Domain.h>
#include <Variational_MultiScale3D.h>
#include <MainUtilities.h>
#include <Output3D.h>

#include <sys/stat.h>

#ifdef _MPI
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

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
    double *u1  = this->systems_.front().solution.block(0);
    double *u2  = this->systems_.front().solution.block(1);
    double *u3  = this->systems_.front().solution.block(2);
    double *p   = this->systems_.front().solution.block(3);
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
void Time_NSE3D_Merged::assemble_matrices_rhs(unsigned int it_counter)
{
  if(it_counter == 0)
  {
    // initialize the rhs from the time discretization
    rhs_from_time_disc = this->systems_.front().rhs;
    rhs_from_time_disc.reset();
    System_per_grid& s = this->systems_.front();
    // only assembles the right-hand side
    call_assembling_routine(s, LocalAssembling3D_type::TNSE3D_Rhs);
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
    rhs_[1] = old_rhs_; // old right hand side is needed for the Crank-Nicolson time stepping
    //NOTE: scale the B blocks only at the first iteration
    for(System_per_grid& sys : this->systems_)
      time_stepping_scheme.scale_descale_all_b_blocks(sys.matrix, "scale");
    // prepare the right hand side for the solver
    time_stepping_scheme.prepare_rhs_from_time_disc(s.matrix, s.mass_matrix,
                     rhs_, oldsolutions);
    rhs_from_time_disc=rhs_[0];
    old_rhs_=s.rhs;
    // copy the non-actives
    rhs_from_time_disc.copy_nonactive(s.solution);
  }// endif it_counter == 0
  
  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(this->systems_.size() > 1)
    this->restrict_function();
  // assemble the nonlinear matrices
  for(System_per_grid & s : systems_)
  {
    call_assembling_routine(s, LocalAssembling3D_type::TNSE3D_LinGAL);
    //TODO: for local projection stabilization
    // if(db["disctype"].is("local_projection"))
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
  }
  // reset   DISCTYPE to VMS_PROJECTION to be correct in the next assembling
  if(db_["disctype"].is("smagorinsky_coarse"))
  {
    db_["disctype"]  = "vms_projection";
    TDatabase::ParamDB->DISCTYPE = VMS_PROJECTION;
  }
  // slip boundary modification of matrices
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1 )
  {
    ErrThrow("not yet supported");
  }
  
  // prepare the system matrices for the solver
  for(System_per_grid& s : this->systems_)
  {
    // call the preparing method
    time_stepping_scheme.prepare_system_matrix(s.matrix, s.mass_matrix, it_counter);
    if(db_["disctype"].is("supg") || db_["disctype"].is("residual_based_vms"))
      time_stepping_scheme.scale_nl_b_blocks(s.matrix);
  }
  Output::info<5>("Assemble non linear terms", "End of the assembling of the nonlinear matrix.");
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
        case 1: case 2:
          blocks = s.matrix.get_blocks_uniquely({{0,0},{1,1},{2,2}});
          sqMatrices.resize(1);
          sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
          break;
        case 3: case 4: 
            sqMatrices.resize(3);
            blocks = s.matrix.get_blocks_uniquely({{0,0},{1,1},{2,2}});
            sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
            sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
            sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
            if(db_["disctype"].is("smagorinsky") 
                || db_["disctype"].is("smagorinsky_coarse")
                || db_["disctype"].is("vms_projection") )
            {
              sqMatrices.resize(9);
              sqMatrices[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
              sqMatrices[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
              sqMatrices[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
              sqMatrices[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
              sqMatrices[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
              sqMatrices[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
              if(db_["disctype"].is("vms_projection"))
              {
                rectMatrices.resize(3);
                rectMatrices[0]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(0).get());
                rectMatrices[1]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(1).get());
                rectMatrices[2]=reinterpret_cast<TMatrix3D*>(matrices_for_turb_mod.at(2).get());
              }
            }
            if(db_["disctype"].is("supg"))
            {
              // mass matrix is nonlinear
              sqMatrices.resize(4);
              sqMatrices[3] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
              rectMatrices.resize(3);
              rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get()); //than the standing B blocks
              rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
              rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
              // right-hand side only for the full nonlinear part
              if(!db_["extrapolate_velocity"].is("constant_extrapolate") && !db_["extrapolate_velocity"].is("linear_extrapolate"))
              {
                rhs_array.resize(3);
                rhs_array[0]=s.rhs.block(0);
                rhs_array[1]=s.rhs.block(1);
                rhs_array[2]=s.rhs.block(2);              
              }
            }
          break;
        case 14:
            if(!db_["disctype"].is("supg") && !db_["disctype"].is("residual_based_vms"))
            {
              ErrThrow("NSTYPE 14 only supports SUPG or RBVMS");
            }
            // all A blocks, all B blocks, and the mass matrix needs to be
            // re-assembled during nonlinear loop
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
            
            if(db_["disctype"].is("supg"))
            {
              sqMatrices.resize(11);
              // mass matrix 
              sqMatrices[9] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
              // C-block
              sqMatrices[10] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(15).get());
              
              rectMatrices.resize(6);
              rectMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks.at(12).get()); // standing B blocks
              rectMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks.at(13).get());
              rectMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks.at(14).get());
              rectMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks.at(3).get()); //than the standing B blocks
              rectMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
              rectMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
              // right hand side only needs to be re-assembled if the fully nonlinear 
              // version has been considered
              if(!db_["extrapolate_velocity"].is("constant_extrapolate") && 
                !db_["extrapolate_velocity"].is("linear_extrapolate"))
              {
                rhs_array.resize(4);
                rhs_array[0]=s.rhs.block(0);
                rhs_array[1]=s.rhs.block(1);
                rhs_array[2]=s.rhs.block(2);
                rhs_array[3]=s.rhs.block(3);
              }
            }
            else
            { // for example equal-order and LPS
               // C-block
              sqMatrices[9] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(15).get());
            }
            
            break;
      }// swith over different nstypes
    }// case LocalAssembling3D_type::TNSE3D_NLGAL:
    break;
    case LocalAssembling3D_type::TNSE3D_Rhs: 
    {
      // right hand side 
      rhs_array.resize(3);
      rhs_array[0]= s.rhs.block(0);
      rhs_array[1]= s.rhs.block(1);
      rhs_array[2]= s.rhs.block(2);
      if(ns_type == 14) // TODO remove the case 4: no need the pressure block 
      {
        rhs_array.resize(4);
        rhs_array[3]=s.rhs.block(3);
      }
    }// case LocalAssembling3D_type::TNSE3D_Rhs:
    break;
  }
  // reset matrices and right hand sides 
  s.rhs.reset();
  for(auto mat : sqMatrices)
    mat->reset();
  for(auto remat : rectMatrices)
    remat->reset();
}

/**************************************************************************** */

void Time_NSE3D_Merged::restrict_function()
{
  for( int block = 0; block < 3 ;++block)
  {
    std::vector<const TFESpace3D*> spaces;
    std::vector<double*> u_entries;
    std::vector<size_t> u_ns_dofs;
    for(auto &s : systems_ )
    {
      spaces.push_back(&s.velocity_space);
      u_entries.push_back(s.solution.block(block));
      u_ns_dofs.push_back(s.solution.length(block));
    }
    GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
  }
}

/**************************************************************************** */
bool Time_NSE3D_Merged::stop_it(unsigned int iteration_counter)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  // stores current norm of the residual. They are normed per default in
  // the class Residuals
  const double normOfResidual        = this->get_full_residual();
  // hold the residual from up to 10 iterations ago
  const double veryOldNormOfResidual  = this->old_residual_.front().fullResidual;
  
  /* Some parameters have to be set up at the first nonlinear iteration */
  if(iteration_counter == 0)
  {
    initial_residual_ = normOfResidual;

    // saves the solution from previous time step with nonActive of current step
    //this->old_solution_ = this->systems_.front().solution_;
  }
  
  /** Parameters for stopping criteria (desired precision epsilon, max number
   *  of iteration, convergence rate, IMEX_scheme : if IMEX is used, we only
   *  need to solve the nonlinear iterations for the first time step in order
   *  to obtain both initial solution and the one after. Then, the extrapolated
   *  scheme is solved just once, as a linear system. ) */
  double epsilon    = db_["nonlinloop_epsilon"];
  size_t max_It     = db_["nonlinloop_maxit"];
  double conv_speed = db_["nonlinloop_slowfactor"];
  bool slow_conv    = false;
  System_per_grid& s=this->systems_.front();  
  
  if ( db_["nonlinloop_scale_epsilon_with_size"] )
  {
    epsilon *= sqrt(this->get_size());
    if (my_rank==0)
      Output::print("stopping tolerance for nonlinear iteration ", epsilon);
  }

  if ( normOfResidual >= conv_speed*veryOldNormOfResidual )
  {
    slow_conv = true;
  }

  // Stopping criteria
  if ( (normOfResidual <= epsilon) || (iteration_counter == max_It)
      || (slow_conv) )
  {
    for(System_per_grid& s: this->systems_)
    {
      s.solution_m2 = s.solution_m1;
      s.solution_m1 = s.solution;
    }
    this->old_solution_ = s.solution;
     
    if(slow_conv && my_rank==0)
      Output::print<1>(" SLOW !!! ", normOfResidual/veryOldNormOfResidual);

    if(my_rank==0)
    {
      Output::print("Last nonlinear iteration : ", setw(3), iteration_counter,
                    "\t\t", "Residual: ", normOfResidual,
                    "\t\t", "Reduction: ", normOfResidual/initial_residual_);
    }
    // descale the matrices, since only the diagonal A block will
    // be reassembled in the next time step
    if (this->imex_scheme(0) && iteration_counter>0)
      return true; // in these conditions, the matrix are already descaled
    else
    {
      for(System_per_grid & s : this->systems_)
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
void Time_NSE3D_Merged::compute_residuals()
{
  System_per_grid& s = this->systems_.front();
  unsigned int number_u_Dof = s.solution.length(0);
  unsigned int number_p_Dof = s.solution.length(3);

#ifdef _MPI
    //MPI: solution in consistency level 3 (TODO: maybe this is superfluous here
    // (because solution might be in level 3 consistency already)!)
    auto comms = s.matrix.get_communicators();
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(s.solution.block(bl), 3);
    }
#endif

  // copy rhs to defect and compute defect
  this->defect_ = s.rhs;
  s.matrix.apply_scaled_add(s.solution, defect_,-1.);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    IntoL20Vector3D(&defect_[3*number_u_Dof], number_p_Dof,
                      TDatabase::ParamDB->PRESSURE_SPACE);
  }

  // This is the calculation of the residual, given the defect.
  BlockVector defect_impuls({number_u_Dof,number_u_Dof,number_u_Dof});
  BlockVector defect_mass({number_p_Dof});
  //copy the entries (BlockVector offers no functionality to do this more nicely)
  for(size_t i = 0; i<3*number_u_Dof ;++i)
    defect_impuls.get_entries()[i] = defect_.get_entries()[i];
  for(size_t i =0 ; i<number_p_Dof ; ++i)
    defect_mass.get_entries()[i] = defect_.get_entries()[3*number_u_Dof + i];

#ifdef _MPI
  double impuls_residual_square = defect_impuls.norm_global({comms[0],comms[1],comms[2]});
  impuls_residual_square *= impuls_residual_square;
  double mass_residual_square = defect_mass.norm_global({comms[3]});
  mass_residual_square *= mass_residual_square;
#else
  double impuls_residual_square = defect_impuls.norm();
  impuls_residual_square *= impuls_residual_square;
  double mass_residual_square = defect_mass.norm();
  mass_residual_square *= mass_residual_square;
#endif

  Residuals current_residual(impuls_residual_square, mass_residual_square);
  old_residual_.add(current_residual);
}

/**************************************************************************** */
void Time_NSE3D_Merged::solve()
{
  System_per_grid& s = systems_.front();

#ifndef _MPI
  solver_.solve(s.matrix, s.rhs, s.solution);
#endif
  
#ifdef _MPI
  if(solver_.get_db()["solver_type"].is("direct"))
  {
    //set up a MUMPS wrapper
    MumpsWrapper mumps_wrapper(s.matrix);
    //kick off the solving process
    mumps_wrapper.solve(s.rhs, s.solution);
  }
  else
    this->solver_.solve(s.matrix, s.rhs, s.solution);
#endif
  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11, A22 and A33 matrices are
  // reset and assembled again but the non-diagonal blocks are scaled, so
  // for the next iteration we have to descale, see assemble_system()
  for(System_per_grid & s : this->systems_)
    time_stepping_scheme.reset_linear_matrices(s.matrix, s.mass_matrix);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       s.p.project_into_L20();
}

/**************************************************************************** */
void Time_NSE3D_Merged::output(int m, int& image)
{
  #ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    bool no_output = !db_["output_write_vtk"] && !db_["output_compute_errors"];
    if(no_output)
      return;

  System_per_grid& s = this->systems_.front();

  TFEFunction3D* u1 = s.u.GetComponent(0);
  TFEFunction3D* u2 = s.u.GetComponent(1);
  TFEFunction3D* u3 = s.u.GetComponent(2);

  if((size_t)db_["verbosity"]> 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    u3->PrintMinMax();
    s.p.PrintMinMax();
  }

  if((m==0) || (m/TDatabase::TimeDB->STEPS_PER_IMAGE) )
  {
    if(db_["output_write_vtk"])
    {
      // last argument in the following is domain but is never used in this class
      TOutput3D output(5, 5, 2, 1, NULL);
      output.AddFEFunction(&s.p);
      output.AddFEVectFunct(&s.u);
#ifdef _MPI
      char SubID[] = "";
      if(my_rank == 0)
        mkdir(db_["output_vtk_directory"], 0777);
      std::string dir = db_["output_vtk_directory"];
      std::string base = db_["output_basename"];
      output.Write_ParVTK(MPI_COMM_WORLD, image, SubID, dir, base);
      image++;
#else
    // Create output directory, if not already existing.
    mkdir(db_["output_vtk_directory"], 0777);
    std::string filename = db_["output_vtk_directory"];
    filename += "/" + db_["output_basename"].value_as_string();

      if(image<10) filename += ".0000";
      else if(image<100) filename += ".000";
      else if(image<1000) filename += ".00";
      else if(image<10000) filename += ".0";
      else filename += ".";
      filename += std::to_string(image) + ".vtk";
      output.WriteVtk(filename.c_str());
      image++;
#endif
    }
  }

  // Measure errors to known solution
  // if an exact solution is not known, it is usually set to be zero, so that
  // in such a case, here only integrals of the solution are computed.
  if(db_["output_compute_errors"])
  {
    double err_u1[4];  // FIXME? Of these arrays only the 2 first entries are
    double err_u2[4];  // used. But the evil GetErrors() will corrupt memory if
    double err_u3[4];  // these have not at least size 4.
    double err_p[4];

    TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    MultiIndex3D allderiv[4]= {D000, D100, D010, D001};
    const TFESpace3D *v_space = &this->get_velocity_space();
    const TFESpace3D *p_space = &this->get_pressure_space();

   double tau = TDatabase::TimeDB->TIMESTEPLENGTH;

    // Errors in velocity components and pressure
    u1 ->GetErrors(example_.get_exact(0), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &v_space, err_u1);
    u2 ->GetErrors(example_.get_exact(1), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &v_space, err_u2);
    u3 ->GetErrors(example_.get_exact(2), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &v_space, err_u3);
    s.p.GetErrors(example_.get_exact(3), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &p_space, err_p);

#ifdef _MPI
    double err_red[8]; //memory for global (across all processes) error
    double err_send[8]; //fill send buffer
    err_send[0]=err_u1[0];
    err_send[1]=err_u1[1];
    err_send[2]=err_u2[0];
    err_send[3]=err_u2[1];
    err_send[4]=err_u3[0];
    err_send[5]=err_u3[1];
    err_send[6]=err_p[0];
    err_send[7]=err_p[1];

    MPI_Allreduce(err_send, err_red, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    int j;
    for(j=0;j<8;j++)
    {//MPI: sqrt was skipped in GetErrors function - do it here globally!
      err_red[j] = sqrt(err_red[j]);
    }
    //fill the reduced errors back where they belong
    err_u1[0] = err_red[0];
    err_u1[1] = err_red[1];
    err_u2[0] = err_red[2];
    err_u2[1] = err_red[3];
    err_u3[0] = err_red[4];
    err_u3[1] = err_red[5];
    err_p[0] = err_red[6];
    err_p[1] = err_red[7];
#else
    int my_rank = 0;
#endif

    errors_.at(0) = err_u1[0]*err_u1[0] + err_u2[0]*err_u2[0] +
                    err_u3[0]*err_u3[0];  // (L2-norm)^2 for u
    errors_.at(1) = err_u1[1]*err_u1[1] + err_u2[1]*err_u2[1] +
                    err_u3[1]*err_u3[1];  // (H1-semi)^2 for u
    errors_.at(2) = err_p[0]*err_p[0];  // (L2-norm)^2 for p
    errors_.at(3) = err_p[1]*err_p[1];  // (H1-norm)^2 for p

    errors_.at(4) += (errors_[0] + this->errors_.at(5))*tau*0.5; //l2(0,t,l2)(u)
    errors_.at(5) = errors_.at(0);
    errors_.at(6) += (errors_.at(1) + errors_.at(7))*tau*0.5; // l2(0,t,h1)(u)
    errors_.at(7) = errors_.at(1);
    
    errors_.at(8) += (errors_.at(2)+errors_.at(9))*tau*0.5; //l2(0,t,l2) (p)
    errors_.at(9) = errors_.at(2);
    errors_.at(10) += (errors_.at(3) + errors_.at(11))*tau*0.5;//l2(0,t,h1)(p)
    errors_.at(11) = errors_.at(3);
    double t = TDatabase::TimeDB->CURRENTTIME;
    // print errors
    if (my_rank == 0 )
    {
      Output::print<1>("time: ", t, " L2(u)         : ", setprecision(10), sqrt(errors_.at(0)));
      Output::print<1>("time: ", t, " H1-semi(u)    : ", setprecision(10), sqrt(errors_.at(1)));
      
      Output::print<1>("time: ", t, " L2(0,t,L2(u)) : ", setprecision(10), sqrt(errors_.at(4)));
      Output::print<1>("time: ", t, " L2(0,t,H1-semi(u)) : ", 
                       setprecision(10), sqrt(errors_.at(6)));
      Output::print<1>("time: ", t, " L2(p)      : ", setprecision(10), sqrt(errors_.at(2)));
      Output::print<1>("time: ", t, " H1-semi(p)): ", setprecision(10), sqrt(errors_.at(3)));
      
      Output::print<1>("time: ", t, " L2(0,t,L2(p)) : ", setprecision(10), sqrt(errors_.at(8)) );
      Output::print<1>("time: ", t, " L2(0,t,H1-semi(p)) : ", setprecision(10), sqrt(errors_.at(10)) );
    }
  }
   delete u1;
   delete u2;
   delete u3;

   // do post-processing step depending on what the example implements, if needed
   // NOTE: TODO: change in the example class 
   // example_.do_post_processing(*this);
}

/**************************************************************************** */
bool Time_NSE3D_Merged::imex_scheme(bool print_info)
{
    //IMEX-scheme needs to get out of the iteration directly after the 1st solve()
  bool interruption_condition  = (db_["extrapolate_velocity"].is("linear_extrapolate"))*
                      (time_stepping_scheme.current_step_>=3);

  // change maximum number of nonlin_iterations to 1 in IMEX case
  if (interruption_condition)
  {
    db_["nonlinloop_maxit"] = 1;
    if(print_info) // condition is here just to print it once
      Output::info<1>("Nonlinear Loop MaxIteration",
                    "The parameter 'nonlinloop_maxit' was changed to 1."
                    " Only one non-linear iteration is done, because the IMEX scheme was chosen.\n");
  }
  return interruption_condition;
}

/**************************************************************************** */
const Residuals& Time_NSE3D_Merged::get_residuals() const
{
  return old_residual_.back();
}

/**************************************************************************** */
double Time_NSE3D_Merged::get_impulse_residual() const
{
  return old_residual_.back().impulsResidual;
}

/**************************************************************************** */
double Time_NSE3D_Merged::get_mass_residual() const
{
  return old_residual_.back().massResidual;
}

/**************************************************************************** */
double Time_NSE3D_Merged::get_full_residual() const
{
  return old_residual_.back().fullResidual;
}

/**************************************************************************** */


std::array< double, int(6) > Time_NSE3D_Merged::get_errors() const
{
  std::array<double, int(6)> error_at_time_points;
  error_at_time_points[0] = sqrt(this->errors_[0]); // L2 velocity error
  error_at_time_points[1] = sqrt(this->errors_[1]); // H1 velocity error
  error_at_time_points[2] = sqrt(this->errors_[2]); // L2 pressure error
  error_at_time_points[3] = sqrt(this->errors_[3]); // H1 pressure error

  return error_at_time_points;
}

/**************************************************************************** */
void Time_NSE3D_Merged::output_problem_size_info() const
{
  double h_min, h_max;
  TCollection * coll = this->get_velocity_space().GetCollection();
  coll->GetHminHmax(&h_min, &h_max);
#ifndef _MPI
  // print out some information about number of DoFs and mesh size
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 3 * n_u + n_p; // total number of degrees of freedom
  int nActive = this->get_velocity_space().GetN_ActiveDegrees();

  Output::print("N_Cells     : ", setw(10), coll->GetN_Cells());
  Output::print("h (min,max) : ", setw(10), h_min ," ", setw(12), h_max);
  Output::print("dof Velocity: ", setw(10), 3* n_u);
  Output::print("dof Pressure: ", setw(10), n_p   );
  Output::print("dof all     : ", setw(10), n_dof );
  Output::print("active dof  : ", setw(10), 3*nActive);
#else
  int size, my_rank;
  auto comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &size);
  
  double global_hmin = 0, global_hmax=0;
  MPI_Reduce(&h_min,&global_hmin, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&h_max,&global_hmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  
  std::vector<int> n_hal_cells(size,0);// hallo cells
  std::vector<int> n_own_cells(size,0);
  int l_hal_cells = coll->GetN_HaloCells();
  int l_own_cells = coll->GetN_OwnCells();
  // gather the halo cells to root
  MPI_Gather(&l_hal_cells, 1, MPI_INT, // send
             &n_hal_cells.at(0), 1, MPI_INT,// receive
             0, comm); // control
  // gather the own cells to root
  MPI_Gather(&l_own_cells, 1, MPI_INT, // send
             &n_own_cells.at(0), 1, MPI_INT, // receive
             0, comm); // control

  auto velocity_comm = systems_.front().velocity_space.get_communicator();
  auto pressure_comm = systems_.front().pressure_space.get_communicator();
  int n_loc_master_dof_velo = velocity_comm.GetN_Master();
  std::vector<int> n_master_dof_velo(size,0);
  MPI_Gather(&n_loc_master_dof_velo, 1, MPI_INT, // send
             &n_master_dof_velo.at(0), 1, MPI_INT, // receive
             0, comm); // control
  int n_loc_master_dof_pres = pressure_comm.GetN_Master();
  std::vector<int> n_master_dof_pres(size,0);
  MPI_Gather(&n_loc_master_dof_pres, 1, MPI_INT, // send
             &n_master_dof_pres.at(0), 1, MPI_INT, // receive
             0, comm); // control
  
  if(my_rank == 0)
  {
    Output::stat("Time_NSE3D", "information on the FE space");
    size_t n_total_cells=0; size_t n_dofs_velo = 0; size_t n_dofs_pres = 0;
    for(int i=0; i<size; ++i)
    {
      Output::dash("Process ", i, "\t n_own_cells ", n_own_cells.at(i), 
                   "\t n_halo_cells ", n_hal_cells.at(i), 
                   "\t h_min/h_max ", h_min, "/", h_max,
                   "\t n_master_dof_velo ", n_master_dof_velo.at(i), 
                   "\t n_master_dof_pres ", n_master_dof_pres.at(i));
      n_total_cells += n_own_cells.at(i);
      n_dofs_velo += n_master_dof_velo.at(i);
      n_dofs_pres += n_master_dof_pres.at(i);
    }
    Output::dash("Total number of cells:                       ", n_total_cells);
    Output::dash("Total number of velcoity degrees of freedom: ", 3*n_dofs_velo);
    Output::dash("Total number of pressure degrees of freedom: ", n_dofs_pres);
    Output::dash("Total number of degrees of freedom:          ", 3*n_dofs_velo+n_dofs_pres);
  }
#endif
}

