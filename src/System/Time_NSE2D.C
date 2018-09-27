#include <Time_NSE2D.h>
#include <Database.h>
#include <Assemble2D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <BoundaryAssembling2D.h>
#include <GridTransfer.h>
#include <Domain.h>
#include <LocalProjection.h>
#include <Hotfixglobal_AssembleNSE.h>
#include <Upwind.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>

/* *************************************************************************** */
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

  // a default solution in out database
  ParameterDatabase in_out_db = ParameterDatabase::default_solution_in_out_database();
  db.merge(in_out_db,true);
  
  // a default local assembling database
  db.merge(LocalAssembling2D::default_local_assembling_database());
  
  return db;
}
/* *************************************************************************** */

Time_NSE2D::System_per_grid::System_per_grid(const Example_TimeNSE2D& example,
                  TCollection& coll, std::pair< int, int > order,
                  Time_NSE2D::Matrix type)
 : velocity_space(new TFESpace2D(&coll, "u", "velocity space",  example.get_bc(0),
                  order.first)),
   pressure_space(new TFESpace2D(&coll, "p", "pressure space", example.get_bc(2),
                  order.second))
{
  switch(type)
  {
    case Time_NSE2D::Matrix::Type1:
      matrix = BlockFEMatrix::NSE2D_Type1(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type1(*velocity_space, *pressure_space);
      break;
    case Time_NSE2D::Matrix::Type2:
      matrix = BlockFEMatrix::NSE2D_Type2(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type2(*velocity_space, *pressure_space);
      break;
    case Time_NSE2D::Matrix::Type3:
      matrix = BlockFEMatrix::NSE2D_Type3(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type3(*velocity_space, *pressure_space);
      break;
    case Time_NSE2D::Matrix::Type4:
      matrix = BlockFEMatrix::NSE2D_Type4(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type4(*velocity_space, *pressure_space);
      break;
    case Time_NSE2D::Matrix::Type14:
      matrix = BlockFEMatrix::NSE2D_Type14(*velocity_space, *pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type4(*velocity_space, *pressure_space);
      break;
  }
  rhs = BlockVector(matrix, true);
  solution = BlockVector(matrix, false);
  u = TFEVectFunct2D(velocity_space.get(), "u", "u", solution.block(0),
    solution.length(0), 2);
  p = TFEFunction2D(pressure_space.get(),"p","p", this->solution.block(2),
    solution.length(2));
  solution_m1 = BlockVector(matrix, false);
  u_m1 = TFEVectFunct2D(velocity_space.get(), "u", "u", solution_m1.block(0),
       solution_m1.length(0), 2);
  p_m1 = TFEFunction2D(pressure_space.get(),"p","p", this->solution_m1.block(2),
    solution_m1.length(2));
  solution_m2 = BlockVector(matrix, false);
  u_m2 = TFEVectFunct2D(velocity_space.get(),"u","u", solution_m2.block(0),
                        solution_m2.length(0), 2);
  p_m2 = TFEFunction2D(pressure_space.get(),"p","p", this->solution_m2.block(2),
       solution_m2.length(2));
  combined_old_sols = BlockVector(matrix, false);
  comb_old_u = TFEVectFunct2D(velocity_space.get(),"u","u", combined_old_sols.block(0),
       combined_old_sols.length(0), 2);
  extrapolate_sol = BlockVector(matrix, false);
  extrapolate_u = TFEVectFunct2D(velocity_space.get(),"u","u", extrapolate_sol.block(0),
       extrapolate_sol.length(0), 2);
}

/** ************************************************************************ */
Time_NSE2D::System_per_grid::System_per_grid(const System_per_grid& other)
 : velocity_space(other.velocity_space), pressure_space(other.pressure_space),
   matrix(other.matrix), rhs(other.rhs), solution(other.solution)
{
  // the fe functions must be newly created, because copying would mean
  // referencing the BlockVectors in 'other'.
  u = TFEVectFunct2D(velocity_space.get(), "u", "u", solution.block(0),
                     solution.length(0), 2);
  p = TFEFunction2D(pressure_space.get(), "p", "p", solution.block(2),
                    solution.length(2));
}


/**************************************************************************** */
Time_NSE2D::Time_NSE2D(const TDomain& domain,
           const ParameterDatabase& param_db, int reference_id)
: Time_NSE2D(domain, param_db, Example_TimeNSE2D(param_db), reference_id)
{
}

/**************************************************************************** */
Time_NSE2D::Time_NSE2D(const TDomain& domain,
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
  // set the velocity and pressure spaces
  // this function returns a pair which consists of
  // velocity and pressure order
  this->get_velocity_pressure_orders(velo_pres_order);
  // determine NSE TYPE from Database TODO change that handling!
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

  // // post-processing
  this->prepared_postprocessing(coll);
  
  if(usingMultigrid)
  {
    // Construct multigrid object
    auto mg = solver.get_multigrid();
    bool mdml = mg->is_using_mdml();
    if(mdml)
    {
      // change the discretization on the coarse grids to lowest order
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velo_pres_order = {-1, -4711};
      this->get_velocity_pressure_orders(velo_pres_order);
    }
    
    // Determine coarser multigrid levels and construct them.
    int second_grid;
    int coarsest_grid = domain.get_ref_level() - mg->get_n_geometric_levels() + 1;
    if(mdml)
      //the finest grid is taken a second time in mdml
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
    // matrix on finest grid is already constructed
    matrices.push_back(&systems.back().matrix);
    for(int grid_no = second_grid; grid_no >= coarsest_grid; --grid_no)
    {
      TCollection *coll = domain.GetCollection(It_EQ, grid_no, reference_id);
      systems.emplace_back(example, *coll, velo_pres_order, type);
      // prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    // initialize the multigrid object with all the matrices on all levels
    mg->initialize(matrices);
  }
  
  // initial solution on finest grid - read-in or interpolation
  if(db["read_initial_solution"].is(true))
  {//initial solution is given
    std::string file = db["initial_solution_file"];
    Output::info("Initial Solution", "Reading initial solution from file ", file);
    systems.front().solution.read_from_file(file);
  }
  else
  {///interpolate initial condition from the example
    Output::info("Initial Solution", "Interpolating initial solution from example.");
    TFEFunction2D * u1 = this->systems.front().u.GetComponent(0);
    TFEFunction2D * u2 = this->systems.front().u.GetComponent(1);
    u1->Interpolate(example.get_initial_cond(0));
    u2->Interpolate(example.get_initial_cond(1));
  }
  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);
  this->rhs_from_time_disc.copy_structure(this->systems.front().rhs);

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
  if(db["time_discretization"].is("forward_euler"))
  {
    ErrThrow("time discretization: ", db["time_discretization"]," does not supported");
  }

  if(db["imex_scheme_"] && db["extrapolate_velocity"])
  {
    // first two steps are performed with the full nonlinear 
    db["extrapolate_velocity"] = true;
  }
  // set the discretization parameters
  // standard Galerkin
  if(db["space_discretization_type"].is("galerkin"))
  {
    space_disc_global = 1;
    /// set scaling factor for B, BT's block
    time_stepping_scheme.n_scale_block = 4;
    time_stepping_scheme.b_bt_linear_nl = "linear";
  }
  
  if(db["space_discretization_type"].is("supg"))
  {
    // supg: NOTE: only tested with BDF2 so far
    if(!db["time_discretization"].is("bdf_two"))
    {
      ErrThrow("supg method is only implemented for BDF2 time stepping scheme");
    }
    space_disc_global = 2;
    /// set scaling factor for B, BT's block
    // depends on how to deal the nonlinearity in the 
    // test function: fully implicit case
    time_stepping_scheme.b_bt_linear_nl = "nonlinear";
    time_stepping_scheme.n_scale_block = 2;
    if(TDatabase::ParamDB->NSTYPE==14)
      time_stepping_scheme.n_scale_block = 5;
  }
  
  // Smagorinsky
  if(db["space_discretization_type"].is("smagorinsky"))
  {
    space_disc_global = 4;
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->LAPLACETYPE =1;
  }

  if(db["space_discretization_type"].is("local_projection"))
  {
     space_disc_global = 14;
  }
  // the only case where one have to re-assemble the right hand side
  if(db["space_discretization_type"].is("supg") && db["time_discretization"].is("bdf_two"))
  {
    is_rhs_and_mass_matrix_nonlinear = true;
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
        case 100: case 201: case 302: case 403: case 504:          
          pressure_order = -(velocity_order%100 + 10)*10;
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
    case 100: case 201: case 302: case 403: case 504:
      // pressure order is chosen correctly
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
  if(systems.size() > 1) //using  multigrid
  {
    this->restrict_function();
  }
  for(auto &s : this->systems)
  {
    call_assembling_routine(s, LocalAssembling_type::TNSE3D_LinGAL);
    //update matrices for local projection stabilization
    if(db["space_discretization_type"].is("local_projection"))
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
void Time_NSE2D::assemble_matrices_rhs(unsigned int it_counter)
{
  // case 1. standard Galerkin case:: Mass matrix, B's blocks
  // and right hand side is linear 
  // case 2. SUPG method: all matrices and right hand sides are
  // nonlinear (NSTYPE 14) 
  if(it_counter == 0)
  {
    // initialize the rhs from the time discretization
    rhs_from_time_disc = this->systems.front().rhs;
    rhs_from_time_disc.reset();
    System_per_grid& s = this->systems.front();
    // only assembles the right-hand side
    call_assembling_routine(s, LocalAssembling_type::TNSE3D_Rhs);
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
    if((time_stepping_scheme.current_step_ == 1 ) 
      && TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1)
        this->modify_slip_bc(true, true);
    //NOTE: scale the B blocks only at the first iteration
    for(System_per_grid& sys : this->systems)
      time_stepping_scheme.scale_descale_all_b_blocks(sys.matrix, "scale");
   
    // prepare the right hand side for the solver
    time_stepping_scheme.prepare_rhs_from_time_disc(s.matrix, s.mass_matrix,
                     rhs_, oldsolutions);
    rhs_from_time_disc=rhs_[0];
    old_rhs=s.rhs;
    // copy the non-actives
    rhs_from_time_disc.copy_nonactive(s.solution);    
  }
  // assemble the nonlinear matrices
  if(systems.size() > 1)
    this->restrict_function();
  for(System_per_grid & s : systems)
  {
    call_assembling_routine(s, LocalAssembling_type::TNSE3D_NLGAL);
    // update matrices with local projection term
    if(db["space_discretization_type"].is("local_projection"))
      update_matrices_lps(s);
  }

  // slip boundary modification of matrices
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1 )
  {
    if(is_rhs_and_mass_matrix_nonlinear)
      this->modify_slip_bc(true, true);
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
    time_stepping_scheme.prepare_system_matrix(s.matrix, s.mass_matrix);
    if(db["space_discretization_type"].is("supg"))
      time_stepping_scheme.scale_nl_b_blocks(s.matrix);
  }
  Output::print<5>("Assembling of matrices and right hand side is done");
}

/**************************************************************************** */
void Time_NSE2D::call_assembling_routine(Time_NSE2D::System_per_grid& s,
                                         LocalAssembling_type type)
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
    s.velocity_space->get_boundary_condition(),
    s.velocity_space->get_boundary_condition(),
    s.pressure_space->get_boundary_condition()};
  
  // boundary values:
  std::vector<BoundValueFunct2D*>bv(3);
  bv[0]=example.get_bd(0);
  bv[1]=example.get_bd(1);
  bv[2]=example.get_bd(2);

  // local assembling settings
  LocalAssembling2D la(this->db, type, fefunctios.data(),
                       this->example.get_coeffs(),
                       this->get_space_disc_global());

  // find out if we have to do upwinding
  bool do_upwinding = false;
  {
    bool mdml =  this->solver.is_using_multigrid()
                && this->solver.get_multigrid()->is_using_mdml();
    bool on_finest_grid = &systems.front() == &s;
    do_upwinding = (db["space_discretization_type"].is("upwind")
                    || (mdml && !on_finest_grid));
    }

    if(do_upwinding)  //HOTFIX: Check the documentation!
      assemble_nse = Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION;
    else
      assemble_nse = Hotfixglobal_AssembleNSE::WITH_CONVECTION;
    
  // assemble all the matrices and ru_coight hand side
  Assemble2D(spaces_mat.size(), spaces_mat.data(),
               sqMatrices.size(), sqMatrices.data(),
               rectMatrices.size(), rectMatrices.data(),
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
               bc, bv.data(), la);

  if(do_upwinding)
  {
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        // do upwinding with one matrix
        UpwindForNavierStokes(la.GetCoeffFct(), sqMatrices[0],
                              la.get_fe_function(0),
                              la.get_fe_function(1));
        Output::print<3>("UPWINDING DONE");
        break;
      case 3:
      case 4:
      case 14:
        // do upwinding with two matrices
        Output::print<3>("UPWINDING DONE");
        UpwindForNavierStokes(la.GetCoeffFct(), sqMatrices[0],
                              la.get_fe_function(0),
                              la.get_fe_function(1));
        UpwindForNavierStokes(la.GetCoeffFct(), sqMatrices[1],
                              la.get_fe_function(0),
                              la.get_fe_function(1));
        break;
    } // endswitch
  }
  // we are assembling only one mass matrix M11, but in general we need the
  // diagonal block M22 to be the same
  // For SUPG method, mass matrix is non-linear and will be changed during 
  // the nonlinear assembling, therefor also copy the M11 to M22.
  //if(db["space_discretization_type"].is("supg"))
  //{
  //if(time_stepping_scheme.current_step_ == 1)
  //  s.mass_matrix.replace_blocks(*s.mass_matrix.get_blocks().at(0).get(), {{1,1}}, {false});    
    
  //}
}

/**************************************************************************** */
void Time_NSE2D::set_matrices_rhs(Time_NSE2D::System_per_grid& s,
                                  LocalAssembling_type type,
                                  std::vector< TSquareMatrix2D* >& sqMat,
                                  std::vector< TMatrix2D* >& reMat,
                                  std::vector< double* >& rhs_array)
{
  sqMat.resize(0);
  reMat.resize(0);
  // right hand side: for NSTYPE: 1,2 and 3, size is 2
  rhs_array.resize(3, nullptr);
  

  std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix.get_blocks_uniquely();
        // get the blocks of the mass matrix
  //std::vector<std::shared_ptr<FEMatrix>> mass_blocks
  //      = s.mass_matrix.get_blocks_uniquely();
  std::vector<std::shared_ptr<FEMatrix>> mass_blocks
              = s.mass_matrix.get_blocks_uniquely(true);
  switch(type)
  {
    case LocalAssembling_type::TNSE3D_LinGAL:
    {
      rhs_array[0] = s.rhs.block(0);
      rhs_array[1] = s.rhs.block(1);
      rhs_array[2] = s.rhs.block(2);
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
          sqMat.resize(6);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
          sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
          sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
          // mass matrix
          sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
	  sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());
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
          sqMat.resize(6);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
          sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
          sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
          
          sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
          sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());
          if(TDatabase::ParamDB->NSTYPE == 14)
          {// C block
            sqMat.resize(7);
            sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
          }
          // rectangular matrices
          reMat.resize(4);
          reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
          reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
          reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
          reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
          break;
      }
      // right hand sides are assembled for the initial time step
      // for the remaining time steps, they are assembled in
      // another function. so reset to zero here
      s.rhs.reset();
      break;
    }// case TNSE2D
// case TNSE2D       
    case LocalAssembling_type::TNSE3D_NLGAL:
    {
      // no right-hand side needs to be assembled here (with a few exceptions)
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
          if(db["space_discretization_type"].is("smagorinsky"))
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
          if(db["space_discretization_type"].is("supg"))
          {
            sqMat.resize(3);
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
            reMat.resize(2); 
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); 
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
            rhs_array[0] = s.rhs.block(0);
            rhs_array[1] = s.rhs.block(1);
            s.rhs.reset(); // reset to zero
          }          
          break;
        case 14:
          if(!db["space_discretization_type"].is("supg") || !db["space_discretization_type"].is("local_projection"))
          {
            ErrThrow("NSTYPE 14 only supports SUPG and Local Projection ");
          }
          // we need to re-assemble all the matrices due to the solution
          // dependency of the stabilization parameters
          sqMat.resize(4);
          sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
          sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
          sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
          sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
          
          if(db["space_discretization_type"].is("supg") )
          {
            sqMat.resize(6);
            sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
            // pressure-pressure block
            sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
          }
          
          reMat.resize(4);
          reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
          reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
          reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
          reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
          
          rhs_array[0] = s.rhs.block(0);
          rhs_array[1] = s.rhs.block(1);
          rhs_array[2] = s.rhs.block(2);
          s.rhs.reset();
          break;
      }// endswitch NSTYPE
      break;
    }// endswitch TNSE2D_NL
    break;
//---------------------------    
    case LocalAssembling_type::TNSE3D_Rhs:
    {
      // no matrices to be assembled
      sqMat.resize(0);
      reMat.resize(0);
      
      rhs_array[0] = s.rhs.block(0);
      rhs_array[1] = s.rhs.block(1);
      rhs_array[2] = s.rhs.block(2);
      // reset rhs to zero
      s.rhs.reset();
      break;
    }
    default:
      ErrThrow("The assembling type ", type, " is unknown to Time_NSE2D.");
      
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
  spaces_rhs.resize(3);

  spaces[0] = s.velocity_space.get();
  spaces[1] = s.pressure_space.get();

      spaces_rhs[0] = s.velocity_space.get();
      spaces_rhs[1] = s.velocity_space.get();
      spaces_rhs[2] = s.pressure_space.get();
  // standard for all methods.
  functions.resize(3);  
  functions[0] = s.u.GetComponent(0);
  functions[1] = s.u.GetComponent(1);
  functions[2] = &s.p;
  bool is_imex = imex_scheme(0);
  if(is_imex && db["extrapolate_velocity"])
  {
    if(db["space_discretization_type"].is("galerkin") || db["space_discretization_type"].is("local_projection"))
    {
      if(db["extrapolation_type"].is("constant_extrapolate"))
      {
        s.extrapolate_sol.reset();
        s.extrapolate_sol = s.solution_m1;
        
        functions[0] = s.extrapolate_u.GetComponent(0);
        functions[1] = s.extrapolate_u.GetComponent(1);
      }
      else if(db["extrapolation_type"].is("linear_extrapolate"))
      {
        s.extrapolate_sol.reset();
        s.extrapolate_sol = s.solution_m1;
        s.extrapolate_sol.scale(2.);
        s.extrapolate_sol.add_scaled(s.solution_m2, -1.);
       
        functions[0] = s.extrapolate_u.GetComponent(0);
        functions[1] = s.extrapolate_u.GetComponent(1);
      }
      else
      {
        ErrThrow("Only constant or linear extrapolation of velocity are used ", db["extrapolation_type"]);
      }
      functions[2] = &s.p;
    }
    // supg: NOTE: only tested with BDF2 so far
    if(db["space_discretization_type"].is("supg") && !db["time_discretization"].is("bdf_two"))
    {
      ErrThrow("supg method is only implemented for BDF2 time stepping scheme");
    }
    
    if(db["space_discretization_type"].is("supg"))
    {
      if(TDatabase::ParamDB->NSTYPE < 4)
      {
        ErrThrow("SUPG is for equal-order implemented so far !!!");        
      }
      functions.resize(4);
      
      s.extrapolate_sol.reset();
      s.extrapolate_sol = s.solution_m1;
      s.extrapolate_sol.scale(2.);
      s.extrapolate_sol.add_scaled(s.solution_m2, -1.);
      functions[0] = s.extrapolate_u.GetComponent(0);
      functions[1] = s.extrapolate_u.GetComponent(1);
      // combination of previous time solutions for assembling the right-hand
      // side, this is used for the pressure part.
      s.combined_old_sols.reset();
      // copy and scale the solution at previous time step with factor 2
      s.combined_old_sols = s.solution_m1;
      s.combined_old_sols.scale(2.);
      // subtract with right factor the solution at pre-previous solution
      s.combined_old_sols.add_scaled(s.solution_m2, -1./2.);
      
      functions[2] = s.comb_old_u.GetComponent(0);
      functions[3] = s.comb_old_u.GetComponent(1);
    }
  }
  else 
  {
    // For the standard methods or symmetric stabilization schemes
    // there is no time derivative involved in combination with the 
    // pressure or velocity as in the supg case. So nothing to do 
    // for those schemes.
    
    // 
    if(db["space_discretization_type"].is("supg"))
    {
      if(time_stepping_scheme.pre_stage_bdf)
      {
        functions.resize(4);
        functions[2] = s.u_m1.GetComponent(0);
        functions[3] = s.u_m1.GetComponent(1);
      }
      else
      {
        s.combined_old_sols.reset();
        // copy and scale the solution at previous time step with factor 2
        s.combined_old_sols = s.solution_m1;
        s.combined_old_sols.scale(2.);
        // subtract with right factor the solution at pre-previous solution
        s.combined_old_sols.add_scaled(s.solution_m2, -1./2.);
        
        functions.resize(4);
        functions[2] = s.comb_old_u.GetComponent(0);
        functions[3] = s.comb_old_u.GetComponent(1);
      }
    }
  }
}

/**************************************************************************** */
void Time_NSE2D::restrict_function()
{
  //assembling requires an approximate velocity solution on every grid
  for( int block = 0; block < 2 ;++block)
  {
    std::vector<const TFESpace2D*> spaces;
    std::vector<double*> u_entries;
    std::vector<size_t> u_ns_dofs;
    for(auto &s : systems )
    {
      spaces.push_back(s.velocity_space.get());
      u_entries.push_back(s.solution.block(block));
      u_ns_dofs.push_back(s.solution.length(block));
    }
    GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
  }
}

/**************************************************************************** */
bool Time_NSE2D::imex_scheme(bool print_info)
{
  if((db["imex_scheme_"] && time_stepping_scheme.current_step_ >= 3))
  {
    db["nonlinloop_maxit"] = 1;
    db["extrapolate_velocity"] = true;
    
    if(print_info) // condition is here just to print it once
      Output::info<1>("Nonlinear Loop MaxIteration",
                      "The parameter 'nonlinloop_maxit' was changed to 1."
                      " Only one non-linear iteration is done, because the IMEX scheme was chosen.\n");
    // this is typical for the problem with slip type boundary condition
    // NOTE: only tested for the mixing layer problem for the moment
    if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1)
      is_rhs_and_mass_matrix_nonlinear = true;
    if(db["space_discretization_type"].is("supg"))
      space_disc_global = -2;
    return true;
  }
  else
    return false;
}

/**************************************************************************** */
void Time_NSE2D::modify_slip_bc(bool BT_Mass, bool slip_A_nl)
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
    ErrThrow("Slip with friction b.c. is only implemented for NSTYPE 4 and 14");
  }
  std::vector<const TFESpace2D*> spaces_mat(1);
  std::vector<double*> rhs_array(2);
  std::vector<const TFESpace2D*> rhs_space(2);

  for(System_per_grid& s: this->systems)
  {
    spaces_mat[0] = s.velocity_space.get();
    rhs_space[0] = spaces_mat[0];
    rhs_space[1] = spaces_mat[0];

    rhs_array[0] = s.rhs.block(0);
    rhs_array[1] = s.rhs.block(1);

    std::vector<std::shared_ptr<FEMatrix>> blocks;
    blocks = s.matrix.get_blocks_uniquely();

    std::vector<std::shared_ptr<FEMatrix>> mass_blocks;
    mass_blocks = s.mass_matrix.get_blocks_uniquely(true);
    
    BoundCondFunct2D* bc[3] = {
    s.velocity_space->get_boundary_condition(),
    s.velocity_space->get_boundary_condition(),
    s.pressure_space->get_boundary_condition()};
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
      sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());//m11
      sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());//m22
      sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(1).get());//m12
      sqMat[7] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(3).get());//m21
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
void Time_NSE2D::prepared_postprocessing(TCollection *coll)
{
  stream_function_space
     = std::make_shared<TFESpace2D>(coll, (char*)"stream function space",
                         (char*)"stream function space", example.get_bc(0), 1);
  n_psi = stream_function_space->GetN_DegreesOfFreedom();
  psi.resize(n_psi, 0.);
  stream_function
     = std::make_shared<TFEFunction2D>(stream_function_space.get(), (char*)"streamfunction",
                          (char*)"streamfunction", psi.data(), n_psi);
  // add to the wrapper
  outputWriter.add_fe_function(stream_function.get());

  if(db["example"].is(3) || db["example"].is(6))
  {
    zero_vorticity = -4711;
    vorticity_space
      =std::make_shared<TFESpace2D>(coll, (char*)"vorticity space", (char*)"vorticity space",
                                example.get_bc(0), ContP_USpace, 1);
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

/**************************************************************************** */
void Time_NSE2D::assemble_rhs_nonlinear()
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
}

/**************************************************************************** */
void Time_NSE2D::update_matrices_lps(System_per_grid &s)
{
  std::vector<std::shared_ptr<FEMatrix>> blocks;
  blocks = s.matrix.get_blocks_uniquely();
  if(TDatabase::ParamDB->NSTYPE==3 || TDatabase::ParamDB->NSTYPE==4)
  {
    //update matrices for local projection stabilization
    std::vector< TSquareMatrix2D* > sqMat(2);
    sqMat[0]=reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    sqMat[1]=reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
    UltraLocalProjection(sqMat[0], false);
    UltraLocalProjection(sqMat[1], false);
  }
  else
  {
    std::vector< TSquareMatrix2D* > sqMat(1);
    sqMat[0]=reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    UltraLocalProjection(sqMat[0], false);
  }
  
}

/**************************************************************************** */
bool Time_NSE2D::stopIte(unsigned int it_counter)
{
  //TODO This has no "slow convergence criterion yet!"
  System_per_grid& s = this->systems.front();
  unsigned int nuDof = s.solution.length(0);
  unsigned int npDof = s.solution.length(2);
  // unsigned int sc_minit = db["nonlinloop_minit"];
 
  this->defect = rhs_from_time_disc;
  s.matrix.apply_scaled_add(s.solution, defect,-1.);  
  
  if(s.matrix.pressure_projection_enabled())
    IntoL20FEFunction(&defect[2*nuDof], npDof, &this->get_pressure_space(),
                      TDatabase::ParamDB->VELOCITY_SPACE,
                      TDatabase::ParamDB->PRESSURE_SPACE);
  double residual =  Ddot(2*nuDof+npDof, &this->defect[0], &this->defect[0]);
  double impulse_residual = Ddot(2*nuDof, &this->defect[0],
         &this->defect[0]);
  double mass_residual    = Ddot(npDof,&this->defect[2*nuDof],
         &this->defect[2*nuDof]);
  
  Output::print<3>("nonlinear step  :  " , setw(3), it_counter);
  Output::print<3>("impulse_residual:  " , setw(3), impulse_residual);
  Output::print<3>("mass_residual   :  " , setw(3), mass_residual);
  Output::print<3>("residual        :  " , setw(3), sqrt(residual));

  if (it_counter>0)
  {
    Output::print<3>("rate:           :  " , setw(3), sqrt(residual)/oldResidual);
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

     Output::print<3>("ITE : ", setw(3), it_counter, "  RES : ", sqrt(residual),
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
void Time_NSE2D::solve()
{
  System_per_grid& s = this->systems.front();
  
  solver.solve(s.matrix, rhs_from_time_disc, s.solution);   
  
  double damping = this->db["nonlinloop_damping_factor"];
  if(damping != 1.0)
  {
    s.solution.scale(damping);
    s.solution.add_scaled(this->old_solution, 1-damping);
  }

  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11 and A22 matrices are
  // reset and assembled again but the A12 and A21 are scaled, so
  // for the next iteration we have to descale, see assemble_system()
  for(System_per_grid & s : this->systems)
    time_stepping_scheme.reset_linear_matrices(s.matrix, s.mass_matrix);
  
  if(s.matrix.pressure_projection_enabled())
       s.p.project_into_L20();
}

/**************************************************************************** */
void Time_NSE2D::output(int m)
{
  System_per_grid& s = this->systems.front();

  TFEFunction2D * u1 = s.u.GetComponent(0);
  TFEFunction2D * u2 = s.u.GetComponent(1);

  if((size_t)db["verbosity"] > 1)
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
    if(db["example"].is(6))
    {
      u1->GetErrors(ExactNull,3, allderiv, 2, L2H1Errors, nullptr,
                  &aux,1, &v_sp,locerr);
      
      u2->GetErrors(ExactNull,3, allderiv, 2, L2H1Errors, nullptr,
                  &aux,1, &v_sp,locerr+2);
      
      Output::print( t, " kinetic energy ", (locerr[0]*locerr[0] + locerr[2]*locerr[2])/2 );
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
      
      Output::print<1>(t, " L2(0,t,L2(u)) : ", setprecision(10), sqrt(this->errors[0]));
      Output::print<1>(t, " L2(0,t,H1-semi(u)) : ", setprecision(10), sqrt(this->errors[2]));
      
      s.p.GetErrors(example.get_exact(2), 3, allderiv, 2, L2H1Errors,
                    nullptr, &aux, 1, &p_sp, locerr);
      
      Output::print<1>(t, " L2(p) : ", setprecision(10), locerr[0]);
      Output::print<1>(t, " H1-semi(p)) : " , setprecision(10), locerr[1] );
      
      errors[4] += (locerr[0]*locerr[0] + this->errors[5])*tau*0.5;
      errors[5] = locerr[0]*locerr[0];
      Output::print<1>(t, " L2(0,t,L2(p)) : ", setprecision(10), sqrt(errors[4]) );
      
      errors[6] += (locerr[1]*locerr[1] + this->errors[7])*tau*0.5;
      errors[7] = locerr[1]*locerr[1];
      Output::print<1>(t, " L2(0,t,H1-semi(p)) : ", setprecision(10), sqrt(errors[6]) );
    }
  }

  
  int n= s.solution.length(0);
  double *sol = s.solution.get_entries();
  StreamFunction(s.velocity_space.get(), sol,sol+n,
                    stream_function_space.get(), psi.data());
  if(db["example"].is(6))// mixing layer example
  {
    ComputeVorticityDivergence(s.velocity_space.get(),u1, u2, vorticity_space.get(),
                              vorticity_funct->GetValues(), divergence->GetValues());
    example.do_post_processing(*this, zero_vorticity);
  }
  else
  {
    double dummy = 0;
    example.do_post_processing(*this, dummy);
  }

  delete u1;
  delete u2;

  outputWriter.write(TDatabase::TimeDB->CURRENTTIME);

  if(db["write_solution_binary"].is(true))
  {
    size_t interval = db["write_solution_binary_all_n_steps"];
    if(m % interval == 0)
    {//write solution to a binary file
      std::string file = db["write_solution_binary_file"];
      if(!db["overwrite_solution_binary"]) //create a new file every time
      {
        file += ".";
        file += std::to_string(TDatabase::TimeDB->CURRENTTIME);
      }
      Output::info("output", "Writing current solution to file ", file);
      systems.front().solution.write_to_file(file);
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


double Time_NSE2D::getFullResidual() const
{
  return this->oldResidual;
}
