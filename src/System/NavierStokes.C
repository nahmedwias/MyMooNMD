#include "NavierStokes.h"
#include "LocalAssembling.h"
#include "Database.h"
#include "Multigrid.h"
#include "GridTransfer.h"
#include "NSE_local_assembling_routines.h"
#include "LPS_scott_zhang.h"
#ifdef __3D__
#include "Assemble3D.h"
#include "SquareMatrix3D.h"
#include "Upwind3D.h"
#include "AuxParam3D.h"
#include "BoundaryAssembling3D.h"
#else
#include "Assemble2D.h"
#include "SquareMatrix2D.h"
#include "Upwind.h"
#include "AuxParam2D.h"
#include "BoundaryAssembling2D.h"
#endif
#ifdef _MPI
#include "ParFECommunicator3D.h"
#include "MumpsWrapper.h"
#endif

template <int d>
ParameterDatabase NavierStokes<d>::default_nse_database()
{
  Output::print<5>("creating a default NSE parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("NSE parameter database");
  
  //NSE3D requires a nonlinear iteration, set up a nonlinit_database and merge
  db.merge(ParameterDatabase::default_nonlinit_database(), true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);
  
  db.merge(LocalAssembling<d>::default_local_assembling_database(), true);
  db.merge(Solver<>::default_solver_database(), true);
  
  //stokes case - reduce no nonlin its TODO remove global database dependency
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 3)
  {
     if (TDatabase::ParamDB->PRESSURE_SEPARATION==1)
     {
        db["nonlinloop_maxit"] = 1;
     }
     else
     {
       db["nonlinloop_maxit"] = 1;
     }
  }

  return db;
}

/* ************************************************************************* */
template <int d>
NavierStokes<d>::System_per_grid::System_per_grid (
  const Example_NSE& example, TCollection& coll,
  std::pair<int,int> velocity_pressure_orders)
 : velocity_space(new FESpace(&coll, "u", "Navier--Stokes velocity",
                              example.get_bc(0), 
                              velocity_pressure_orders.first)),
   pressure_space(new FESpace(&coll, "p", "Navier--Stokes pressure",
                              example.get_bc(d),
                              velocity_pressure_orders.second))
{
  // build the matrix due to NSE type
  switch (TDatabase::ParamDB->NSTYPE)
    {
#ifdef __2D__
    case 1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
      break;
    case 2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
      break;
    case 3:
      matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
      break;
    case 4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      break;
    case 14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
      break;
#else // 3D
    case 1:
      matrix = BlockFEMatrix::NSE3D_Type1(velocity_space, pressure_space);
      break;                                               
    case 2:                                                
      matrix = BlockFEMatrix::NSE3D_Type2(velocity_space, pressure_space);
      break;                                               
    case 3:                                                
      matrix = BlockFEMatrix::NSE3D_Type3(velocity_space, pressure_space);
      break;                                               
    case 4:                                                
      matrix = BlockFEMatrix::NSE3D_Type4(velocity_space, pressure_space);
      break;                                               
    case 14:                                               
      matrix = BlockFEMatrix::NSE3D_Type14(velocity_space, pressure_space);
      break;
#endif
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }
  rhs = BlockVector(matrix, true);
  solution = BlockVector(matrix, false);

  u = FEVectFunct(velocity_space, "u", "u", solution.block(0),
                  solution.length(0), d);
  p = FEFunction(pressure_space, "p", "p", solution.block(d),
                 solution.length(d));
  
#ifdef _MPI
  //print some information
  velocity_space->get_communicator().print_info();
  pressure_space->get_communicator().print_info();
#endif
}

/* ************************************************************************* */
template <int d>
NavierStokes<d>::System_per_grid::System_per_grid(const System_per_grid& other)
 : velocity_space(other.velocity_space), pressure_space(other.pressure_space),
   matrix(other.matrix), rhs(other.rhs), solution(other.solution)
{
  // the fe functions must be newly created, because copying would mean 
  // referencing the BlockVectors in 'other'.
  u = FEVectFunct(velocity_space, "u", "u", solution.block(0),
                  solution.length(0), d);
  p = FEFunction(pressure_space, "p", "p", solution.block(d),
                 solution.length(d));
}

/* ************************************************************************* */
template <int d>
NavierStokes<d>::NavierStokes(const TDomain& domain,
                              const ParameterDatabase& param_db)
 : NavierStokes<d>(domain, param_db, Example_NSE(param_db))
{
}

/* ************************************************************************* */
template <int d>
NavierStokes<d>::NavierStokes(const TDomain& domain,
                              const ParameterDatabase& param_db,
                              Example_NSE example_in)
  : systems(), example(example_in), db(default_nse_database()),
    outputWriter(param_db), solver(param_db), defect(), old_residuals(),
    initial_residual(1e10), errors()
{
  this->db.merge(param_db, false);
  this->check_parameters();

  std::pair <int,int> 
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and pressure spaces
  // this function returns a pair which consists of 
  // velocity and pressure order
  this->get_velocity_pressure_orders(velocity_pressure_orders);

  bool usingMultigrid = solver.is_using_multigrid();
  auto collections = domain.get_grid_collections();
  TCollection *coll = collections.front(); //the finest grid collection
  // create finite element space and function, a matrix, rhs, and solution
  systems.emplace_back(example, *coll, velocity_pressure_orders);
  
  if(usingMultigrid)
  {
    // Construct multigrid object
    auto mg = solver.get_multigrid();
    bool mdml = mg->is_using_mdml();

    //Check whether number of given grids is alright
    size_t n_geo_multigrid_levels = mg->get_n_geometric_levels();
    size_t n_grids = collections.size();
    if(n_geo_multigrid_levels > n_grids )
      ErrThrow("Wrong number of grids for multigrid! I was expecting ",
               n_geo_multigrid_levels, " geometric grids but only got ", n_grids,".");
    // remove not needed coarser grid from list of collections
    for(size_t i = n_geo_multigrid_levels; i < n_grids; ++i)
    {
      collections.pop_back();
    }

    if(mdml)
    {
      // change the discretization on the coarse grids to lowest order 
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velocity_pressure_orders = {-1, -4711};
      this->get_velocity_pressure_orders(velocity_pressure_orders);
    }
    else
    {
      // for standard multigrid, pop the finest collection - it was already
      // used to construct a space before the "if(usingMultigrid)" clause
      // and will not (as in mdml) be used a second time with a different discretization
      collections.pop_front();
    }
    
     // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    // matrix on finest grid is already constructed
    matrices.push_back(&systems.back().matrix);

    for(auto coll : collections) // initialize the coarse grid space hierarchy
    {
      systems.emplace_back(example, *coll, velocity_pressure_orders);
      // prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    // initialize the multigrid object with all the matrices on all levels
    mg->initialize(matrices);
  }
  
  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());

  if(db["output_write_exact_solution"])
  {
    // initialize variables for storing exact solution (FE functions)
    solution_exact = BlockVector(this->get_solution());
    u_exact = FEVectFunct(this->systems.front().velocity_space, "u_exact",
                          "u_exact", solution_exact.block(0),
                          solution_exact.length(0), d);
    p_exact = FEFunction(this->systems.front().pressure_space, "p_exact",
                         "p_exact", solution_exact.block(d),
                         solution_exact.length(d));
    // Interpolating the exact solution
    for(int i = 0; i < d; ++i)
    {
      auto ui = u_exact.GetComponent(i);
      ui->Interpolate(example.get_exact(i));
      delete ui;
    }
    p_exact.Interpolate(example.get_exact(d));
    p_exact.PrintMinMax(std::string("p_exact"));
    outputWriter.add_fe_vector_function(&u_exact);
    outputWriter.add_fe_function(&p_exact);
  }
  
  output_problem_size_info();
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::check_parameters()
{
  if(!db["problem_type"].is(3) && !db["problem_type"].is(5) &&
     !db["problem_type"].is(7))
  {
    Output::warn<2>("The parameter problem_type doesn't correspond neither to "
        "NSE nor to Stokes or Brinkman. It is now reset to the default value "
        "for NSE (=5).");
    db["problem_type"] = 5;
  }
  
  // more check needed
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::get_velocity_pressure_orders(
  std::pair<int, int>& velocity_pressure_orders)
{
  int velocity_order = velocity_pressure_orders.first;
  int pressure_order = velocity_pressure_orders.second;
  int order = 0;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5: // P_k/Q_k
      order = velocity_order;
      break;
    case 12: case 13: case 14: case 15:
      if(d == 3 && velocity_order == 12)
      {
        ErrThrow("Scott-Vogelius P2/P1disc is not stable in 3D even if the final "
                 "refinement step is barycentric.");
      }
      // P2/P1disc and P3/P2disc elements are in general not inf-sup stable on 
      // triangles. If the last refinement step was barycentric refinement,
      // then these elements are inf-sup stable.  
      Output::warn("You chose to use Scott-Vogelius finite elements, make sure "
                   "that the final refinement step is barycentric, see the "
                   "parameter 'refinement_final_step_barycentric'.");
      order = velocity_order-10;
      break;
    case -1: case -2: case -3: case -4: case -5:
    case -101:
      order = velocity_order;
      break;
    // conforming fe spaces with bubbles on triangles
    case 22: case 23: case 24: case 25:
      order = velocity_order;
      if(d == 3 && velocity_order == 25)
      {
        ErrThrow("Velocity space 25 not supported in 3D");
      }
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
          Output::warn("NSE3D", "The P1/P0 element pair (Q1/Q0 on hexa) is "
                       " not stable. Make sure to use stabilization!");
          break;
        case 2: case 3: case 4: case 5:
        // standard conforming velo and continuous pressure
          pressure_order = velocity_order-1;
          break;
          // Scott-Vogelius: discontinuous pressure spaces with standard 
          // conforming velocity space. This is not stable on general triangles,
          // be sure to use barycentric refinement.
        case 12: case 13: case 14: case 15:
          pressure_order = -velocity_order+1;
          break;
        case 22: case 23: case 24: case 25:
          pressure_order = -(velocity_order-11)*10;
          break;
      }
      break;
    // continuous pressure spaces
    case 1: case 2: case 3: case 4: case 5:
      // nothing to do
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
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::output_problem_size_info() const
{
  int my_rank = 0;
#ifndef _MPI
  auto & velocity_space = *this->systems.front().velocity_space;
  auto & pressure_space = *this->systems.front().pressure_space;

  size_t nDofu  = velocity_space.GetN_DegreesOfFreedom();
  size_t nDofp  = pressure_space.GetN_DegreesOfFreedom();
  size_t nTotal = d*nDofu + nDofp;
  size_t nActive= d*velocity_space.GetActiveBound();

  auto coll = velocity_space.GetCollection();

  double hmin, hmax;
  coll->GetHminHmax(&hmin, &hmax);
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  auto velocity_comm = systems.front().velocity_space->get_communicator();
  auto pressure_comm = systems.front().pressure_space->get_communicator();
  int nDofu  = velocity_comm.get_n_global_dof();
  int nDofp  = pressure_comm.get_n_global_dof();
  int nTotal = d*nDofu + nDofp;
#endif
  if(my_rank ==0)
  {
    Output::stat("NavierStokes", "Mesh data and problem size on finest grid");
#ifndef _MPI
    Output::dash("N_Cells            :  ", setw(10), coll->GetN_Cells());
    Output::dash("h(min, max)        :  ", setw(10), hmin, setw(10), " ", hmax);
#endif
    Output::dash("dof velocity       :  ", setw(10), d*nDofu);
#ifndef _MPI
    Output::dash("dof velocity active:  ", setw(10), nActive);
#endif
    Output::dash("dof pressure       :  ", setw(10), nDofp);
    Output::dash("dof total          :  ", setw(10), nTotal);
  }
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::assemble_linear_terms()
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  int n_fe_spaces = 2; // spaces used for assembling matrices
  int n_square_matrices=d*d + 1; // maximum no of square matrices (type 14)
  std::vector<SquareMatrixD*> sq_matrices(n_square_matrices, nullptr);
  int n_rectangular_matrices = 2*d; // maximum no of rectangular matrices
  std::vector<MatrixD*> re_matrices(n_rectangular_matrices, nullptr);
  constexpr int n_rhs = d+1; // maximum number of right hand sides
  std::array<double*, n_rhs> rhs_array; // right hand side 
  // finite element function used for nonlinear term
  std::array<FEFunction*, d+1> feFunction;
  // boundary conditions and boundary values
  std::array<BoundaryConditionFunction*, d+1> boundCondition;
  std::array<BoundaryValuesFunction*, d+1> boundValues;
  for(int i = 0; i < d+1; ++i)
    boundValues[i] = example.get_bd()[i];
  
  std::array<const FESpace*, d+1> rhs_spaces;
  
  for(auto &s : this->systems)
  {
    s.rhs.reset(); // right hand side reset (is that necessary?)
    s.matrix.reset(); // reset matrix (needed for mdml where this is called)
    
    const FESpace *v_space = s.velocity_space.get();
    const FESpace *p_space = s.pressure_space.get();

    // spaces for matrices
    const FESpace *spaces[2] = {v_space, p_space};
    for(int i = 0; i < d; ++i)
      rhs_spaces[i] = v_space;
    rhs_spaces[d] = p_space;

    // spaces for right hand side
    for(int i = 0; i <= d; ++i)
      rhs_array[i] = s.rhs.block(i);
    
    std::vector<std::shared_ptr<FEMatrix>> blocks 
      = s.matrix.get_blocks_uniquely();
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[1+i].get());
        break;
      case 2:
        sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[d+1+i].get());
        for(int i = 0; i < d; ++i)
          re_matrices[i+d] = reinterpret_cast<MatrixD*>(blocks[1+i].get());
        break;
      case 3:
        for(int i = 0, j = 0; i < d*d; ++i, ++j)
        {
          if(i%d == 0 && i > 0)
            j++;
          sq_matrices[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
        }
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[d+(d+1)*i].get());
        break;
      case 4:
        for(int i = 0, j = 0; i < d*d; ++i, ++j)
        {
          if(i%d == 0 && i > 0)
            j++;
          sq_matrices[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
        }
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[d*(d+1)+i].get());
        for(int i = 0; i < d; ++i)
          re_matrices[d+i] = reinterpret_cast<MatrixD*>(blocks[d+(d+1)*i].get());
        break;
      case 14:        
        for(int i = 0, j = 0; i < d*d; ++i, ++j)
        {
          if(i%d == 0 && i > 0)
            j++;
          sq_matrices[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
        }
        sq_matrices[d*d] = reinterpret_cast<SquareMatrixD*>(blocks[(d+1)*(d+1)-1].get());
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[d*(d+1)+i].get());
        for(int i = 0; i < d; ++i)
          re_matrices[d+i] = reinterpret_cast<MatrixD*>(blocks[d+(d+1)*i].get());
        break;
    }// endswitch nstype

    for(int i = 0; i < d; ++i)
      boundCondition[i] = spaces[0]->get_boundary_condition();
    boundCondition[d] = spaces[1]->get_boundary_condition();

    // finite element functions
    for(int i = 0; i < d; ++i)
      feFunction[i] = s.u.GetComponent(i);
    feFunction[d] = &s.p;

    // local assembling object    
    LocalAssembling<d> la(this->db, LocalAssembling_type::NSE3D_Linear, 
                          feFunction.data(), example.get_coeffs());
    
    // assemble now the matrices and right hand side 
#ifdef __3D__
    Assemble3D(
#else
    Assemble2D(
#endif
               n_fe_spaces, spaces, 
               n_square_matrices, sq_matrices.data(),
               n_rectangular_matrices, re_matrices.data(), 
               n_rhs, rhs_array.data(), rhs_spaces.data(),
               boundCondition.data(), boundValues.data(), la);

    // assemble on the boundary if needed
    assemble_boundary_terms();
    
    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);
    
    //delete the temorary feFunctions gained by GetComponent
    for(int i = 0; i < d; ++i)
      delete feFunction[i];
    
    if(db["space_discretization_type"].is("local_projection"))
    {
      if(d == 3)
        ErrThrow("local_projection stabilization is not implemented in 3D");
      LPS_parameter_set lp{db["lps_coeff_type"], db["lps_delta0"],
                           db["lps_delta1"]};
      auto C = LPS_for_pressure_Scott_Zhang(blocks.at(8), false,
                                            this->example.get_nu(), lp);
      s.matrix.replace_blocks(*C, {{d,d}}, {false}); // creates a copy
    }
  }// endfor auto grid

/** When we call copy_nonactive in MPI-case, we have to remember the following:
   * it can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active.
   * That's why we have to update the values so that the vector becomes consistent again.
   * This is done here.
   */
#ifdef _MPI
  auto s = this->systems.front();
  for(int i = 0; i < d; ++i)
  {
    auto ui = s.solution.block(i);
    s.velocity_space->get_communicator().consistency_update(ui, 3);
  }
  auto p = s.solution.block(3);
  s.pressure_space->get_communicator().consistency_update(p, 3);
#endif
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::assemble_nonlinear_term()
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  size_t n_fe_spaces = 2; // space needed for assembling matrices
  size_t n_square_matrices = d*d + 1; // no of square matrices
  std::vector<SquareMatrixD*> sq_matrices(n_square_matrices, nullptr);
  int n_rectangular_matrices = 2*d; // maximum no of rectangular matrices
  std::vector<MatrixD*> re_matrices(n_rectangular_matrices, nullptr);
  int n_rhs = d+1; // maximum number of right hand sides
  std::vector<double*> rhs_array(n_rhs, nullptr); // right hand side 
  
  std::array<FEFunction*, d+1> feFunction;
  // boundary conditions and boundary values
  std::array<BoundaryConditionFunction*, d+1> boundCondition;
  std::array<BoundaryValuesFunction*, d+1> boundValues;
  for(int i = 0; i < d+1; ++i)
    boundValues[i] = example.get_bd()[i];
  
  std::array<const FESpace*, d+1> rhs_spaces;
  
  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(systems.size() > 1)
  {
    for( int block = 0; block < d ;++block)
    {
      std::vector<const FESpace*> spaces;
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
  
  bool mdml =  this->solver.is_using_multigrid() 
            && this->solver.get_multigrid()->is_using_mdml();
  bool is_stokes = this->db["problem_type"].is(3); // otherwise Navier-Stokes
  if ((mdml && !is_stokes)|| db["space_discretization_type"].is("upwind"))
  {
    // in case of upwinding we only assemble the linear terms. The nonlinear
    // term is not assembled but replaced by a call to the upwind method.
    // Note that we assemble the same terms over and over again here. Not 
    // nice, but otherwise we would have to store the linear parts in a 
    // separate BlockFEMatrix.
    this->assemble_linear_terms();
  }

  for(auto &s : this->systems)
  {
    //hold the velocity space, we'll need it...
    const FESpace * v_space = s.velocity_space.get();
    const FESpace * p_space = s.pressure_space.get();
    
#ifdef _MPI
    //MPI: solution in consistency level 3 (TODO: this might be superfluous here)
    for (size_t bl = 0; bl < s.solution.n_blocks(); ++bl)
    {
      s.matrix.get_communicators()[bl]->consistency_update(s.solution.block(bl),
                                                           3);
    }
#endif

    // spaces for matrices
    const FESpace *spaces[2] = {v_space, p_space};
    for(int i = 0; i < d; ++i)
      rhs_spaces[i] = v_space;
    rhs_spaces[d] = p_space;
    std::vector<std::vector<size_t>> cells;
    for(size_t i = 0; i < d; ++i)
    {
      for(size_t j = 0; j < d; ++j)
        cells.push_back({{i, j}});
    }
    auto blocks = s.matrix.get_blocks_uniquely(cells);
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());
        break;
      case 3: 
      case 4:
      case 14:
        for(int i = 0; i < d*d; ++i)
          sq_matrices[i] = reinterpret_cast<SquareMatrixD*>(blocks[i].get());
        break;
    }// endswitch nstype

    for(int i = 0; i < d; ++i)
      boundCondition[i] = spaces[0]->get_boundary_condition();
    boundCondition[d] = spaces[1]->get_boundary_condition();
    
    // finite element functions
    for(int i = 0; i < d; ++i)
      feFunction[i] = s.u.GetComponent(i);
    feFunction[d] = &s.p;
    

    //decide wether to assemble by upwinding or not
    bool finest_grid = (&s == &systems.at(0));
    bool do_upwinding = (db["space_discretization_type"].is("upwind")
                        || (mdml && !finest_grid))
                        && !is_stokes;
    // local assembling object    
    LocalAssembling<d> la(this->db, LocalAssembling_type::NSE3D_NonLinear, 
                          feFunction.data(), example.get_coeffs());
    if(!do_upwinding)
    {
      for(auto* mat : sq_matrices)
      {
        if(mat != nullptr)
          mat->reset();
      }

      // assemble now the matrices and right hand side 
#ifdef __3D__
      Assemble3D(
#else
      Assemble2D(
#endif
                 n_fe_spaces, spaces, n_square_matrices, sq_matrices.data(),
                 n_rectangular_matrices, re_matrices.data(), n_rhs,
                 rhs_array.data(), rhs_spaces.data(), boundCondition.data(),
                 boundValues.data(), la);
    }
    else
    {
#ifdef __3D__
      // the inverse of the example's diffusion coefficient
      double one_over_nu = 1./example.get_nu();
      for(auto mat : sq_matrices)
      {
        UpwindForNavierStokes3D(mat, feFunction[0], feFunction[1],
                                feFunction[2], one_over_nu);
      }
#else // 3D -> 2D
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                                la.get_fe_function(0), la.get_fe_function(1));
          break;

        case 3:
        case 4:
        case 14:
          // do upwinding with two matrices
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                                la.get_fe_function(0), la.get_fe_function(1));
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[1],
                                la.get_fe_function(0), la.get_fe_function(1));
          break;
      } // endswitch
#endif
    }
    

    //delete the temporary feFunctions gained by GetComponent
    for(int i = 0; i < d; ++i)
      delete feFunction[i];

    //TODO: Copying non-actives??

  }// endfor auto grid
}

/* ************************************************************************* */
template <int d>
bool NavierStokes<d>::stop_it(unsigned int iteration_counter)
{
  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  //compute and update defect and residuals
  compute_residuals();
  
  // the current norm of the residual
  const double normOfResidual = this->get_full_residual();
  // store initial residual, so later we can print the overall reduction
  if(iteration_counter == 0)
  {
    initial_residual = normOfResidual;
#ifdef _MPI
    auto comms = this->get_matrix().get_communicators();
    initial_rhs_norm = this->systems.front().rhs.norm(comms);
#else
    initial_rhs_norm = this->systems.front().rhs.norm();
#endif
  }
  
  // check if minimum number of iterations was performed already
  size_t min_it = db["nonlinloop_minit"];
  if(iteration_counter < min_it)
	  return false;

  // hold the residual from 10 iterations ago
  const double oldNormOfResidual = this->old_residuals.front().fullResidual;

  size_t max_it = db["nonlinloop_maxit"];
  double convergence_speed = db["nonlinloop_slowfactor"];
  bool slow_conv = false;

  if(normOfResidual >= convergence_speed*oldNormOfResidual)
    slow_conv = true;

  double limit = db["nonlinloop_epsilon"];
  if(db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= sqrt(this->get_size());
    if(my_rank==0)
     Output::print<1>("stopping tolerance for nonlinear iteration ", limit);
  }
  
  //check residual relative to initial right hand side
  if(db["nonlinloop_residual_relative_to_rhs"])
    limit *= initial_rhs_norm;

  // check if the iteration has converged, or reached the maximum number of
  // iterations or if convergence is too slow. Then return true otherwise false
  if((normOfResidual <= limit) || (iteration_counter == max_it) || (slow_conv))
  {
    if(slow_conv && my_rank==0)
      Output::print<1>(" SLOW !!! ", normOfResidual/oldNormOfResidual);

    // stop iteration
    adjust_pressure();
    return true;
  }
  else
    return false;
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::compute_residuals()
{
  System_per_grid& s = this->systems.front();
  unsigned int n_u_dof = s.solution.length(0);
  unsigned int n_p_dof = s.solution.length(d);

  // copy rhs to defect and compute defect
#ifdef _MPI
    //MPI: solution in consistency level 3 (TODO: maybe this is superfluous here
    // (because solution might be in level 3 consistency already)!)
    auto comms = s.matrix.get_communicators();
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(s.solution.block(bl), 3);
    }
#endif

  defect = s.rhs;
  s.matrix.apply_scaled_add(s.solution, defect, -1.);

  if(s.matrix.pressure_projection_enabled())
  {
    FEFunction defect_fctn(s.pressure_space, "p_defect",
                           "pressure defect function", &defect[d*n_u_dof],
                           n_p_dof);
    defect_fctn.project_into_L20();
  }
  
  std::vector<unsigned int> velocity_blocks(d, 0);
  std::iota(std::begin(velocity_blocks), std::end(velocity_blocks), 0);
#ifdef _MPI
  std::vector<const TParFECommunicator3D*> velocity_comms(d, nullptr);
  for(int i = 0; i < d; ++i)
    velocity_comms[i] = comms[i];
  double impuls_residual_square = defect.norm(velocity_blocks, velocity_comms);
  double mass_residual_square = defect.norm({d}, {comms[d]});
#else
  double impuls_residual_square = defect.norm(velocity_blocks);
  double mass_residual_square = defect.norm({d});
#endif

  // the struct 'Residuals' takes squares:
  impuls_residual_square *= impuls_residual_square;
  mass_residual_square *= mass_residual_square;
  Residuals current_residuals(impuls_residual_square, mass_residual_square);
  old_residuals.add(current_residuals);
  /*
  for(int i = 0; i < defect.length(d); ++i)
  {
    Output::print(i, "  ", defect.block(d)[i]);
  }
  Output::print("residual: ", current_residuals);
  defect.write("defect.txt");
  exit(0);*/
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::solve()
{
  System_per_grid& s = this->systems.front();
  double damping = this->db["nonlinloop_damping_factor"];
  // store previous solution for damping, it is a pointer so that we can avoid
  // the copy in case of no damping
  std::shared_ptr<BlockVector> old_solution(nullptr);
  if(damping != 1.0)
    old_solution = std::make_shared<BlockVector>(s.solution);  
  // solving:
#ifdef _MPI
  if(this->solver.get_db()["solver_type"].is("direct"))
  {
    if(damping != 1.0)
      Output::warn("NSE3D::solve", "damping in an MPI context is not tested");

    //set up a MUMPS wrapper
    MumpsWrapper mumps_wrapper(s.matrix);

    //kick off the solving process
    mumps_wrapper.solve(s.rhs, s.solution);
  }
  else
#endif
  this->solver.solve(s.matrix, s.rhs, s.solution);
  
  if(damping != 1.0)
  {
    s.solution.scale(damping);
    s.solution.add_scaled(*old_solution, 1-damping);
  }

  // project pressure if necessary
  if(s.matrix.pressure_projection_enabled())
    s.p.project_into_L20();
}

/* ************************************************************************* */
// this is a bad implementation, we need to find ways to pass parameters to the
// TFEFunction2D::GetErrors method
double delta0 = 0.;
template <int d>
void natural_error_norm_infsup_stabilizations(
  int N_Points, std::array<double*, d> xyz, double *AbsDetjk,
  const double *Weights, double hK, double **Der, double **Exact,
  double **coeffs, double *LocError);
template <int d>
void parameter_function_for_errors(const double *in, double *out)
{
  // d=2: u1, u2, u1x u1y, u2x, u2y
  // d=3: u1, u2, u3, u1x u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z
  for(int i = 0; i < d*(d+1); ++i)
    out[i] = in[2+i];
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::output(int i)
{
  int my_rank =0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  System_per_grid& s = this->systems.front();
  std::array<FEFunction*, d> velocity_components;
  for(int i = 0; i < d; ++i)
    velocity_components[i] = s.u.GetComponent(i);
  
  if((size_t)db["verbosity"] > 1)
  {
    for(int i = 0; i < d; ++i)
      velocity_components[i]->PrintMinMax("u" + std::to_string(i));
    s.p.PrintMinMax(std::string("p"));
  }
  
  if(i < 0)
    outputWriter.write();
  else
    outputWriter.write(i);
  
#ifdef __2D__
 // velocity_over_line({0.5, 0.}, {0.5, 1.}, 81, velocity_components);

  // geothermal tests flow rate
 /*velocity_over_line({0.15, 0.4}, {0.15, 0.6}, 30, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.15.txt");
 velocity_over_line({0.3, 0.4}, {0.3, 0.6}, 30, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.3.txt");
 velocity_over_line({0.45, 0.4}, {0.45, 0.6}, 3, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.45.txt");
 velocity_over_line({0.6, 0.4}, {0.6, 0.6}, 30, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.6.txt");
 velocity_over_line({0.75, 0.4}, {0.75, 0.4}, 30, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.75.txt");
*/
#endif

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    std::vector<std::array<double, 4>> computed_errors(d+1);
#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D nsAllDerivs[d+1] = {D000, D100, D010, D001};
#else
    TAuxParam2D aux, aux2;
    MultiIndex2D nsAllDerivs[d+1] = {D00, D10, D01};
#endif
    const FESpace *velocity_space = this->systems.front().velocity_space.get();
    const FESpace *pressure_space = this->systems.front().pressure_space.get();
    
    // errors in the velocity components
    for(int i = 0; i < d; ++i)
    {
      auto ui = velocity_components[i];
      ui->GetErrors(example.get_exact(i), d+1, nsAllDerivs, d, L2H1Errors,
                    nullptr, &aux, 1, &velocity_space,
                    computed_errors[i].data());
    }
    // error in divergence
    double div_error = s.u.GetL2NormDivergenceError(example.get_exact(0),
                                                    example.get_exact(1)
#ifdef __3D__
                                                    , example.get_exact(2)
#endif
                                                   );
    // errors in pressure
    s.p.GetErrors(example.get_exact(d), d+1, nsAllDerivs, d, L2H1Errors,
                  nullptr, &aux, 1, &pressure_space, computed_errors[d].data());


#ifdef __2D__
    int boundary_component_id;
    double un_boundary_error = 0.;
    double boundary_error_l2_squared = 0.;
    double boundary_error_l2 = 0;
    const ParameterDatabase e_db = example.get_database();
    int n_nitsche_bd = e_db["n_nitsche_bd"];
  if (n_nitsche_bd)
  {
    std::vector<size_t> nitsche_id = e_db["nitsche_id"];
      std::vector<double> nitsche_penalty = e_db["nitsche_penalty"];
      for(unsigned int k = 0; k < nitsche_id.size(); k++)
      {
        boundary_component_id = nitsche_id[k];

        double boundary_error_l2_u0[1], boundary_error_l2_u1[1];

        velocity_components[0]->GetL2BoundaryError(example.get_bd(0),
                &aux2, 1, &velocity_space,
                boundary_error_l2_u0, boundary_component_id);
        velocity_components[1]->GetL2BoundaryError(example.get_bd(1),
                &aux2, 1, &velocity_space,
                boundary_error_l2_u1, boundary_component_id);

        boundary_error_l2_squared += boundary_error_l2_u0[0] + boundary_error_l2_u1[0];

        // compute the L2-norm of the normal velocity error at the Nitsche boundaries
        un_boundary_error += s.u.GetL2NormNormalComponentError(example.get_bd(0),
                example.get_bd(1), boundary_component_id);
      }
      un_boundary_error = sqrt(un_boundary_error);
      boundary_error_l2 = sqrt(boundary_error_l2_squared);
  }
#endif


#ifdef _MPI
    int n_send = 2*(d+1)+1;
    double err_red[n_send]; //memory for global (across all processes) error
    double err_send[n_send]; //fill send buffer
    err_send[0] = computed_errors[0][0];
    err_send[1] = computed_errors[0][1];
    err_send[2] = computed_errors[1][0];
    err_send[3] = computed_errors[1][1];
    err_send[4] = computed_errors[2][0];
    err_send[5] = computed_errors[2][1];
    err_send[6] = div_error*div_error;
    err_send[7] = computed_errors[d][0];
    err_send[8] = computed_errors[d][1];

    MPI_Allreduce(err_send, err_red, n_send, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(i = 0; i < n_send; i++)
    { //MPI: sqrt was skipped in GetErrors function - do it here globally!
      err_red[i] = std::sqrt(err_red[i]);
    }
    //fill the reduced errors back where they belong
    computed_errors[0][0] = err_red[0];
    computed_errors[0][1] = err_red[1];
    computed_errors[1][0] = err_red[2];
    computed_errors[1][1] = err_red[3];
    computed_errors[2][0] = err_red[4];
    computed_errors[2][1] = err_red[5];
    div_error = err_red[6];
    computed_errors[d][0] = err_red[7];
    computed_errors[d][1] = err_red[8];
#endif

    double l2_error = 0;
    double h1_error = 0;
    for(int i = 0; i < d; ++i)
    {
      l2_error += computed_errors[i][0] * computed_errors[i][0];
      h1_error += computed_errors[i][1] * computed_errors[i][1];
    }

    errors.at(0) = std::sqrt(l2_error);
    errors.at(1) = div_error;
    errors.at(2) = std::sqrt(h1_error);
    errors.at(3) = computed_errors[d][0];
    errors.at(4) = computed_errors[d][1];
#ifdef __2D__
    errors.at(6) = boundary_error_l2;
    errors.at(7) = un_boundary_error;
#endif

    //print errors
    if(my_rank == 0)
    {
      Output::print("--------------------------------------------");
      Output::stat("NavierStokes", "Measured errors");
      Output::dash("L2(u)     : ", setprecision(14), errors.at(0));
      Output::dash("L2(div(u)): ", setprecision(14), errors.at(1));
      Output::dash("H1-semi(u): ", setprecision(14), errors.at(2));
      Output::dash("L2(p)     : ", setprecision(14), errors.at(3));
      Output::dash("H1-semi(p): ", setprecision(14), errors.at(4));
#ifdef __2D__
      Output::dash("L2(u)_boundary: ", setprecision(14), errors.at(6));
      Output::dash("L2(u.n)_boundary: ", setprecision(14), errors.at(7));
#endif
      Output::print("--------------------------------------------");
    }

    auto sdt(db["space_discretization_type"]);
    bool pspg = sdt.is("pspg");
    bool symm_gls = sdt.is("symm_gls");
    bool nonsymm_gls = sdt.is("nonsymm_gls");
    if(pspg || symm_gls || nonsymm_gls)
    {
      double nu = example.get_nu();
      double error_in_natural_norm = nu * errors[2]*errors[2];
      if(symm_gls)
        error_in_natural_norm += (1./nu) * errors[3]*errors[3];
      delta0 = db["pspg_delta0"];
      int beginparameter = 0;
      ParamFct * parameter_function = &parameter_function_for_errors<d>;
#ifdef __3D__
      int fevalue_fctindex[12] = {0, 1, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2};
      MultiIndex3D fevalue_multiindex[] = {D000, D000, D000, D200, D020, D002,
                                           D200, D020, D002, D200, D020, D002};
      TAuxParam3D NSE_aux(1, 3, 1, d*(d+1), nullptr, velocity_components.data(),
                          &parameter_function, fevalue_fctindex,
                          fevalue_multiindex, d*(d+1), &beginparameter);
#else
      int fevalue_fctindex[6] = {0, 1, 0, 0, 1, 1};
      MultiIndex2D fevalue_multiindex[6] = {D00, D00, D20, D02, D20, D02};
      TAuxParam2D NSE_aux(1, d*(d+1), velocity_components.data(),
                          &parameter_function, fevalue_fctindex,
                          fevalue_multiindex, d*(d+1), &beginparameter);
#endif
      double err[4];
      s.p.GetErrors(example.get_exact(2), d+1, nsAllDerivs, 1,
                    natural_error_norm_infsup_stabilizations<d>,
                    get_example().get_coeffs(), &NSE_aux, 1,
                    &pressure_space, err);
      error_in_natural_norm += err[0]*err[0];
      error_in_natural_norm = std::sqrt(error_in_natural_norm);
      errors[5] = error_in_natural_norm;
      Output::print("Error in natural(", sdt, ") norm: ", std::setprecision(14),
                    errors[5]);
    }
  } // if(this->db["compute_errors"])
  
  for(int i = 0; i < d; ++i)
    delete velocity_components[i];

  bool write_matrix = false;
  if (write_matrix)
  {
    // create an output file containing the whole FE matrix. This can be read into Matlab using the Matlab function mmread.m
    std::stringstream matrix_name;
    matrix_name << "Coeff_Matrix_matrixmarket";
    s.matrix.get_combined_matrix()->write(matrix_name.str());
  }

  //do postprocessing step depending on what the example implements
  example.do_post_processing(*this);
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::adjust_pressure()
{
  System_per_grid& s = this->systems.front();
  if(db["problem_type"].is(5)) // Navier--Stokes
  {
    int sign = 0;
    if(db["nse_nonlinear_form"].is("rotational"))
      sign = -1;
    if(db["nse_nonlinear_form"].is("emac"))
      sign = 1;
    if(sign)
    {
      std::array<FEFunction*, d> velocity_components;
      for(int i = 0; i < d; ++i)
        velocity_components[i] = s.u.GetComponent(i);
      typename FEFunction::AnalyticFunction 
      f = [&velocity_components, sign](const TBaseCell* cell, int i, 
                                       std::array<double, d> xyz)
          {
            double val = 0;
            for(int c = 0; c < d; ++c)
            {
              double val_ui;
#ifdef __3D__
              velocity_components[c]->FindValueLocal(cell, i, xyz[0], xyz[1],
                                                     xyz[2], &val_ui);
#else
              velocity_components[c]->FindValueLocal(cell, i, xyz[0], xyz[1],
                                                     &val_ui);
#endif
              val += val_ui*val_ui;
            }
            return sign * 0.5 * val;
          };
      s.p.add(f);
      for(int i = 0; i < d; ++i)
        delete velocity_components[i];
      // project pressure if necessary
      if(s.matrix.pressure_projection_enabled())
        s.p.project_into_L20();
    }
  }
}

/* ************************************************************************* */
template <int d>
typename NavierStokes<d>::FEFunction* NavierStokes<d>::get_velocity_component(int i)
{
  if(i == 0)
    return this->systems.front().u.GetComponent(0);
  else if(i == 1)
    return this->systems.front().u.GetComponent(1);
  else if(i == 2 && d == 3)
    return this->systems.front().u.GetComponent(2);
  else
    ErrThrow("There are only ", d, " velocity components!");
}

/* ************************************************************************* */
template <int d>
const Residuals& NavierStokes<d>::get_residuals() const
{
  return old_residuals.back();
}

/* ************************************************************************* */
template <int d>
double NavierStokes<d>::get_impuls_residual() const
{
  return old_residuals.back().impulsResidual;
}

/* ************************************************************************* */
template <int d>
double NavierStokes<d>::get_mass_residual() const
{
  return old_residuals.back().massResidual;
}

/* ************************************************************************* */
template <int d>
double NavierStokes<d>::get_full_residual() const
{
  return old_residuals.back().fullResidual;
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::reset_residuals()
{
  this->old_residuals = FixedSizeQueue<10, Residuals>();
}

/* ************************************************************************* */
template<int d> std::array<double, 8> NavierStokes<d>::get_errors() const
{
  return errors;
}

/* ************************************************************************* */
template <int d>
void natural_error_norm_infsup_stabilizations(int N_Points,
                                              std::array<double*, d>,
                                              double *AbsDetjk,
                                              const double *Weights, double hK,
                                              double **Der, double **Exact,
                                              double **coeffs, double *LocError)
{
  LocError[0] = 0.0;
  for(int i=0;i<N_Points;i++)
  {
    double nu = coeffs[i][0];
    double delta = compute_PSPG_delta(delta0, hK, nu);
    double *deriv = Der[i];
    double *exactval = Exact[i];
    double w = delta*Weights[i]*AbsDetjk[i];

    double t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;
      
    t = deriv[2]-exactval[2];
    LocError[0] += w*t*t;
    if(d == 3)
    {
      t = deriv[3]-exactval[3];
      LocError[0] += w*t*t;
    }
  }
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::assemble_boundary_terms()
{
  const ParameterDatabase e_db = example.get_database();
  int n_neumann_bd = e_db["n_neumann_bd"];
  int n_nitsche_bd = e_db["n_nitsche_bd"];

  for(System_per_grid& s : this->systems)
  {

    ///@todo this part of the code needs still to be implemented dimension-independent
#ifdef __2D__

    if (n_neumann_bd)
    {
      // Neumann BC
      std::vector<size_t> neumann_id = e_db["neumann_id"];
      std::vector<double> neumann_value = e_db["neumann_value"];

      for(unsigned int k = 0; k < neumann_id.size(); k++)
      {
        Output::print<2>(" Neumann BC on boundary: ", neumann_id[k],
                ", value = ",neumann_value[k] );
        const FESpace* v_space = s.velocity_space.get();

        BoundaryAssembling2D::rhs_g_v_n(s.rhs, v_space,
                nullptr,neumann_id[k],
                -1.*neumann_value[k]);
      }
    }

    if (n_nitsche_bd)
    {
      // Nitsche penalty for weak essential BC
      std::vector<size_t> nitsche_id = e_db["nitsche_id"];
      std::vector<double> nitsche_penalty = e_db["nitsche_penalty"];
      double effective_viscosity = this->example.get_nu();

      for (unsigned int k = 0; k < nitsche_id.size(); k++)
      {
        const FESpace * v_space = s.velocity_space.get();
        const FESpace * p_space = s.pressure_space.get();
        Output::print<2>(" Nitsche BC on boundary: ", nitsche_id[k], ", nitsche penalty: ", nitsche_penalty[k]);
        int sym_u = e_db["symmetric_nitsche_u"];
        int sym_p = e_db["symmetric_nitsche_p"];
        double sigma = this->example.get_inverse_permeability();
        double L_0 = db["L_0"];

        BoundaryAssembling2D::nitsche_bc(s.matrix, s.rhs,
                v_space, p_space,
                this->example.get_bd(0), this->example.get_bd(1),
                nitsche_id[k], nitsche_penalty[k],
                effective_viscosity, sigma, L_0,
                sym_u, sym_p);
      }

      double corner_stab = e_db["corner_stab"];
      if (corner_stab)
      {
        const TFESpace2D * v_space = s.velocity_space.get();
        Output::print<2>(" Corner stabilization is applied, corner_stab = ",corner_stab);
        double sigma = this->example.get_inverse_permeability();
        double L_0 = db["L_0"];
        corner_stab = corner_stab * (effective_viscosity + sigma * L_0* L_0);

        BoundaryAssembling2D::matrix_and_rhs_corner_stabilization(s.matrix, s.rhs, v_space,
                this->example.get_bd(0),
                this->example.get_bd(1),
                nitsche_id,
                corner_stab);
      }
    }

#else

    auto coll = s.velocity_space.get()->GetCollection();
    BoundaryAssembling3D ba;
    if (n_neumann_bd)
    {
      // Neumann BC
      std::vector<TBoundFace*> boundaryFaceList;
      boundaryFaceList.clear();
      std::vector<size_t> neumann_id = e_db["neumann_id"];
      std::vector<double> neumann_value = e_db["neumann_value"];

      std::vector<TBaseCell*> dummy;
      for (size_t k = 0; k < neumann_id.size(); k++)
      {
        Output::print<2>(" Neumann BC on boundary: ", neumann_id[k]);
        coll->get_face_list_on_component(neumann_id[k], boundaryFaceList);
        const TFESpace3D * v_space = s.velocity_space.get();
        ba.rhs_g_v_n(s.rhs, v_space, nullptr, boundaryFaceList,
                (int) neumann_id[k], -1.*neumann_value[k]);
      }
    }
    if (n_nitsche_bd)
    {
      // Nitsche penalty for weak essential BC
      std::vector<TBoundFace*> boundaryFaceList;
      boundaryFaceList.clear();
      std::vector<size_t> nitsche_id = e_db["nitsche_id"];
      std::vector<double> nitsche_penalty = e_db["nitsche_penalty"];

      for (size_t k = 0; k < nitsche_id.size(); k++)
      {
        Output::print<2>(" Nitsche BC on boundary: ", nitsche_id[k], ", nitsche penalty: ", nitsche_penalty[k]);
        coll->get_face_list_on_component(nitsche_id[k], boundaryFaceList);
        Output::print<5>("boundaryFaceList.size(): ", boundaryFaceList.size() );
        const TFESpace3D * v_space = s.velocity_space.get();
        const TFESpace3D * p_space = s.pressure_space.get();

        double effective_viscosity = this->example.get_nu();
        int sym_u = e_db["symmetric_nitsche_u"];
        int sym_p = e_db["symmetric_nitsche_p"];
        //double sigma = this->example.get_inverse_permeability();
        //double L_0 = db["L_0"];

        ba.nitsche_bc(s.matrix, s.rhs, v_space, p_space,
                nullptr, nullptr, nullptr,
                boundaryFaceList,
                nitsche_id[k], nitsche_penalty[k],
                effective_viscosity,// sigma, L_0,
                sym_u, sym_p);
      }
    }

#endif

  }
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::velocity_over_line(
  const std::vector<double>& start_point, const std::vector<double>& end_point,
  size_t number_of_points, std::array<FEFunction*, d> velocity_components, std::string name_of_file)
{
  //----------------------------------------------------------------------------------
  // The following output is made for the geometry channel.mesh or channel_simple.mesh
  Output::print("The values the solution u1, u2 (, u3) takes (and the velocity magnitude) at the line between [x1, y1] and [x2, y2] are saved in", name_of_file);
  std::ostringstream oss;
  oss << name_of_file; //"u_values_over_line.txt";

  std::string var = oss.str();
  std::ofstream velfile(var);

#ifdef __2D__
  double values_u1[3];
  double values_u2[3];

  for (unsigned int k = 0; k < number_of_points; k++)
  {
    double Y = start_point[1] + k * (end_point[1]-start_point[1])/(number_of_points-1);
    double X = start_point[0] + k * (end_point[0]-start_point[0])/(number_of_points-1); //x[0];

    // for (int i = 0; i < d; i++)
    //{
    velocity_components[0]->FindGradient(X, Y, values_u1);
    velocity_components[1]->FindGradient(X, Y, values_u2);
    //}
    velfile << "(X, Y) = (" << X << ", " << Y << "); " << " " << " " << " " <<  " (u0, u1, |u|) = ( " << values_u1[0] << ", " << values_u2[0] << ", " << sqrt( values_u1[0] * values_u1[0] + values_u2[0] * values_u2[0]) << ")" << endl;
  }
#endif
  velfile.close();
}


/* ************************************************************************* */
#ifdef __3D__
template class NavierStokes<3>;
#else
template class NavierStokes<2>;
#endif
