#include <NSE2D.h>
#include <Database.h>
#include <LinAlg.h> // Ddot, IntoL20FEFunction
#include <Upwind.h>
#include <GridTransfer.h>
#include <Multigrid.h>
#include <Assemble2D.h>
#include <NSE_local_assembling_routines.h>
#include <BoundaryAssembling2D.h>

#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!
#include <AuxParam2D.h>

ParameterDatabase NSE2D::default_NSE_database()
{
  Output::print<5>("creating a default NSE2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("NSE2D parameter database");
  
  //NSE2D requires a nonlinear iteration, set up a nonlinit_database and merge
  db.merge(ParameterDatabase::default_nonlinit_database());

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);
  
  db.merge(LocalAssembling2D::default_local_assembling_database(), true);
  db.merge(Example2D::default_example_database(), true);
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

/** ************************************************************************ */
NSE2D::System_per_grid::System_per_grid (const Example_NSE2D& example,
               TCollection& coll, std::pair<int,int> velocity_pressure_orders,
               NSE2D::Matrix type)
 : velocity_space(new TFESpace2D(&coll, "u", "Navier--Stokes velocity", example.get_bc(0),
                  velocity_pressure_orders.first)),
   pressure_space(new TFESpace2D(&coll, "p", "Navier--Stokes pressure", example.get_bc(2),
                  velocity_pressure_orders.second))
{
  // build the matrix due to NSE type
  switch (type)
  {
    case NSE2D::Matrix::Type1:
      matrix = BlockFEMatrix::NSE2D_Type1(*velocity_space, *pressure_space);
    break;
    case NSE2D::Matrix::Type2:
      matrix = BlockFEMatrix::NSE2D_Type2(*velocity_space, *pressure_space);
    break;
    case NSE2D::Matrix::Type3:
      matrix = BlockFEMatrix::NSE2D_Type3(*velocity_space, *pressure_space);
    break;
    case NSE2D::Matrix::Type4:
      matrix = BlockFEMatrix::NSE2D_Type4(*velocity_space, *pressure_space);
    break;
    case NSE2D::Matrix::Type14:
      matrix = BlockFEMatrix::NSE2D_Type14(*velocity_space, *pressure_space);
      break;
    default:
      ErrThrow("Unknown NSE type given to constructor of NSE2D::System_per_grid.");
  }

  rhs = BlockVector(matrix, true);
  solution = BlockVector(matrix, false);

  u = TFEVectFunct2D(velocity_space.get(), "u", "u", solution.block(0),
                     solution.length(0), 2);
  p = TFEFunction2D(pressure_space.get(), "p", "p", solution.block(2),
                    solution.length(2));

}

/** ************************************************************************ */
NSE2D::System_per_grid::System_per_grid(const System_per_grid& other)
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

/** ************************************************************************ */
NSE2D::NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
             int reference_id)
 : NSE2D(domain, param_db, Example_NSE2D(param_db), reference_id)
{
  // note that the way we construct the example above will produce a memory 
  // leak, but that class is small.
  // FIXME Find a workaround - we do not like memory leaks at all,
  // because they pollute our valgrind tests.
}

/** ************************************************************************ */
NSE2D::NSE2D(const TDomain & domain, const ParameterDatabase& param_db,
             const Example_NSE2D e, unsigned int reference_id)
    : systems(), example(e), db(default_NSE_database()), 
      outputWriter(param_db), solver(param_db), defect(), oldResiduals(), 
      initial_residual(1e10), errors()
{
  this->db.merge(param_db, false);
  this->set_parameters();
  
  std::pair <int,int> 
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and pressure spaces
  // this function returns a pair which consists of 
  // velocity and pressure order
  this->get_velocity_pressure_orders(velocity_pressure_orders);

  // determine NSE TYPE from Database TODO change that handling!
  NSE2D::Matrix type;
  switch (TDatabase::ParamDB->NSTYPE)
  {
    case  1: type = Matrix::Type1; break;
    case  2: type = Matrix::Type2; break;
    case  3: type = Matrix::Type3; break;
    case  4: type = Matrix::Type4; break;
    case 14: type = Matrix::Type14; break;
    default:
      ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
               " That NSE Block Matrix Type is unknown to class NSE2D.");
  }

  bool usingMultigrid = this->solver.is_using_multigrid();
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  // create finite element space and function, a matrix, rhs, and solution
  systems.emplace_back(example, *coll, velocity_pressure_orders, type);

  if(usingMultigrid)
  {
    // Construct multigrid object
    auto mg = solver.get_multigrid();
    bool mdml = mg->is_using_mdml();
    if(mdml)
    {
      // change the discretization on the coarse grids to lowest order 
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velocity_pressure_orders = {-1, -4711};
      this->get_velocity_pressure_orders(velocity_pressure_orders);
    }

    //determine multigrid levels
    int second_grid;
    int coarsest_grid = domain.get_ref_level() - mg->get_n_geometric_levels() + 1;
    if(mdml)
      second_grid = domain.get_ref_level(); //the finest grid is taken a second time in mdml
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
      systems.emplace_back(example, *coll, velocity_pressure_orders, type);
      // prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    // initialize the multigrid object with all the matrices on all levels
    mg->initialize(matrices);
  }
  
  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  
  // print out some information  
  this->output_problem_size_info();
}

/** ************************************************************************ */
void NSE2D::get_velocity_pressure_orders(
  std::pair<int,int>& velocity_pressure_orders)
{
  int velocity_order = velocity_pressure_orders.first;
  int pressure_order = velocity_pressure_orders.second;
  int order = 0;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5:
      order = velocity_order;
      break;
    case 12: case 13: case 14: case 15:
      // P2/P1disc and P3/P2disc elements are in general not inf-sup stable on 
      // triangles. If the last refinement step was barycentric refinement, then
      // these elements are inf-sup stable.  
      Output::warn("You chose to use Scott-Vogelius finite elements, make sure "
                   "that the final refinement step is barycentric, see the "
                   "parameter 'refinement_final_step_barycentric'.");
      order = velocity_order-10;
      break;
    case -1: case -2: case -3: case -4: case -5:
    case -101:
      order = velocity_order;
      break;
    // conforming fe spaces, with bubbles on triangles, regular Q_k on quads
    // (25, ie P5 with bubbles, is not implemented on Triangles)
    case 22: case 23: case 24: case 25:
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
      //pressure_order = 1;
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

void NSE2D::set_parameters()
{
  if(!db["problem_type"].is(3) && !db["problem_type"].is(5) &&
     !db["problem_type"].is(7) )
  {
    Output::warn<2>("The parameter problem_type doesn't correspond neither to NSE "
        "nor to Stokes. It is now reset to the default value for NSE (=5).");
    db["problem_type"] = 5;
  }

  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
  {
    //Assemble2DSlipBC does not work, and is not implemented yet
    ErrThrow("Set INTERNAL_SLIP_WITH_FRICTION to 0, this feature is not yet included.");
  }
}

/** ************************************************************************ */
void NSE2D::assemble()
{
  for(System_per_grid& s : this->systems)
  {
    s.rhs.reset(); //right hand side reset (TODO: is that necessary?)
    s.matrix.reset(); // reset matrix (needed for mdml where this is called)

    const TFESpace2D * v_space = s.velocity_space.get();
    const TFESpace2D * p_space = s.pressure_space.get();

    // declare the variables which Assemble2D needs and each nstype has to fill
    size_t N_FESpaces = 2;

    const TFESpace2D *fespmat[2] = {v_space, p_space};
    size_t n_sq_mat;
    TSquareMatrix2D *sq_matrices[5]{nullptr};//it's five pointers maximum (Type14)

    size_t n_rect_mat;
    TMatrix2D *rect_matrices[4]{nullptr};//it's four pointers maximum (Types 2, 4, 14)

    size_t N_Rhs = 3;
    double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), s.rhs.block(2)};
    const TFESpace2D *fesprhs[3] = {v_space, v_space, p_space};

    BoundCondFunct2D * boundary_conditions[3] = {
      v_space->get_boundary_condition(), v_space->get_boundary_condition(),
      p_space->get_boundary_condition() };

    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    //same for all: the local asembling object
    TFEFunction2D *fe_functions[3] =
      { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    LocalAssembling2D la(this->db, LocalAssembling_type::NSE3D_Linear,
                         fe_functions, example.get_coeffs());

    std::vector<std::shared_ptr<FEMatrix>> blocks =
        s.matrix.get_blocks_uniquely();
    n_sq_mat = 5;
    n_rect_mat = 4;
    switch(TDatabase::ParamDB->NSTYPE)
    {// switch over known Block Matrix types, treat each one individually,
      // using a priori knowledge about the structure and the way it fits
      // to the LocalAssembling2D object
      // TODO remove all reinterpret_casts as soon as Assembling process takes only FEMatrices
      // we have to use reinterpret_casts because dynamic downcasting won't work here
      // FIXME replace global switch by local checking of blockmatrix type!
      case 1:
        sq_matrices[0] =  reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get());
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
        break;
      case 2:
        sq_matrices[0] =  reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(3).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(4).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
        break;
      case 3:
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
        break;
      case 4:
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
        break;
      case 14:
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        sq_matrices[4] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
        break;
      default:
        ErrThrow("Sorry, the structure of that BlockMatrix is unknown to class NSE2D. "
            "I don't know how to pass its blocks to Assemble2D.");
    }

    //HOTFIX: Check the documentation!
    assemble_nse = Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION;

    // call the assemble method with the information that has been patched together
    Assemble2D(N_FESpaces, fespmat, n_sq_mat, sq_matrices,
               n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
               boundary_conditions, non_const_bound_values.data(), la);

    // assemble on the boundary
    assemble_boundary_terms();

    
    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);

    // TODO Maybe we have to explicitely set non-actives in non-diagonal blocks
    // to zero here, that was done in former code, but maybe we can move it to the solver part

    //tidy up
    delete fe_functions[0];
    delete fe_functions[1];
  }
}

/** ************************************************************************ */
void NSE2D::assemble_nonlinear_term()
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
        spaces.push_back(s.velocity_space.get());
        u_entries.push_back(s.solution.block(block));
        u_ns_dofs.push_back(s.solution.length(block));
      }
      GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
    }
  }
  
  bool mdml =  this->solver.is_using_multigrid() 
            && this->solver.get_multigrid()->is_using_mdml();
  bool is_stokes = this->db["problem_type"].is(3) ||
    this->db["problem_type"].is(7); // otherwise Navier-Stokes
  
  if ((mdml && !is_stokes)|| db["space_discretization_type"].is("upwind"))
  {
    // in case of upwinding we only assemble the linear terms. The nonlinear
    // term is not assembled but replaced by a call to the upwind method.
    // Note that we assemble the same terms over and over again here. Not 
    // nice, but otherwise we would have to store the linear parts in a 
    // separate BlockFEMatrix.
    this->assemble();
  }

  for(System_per_grid& s : this->systems)
  {
    //hold the velocity space, we'll need it...
    const TFESpace2D * v_space = s.velocity_space.get();
    const TFESpace2D * p_space = s.pressure_space.get();

    //the variables we will have to fill for the call to Assemble2D
    size_t n_fe_spaces = 2;
    const TFESpace2D* fe_spaces[2]{v_space, p_space};

    size_t n_sq_mat = 5;
    TSquareMatrix2D* sq_mat[5]{nullptr};//two pointers maximum

    size_t n_rect_mat = 4;
    TMatrix2D* rect_mat[4]{nullptr};

    size_t n_rhs = 3;
    double *rhs[3] = {nullptr, nullptr, nullptr};
    const TFESpace2D* fe_rhs[3] = {v_space, v_space, p_space};

    BoundCondFunct2D * boundary_conditions[2]
      = { v_space->get_boundary_condition(), p_space->get_boundary_condition()};
    std::array<BoundValueFunct2D*, 4> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];
    non_const_bound_values[3] = example.get_bd()[3];

    TFEFunction2D *fe_functions[3] = 
    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    LocalAssembling2D la_nonlinear(this->db,
                                   LocalAssembling_type::NSE3D_NonLinear,
                                   fe_functions, this->example.get_coeffs());

    //fetch us (a) pointer(s) to the diagonal A block(s)
    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely({{0,0},{1,1}});

    switch(TDatabase::ParamDB->NSTYPE)
    {// switch over known Block Matrix types, treat each one individually,
      // using a priori knowledge about the structure and the way it fits
      // to the LocalAssembling2D object
      // TODO remove all reinterpret casts as soon as Assembling process takes only FEMatrices
      // FIXME replace global switch by local checking of blockmatrix type!
      case 1:
      case 2:
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        break;
      case 3:
      case 4:
      case 14:
        blocks = s.matrix.get_blocks_uniquely({{0,0},{0,1},{1,0},{1,1}});
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_mat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_mat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(2).get());
        sq_mat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        break;
      default:
        ErrThrow("Sorry, the structure of that BlockMatrix is unknown to class NSE2D. "
            "I don't know how to pass its blocks to Assemble2D.");
    }
    
    bool on_finest_grid = &systems.front() == &s;
    bool do_upwinding = (db["space_discretization_type"].is("upwind")
                         || (mdml && !on_finest_grid))
                        && !is_stokes;

    // assemble the nonlinear part of NSE
    if(!do_upwinding)
    {
      for(size_t i =0; i < n_sq_mat; ++i)
      {
        //reset the matrices, linear part is assembled anew
        if(sq_mat[i] != nullptr)
          sq_mat[i]->reset();
      }

      //HOTFIX: Check the documentation!
      assemble_nse = Hotfixglobal_AssembleNSE::WITH_CONVECTION;

      //do the actual assembling
      Assemble2D(n_fe_spaces, fe_spaces, n_sq_mat, sq_mat, n_rect_mat, rect_mat,
                 n_rhs, rhs, fe_rhs, boundary_conditions,
                 non_const_bound_values.data(), la_nonlinear);
    }
    else
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes(la_nonlinear.GetCoeffFct(), sq_mat[0],
                                la_nonlinear.get_fe_function(0),
                                la_nonlinear.get_fe_function(1));
          Output::print<3>("UPWINDING DONE : level ");
          break;

        case 3:
        case 4:
        case 14:
          // do upwinding with two matrices
          Output::print<3>("UPWINDING DONE : level ");
          UpwindForNavierStokes(la_nonlinear.GetCoeffFct(), sq_mat[0],
                                la_nonlinear.get_fe_function(0),
                                la_nonlinear.get_fe_function(1));
          UpwindForNavierStokes(la_nonlinear.GetCoeffFct(), sq_mat[1],
                                la_nonlinear.get_fe_function(0),
                                la_nonlinear.get_fe_function(1));
          break;
      } // endswitch
    } // endif

    //end nonlinear assembling

    // copy Dirichlet values from right hand side into solution
    //(is this necessary here? solution has not been touched!)
    s.solution.copy_nonactive(s.rhs);

    //tidy up
    delete fe_functions[0];
    delete fe_functions[1];

  }


}

/** ************************************************************************ */
bool NSE2D::stopIt(unsigned int iteration_counter)
{
  // compute the residuals with the current matrix and solution
  this->computeNormsOfResiduals();

  // check if minimum number of iterations was performed already
  size_t min_it = db["nonlinloop_minit"];
  if(iteration_counter < min_it)
	  return false;

  // the current norm of the residual
  double normOfResidual = this->getFullResidual();
  // store initial residual, so later we can print the overall reduction
  if(iteration_counter == 0)
  {
    initial_residual = normOfResidual;
    initial_rhs_norm = this->systems.front().rhs.norm();
    Output::print("Initial rhs norm ", initial_rhs_norm);
  }
  // the residual from 10 iterations ago
  double oldNormOfResidual = this->oldResiduals.front().fullResidual;
  
  size_t Max_It = db["nonlinloop_maxit"];
  double convergence_speed = db["nonlinloop_slowfactor"];
  bool slow_conv = false;
  
  
  if(normOfResidual >= convergence_speed*oldNormOfResidual)
    slow_conv = true;
  
  double limit = db["nonlinloop_epsilon"];
  if (db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= sqrt(this->get_size());
    Output::print<1>("stopping tolerance for nonlinear iteration ", limit);
  }
  //check residual relative to initial right hand side
  if(db["nonlinloop_residual_relative_to_rhs"])
    limit *= initial_rhs_norm;


  // check if the iteration has converged, or reached the maximum number of
  // iterations or if convergence is too slow. Then return true otherwise false
  if( (normOfResidual<=limit) || (iteration_counter==Max_It) || (slow_conv) )
  {
    if(slow_conv)
      Output::print<1>(" SLOW !!! ", normOfResidual/oldNormOfResidual);
    // stop iteration
    return true;
  }
  else
    return false;
}

/** ************************************************************************ */
void NSE2D::computeNormsOfResiduals()
{
  System_per_grid& s = this->systems.front();
  unsigned int n_u_dof = s.solution.length(0);
  unsigned int n_p_dof = s.solution.length(2);
  
  // copy rhs to defect
  this->defect = s.rhs;
  s.matrix.apply_scaled_add(s.solution, defect,-1.);

  if(s.matrix.pressure_projection_enabled())
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
  double damping = this->db["nonlinloop_damping_factor"];
  // store previous solution for damping, it is a pointer so that we can avoid
  // the copy in case of no damping
  std::shared_ptr<BlockVector> old_solution(nullptr);
  if(damping != 1.0)
    old_solution = std::make_shared<BlockVector>(s.solution);
  // solve linear problem 
  solver.solve(s.matrix,s.rhs, s.solution);
  // apply damping if prescribed
  if(damping != 1.0)
  {
    s.solution.scale(damping);
    s.solution.add_scaled(*old_solution, 1-damping);
  }
  // project pressure if necessary
  if(s.matrix.pressure_projection_enabled())
    s.p.project_into_L20();
}


// this is a bad implementation, we need to find ways to pass parameters to the
// TFEFunction2D::GetErrors method
double delta0 = 0.;
void natural_error_norm_infsup_stabilizations(
  int N_Points, double *X, double *Y, double *AbsDetjk, double *Weights,
  double hK, double **Der, double **Exact, double **coeffs, double *LocError);
void parameter_function_for_errors(double *in, double *out)
{
  out[0] = in[2]; // u1 (for conformity within the coefficient function)
  out[1] = in[3]; // u2 (for conformity within the coefficient function)
  out[2] = in[4]; // u1_xx
  out[3] = in[5]; // u1_yy
  out[4] = in[6]; // u2_xx
  out[5] = in[7]; // u2_yy
}
/** ************************************************************************ */
void NSE2D::output(int i)
{
  System_per_grid& s = this->systems.front();
  TFEFunction2D* u1 = s.u.GetComponent(0);
  TFEFunction2D* u2 = s.u.GetComponent(1);
  
  // print the value of the largest and smallest entry in the finite element 
  // vector
  if((size_t)db["verbosity"]> 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    s.p.PrintMinMax();
  }

  if(i < 0)
    outputWriter.write();
  else
    outputWriter.write(i);
    
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
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
    errors[1] = s.u.GetL2NormDivergenceError(example.get_exact(0), 
                                             example.get_exact(1));
    
    errors.at(0) = sqrt(err[0]*err[0] + err[2]*err[2]);
    errors.at(2) = sqrt(err[1]*err[1] + err[3]*err[3]);    
    Output::print<1>("L2(u)     : ", setprecision(14), errors[0]);
    Output::print<1>("L2(div(u)): ", setprecision(14), errors[1]);
    Output::print<1>("H1-semi(u): ", setprecision(14), errors[2]);
    // errors in pressure
    s.p.GetErrors(example.get_exact(2), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &pressure_space, err);
    
    errors.at(3) = err[0];
    errors.at(4) = err[1];    
    Output::print<1>("L2(p)     : ", setprecision(14), errors[3]);
    Output::print<1>("H1-semi(p): ", setprecision(14), errors[4]);    
    
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
      TFEFunction2D* fe_functions_u[2] = {u1, u2};
      int fevalue_fctindex[6] = {0, 1, 0, 0, 1, 1};
      MultiIndex2D fevalue_multiindex[6] = {D00, D00, D20, D02, D20, D02};
      int beginparameter = 0;
      ParamFct * parameter_function = &parameter_function_for_errors;
      TAuxParam2D NSE_aux(1, 6, fe_functions_u, &parameter_function, 
                          fevalue_fctindex, fevalue_multiindex, 6,
                          &beginparameter);
      s.p.GetErrors(example.get_exact(2), 3, NSAllDerivatives, 1,
                    natural_error_norm_infsup_stabilizations,
                    get_example().get_coeffs(), &NSEaux_error, 1,
                    &pressure_space, err);
      error_in_natural_norm += err[0]*err[0];
      error_in_natural_norm = std::sqrt(error_in_natural_norm);
      errors[5] = error_in_natural_norm;
      Output::print("Error in natural(", sdt, ") norm: ", std::setprecision(14),
                    errors[5]);
    }
  } // if(this->db["compute_errors"])
  delete u1;
  delete u2;

  //do postprocessing step depending on what the example implements
  example.do_post_processing(*this);
}

/** ************************************************************************ */
void NSE2D::output_problem_size_info() const
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

/** ************************************************************************ */
std::array< double, int(6) > NSE2D::get_errors() const
{
  return errors;
}

/** ************************************************************************ */
const Residuals& NSE2D::getResiduals() const
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
void NSE2D::reset_residuals()
{
  this->oldResiduals = FixedSizeQueue<10, Residuals>();
}


void natural_error_norm_infsup_stabilizations(int N_Points, double *X,
                                              double *Y, double *AbsDetjk,
                                              double *Weights, double hK,
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
  }
}


void NSE2D::assemble_boundary_terms()
{
  for(System_per_grid& s : this->systems)
  {
    // Neumann BC
    for (int k = 0; k < TDatabase::ParamDB->n_neumann_boundary; k++)
    {
      // add the term int_Gamma -p v.n (with constant p)
      double bd_comp = TDatabase::ParamDB->neumann_boundary_id[k];
      double bd_pressure = TDatabase::ParamDB->neumann_boundary_value[k];
      Output::print<1>(" Neumann BC on boundary: ", bd_comp);
      const TFESpace2D * v_space = s.velocity_space.get();

      BoundaryAssembling2D::rhs_g_v_n(s.rhs, v_space,
          nullptr, bd_comp, -1. * bd_pressure);
    }

    // Nitsche penalty for weak essential BC
    for (int k = 0; k < TDatabase::ParamDB->n_nitsche_boundary; k++)
    {
      const TFESpace2D * v_space = s.velocity_space.get();
      const TFESpace2D * p_space = s.pressure_space.get();
      int bd_comp = TDatabase::ParamDB->nitsche_boundary_id[k];
      Output::print<1>(" Nitsche BC on boundary: ", bd_comp);
      double effective_viscosity = this->example.get_nu();
      double nitsche_gamma = TDatabase::ParamDB->nitsche_penalty[k];

      // gamma/h (u,v)
      BoundaryAssembling2D::matrix_u_v(s.matrix, v_space, bd_comp,
          nitsche_gamma * effective_viscosity,
          true);  // rescale local integral by edge values

      // gamma/h (uD,v) [rhs]
      BoundaryAssembling2D::rhs_uD_v(s.rhs, v_space,
          this->example.get_bd(0), this->example.get_bd(1), bd_comp,
          nitsche_gamma * effective_viscosity,
          true);   // rescale local integral by edge values

      // - (mu grad(u)n,v)
      BoundaryAssembling2D::matrix_gradu_n_v(s.matrix, v_space, bd_comp,
          -1. * effective_viscosity);

      // - sign_u * (u,mu grad(v)n) [sign_u=1: symmetrix, -1: skew-symmetric]
      BoundaryAssembling2D::matrix_gradv_n_u(s.matrix, v_space, bd_comp,
          -1. * TDatabase::ParamDB->s1* effective_viscosity);

      // - sign_u * (uD,mu grad(v)n) [rhs]
      BoundaryAssembling2D::rhs_gradv_n_uD(s.rhs, v_space,
          this->example.get_bd(0),this->example.get_bd(1), bd_comp,
          -1. * TDatabase::ParamDB->s1 * effective_viscosity);

      // (pn,v)
      BoundaryAssembling2D::matrix_p_v_n(s.matrix, v_space, p_space, bd_comp,
          1.);

      // sign_div * (u,qn)
      BoundaryAssembling2D::matrix_q_u_n(s.matrix, v_space, p_space, bd_comp,
          1. * TDatabase::ParamDB->s2);

      // sign_div * (uD,qn) [rhs]
      BoundaryAssembling2D::rhs_q_uD_n(s.rhs, v_space, p_space,
          this->example.get_bd(0), this->example.get_bd(1), bd_comp,
          1. * TDatabase::ParamDB->s2);
    }

    double corner_stab = db["corner_stab"];
    if (corner_stab)
    {
      const TFESpace2D * v_space = s.velocity_space.get();
      Output::print<1>(" Corner stabilization is applied. ");
      double effective_viscosity = this->example.get_nu();
      double sigma = this->example.get_inverse_permeability();
      double L_0 = TDatabase::ParamDB->L_0; //db["L_0"];
      BoundaryAssembling2D::matrix_cornerjump_u_n_cornerjump_v_n(s.matrix, v_space, // TDatabase::ParamDB->nitsche_boundary_id[k],
          1, // nBoundaryParts
          corner_stab * effective_viscosity + sigma * L_0 * L_0 );
    }
  }
}
