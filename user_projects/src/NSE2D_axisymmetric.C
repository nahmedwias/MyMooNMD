/**
 * Implementation of class NSE2D_axisymmetric.
 * Compared to NSE2D this
 * 	- holds two spaces for the velocity (radial 'r' and axial 'z')
 * 	- holds two fe functions for the velocity (radial 'r' and axial 'z')
 * 	- supports only 'NSType4' (though implicetely)
 *
 * Plus I did some changes which I found more suitable, like
 * 	- renaming of some functions towards_this_convention
 * 	- passing a grid hierarchy instead of a domain
 */
#include <NSE2D_axisymmetric.h>
#include <Database.h>
#include <LinAlg.h> // Ddot, IntoL20FEFunction
#include <Upwind.h>
#include <GridTransfer.h>
#include <Multigrid.h>
#include <Assemble2D.h>

#include <memory>

#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!

ParameterDatabase get_default_NSE2D_axisymmetric_parameters()
{
  Output::print<5>("creating a default NSE2D_axisymmetric parameter database");

  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("NSE2D parameter database");

  //NSE2D requires a nonlinear iteration, set up a nonlinit_database and merge
  db.merge(ParameterDatabase::default_nonlinit_database());

  // a default output database - needed here as long as there's no class handling the output
  db.merge(ParameterDatabase::default_output_database());

  //TODO this might become necessary for Stokes.
  //db["nonlinloop_maxit"] = 1;

  return db;
}

void AxialSymmetricNSE2D_AssembleFunction(double Mult, double *coeff,
                          double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs);

void AxialSymmetricNSE2D_ParamFunction(double *in, double *out);

/** ************************************************************************ */
NSE2D_axisymmetric::System_per_grid::System_per_grid (const Example_NSE2D& example,
               TCollection& coll, std::pair<int,int> velocity_pressure_orders)
: velocity_space_z(&coll, "u_z_space", "axial velocity", example.get_bc(0),
	velocity_pressure_orders.first, nullptr),
  velocity_space_r(&coll, "u_r_space", "radial velocity", example.get_bc(1),
	velocity_pressure_orders.first, nullptr),
  pressure_space(&coll, "p", "pressure", example.get_bc(2),
	velocity_pressure_orders.second, nullptr),
  matrix({&velocity_space_z, &velocity_space_r, &pressure_space}),
  rhs(matrix, true),
  solution(matrix, false),
  u_z(&velocity_space_z, "u_z", "axial velocity component", solution.block(0),
     solution.length(0)),
  u_r(&velocity_space_r, "u_r", "radial velocity component", solution.block(1),
     solution.length(1)),
  p(&pressure_space, "p", "p", solution.block(2),
     solution.length(2))
{
	// TODO Check whether the boundary conditions fit.
	// The example should inform us, which bdry is the rotation axis.
	// Then it can be checked whether the conditions there are correct
	// (Dirichlet 0 for radial, Neumann 0 for axial component).

	  //create new blocks with correct structures filled with 0
	  FEMatrix velo_velo_0_0(&velocity_space_z, &velocity_space_z); //A blocks
	  FEMatrix velo_velo_0_1(&velocity_space_z, &velocity_space_r);
	  FEMatrix velo_velo_1_0(&velocity_space_r, &velocity_space_z);
	  FEMatrix velo_velo_1_1(&velocity_space_r, &velocity_space_r);

	  FEMatrix pressure_velo_1(&pressure_space, &velocity_space_z);
	  FEMatrix pressure_velo_2(&pressure_space, &velocity_space_r);

	  FEMatrix velo_pressure_1(&velocity_space_z, &pressure_space);
	  FEMatrix velo_pressure_2(&velocity_space_r, &pressure_space);

	  // fill in the velo-velo blocks
	  matrix.replace_blocks(velo_velo_0_0, {{0,0}}, {false});
	  matrix.replace_blocks(velo_velo_0_1, {{0,1}}, {false});
	  matrix.replace_blocks(velo_velo_1_0, {{1,0}}, {false});
	  matrix.replace_blocks(velo_velo_1_1, {{1,1}}, {false});

	  // fill in the pressure_velo blocks B_1 and B_1^T
	  matrix.replace_blocks(pressure_velo_1, {{2,0}}, {false});
	  matrix.replace_blocks(velo_pressure_1, {{0,2}}, {false});

	  // fill in the pressure_velo blocks B_2 and B_2^T
	  matrix.replace_blocks(pressure_velo_2, {{2,1}}, {false});
	  matrix.replace_blocks(velo_pressure_2, {{1,2}}, {false});
}

/** ************************************************************************ */
NSE2D_axisymmetric::NSE2D_axisymmetric(
		std::list<TCollection*> grids,
		const ParameterDatabase& param_db,
        const Example_NSE2D ex)
    : systems(), example(ex), db(get_default_NSE2D_axisymmetric_parameters()),
      outputWriter(param_db), solver(param_db), defect(), oldResiduals(),
      initial_residual(1e10), initial_rhs_norm(0), errors()
{

  // read in additional parameters
  db.merge(param_db, false);

  std::pair <int,int>
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and pressure spaces
  // this function returns a pair which consists of
  // velocity and pressure order
  get_velocity_pressure_orders(velocity_pressure_orders);

  bool usingMultigrid = solver.is_using_multigrid();
  TCollection *coll = grids.front(); //the finest grid collection
  // create finite element space and function, a matrix, rhs, and solution
  systems.emplace_back(example, *coll, velocity_pressure_orders);

  if(usingMultigrid)
  {
    // Construct multigrid object
    std::shared_ptr<Multigrid> mg = solver.get_multigrid();
    bool mdml = mg->is_using_mdml();

    //Check whether number of given grids is alright
    size_t n_geo_multigrid_levels = mg->get_n_geometric_levels();
    size_t n_grids = grids.size();
    if(n_geo_multigrid_levels != n_grids )
      ErrThrow("Wrong number of grids for multigrid! I was expecting ",
               n_geo_multigrid_levels, " geometric grids but only got ", n_grids,".");

    if(mdml)
    {
      // change the discretization on the coarse grids to lowest order
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velocity_pressure_orders = {-1, -4711};
      get_velocity_pressure_orders(velocity_pressure_orders);
    }
    else
    {
      // for standard multigrid, pop the finest collection - it was already
      // used to construct a space before the "if(usingMultigrid)" clause
      // and will not (as in mdml) be used a second time with a different discretization
      grids.pop_front();
    }

     // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    // matrix on finest grid is already constructed
    matrices.push_back(&systems.back().matrix);

    for(auto coll : grids) // initialize the coarse grid space hierarchy
    {
      systems.emplace_back(example, *coll, velocity_pressure_orders);
      // prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    // initialize the multigrid object with all the matrices on all levels
    mg->initialize(matrices);
  }

  outputWriter.add_fe_function(&get_axial_velocity());
  outputWriter.add_fe_function(&get_radial_velocity());
  outputWriter.add_fe_function(&get_pressure());

  // print out some information
  this->output_problem_size_info();
}

/** ************************************************************************ */
//Absolute copy-paste from NSE2D. Nothing changed here.
void NSE2D_axisymmetric::get_velocity_pressure_orders(
  std::pair<int,int>& velocity_pressure_orders)
{
  int velocity_order = velocity_pressure_orders.first;
  int pressure_order = velocity_pressure_orders.second;
  int order = 0;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5:
    case 14: case 15:
      if(velocity_order > 10)
        order = velocity_order-10;
      else
        order = velocity_order;
      break;
    // P2/P1disc and P3/P2disc elements are not stable on triangles and not allowed
    // They are transformed into P2bubbles/P1disc and P3bubbles/P2disc
    // On quads, it will be Q2/P1disc, Q3/P2disc
    case 12: case 13:
      velocity_order += 10;
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
}

/** ************************************************************************ */
void NSE2D_axisymmetric::assemble_linear_terms()
{
  for(System_per_grid& s : this->systems)
  {
    s.rhs.reset(); //right hand side reset
    s.matrix.reset(); // reset matrix (needed for mdml where this is called)

    const TFESpace2D * v_z_space = &s.velocity_space_z;
    const TFESpace2D * v_r_space = &s.velocity_space_r;
    const TFESpace2D * p_space = &s.pressure_space;

    // declare the variables which Assemble2D needs
    size_t N_FESpaces = 3;

    const TFESpace2D *fespmat[3] = {v_z_space, v_r_space, p_space};

    size_t n_sq_mat = 2;
    TSquareMatrix2D *sq_matrices[n_sq_mat];//it's two pointers here (only A11 and A22)

    size_t n_rect_mat = 6;
    TMatrix2D *rect_matrices[n_rect_mat];//it's six pointers here (all other blocks, no C block)

    size_t N_Rhs = 3;
    double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), s.rhs.block(2)};
    const TFESpace2D *fesprhs[3] = {v_z_space, v_r_space, p_space};

    BoundCondFunct2D * boundary_conditions[3] = {
      v_z_space->GetBoundCondition(),
	  v_r_space->GetBoundCondition(),
      p_space->GetBoundCondition() };

    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    TFEFunction2D *fe_functions[3] = { &s.u_z, &s.u_r, &s.p };

    std::vector<std::shared_ptr<FEMatrix>> blocks =
        s.matrix.get_blocks_uniquely();

    // The following parameters do all mimic 'NSTYPE 4', which is the only
    // NSType for axisymmetric problems for which we have assemble routines ready.
    sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());

    rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get()); //two off-diagonal A blocks,
    rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(3).get());
    rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //then the lying B blocks
    rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
    rect_matrices[4] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //then the standing B blocks
    rect_matrices[5] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

    //HOTFIX: Check the documentation!
    assemble_nse = Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION;

     //Set up a customized local assembling object.
	 int myN_Terms = 7; // we need u_z, du_z/dz, du_z/dr, u_r, du_r/d_z, du_r/dr and p.
	 std::vector<MultiIndex2D> myDerivatives = { D00, D10, D01, D00, D10, D01, D00};
	 std::vector<int> myFESpaceNumber = {0, 0, 0, 1, 1, 1, 2};

	 //the following depends on the ordering of matrices. We must reorder them somewhat,
	 // because A21 and A12 are "non-square" (in the sense of ansatz != test)
	 //	 	S_z		S_r		S_p
	 //	S_z	(A11	A12		B1T )				(0 2 6)
	 // S_r	(A21	A22		B2T )	ordered as	(3 1 7)
	 // S_p	(B1		B2		0	)				(4 5 -)
	 //
	 //
	 std::vector<int> myRowSpace = {0, 1, 0, 1, 2, 2, 0, 1};
	 std::vector<int> myColumnSpace = {0, 1, 1, 0, 0, 1, 2, 2};
	 std::vector<int> myRhsSpace = {0, 1, 2};
	 CoeffFct2D* myCoeffs = example.get_coeffs();
	 AssembleFctParam2D* myAssembleParam = AxialSymmetricNSE2D_AssembleFunction;
	 ManipulateFct2D* myManipulate = NULL;
	 int myN_Matrices = 8;
	 int myN_Rhs = 3;
	 int myN_ParamFct = 1;
	 std::vector<ParamFct*> myParameterFct = {AxialSymmetricNSE2D_ParamFunction};
	 std::vector<int> myBeginParameter = { 0 };
	 int myN_Parameters = 3;
	 TFEFunction2D **myFEFunctions2D = fe_functions;
	 int myN_FEValues = 2;
	 std::vector<int> myFEValue_FctIndex = { 0, 1 };
	 std::vector<MultiIndex2D> myFEValue_MultiIndex = { D00, D00 };

	// 'type' is per default LocalAssembling2D_type::Custom
    // and 'discretization_type' is GALERKIN
    LocalAssembling2D la(
    		 myN_Terms,
    		 myDerivatives,
			 myFESpaceNumber,
			 myRowSpace,
			 myColumnSpace,
			 myRhsSpace,
			 myCoeffs,
			 myAssembleParam,
			 myManipulate,
			 myN_Matrices,
			 myN_Rhs,
			 myN_ParamFct,
			 myParameterFct,
			 myBeginParameter,
			 myN_Parameters,
			 myFEFunctions2D,
			 myN_FEValues,
			 myFEValue_FctIndex,
			 myFEValue_MultiIndex);

    // call the assemble method with the information that has been patched together
    Assemble2D(N_FESpaces, fespmat, n_sq_mat, sq_matrices,
               n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
               boundary_conditions, non_const_bound_values.data(), la);

    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);

  }
}
//
///** ************************************************************************ */
//void NSE2D_axisymmetric::assemble_nonlinear_term()
//{
// //Nonlinear assembling requires an approximate velocity solution on every grid!
//  if(systems.size() > 1)
//  {
//    for( int block = 0; block < 2 ;++block)
//    {
//      std::vector<const TFESpace2D*> spaces;
//      std::vector<double*> u_entries;
//      std::vector<size_t> u_ns_dofs;
//      for(auto &s : systems )
//      {
//        spaces.push_back(&s.velocity_space);
//        u_entries.push_back(s.solution.block(block));
//        u_ns_dofs.push_back(s.solution.length(block));
//      }
//      GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
//    }
//  }
//
//  bool mdml =  this->solver.is_using_multigrid()
//            && this->solver.get_multigrid()->is_using_mdml();
//  bool is_stokes = this->db["problem_type"].is(3); // otherwise Navier-Stokes
//
//  if ((mdml && !is_stokes)|| db["space_discretization_type"].is("upwind"))
//  {
//    // in case of upwinding we only assemble the linear terms. The nonlinear
//    // term is not assembled but replaced by a call to the upwind method.
//    // Note that we assemble the same terms over and over again here. Not
//    // nice, but otherwise we would have to store the linear parts in a
//    // separate BlockFEMatrix.
//    this->assemble();
//  }
//
//  for(System_per_grid& s : this->systems)
//  {
//    //hold the velocity space, we'll need it...
//    const TFESpace2D * v_space = &s.velocity_space;
//
//    //the variables we will have to fill for the call to Assemble2D
//    size_t n_fe_spaces = 1;
//    const TFESpace2D* fe_spaces[1]{v_space};
//
//    size_t n_sq_mat;
//    TSquareMatrix2D* sq_mat[2]{nullptr};//two pointers maximum
//
//    size_t n_rect_mat = 0;
//    TMatrix2D** rect_mat = nullptr;
//
//    size_t n_rhs = 0;
//    double** rhs = nullptr;
//    const TFESpace2D** fe_rhs = nullptr;
//
//    BoundCondFunct2D * boundary_conditions[1] = { v_space->GetBoundCondition() };
//    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
//    non_const_bound_values[0] = example.get_bd()[0];
//    non_const_bound_values[1] = example.get_bd()[1];
//    non_const_bound_values[2] = example.get_bd()[2];
//
//    TFEFunction2D *fe_functions[3] =
//    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
//    LocalAssembling2D la_nonlinear(NSE2D_NL, fe_functions,
//                                   this->example.get_coeffs());
//
//    //fetch us (a) pointer(s) to the diagonal A block(s)
//    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely({{0,0},{1,1}});
//
//    switch(TDatabase::ParamDB->NSTYPE)
//    {// switch over known Block Matrix types, treat each one individually,
//      // using a priori knowledge about the structure and the way it fits
//      // to the LocalAssembling2D object
//      // TODO remove all reinterpret casts as soon as Assembling process takes only FEMatrices
//      // FIXME replace global switch by local checking of blockmatrix type!
//      case 1:
//      case 2:
//        n_sq_mat = 1;
//        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
//        break;
//      case 3:
//      case 4:
//      case 14:
//        n_sq_mat = 2;
//        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
//        sq_mat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
//        break;
//      default:
//        ErrThrow("Sorry, the structure of that BlockMatrix is unknown to class NSE2D. "
//            "I don't know how to pass its blocks to Assemble2D.");
//    }
//
//    bool on_finest_grid = &systems.front() == &s;
//    bool do_upwinding = (db["space_discretization_type"].is("upwind")
//                         || (mdml && !on_finest_grid))
//                        && !is_stokes;
//
//    // assemble the nonlinear part of NSE
//    if(!do_upwinding)
//    {
//      for(size_t i =0; i < n_sq_mat; ++i)
//      {
//        //reset the matrices, linear part is assembled anew
//        sq_mat[i]->reset();
//      }
//
//      //HOTFIX: Check the documentation!
//      assemble_nse = Hotfixglobal_AssembleNSE::WITH_CONVECTION;
//
//      //do the actual assembling
//      Assemble2D(n_fe_spaces, fe_spaces, n_sq_mat, sq_mat, n_rect_mat, rect_mat,
//                 n_rhs, rhs, fe_rhs, boundary_conditions,
//                 non_const_bound_values.data(), la_nonlinear);
//    }
//    else
//    {
//      switch(TDatabase::ParamDB->NSTYPE)
//      {
//        case 1:
//        case 2:
//          // do upwinding with one matrix
//          UpwindForNavierStokes(la_nonlinear.GetCoeffFct(), sq_mat[0],
//                                la_nonlinear.get_fe_function(0),
//                                la_nonlinear.get_fe_function(1));
//          Output::print<3>("UPWINDING DONE : level ");
//          break;
//
//        case 3:
//        case 4:
//        case 14:
//          // do upwinding with two matrices
//          Output::print<3>("UPWINDING DONE : level ");
//          UpwindForNavierStokes(la_nonlinear.GetCoeffFct(), sq_mat[0],
//                                la_nonlinear.get_fe_function(0),
//                                la_nonlinear.get_fe_function(1));
//          UpwindForNavierStokes(la_nonlinear.GetCoeffFct(), sq_mat[1],
//                                la_nonlinear.get_fe_function(0),
//                                la_nonlinear.get_fe_function(1));
//          break;
//      } // endswitch
//    } // endif
//
//    //end nonlinear assembling
//
//    // copy Dirichlet values from right hand side into solution
//    //(is this necessary here? solution has not been touched!)
//    s.solution.copy_nonactive(s.rhs);
//
//    //tidy up
//    delete fe_functions[0];
//    delete fe_functions[1];
//
//  }
//
//
//}
//
/** ************************************************************************ */
bool NSE2D_axisymmetric::stop_iteration(unsigned int iteration_counter)
{
  // compute the residuals with the current matrix and solution
  compute_norms_of_residuals();
  // the current norm of the residual
  double normOfResidual = getFullResidual();
  // store initial residual, so later we can print the overall reduction
  if(iteration_counter == 0)
  {
    initial_residual = normOfResidual;
    initial_rhs_norm = systems.front().rhs.norm();
    Output::print("Initial rhs norm ", initial_rhs_norm);
  }
  // the residual from 10 iterations ago
  double oldNormOfResidual = oldResiduals.front().fullResidual;

  size_t Max_It = db["nonlinloop_maxit"];
  double convergence_speed = db["nonlinloop_slowfactor"];
  bool slow_conv = false;


  if(normOfResidual >= convergence_speed*oldNormOfResidual)
    slow_conv = true;

  double limit = db["nonlinloop_epsilon"];
  if (db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= sqrt(get_size());
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
void NSE2D_axisymmetric::compute_norms_of_residuals()
{
  System_per_grid& s = systems.front();
  size_t n_uz_dof = s.solution.length(0);
  size_t n_ur_dof = s.solution.length(0);
  size_t n_p_dof = s.solution.length(2);

  // copy rhs to defect
  defect = s.rhs;
  s.matrix.apply_scaled_add(s.solution, defect,-1.);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    IntoL20FEFunction(&defect[n_uz_dof+n_ur_dof], n_p_dof, &get_pressure_space(),
                      TDatabase::ParamDB->VELOCITY_SPACE,
                      TDatabase::ParamDB->PRESSURE_SPACE);
  }

  // square of the norms of the residual components
  double impuls_Residual = Ddot(n_uz_dof+n_ur_dof, &defect[0],&defect[0]);
  double mass_residual = Ddot(n_p_dof, &defect[n_uz_dof+n_ur_dof],
                              &defect[n_uz_dof+n_ur_dof]);

  Residuals currentResiduals(impuls_Residual, mass_residual);
  oldResiduals.add(currentResiduals);
}

/** ************************************************************************ */
void NSE2D_axisymmetric::solve()
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
    s.solution.add_scaled(*old_solution, damping);
  }
  // project pressure if necessary
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    s.p.project_into_L20();
}

/** ************************************************************************ */
void NSE2D_axisymmetric::output(int i)
{
  bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
  if(no_output)
	return;

  System_per_grid& s = this->systems.front();

  // print the value of the largest and smallest entry in the finite element
  // vector
  if((size_t) db["verbosity"] > 1 )
  {
    s.u_z.PrintMinMax();
    s.u_r.PrintMinMax();
    s.p.PrintMinMax();
  }

  outputWriter.add_fe_function(&s.u_z);
  outputWriter.add_fe_function(&s.u_r);
  outputWriter.add_fe_function(&s.p);
  outputWriter.write();

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    double err[4];
    TAuxParam2D NSEaux_error;
    MultiIndex2D NSAllDerivatives[3] = {D00, D10, D01};

    auto radial_velo_space = &get_radial_velocity_space();
    auto axial_velo_space = &get_axial_velocity_space();
    auto pressure_space = &get_pressure_space();

    // errors in first velocity component
    s.u_z.GetErrors(example.get_exact(0), 3, NSAllDerivatives, 2, L2H1Errors,
                  nullptr, &NSEaux_error, 1, &axial_velo_space, err);
    // errors in second velocity component
    s.u_r.GetErrors(example.get_exact(1), 3, NSAllDerivatives, 2, L2H1Errors,
                  nullptr, &NSEaux_error, 1, &radial_velo_space, err + 2);

    errors.at(0) = sqrt(err[0]*err[0] + err[2]*err[2]);
    errors.at(1) = sqrt(err[1]*err[1] + err[3]*err[3]);
    Output::print<1>("L2(u)     : ", setprecision(10), errors[0]);
    Output::print<1>("H1-semi(u): ", setprecision(10),errors[1]);
    // errors in pressure
    s.p.GetErrors(example.get_exact(2), 3, NSAllDerivatives, 2, L2H1Errors,
                  nullptr, &NSEaux_error, 1, &pressure_space, err);

    errors.at(2) = err[0];
    errors.at(3) = err[1];
    Output::print<1>("L2(p)     : ", setprecision(10), errors[2]);
    Output::print<1>("H1-semi(p): ", setprecision(10), errors[3]);

  } // if(this->db["compute_errors"])

//  //do postprocessing step depending on what the example implements
//  example.do_post_processing(*this);
}

/** ************************************************************************ */
void NSE2D_axisymmetric::output_problem_size_info() const
{
  int n_u = get_radial_velocity_space().GetN_DegreesOfFreedom();
  int n_u_active = get_radial_velocity_space().GetN_ActiveDegrees();
  int n_p = get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 2 * n_u + n_p; // total number of degrees of freedom

  //double h_min, h_max;
  //TCollection * coll = this->grids.front();
  //coll->GetHminHmax(&h_min, &h_max);
  Output::stat("NSE2D", "Mesh data and problem size");
  //Output::dash("cells              :  ", setw(10), coll->GetN_Cells());
  //Output::dash("h (min, max)       :  ", setw(10), h_min, setw(10), " ", h_max);
  Output::dash("dof velocity       :  ", setw(10), 2*n_u );
  Output::dash("dof velocity active:  ", setw(10), n_u_active + n_u);
  Output::dash("dof pressure       :  ", setw(10), n_p);
  Output::dash("dof all            :  ", setw(10), n_dof);
}

/** ************************************************************************ */
std::array< double, int(4) > NSE2D_axisymmetric::get_errors() const
{
  return errors;
}

/** ************************************************************************ */
const Residuals& NSE2D_axisymmetric::getResiduals() const
{
  return oldResiduals.back();
}

/** ************************************************************************ */
double NSE2D_axisymmetric::getImpulsResidual() const
{
  return oldResiduals.back().impulsResidual;
}

/** ************************************************************************ */
double NSE2D_axisymmetric::getMassResidual() const
{
  return oldResiduals.back().massResidual;
}

/** ************************************************************************ */
double NSE2D_axisymmetric::getFullResidual() const
{
  return oldResiduals.back().fullResidual;
}

void AxialSymmetricNSE2D_ParamFunction(double *in, double *out)
{
  out[0] = in[1]; // y value - is radius r
  out[1] = in[2]; // u1old
  out[2] = in[3]; // u2old
}

void AxialSymmetricNSE2D_AssembleFunction(double Mult, double *coeff,
                          double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs)
{
  double val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;

  double ** MatrixA11 = LocMatrices[0]; 	//square matrices
  double ** MatrixA22 = LocMatrices[1];
//double ** MatrixA12 = LocMatrices[2]; //rect matrices
//double ** MatrixA21 = LocMatrices[3];
  double ** MatrixB1 = LocMatrices[4];
  double ** MatrixB2 = LocMatrices[5];
  double ** MatrixB1T = LocMatrices[6];
  double ** MatrixB2T = LocMatrices[7];

  double* Rhs1 = LocRhs[0];
  double* Rhs2 = LocRhs[1];
  //no pressure rhs;

  int N_UZ = N_BaseFuncts[0];
  int N_UR = N_BaseFuncts[1];
  int N_P  = N_BaseFuncts[2];

  double* uz   = OrigValues[0];
  double* uz_z = OrigValues[1];
  double* uz_r = OrigValues[2];
  double* ur   = OrigValues[3];
  double* ur_z = OrigValues[4];
  double* ur_r = OrigValues[5];
  double* p    = OrigValues[6];

  double nu = coeff[0];
  double f1 = coeff[1];
  double f2 = coeff[2];

  double r = fabs(param[0]); // the radius (param[0] must be y coordinate)
  double uz_old = param[1]; // u1old
  double ur_old = param[2]; // u2old


  int sign = 1;
  if (param[0]>=0)
     sign = 1;
  else
     sign = -1;

  //go through the first block row (test space is space of uz)
  for(int i=0;i<N_UZ;i++)
  {
    double* MatrixA11Row = MatrixA11[i];

    test00 = uz[i];
    test10 = uz_z[i];
    test01 = uz_r[i]*sign;

    Rhs1[i] += Mult*test00*f1*r;

    for(int j=0;j<N_UZ;j++)
    {
      ansatz00 = uz[j];
      ansatz10 = uz_z[j];
      ansatz01 = uz_r[j]*sign;

      val  = nu*r*(test10*ansatz10+test01*ansatz01);

      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += r*(uz_old*ansatz10+ur_old*ansatz01)*test00;

      MatrixA11Row[j] += Mult * val;
    } // endfor j

    double* MatrixB1TRow = MatrixB1T[i];
    for(int j=0;j<N_P;j++)
    {
      ansatz00 = p[j];

      val = -Mult*ansatz00*test10*r;
      MatrixB1TRow[j] += val;
    }
  } // endfor i

  //go through the second block row (test space is space of uR)
  for(int i=0;i<N_UR;i++)
  {
    double* MatrixA22Row = MatrixA22[i];
    test00 = ur[i];
    test10 = ur_z[i];
    test01 = ur_r[i]*sign;

    Rhs2[i] += Mult*test00*f2*r;

    for(int j=0;j<N_UR;j++)
    {
      ansatz00 = ur[j];
      ansatz10 = ur_z[j];
      ansatz01 = ur_r[j]*sign;

      val  = nu*(r*(test10*ansatz10+test01*ansatz01)
		 -(ansatz10+ansatz01)*test00);

      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += r*(ur_old*ansatz10+uz_old*ansatz01)*test00;

      MatrixA22Row[j] += Mult * val;

    } // endfor j

    double* MatrixB2TRow = MatrixB2T[i];
    for(int j=0;j<N_P;j++)
    {
      ansatz00 = p[j];

      val = -Mult*ansatz00*(r*test01+test00);
      MatrixB2TRow[j] += val;
    }
  } // endfor i

  //go through third block row, where the test space is the pressure space
  for(int i=0;i<N_P;i++)
  {
	double* MatrixB1Row = MatrixB1[i];
	double* MatrixB2Row = MatrixB2[i];

    test00 = p[i];

    //fill B1 (ansatz space is velocity r space)
    for(int j=0;j<N_UZ;j++)
    {
      ansatz00 = uz[j];
      ansatz10 = uz_z[j];
      ansatz01 = uz_r[j]*sign;

      val = -Mult*test00*ansatz10*r;
      MatrixB1Row[j] += val;

    } // endfor j

    //fill B2 (ansatz space is velocity z space)
    for(int j=0;j<N_UR;j++)
    {
      ansatz00 = ur[j];
      ansatz10 = ur_z[j];
      ansatz01 = ur_r[j]*sign;

      val = -Mult*(r*ansatz01+ansatz00)*test00;
      MatrixB2Row[j] += val;
    } // endfor j
  } // endfor i
}

//void AxialSymmetricNSE2D_AssembleFunction(double Mult, double *coeff,
//                          double *param, double hK,
//                          double **OrigValues, int *N_BaseFuncts,
//                          double ***LocMatrices, double **LocRhs)
//{
//  double val;
//  double ansatz00, ansatz10, ansatz01;
//  double test00, test10, test01;
//
//  double ** MatrixA11 = LocMatrices[0]; 	//square matrices
//  double ** MatrixA22 = LocMatrices[1];
//  double ** MatrixA12 = LocMatrices[2]; //rect matrices
//  double ** MatrixA21 = LocMatrices[3];
//  double ** MatrixB1 = LocMatrices[4];
//  double ** MatrixB2 = LocMatrices[5];
//  double ** MatrixB1T = LocMatrices[6];
//  double ** MatrixB2T = LocMatrices[7];
//
//  double* Rhs1 = LocRhs[0];
//  double* Rhs2 = LocRhs[1];
//  //no pressure rhs;
//
//  int N_UR = N_BaseFuncts[0];
//  int N_UZ = N_BaseFuncts[1];
//  int N_P  = N_BaseFuncts[2];
//
//  double* ur   = OrigValues[0];
//  double* ur_r = OrigValues[1];
//  double* ur_z = OrigValues[2];
//  double* uz   = OrigValues[3];
//  double* uz_r = OrigValues[4];
//  double* uz_z = OrigValues[5];
//  double* p    = OrigValues[6];
//
//  double nu = coeff[0];
//  double f1 = coeff[1];
//  double f2 = coeff[2];
//  double r = fabs(param[0]); // the radius (param[0] must be y coordinate)
////  double ur_old = param[1]; // u1old
////  double uz_old = param[2]; // u2old
//
//  //go through the first block row (test space is space of ur)
//  for(int i=0;i<N_UR;i++)
//  {
//    double* MatrixA11Row = MatrixA11[i];
//
//    test00 = ur[i];
//    test10 = ur_r[i];
//    test01 = ur_z[i];
//
//    Rhs1[i] += Mult*test00*f1*r;
//
//    for(int j=0;j<N_UR;j++)
//    {
//      ansatz00 = ur[j];
//      ansatz10 = ur_r[j];
//      ansatz01 = ur_z[j];
//
//      val  = nu * (ansatz10*test10*r + (ansatz00*test00) / r + 0.5 * ansatz01*test01*r);
//
////      //HOTFIX: Check the documentation!
////      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
////        val += r*(ur_old*ansatz10+uz_old*ansatz01)*test00;
//
//      MatrixA11Row[j] += Mult * val;
//    } // endfor j
//
//    double* MatrixA12Row = MatrixA12[i];
//    for(int j=0;j<N_UZ;j++)
//    {
//        ansatz00 = uz[j];
//        ansatz10 = uz_r[j];
//        ansatz01 = uz_z[j];
//
//        val = nu*(0.5*ansatz10*test01*r);
//
//        MatrixA12Row[j] += Mult*val;
//    }
//
//    double* MatrixB1TRow = MatrixB1T[i];
//    for(int j=0;j<N_P;j++)
//    {
//      ansatz00 = p[j];
//
//      val = ansatz00*test10*r + test00;
//      MatrixB1TRow[j] += - Mult * val;
//
//    }
//  } // endfor i
//
//  //go through the second block row (test space is space of uz)
//  for(int i=0;i<N_UZ;i++)
//  {
//
//    test00 = uz[i];
//    test10 = uz_r[i];
//    test01 = uz_z[i];
//
//    Rhs2[i] += Mult*test00*f2*r;
//
//    double* MatrixA21Row = MatrixA21[i];
//    for(int j=0;j<N_UR;j++)
//    {
//        ansatz00 = ur[j];
//        ansatz10 = ur_r[j];
//        ansatz01 = ur_z[j];
//
//        val = nu*(0.5*ansatz01*test10*r);
//
//        MatrixA21Row[j] += Mult*val;
//    }
//
//    double* MatrixA22Row = MatrixA22[i];
//    for(int j=0;j<N_UZ;j++)
//    {
//      ansatz00 = uz[j];
//      ansatz10 = uz_r[j];
//      ansatz01 = uz_z[j];
//
//      val = nu*(0.5*ansatz10*test10*r+ansatz01*test01*r);
//
////      //HOTFIX: Check the documentation!
////      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
////        val += r*(ur_old*ansatz10+uz_old*ansatz01)*test00;
//
//      MatrixA22Row[j] += Mult * val;
//
//    } // endfor j
//
//    double* MatrixB2TRow = MatrixB2T[i];
//    for(int j=0;j<N_P;j++)
//    {
//      ansatz00 = p[j];
//
//      val = ansatz00*test01*r;
//      MatrixB2TRow[j] += - Mult*val;
//    }
//  } // endfor i
//
//  //go through third block row, where the test space is the pressure space
//  for(int i=0;i<N_P;i++)
//  {
//    test00 = p[i];
//
//    //fill B1 (ansatz space is velocity r space)
//    double* MatrixB1Row = MatrixB1[i];
//    for(int j=0;j<N_UR;j++)
//    {
//      ansatz00 = ur[j];
//      ansatz10 = ur_r[j];
//
//      val = test00*ansatz10*r + ansatz00;
//      MatrixB1Row[j] += -Mult * val;
//
//    } // endfor j
//
//    //fill B2 (ansatz space is velocity z space)
//    double* MatrixB2Row = MatrixB2[i];
//    for(int j=0;j<N_UZ;j++)
//    {
//      ansatz01 = uz_z[j];
//
//      val = test00*ansatz01*r;
//      MatrixB2Row[j] += -Mult * val;
//    } // endfor j
//  } // endfor i
//}

void NSType3_4NLGalerkinAxialSymm3D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01; // ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U, sign;
  double c0;
  double r, u1, u2, y;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu

  y  = param[0]; // y
  u1 = param[1]; // u1old
  u2 = param[2]; // u2old
  r = fabs(y);

  if (y>=0)
     sign = 1;
  else
     sign = -1;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i]*sign;
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j]*sign;
//      ansatz00 = Orig2[j];

      val  = c0*r*(test10*ansatz10+test01*ansatz01);
      val += r*(u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

       val  = c0*(r*(test10*ansatz10+test01*ansatz01)
		 -(ansatz10+ansatz01)*test00);
                  //+ansatz00*test01+ansatz01*test00+ansatz00*test00/r);
		  // old +ansatz00*test01-ansatz10*test00);
       val += r*(u1*ansatz10+u2*ansatz01)*test00;
       Matrix22Row[j] += Mult * val;
    } // endfor j
  } // endfor i
}
