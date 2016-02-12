#include <NSE2D.h>
#include <MainUtilities.h> // GetVelocityAndPressureSpace
#include <Database.h>
#include <Output2D.h>
#include <LinAlg.h> // DDot
#include <MultiGridIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <DirectSolver.h>
#include <Upwind.h>

#include<Assemble2D.h>


/** ************************************************************************ */
NSE2D::System_per_grid::System_per_grid (const Example_NSE2D& example,
               TCollection& coll, std::pair<int,int> velocity_pressure_orders,
               NSE2D::Matrix type)
 : velocity_space(&coll, (char*)"u", (char*)"Navier--Stokes velocity", example.get_bc(0),
                  velocity_pressure_orders.first, nullptr),
   pressure_space(&coll, (char*)"p", (char*)"Navier--Stokes pressure", example.get_bc(2),
                  velocity_pressure_orders.second, nullptr)
// TODO CB: Building the matrix here and rebuilding later is due to the
// highly non-functional class TFEVectFunction2D (and TFEFunction2D,
// which do neither provide default constructors nor working copy assignments.)
   ,matrix({&velocity_space, &velocity_space, &pressure_space}),
   rhs(matrix, true),
   solution(matrix, false),
   u(&velocity_space, (char*)"u", (char*)"u", solution.block(0),
     solution.length(0), 2),
   p(&pressure_space, (char*)"p", (char*)"p", solution.block(2),
     solution.length(2))
{
// rebuild the matrix due to NSE type. We must be sure, that the rhs and solution
// vector which we built above do fit the new matrix, too!
  switch (type)
  {
    case NSE2D::Matrix::Type1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
    break;
    case NSE2D::Matrix::Type2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
    break;
    case NSE2D::Matrix::Type3:
      matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
    break;
    case NSE2D::Matrix::Type4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
    break;
    case NSE2D::Matrix::Type14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
      break;
    default:
      ErrThrow("Unknown NSE type given to constructor of NSE2D::System_per_grid.");
  }

}

/** ************************************************************************ */
NSE2D::NSE2D(const TDomain& domain, int reference_id)
 : NSE2D(domain, *(new Example_NSE2D()), reference_id)
{
  // note that the way we construct the example above will produce a memory 
  // leak, but that class is small.
  // FIXME Find a workaround - we do not like memory leaks at all,
  // because they pollute our valgrind tests.
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
  this->systems.emplace_back(example, *coll, velocity_pressure_orders, type);
  
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
  
  // done with the constructor in case we're not using multigrid
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
    this->systems.emplace_back(example, *coll, velocity_pressure_orders, type);
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
    case 2: // BDM2
      TDatabase::ParamDB->PROJECTION_SPACE = 1012;
      break;
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

void NSE2D::set_parameters()
{
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

    const TFESpace2D * v_space = &s.velocity_space;
    const TFESpace2D * p_space = &s.pressure_space;

    // declare the variables which Assemble2D needs and each nstype has to fill
    size_t N_FESpaces = 2;

    const TFESpace2D *fespmat[2] = {v_space, p_space};
    size_t n_sq_mat;
    TSquareMatrix2D *sq_matrices[5]{nullptr};//it's five pointers maximum (Type14)

    size_t n_rect_mat;
    TMatrix2D *rect_matrices[4]{nullptr};//it's four pointers maximum (Types 2, 4, 14)

    size_t N_Rhs = 2; //is 3 if NSE type is 4 or 14
    double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), nullptr}; //third place gets only filled
    const TFESpace2D *fesprhs[3] = {v_space, v_space, nullptr};  // if NSE type is 4 or 14

    BoundCondFunct2D * boundary_conditions[3] = {
      v_space->GetBoundCondition(), v_space->GetBoundCondition(),
      p_space->GetBoundCondition() };
    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    //same for all: the local asembling object
    TFEFunction2D *fe_functions[3] =
      { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    LocalAssembling2D la(NSE2D_Galerkin, fe_functions,
                     this->example.get_coeffs());

    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();

    switch(TDatabase::ParamDB->NSTYPE)
    {// switch over known Block Matrix types, treat each one individually,
      // using a priori knowledge about the structure and the way it fits
      // to the LocalAssembling2D object
      // TODO remove all reinterpret_casts as soon as Assembling process takes only FEMatrices
      // we have to use reinterpret_casts because dynamic downcasting won't work here
      // FIXME replace global switch by local checking of blockmatrix type!
      case 1:
        //CB DEBUG
        if(blocks.size() != 3)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        //END DEBUG

        n_sq_mat = 1;
        sq_matrices[0] =  reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

        n_rect_mat = 2;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get());
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());

        break;
      case 2:
        //CB DEBUG
        if(blocks.size() != 5)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        //END DEBUG
        n_sq_mat = 1;
        sq_matrices[0] =  reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

        n_rect_mat = 4;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(3).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(4).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());

        break;
      case 3:
        //CB DEBUG
        if(blocks.size() != 6)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        //END DEBUG
        n_sq_mat = 4;
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());

        n_rect_mat = 2;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
        break;

      case 4:
        //CB DEBUG
        if(blocks.size() != 8)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        //END DEBUG

        n_sq_mat = 4;
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());

        n_rect_mat = 4;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

        RHSs[2] = s.rhs.block(2); // NSE type 4 includes pressure rhs
        fesprhs[2]  = p_space;
        N_Rhs = 3;

        break;

      case 14:
        //CB DEBUG
        if(blocks.size() != 9)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        //END DEBUG
        n_sq_mat = 5;
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        sq_matrices[4] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());

        n_rect_mat = 4;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

        RHSs[2] = s.rhs.block(2); // NSE type 14 includes pressure rhs
        fesprhs[2]  = p_space;
        N_Rhs = 3;

        break;
      default:
        ErrThrow("Sorry, the structure of that BlockMatrix is unknown to class NSE2D. "
            "I don't know how to pass its blocks to Assemble2D.");
    }

    // call the assemble method with the information that has been patched together
    Assemble2D(N_FESpaces, fespmat, n_sq_mat, sq_matrices,
               n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
               boundary_conditions, non_const_bound_values.data(), la);

    // do upwinding TODO remove dependency of global values
    if((TDatabase::ParamDB->DISCTYPE == UPWIND)
       && !(TDatabase::ParamDB->PROBLEM_TYPE == 3))
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                                la.get_fe_function(0), la.get_fe_function(1));
          Output::print<3>("UPWINDING DONE : level ");
          break;
        case 3:
        case 4:
        case 14:
          // do upwinding with two matrices
          Output::print<3>("UPWINDING DONE : level ");
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                                la.get_fe_function(0), la.get_fe_function(1));
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[3],
                                la.get_fe_function(0), la.get_fe_function(1));
          break;
      } // endswitch
    } // endif

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
  // the class LocalAssembling2D which we will need next, requires an array of
  // pointers to finite element functions, i.e. TFEFunction2D **.
  for(System_per_grid& s : this->systems)
  {
    //hold the velocity space, we'll need it...
    const TFESpace2D * v_space = &s.velocity_space;

    //the variables we will have to fill for the call to Assemble2D

    size_t n_fe_spaces = 1;
    const TFESpace2D* fe_spaces[1]{v_space};

    size_t n_sq_mat;
    TSquareMatrix2D* sq_mat[2]{nullptr};//two pointers maximum

    size_t n_rect_mat = 0;
    TMatrix2D** rect_mat = nullptr;

    size_t n_rhs = 0;
    double** rhs = nullptr;
    const TFESpace2D** fe_rhs = nullptr;

    BoundCondFunct2D * boundary_conditions[1] = { v_space->GetBoundCondition() };
    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    TFEFunction2D *fe_functions[3] = 
    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    LocalAssembling2D la_nonlinear(NSE2D_Galerkin_Nonlinear, fe_functions,
                                   this->example.get_coeffs());

    //fetch us (a) pointer(s) to the diagonal A block(s)
    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely({{0,0},{1,1}});

    switch(TDatabase::ParamDB->NSTYPE)
    {// switch over known Block Matrix types, treat each one individually,
      // using a priori knowledge about the structure and the way it fits
      // to the LocalAssembling2D object
      // TODO remove all reinterpret casts as soon as Assembling process takes only FEMatrices
      // FIXME replace global switch by local checking of blockmatrix type!
      case 1:
        n_sq_mat = 1;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        break;
      case 2:
        n_sq_mat = 1;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        break;
      case 3:
        n_sq_mat = 2;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_mat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        break;
      case 4:
        n_sq_mat = 2;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_mat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        break;
      case 14:
        n_sq_mat = 2;
        sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_mat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        break;
      default:
        ErrThrow("Sorry, the structure of that BlockMatrix is unknown to class NSE2D. "
            "I don't know how to pass its blocks to Assemble2D.");
    }

    // assemble the nonlinear part of NSE
    for(size_t i =0; i < n_sq_mat; ++i)
    {//reset the matrices, linear part is assembled anew
      sq_mat[i]->reset();
    }
    //do the actual assembling
    Assemble2D(n_fe_spaces, fe_spaces, n_sq_mat, sq_mat, n_rect_mat, rect_mat,
               n_rhs, rhs, fe_rhs, boundary_conditions,
               non_const_bound_values.data(), la_nonlinear);

    // do upwinding TODO remove dependency of global values
    if((TDatabase::ParamDB->DISCTYPE == UPWIND)
        && !(TDatabase::ParamDB->PROBLEM_TYPE == 3))
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

  }
}

/** ************************************************************************ */
bool NSE2D::stopIt(unsigned int iteration_counter)
{
  // compute the residuals with the current matrix and solution
  this->computeNormsOfResiduals();
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
void NSE2D::computeNormsOfResiduals()
{
  System_per_grid& s = this->systems.front();
  unsigned int n_u_dof = s.solution.length(0);
  unsigned int n_p_dof = s.solution.length(2);
  
  // copy rhs to defect
  this->defect = s.rhs;
  s.matrix.apply_scaled_add(s.solution, defect,-1.);

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
    if(TDatabase::ParamDB->SOLVER_TYPE != 2)
      ErrThrow("only the direct solver is supported currently");
    
    /// @todo consider storing an object of DirectSolver in this class
    DirectSolver direct_solver(s.matrix, 
                               DirectSolver::DirectSolverTypes::umfpack);
    direct_solver.solve(s.rhs, s.solution);
  }
  else
  { // multigrid preconditioned iterative solver
    mg_solver();
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
  
  //get all blocks for the solver - here blocks may appear more than once
  // we must make use of non-const get_blocks_TERRIBLY_UNSAFE, because the
  // entire multigrid apparatus expects non-const TMatrix pointers
  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_TERRIBLY_UNSAFE();

  TSquareMatrix2D* A11 =  reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
  TSquareMatrix2D* A12 =  reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
  TMatrix2D* B1T =        reinterpret_cast<TMatrix2D*>( blocks.at(2).get());
  TSquareMatrix2D* A21 =  reinterpret_cast<TSquareMatrix2D*>( blocks.at(3).get());
  TSquareMatrix2D* A22 =  reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
  TMatrix2D* B2T =        reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
  TMatrix2D* B1 =         reinterpret_cast<TMatrix2D*>(blocks.at(6).get());
  TMatrix2D* B2 =         reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
  //TSquareMatrix2D* C = (TSquareMatrix2D*) blocks.at(8).get();


  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      //TODO The files are there, why should that not be supported?
      ErrThrow("NSE2D::mg_levels: NSTYPE 1 is not supported");
      break;
    
    case 2:
      mg_l = new TNSE_MGLevel2(i, A11, B1, B2, B1T, B2T,
                               s.rhs.get_entries(), 
                               s.solution.get_entries(), 
                               n_aux, alpha, v_space_code, p_space_code, 
                               nullptr, nullptr);
      break;
      
    case 3:
      //TODO The files are there, why should that not be supported?
      ErrThrow("NSE2D::mg_levels: NSTYPE 3 is not supported");
      break;
      
    case 4:

       mg_l = new TNSE_MGLevel4(i, A11, A12, A21, A22, B1, B2, B1T, B2T,
                                s.rhs.get_entries(), 
                                s.solution.get_entries(), 
                                n_aux, alpha, v_space_code, p_space_code, 
                                nullptr, nullptr);
       break;
    case 14:
      //TODO The files are there, why should that not be supported?
      ErrThrow("NSE2D::mg_levels: NSTYPE 14 is not supported");

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
  
  //get all blocks for the solver - here blocks may appear more than once
  // we must make use of non-const get_blocks_TERRIBLY_UNSAFE, because the
  // entire multigrid apparatus expects non-const TMatrix pointers
  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_TERRIBLY_UNSAFE();

  TSquareMatrix2D* A11 = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
  TSquareMatrix2D* A12 = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
  TMatrix2D* B1T =  reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
  TSquareMatrix2D* A21 = reinterpret_cast< TSquareMatrix2D*>(blocks.at(3).get());
  TSquareMatrix2D* A22 = reinterpret_cast< TSquareMatrix2D*>(blocks.at(4).get());
  TMatrix2D* B2T = reinterpret_cast< TMatrix2D*>(blocks.at(5).get());
  TMatrix2D* B1 = reinterpret_cast< TMatrix2D*>(blocks.at(6).get());
  TMatrix2D* B2 = reinterpret_cast< TMatrix2D*>(blocks.at(7).get());
  //TSquareMatrix2D* C = reinterpret_cast< TSquareMatrix2D*>(blocks.at(8).get());

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      sqMat[0] = A11;
      recMat[0] = B1;
      recMat[1] = B2;
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
      break;
    case 2:
      sqMat[0] = A11;
      
      recMat[0] = B1;
      recMat[1] = B2;
      recMat[2] = B1T;
      recMat[3] = B2T;
      
      MatVect = MatVect_NSE2;
      Defect = Defect_NSE2;
      break;
    case 3:
      sqMat[0] = A11;
      sqMat[1] = A12;
      sqMat[2] = A21;
      sqMat[3] = A22;
      recMat[0] = B1;
      recMat[1] = B2;
      MatVect = MatVect_NSE3;
      Defect = Defect_NSE3;
      break;
    case 4:
      sqMat[0] = A11;
      sqMat[1] = A12;
      sqMat[2] = A21;
      sqMat[3] = A22;
      recMat[0] = B1;
      recMat[1] = B2;
      recMat[2] = B1T;
      recMat[3] = B2T;

      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
      break;
    case 14:
      //TODO The files are there, why should that not be supported?
      ErrThrow("NSTYPE 14 is not fully supported, take NSTYPE 4");
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

