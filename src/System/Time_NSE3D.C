#include <Time_NSE3D.h>
#include <Database.h>
#include <Assemble3D.h>
#include <LocalAssembling3D.h>
#include <LinAlg.h>
#include <ItMethod.h>
#include <MultiGridIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <Output3D.h>
#include <DirectSolver.h>
#include <MainUtilities.h>

/**************************************************************************** */
Time_NSE3D::System_per_grid::System_per_grid(const Example_NSE3D& example,
                  TCollection& coll, std::pair< int, int > order, 
                  Time_NSE3D::Matrix type)
 : velocitySpace_(&coll, (char*)"u", (char*)"velocity space",  example.get_bc(0),
                  order.first),
   pressureSpace_(&coll, (char*)"p", (char*)"pressure space", example.get_bc(3),
                  order.second),
   matrix_({&velocitySpace_, &velocitySpace_, &velocitySpace_, &pressureSpace_}),
   massMatrix_({&velocitySpace_, &velocitySpace_, &velocitySpace_}),
   rhs_(matrix_, true),
   solution_(matrix_, false),
   u_(&velocitySpace_, (char*)"u", (char*)"u", solution_.block(0),
     solution_.length(0), 3),
   p_(&pressureSpace_, (char*)"p", (char*)"p", solution_.block(3),
     solution_.length(3))
{
  // Mass Matrix
  // Output::increaseVerbosity(5);

  massMatrix_ = BlockFEMatrix::Mass_NSE3D(velocitySpace_);
      
  switch(type)
  {
    case Time_NSE3D::Matrix::Type1:
      matrix_ = BlockFEMatrix::NSE3D_Type1(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type2:
      matrix_ = BlockFEMatrix::NSE3D_Type2(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type3:
      matrix_ = BlockFEMatrix::NSE3D_Type3(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type4:
      matrix_ = BlockFEMatrix::NSE3D_Type4(velocitySpace_, pressureSpace_);
      break;
    case Time_NSE3D::Matrix::Type14:
      matrix_ = BlockFEMatrix::NSE3D_Type14(velocitySpace_, pressureSpace_);
      break;
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }
}

/**************************************************************************** */
Time_NSE3D::Time_NSE3D(const TDomain& domain, int reference_id)
  : Time_NSE3D(domain, *(new Example_NSE3D()), reference_id)
{
  
}

/**************************************************************************** */
Time_NSE3D::Time_NSE3D(const TDomain& domain, const Example_NSE3D& ex,
                       int reference_id)
 : systems_(), example_(ex), multigrid_(), defect_(),
   old_residual_(0), initial_residual_(1e10), errors_(), oldtau_()
{
 // TODO Implement the method "set_parameters" or "Check_parameters". Check
 // if it has to be called here or in the main program (see difference between
 // TNSE2D and NSE3D.
//  this->set_parameters();
  
  std::pair <int,int>
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                               TDatabase::ParamDB->PRESSURE_SPACE);

  // Set the velocity and pressure spaces.
  // This function returns a pair which consists of
  // velocity and pressure order. It takes the input read in
  // ParamDB and return a couple of VALID space orders.
  // TODO In this method there is an important thing to check about
  // pressure space.
  // NOTE: this method has the same purpose as set_parameters or check_parameters
  // In the future, these 3 different functions should be merged in one
  // check function. This will be the case with the implementation of
  // the new database.
  this->get_velocity_pressure_orders(velocity_pressure_orders);


  Time_NSE3D::Matrix type;
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

  bool usingMultigrid = (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5
                        && TDatabase::ParamDB->SOLVER_TYPE == 1 );

  if(!usingMultigrid)
  {
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);

  // create finite element space and function, a matrix, rhs and solution
  // all this by calling constructor of System_Per_Grid
  this->systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);

  // initialize the defect of the system. It has the same structure as
  // the rhs (and as the solution)
  this->defect_.copy_structure(this->systems_.front().rhs_);

  // print out some information about number of DoFs and mesh size
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 3 * n_u + n_p; // total number of degrees of freedom
  int nActive = this->get_velocity_space().GetN_ActiveDegrees();
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);

  Output::print<1>("N_Cells     : ", setw(10), coll->GetN_Cells());
  Output::print<1>("h (min,max) : ", setw(10), h_min ," ", setw(12), h_max);
  Output::print<1>("dof Velocity: ", setw(10), 3* n_u);
  Output::print<1>("dof Pressure: ", setw(10), n_p   );
  Output::print<1>("dof all     : ", setw(10), n_dof );
  Output::print<1>("active dof  : ", setw(10), 3*nActive);

  // Initial velocity = interpolation of initial conditions
  TFEFunction3D *u1 = this->systems_.front().u_.GetComponent(0);
  TFEFunction3D *u2 = this->systems_.front().u_.GetComponent(1);
  TFEFunction3D *u3 = this->systems_.front().u_.GetComponent(2);
  u1->Interpolate(example_.get_initial_cond(0));
  u2->Interpolate(example_.get_initial_cond(1));
  u3->Interpolate(example_.get_initial_cond(2));
  }
  else // multigrid  TODO: Multigrid in TNSE3D is not implemented yet.
    // it has to be constructed here
  {
//  // create spaces, functions, matrices on coarser levels
//  double *param = new double[10];
//  param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
//  param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
//  param[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
//  param[9] = 0;
//  this->multigrid.reset(new TNSE_MultiGrid(1, 2, param));
//  // number of refinement levels for the multigrid
//  int LEVELS = TDatabase::ParamDB->LEVELS;
//  if(LEVELS > domain.get_ref_level() + 1)
//    LEVELS = domain.get_ref_level() + 1;
//
//  // the matrix and rhs side on the finest grid are already constructed
//  // now construct all matrices, rhs, and solutions on coarser grids
//  for(int i = LEVELS - 2; i >= 0; i--)
//  {
//    unsigned int grid = i + domain.get_ref_level() + 1 - LEVELS;
//    TCollection *coll = domain.GetCollection(It_EQ, grid, reference_id);
//    this->systems.emplace_back(example, *coll, velocity_pressure_orders, type);
//  }
//
//  // create multigrid-level-objects, must be coarsest first
//  unsigned int i = 0;
//  for(auto it = this->systems.rbegin(); it != this->systems.rend(); ++it)
//  {
//    ErrThrow("NSE3D-multigrid needs to be checked");
//    this->multigrid->AddLevel(this->mg_levels(i, *it));
//    i++;
//  }
  }
}

//
//
///**************************************************************************** */
//void Time_NSE3D::set_parameters()
//{
//  if(TDatabase::ParamDB->EXAMPLE < 101)
//  {
//    ErrMsg("Example " << TDatabase::ParamDB->EXAMPLE
//    <<" is not supported for time dependent problem");
//    exit(1);
//  }
//
//  if(TDatabase::TimeDB->TIME_DISC == 0)
//  {
//    ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC
//          << " does not supported");
//    throw("TIME_DISC: 0 is not supported");
//  }
//}

/**************************************************************************** */
void Time_NSE3D::get_velocity_pressure_orders(std::pair< int, int > &velocity_pressure_orders)
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
      // TODO CHECK IF THIS IS CORRECT!!! IN NSE3D,pressure_order=1  ?!!
      pressure_order = -(velocity_order-1)*10;
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
void Time_NSE3D::assemble_initial_time()
{
  size_t nFESpace = 2;         // number of spaces for assembling matrices
  size_t nSquareMatrices = 10; // maximum number of square matrices (type 14)
  size_t nRectMatrices   = 6;  // maximum number of rectangular matrices (type 14)
  size_t nRhs     = 4;         // maximum number of right hand sides (type 14)

  std::vector<TSquareMatrix3D*> sqMatrices(nSquareMatrices);
  std::vector<TMatrix3D*>       rectMatrices(nRectMatrices);
  std::vector<double*>          rhsArray(nRhs);

  for(System_per_grid& s : this->systems_) // from back to front (coarse to fine)
  {
    const TFESpace3D *v_space = &s.velocitySpace_;
    const TFESpace3D *p_space = &s.pressureSpace_;

    // spaces for matrices
    const TFESpace3D *spaces[2] = {v_space, p_space};
    const TFESpace3D *rhsSpaces[4] = {v_space, v_space, v_space, p_space}; // NSTYPE 4 or 14

    // spaces for right hand sides
    s.rhs_.reset();
    rhsArray[0] = s.rhs_.block(0);
    rhsArray[1] = s.rhs_.block(1);
    rhsArray[2] = s.rhs_.block(2);
    rhsArray[3] = nullptr; //will be reset for type 4 and 14

    std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix_.get_blocks_uniquely();
    std::vector<std::shared_ptr<FEMatrix>> mass_blocks
         = s.massMatrix_.get_blocks_uniquely();

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        if(blocks.size() != 4)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size(), " instead of 4.");
        }
        nSquareMatrices = 2;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        // mass matrix
        sqMatrices[1] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());

        // rectangular matrices
        nRectMatrices = 3;
        rectMatrices[0] = reinterpret_cast<TMatrix3D*>(blocks.at(1).get());
        rectMatrices[1] = reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
        rectMatrices[2] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());

        nRhs = 3;
        break;
      case 2:
        if(blocks.size() != 7)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size(), " instead of 7.");
        }
        nSquareMatrices = 2;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        // mass matrix
        sqMatrices[1] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());

        // rectangular matrices
        nRectMatrices = 6;
        rectMatrices[0] = reinterpret_cast<TMatrix3D*>(blocks.at(4).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix3D*>(blocks.at(5).get());
        rectMatrices[2] = reinterpret_cast<TMatrix3D*>(blocks.at(6).get());
        rectMatrices[3] = reinterpret_cast<TMatrix3D*>(blocks.at(1).get()); //then the standing B blocks
        rectMatrices[4] = reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
        rectMatrices[5] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());

        nRhs = 3;
        break;
      case 3:
        if(blocks.size() != 12)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size(), " instead of 12.");
        }
        nSquareMatrices = 10;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
        sqMatrices[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
        sqMatrices[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
        sqMatrices[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
        sqMatrices[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
        sqMatrices[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
        sqMatrices[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
        sqMatrices[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
        // mass matrices
        sqMatrices[9] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
        ErrThrow("not tested yet!!!, kindly remove one mass matrix from the LocalAssembling3D routine");

        // rectangular matrices
        nRectMatrices = 3;
        rectMatrices[0] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());  // standing B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
        rectMatrices[1] = reinterpret_cast<TMatrix3D*>(blocks.at(11).get());

        nRhs = 3;
        break;
      case 4: // TODO THIS CASE GIVES A SEGMENTATION FAULT. HAS TO BE CORRECTED
        if(blocks.size() != 15)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size(), " instead of 15.");
        }
        nSquareMatrices = 10;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
        sqMatrices[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
        sqMatrices[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
        sqMatrices[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
        sqMatrices[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
        sqMatrices[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
        sqMatrices[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
        sqMatrices[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
        // mass matrices
        sqMatrices[9] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());

        // rectangular matrices
        nRectMatrices = 6;
        rectMatrices[0] = reinterpret_cast<TMatrix3D*>(blocks.at(12).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix3D*>(blocks.at(13).get());
        rectMatrices[2] = reinterpret_cast<TMatrix3D*>(blocks.at(14).get());
        rectMatrices[3] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());  //than the standing B blocks
        rectMatrices[4] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
        rectMatrices[5] = reinterpret_cast<TMatrix3D*>(blocks.at(11).get());

        // right hand side must be adapted
        rhsArray[3] = s.rhs_.block(3); // NSE type 4 includes pressure rhs_
        rhsSpaces[3] = p_space;
        nRhs = 4;
        break;
      case 14:
        if(blocks.size() != 16)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size(), " instead of 16.");
        }
        nSquareMatrices = 11;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
        sqMatrices[2] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
        sqMatrices[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
        sqMatrices[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
        sqMatrices[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
        sqMatrices[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
        sqMatrices[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
        sqMatrices[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
        // mass matrices
        sqMatrices[9] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get());
        // C block pressure pressure
        sqMatrices[10] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(15).get());

        // rectangular matrices
        nRectMatrices = 6;
        rectMatrices[0] = reinterpret_cast<TMatrix3D*>(blocks.at(12).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix3D*>(blocks.at(13).get());
        rectMatrices[2] = reinterpret_cast<TMatrix3D*>(blocks.at(14).get());
        rectMatrices[3] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());  //than the standing B blocks
        rectMatrices[4] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
        rectMatrices[5] = reinterpret_cast<TMatrix3D*>(blocks.at(11).get());

        // right hand side must be adapted
        rhsArray[3] = s.rhs_.block(3); // NSE type 14 includes pressure rhs_
        rhsSpaces[3]  = p_space;
        nRhs = 4;
        break;
      default:
        ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
               " That NSE Block Matrix Type is unknown to class Time_NSE3D.");
    } // end switch NSType

    for(unsigned int i=0; i<nSquareMatrices; i++)
      sqMatrices[i]->reset();
    for(unsigned int i=0; i<nRectMatrices;i++)
      rectMatrices[i]->reset();

    // Boundary conditions and value
    BoundCondFunct3D * boundary_conditions[4] = {
      v_space->getBoundCondition(), v_space->getBoundCondition(),
      v_space->getBoundCondition(), p_space->getBoundCondition() };

    std::array<BoundValueFunct3D*, 4> boundary_values;
    boundary_values[0] = example_.get_bd(0);
    boundary_values[1] = example_.get_bd(1);
    boundary_values[2] = example_.get_bd(2);
    boundary_values[3] = example_.get_bd(3);

    // Finite element functions for non-linear terms
    // for initial assembling, they correspond to the initial conditions
    TFEFunction3D *fe_functions[4] =
      { s.u_.GetComponent(0),
        s.u_.GetComponent(1),
        s.u_.GetComponent(2),
        &s.p_ };

    // local assembling object - used in Assemble3D
    const LocalAssembling3D
              localAssembling(LocalAssembling3D_type::TNSE3D_LinGAL,
                              fe_functions,this->example_.get_coeffs());

    // assemble all the matrices and right hand side
    Assemble3D(nFESpace, spaces,
               nSquareMatrices, sqMatrices.data(),
               nRectMatrices, rectMatrices.data(),
               nRhs, rhsArray.data(), rhsSpaces,
               boundary_conditions, boundary_values.data(), localAssembling);

  }// end for system per grid - the last system is the finer one (front)

  /** manage dirichlet condition by copying non-actives DoFs
  * from rhs to solution of front grid (=finest grid)
  * Note: this operation can also be done inside the loop, so that
  * the s.solution is corrected on every grid. This is the case in
  * TNSE2D.
  * TODO: CHECK WHAT IS THE DIFFERENCE BETWEEN doing this on every grid
  * and doing it only on the finest grid!
  * **/
  this->systems_.front().solution_.copy_nonactive(systems_.front().rhs_);

  // copy the last right hand side and solution vectors to the old ones
  this->old_rhs_      = this->systems_.front().rhs_;
  this->old_solution_ = this->systems_.front().solution_;
}

/**************************************************************************** */
void Time_NSE3D::assemble_rhs()
{
  // TODO Should it be timesteplength or currenttimesteplength
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  const double theta1 = TDatabase::TimeDB->THETA1;
  const double theta2 = TDatabase::TimeDB->THETA2;
  const double theta3 = TDatabase::TimeDB->THETA3;
  const double theta4 = TDatabase::TimeDB->THETA4;

  // reset the right hand side of the grid of interest (finest)
  System_per_grid& s = this->systems_.front();
  s.rhs_.reset();

  // some definitions necessary for assembling
  TFEFunction3D *fe_functions[4] =
  { s.u_.GetComponent(0),
    s.u_.GetComponent(1),
    s.u_.GetComponent(2),
    &s.p_ };

  int nRhs = 3;  // number of rhs blocks - TODO for NSType 4 and 14, it is 4
  const TFESpace3D *v_space = &this->get_velocity_space();
  const TFESpace3D *p_space = &this->get_pressure_space();

  // TODO Implement the case NSType 4 and 14 where there's a 4th rhs block
  std::vector<double*> rhsArray(nRhs);
  rhsArray[0] = s.rhs_.block(0);
  rhsArray[1] = s.rhs_.block(1);
  rhsArray[2] = s.rhs_.block(2);

  const TFESpace3D *spaces[2] = {v_space, p_space};
  const TFESpace3D *rhsSpaces[4] = {v_space, v_space, v_space, p_space};

  BoundCondFunct3D *boundary_conditions[4] = {
             v_space->getBoundCondition(), v_space->getBoundCondition(),
             v_space->getBoundCondition(), p_space->getBoundCondition() };

   std::array<BoundValueFunct3D*, 4> boundary_values;
   boundary_values[0] = this->example_.get_bd(0);
   boundary_values[1] = this->example_.get_bd(1);
   boundary_values[2] = this->example_.get_bd(2);
   boundary_values[3] = this->example_.get_bd(3);

   // Assembling the right hand side
  LocalAssembling3D
      localAssembling(LocalAssembling3D_type::TNSE3D_Rhs,
                      fe_functions,this->example_.get_coeffs());

  Assemble3D(1, spaces,
             0, nullptr,
             0, nullptr,
             nRhs, rhsArray.data(), rhsSpaces,
             boundary_conditions, boundary_values.data(), localAssembling);

  // copy the non active to the solution vector
  // since the rhs vector will be passed to the solver
  // and is modified with matrix vector multiplication
  // which also uses the non-actives
  s.solution_.copy_nonactive(s.rhs_);

  // now it is this->systems[i].rhs = f^k
  // scale by time step length and theta4 (only active dofs)
  s.rhs_.scaleActive(tau*theta4);

  // add rhs from previous time step
  if(theta3 != 0)
  {
    s.rhs_.addScaledActive((this->old_rhs_), tau*theta3);

    // now it is this->systems[i].rhs = tau*theta3*f^{k-1} + tau*theta4*f^k
    // next we want to set old_rhs to f^k (to be used in the next time step)
    this->old_rhs_.addScaledActive(s.rhs_, -1./(tau*theta3));
    this->old_rhs_.scaleActive(-theta3/theta4);
    this->old_rhs_.copy_nonactive(s.rhs_);
  }

  // FIXME Find other solution than this submatrix method.
  // M u^{k-1}
  s.massMatrix_.apply_scaled_submatrix(old_solution_, s.rhs_, 2, 2, 1.0);
  // -tau*theta2 * A u^{k-1}
  double factor = -tau*theta2;
  s.matrix_.apply_scaled_submatrix(old_solution_, s.rhs_, 2, 2, factor);

  // scale the BT blocks with time step length
  for(System_per_grid& s : this->systems_)
  {
    if(tau != oldtau_)
    {
      // TODO: change the factor to be THETA1*tau;
      factor = /*TDatabase::TimeDB->THETA1**/tau;
      if(this->oldtau_ != 0.0)
      {
        factor /= this->oldtau_;
        Output::print<1>("change in tau", this->oldtau_, "->", tau);
      }
      // scale the BT transposed blocks with the current time step
      const std::vector<std::vector<size_t>> cell_positions = {{0,3},
                                                               {1,3},
                                                               {2,3}};
	s.matrix_.scale_blocks(factor, cell_positions);
      if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
      {
        const std::vector<std::vector<size_t>> cell_positions_t = {{3,0},
                                                                   {3,1},
                                                                   {3,2}};
	s.matrix_.scale_blocks(factor, cell_positions_t);
      }
    }
  }

  this->oldtau_ = tau;
  // copy non active from solution into rhs vector
  s.rhs_.copy_nonactive(s.solution_);

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG
     && TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
     this->multigrid_->RestrictToAllGrids();

  Output::print<5>("assembled the system right hand side ");
}

/**************************************************************************** */
void Time_NSE3D::assemble_nonlinear_term()
{
  size_t nFESpace = 1; // space needed to assemble matrices
  size_t nSquareMatrices = 3;  // maximum 3 square matrices to be assembled
  size_t nRectMatrices = 0;
  size_t nRhs = 0; // no right hand side to be assembled here

  std::vector<TSquareMatrix3D*> sqMatrices(nSquareMatrices);
  std::vector<TMatrix3D*>       rectMatrices{nullptr};
  std::vector<double*>          rhsArray{nullptr};

  const TFESpace3D **rhsSpaces{nullptr};

  for(System_per_grid& s : this->systems_)
  {
    // spaces for matrices
    const TFESpace3D *spaces[1]={ &s.velocitySpace_ };

    std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix_.get_blocks_uniquely({{0,0},{1,1},{2,2}});

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        nSquareMatrices = 1;
        sqMatrices.resize(nSquareMatrices);
        sqMatrices.at(0) = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        break;
      case 3:
      case 4:
      case 14:
        nSquareMatrices = 3;
        sqMatrices.at(0) = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        sqMatrices.at(1) = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
        sqMatrices.at(2) = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
        break;
      default:
        ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
               " That NSE Block Matrix Type is unknown to class Time_NSE3D.");
    } // end switch nstypes

    // reset matrices to zero
    for(auto mat : sqMatrices)
    {
      mat->reset();
    }

    // Prepare info about boundary condition for assembling routines
    BoundCondFunct3D *boundary_conditions[1] = {spaces[0]->getBoundCondition()};

    std::array<BoundValueFunct3D*, 3> boundary_values;
    boundary_values[0] = this->example_.get_bd(0);
    boundary_values[1] = this->example_.get_bd(1);
    boundary_values[2] = this->example_.get_bd(2);

    TFEFunction3D *fe_functions[3] =
      { s.u_.GetComponent(0),
        s.u_.GetComponent(1),
        s.u_.GetComponent(2)};

    // assemble nonlinear matrices
    LocalAssembling3D
        localAssembling(LocalAssembling3D_type::TNSE3D_NLGAL,
                        fe_functions, this->example_.get_coeffs());

    Assemble3D(nFESpace, spaces,
               nSquareMatrices, sqMatrices.data(),
               nRectMatrices, rectMatrices.data(),
               nRhs, rhsArray.data(), rhsSpaces,
               boundary_conditions, boundary_values.data(), localAssembling);
  }

  Output::print<5>("Assembled the nonlinear matrix only ");
}

/**************************************************************************** */
void Time_NSE3D::assemble_system()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;

  for(System_per_grid& s : this->systems_)
  {
    const std::vector<std::vector<size_t>>
      cell_positions = {{0,0}, {0,1}, {0,2},
                        {1,0}, {1,1}, {1,2},
                        {2,0}, {2,1}, {2,2}};

    // note: declaring the auxiliary cell_positions is needed by the compiler
    // to sort out the overriding of the function scale_blocks_actives(...,...)
    s.matrix_.scale_blocks_actives(factor, cell_positions);

    const FEMatrix& mass_blocks =
        *s.massMatrix_.get_blocks().at(0).get();

    s.matrix_.add_matrix_actives(mass_blocks, 1.0,
                                 {{0,0}, {1,1}, {2,2}},
                                 {false, false, false});
  }
  cout << "test" << endl;
  Output::print<5>("Assembled the system matrix which will be passed to the ",
                   "solver");
}

///**************************************************************************** */
//bool Time_NSE3D::stop_it(unsigned int iteration_counter)
//{
//  System_per_grid& s = this->systems_.front();
//  unsigned int nuDof = s.solution_.length(0);
//  unsigned int npDof = s.solution_.length(2);
//
//  this->defect_ = s.rhs_;
//  s.matrix_.apply_scaled_add(s.solution_, defect_,-1.);
//  //
//  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
//    IntoL20FEFunction(&defect_[2*nuDof], npDof, &this->get_pressure_space(),
//                      TDatabase::ParamDB->VELOCITY_SPACE,
//                      TDatabase::ParamDB->PRESSURE_SPACE);
//  double residual =  Ddot(2*nuDof+npDof, &this->defect_[0], &this->defect_[0]);
//  double impulse_residual = Ddot(2*nuDof, &this->defect_[0],
//         &this->defect_[0]);
//  double mass_residual    = Ddot(npDof,&this->defect_[2*nuDof],
//         &this->defect_[2*nuDof]);
//
//  Output::print("nonlinear step  :  " , setw(3), iteration_counter);
//  Output::print("impulse_residual:  " , setw(3), impulse_residual);
//  Output::print("mass_residual   :  " , setw(3), mass_residual);
//  Output::print("residual        :  " , setw(3), sqrt(residual));
//
//  if (iteration_counter>0)
//  {
//  Output::print("rate:           :  " , setw(3), sqrt(residual)/old_residual_);
//  }
//
//  old_residual_ = sqrt(residual);
//  if(iteration_counter == 0)
//    initial_residual_ = sqrt(residual);
//
//  int Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
//  double limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
//  if (TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE)
//  {
//    limit *= sqrt(this->get_size());
//    Output::print("stopping tolerance for nonlinear iteration ", limit);
//  }
//
//  if ((((sqrt(residual)<=limit)||(it_counter==Max_It)))
//   && (iteration_counter>=TDatabase::ParamDB->SC_MINIT))
//   {
//     Output::print("ITE : ", setw(3), iteration_counter, "  RES : ", sqrt(residual),
//                   " Reduction : ",  sqrt(residual)/initial_residual_);
//
//     // descale the matrices, since only the diagonal A block will
//     // be reassembled in the next time step
//     this->deScaleMatrices();
//     return true;
//   }
//   else
//     return false;
//}

///**************************************************************************** */
//void Time_NSE3D::solve()
//{
//  System_per_grid& s = this->systems.front();
//
//  if((TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE !=5)
//    || (TDatabase::ParamDB->SOLVER_TYPE !=1 ))
//  {
//    if(TDatabase::ParamDB->SOLVER_TYPE != 2)
//      ErrThrow("only the direct solver is supported currently");
//
//    /// @todo consider storing an object of DirectSolver in this class
//    DirectSolver direct_solver(s.matrix,
//                               DirectSolver::DirectSolverTypes::umfpack);
//    direct_solver.solve(s.rhs, s.solution);
//  }
//  else
//    this->mg_solver();
//  // Important: We have to descale the matrices, since they are scaled
//  // before the solving process. Only A11 and A22 matrices are
//  // reset and assembled again but the A12 and A21 are scaled, so
//  // for the next iteration we have to descale, see assemble_system()
//    this->deScaleMatrices();
//
//this->old_solution = s.solution;
//  Output::print<5>("solver done");
//}
//
///**************************************************************************** */
//void Time_NSE3D::deScaleMatrices()
//{
//  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
//  double factor = tau*TDatabase::TimeDB->THETA1;
//  for(System_per_grid& s : this->systems)
//  {
//    const FEMatrix& mass_bloks = *s.Mass_Matrix.get_blocks().at(0).get();
//    s.matrix.add_matrix_actives(mass_bloks, -1.0, {{0,0}, {1,1}}, {false, false});
//    const std::vector<std::vector<size_t>>
//      cell_positions = {{0,0}, {0,1}, {1, 0}, {1, 1}};
//    // note: declaring the auxiliary cell_positions is needed by the compiler
//    // to sort out the overriding of the function scale_blocks_actives(...,...)
//    s.matrix.scale_blocks_actives(1./factor, cell_positions);
//  }
//}
//
///**************************************************************************** */
//TNSE_MGLevel* Time_NSE3D::mg_levels(int i, Time_NSE3D::System_per_grid& s)
//{
//  TNSE_MGLevel *mg_l;
//  int n_aux;
//  double alpha[2];
//
//  int v_space_code = TDatabase::ParamDB->VELOCITY_SPACE;
//  int p_space_code = TDatabase::ParamDB->PRESSURE_SPACE;
//
//  if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
//        || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE))
//     n_aux=4;
//  else
//     n_aux=2;
//
//  if (i==0)
//  {
//    alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
//    alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
//  }
//  else
//  {
//    alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
//    alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
//  }
//
//  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_TERRIBLY_UNSAFE();
//  switch(TDatabase::ParamDB->NSTYPE)
//  {
//    case 1:
//      ErrThrow("NSE3D::mg_levels: NSTYPE 1 is not supported");
//      break;
//
//    case 2:
//      mg_l = new TNSE_MGLevel2(i,
//                               reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get()),
//                               reinterpret_cast<TMatrix3D*>(blocks.at(3).get()), // B blocks
//                               reinterpret_cast<TMatrix3D*>(blocks.at(4).get()),
//                               reinterpret_cast<TMatrix3D*>(blocks.at(1).get()), // transposed B blocks
//                               reinterpret_cast<TMatrix3D*>(blocks.at(2).get()),
//                               s.rhs.get_entries(),
//                               s.solution.get_entries(),
//                               n_aux, alpha, v_space_code, p_space_code,
//                               nullptr, nullptr);
//      break;
//
//    case 3:
//      ErrThrow("NSE3D::mg_levels: NSTYPE 3 is not supported");
//      break;
//
//    case 4:
//       mg_l = new TNSE_MGLevel4(i,
//                                reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get()),
//                                reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get()),
//                                reinterpret_cast<TSquareMatrix3D*>(blocks.at(3).get()),
//                                reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get()),
//                                reinterpret_cast<TMatrix3D*>(blocks.at(6).get()),   // B blocks
//                                reinterpret_cast<TMatrix3D*>(blocks.at(7).get()),
//                                reinterpret_cast<TMatrix3D*>(blocks.at(2).get()),  // transposed B-blocks
//                                reinterpret_cast<TMatrix3D*>(blocks.at(5).get()),
//                                s.rhs.get_entries(),
//                                s.solution.get_entries(),
//                                n_aux, alpha, v_space_code, p_space_code,
//                                nullptr, nullptr);
//    break;
//    case 14:
//      mg_l = new TNSE_MGLevel14(i,
//                                reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get()),
//                                reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get()),
//                                reinterpret_cast<TSquareMatrix3D*>(blocks.at(3).get()),
//                                reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get()),
//                                reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get()),
//                                reinterpret_cast<TMatrix3D*>(blocks.at(6).get()),
//                                reinterpret_cast<TMatrix3D*>(blocks.at(7).get()),
//                                reinterpret_cast<TMatrix3D*>(blocks.at(2).get()),
//                                reinterpret_cast<TMatrix3D*>(blocks.at(5).get()),
//                                s.rhs.get_entries(),
//                                s.solution.get_entries(),
//                                n_aux, alpha, v_space_code, p_space_code,
//                                nullptr, nullptr);
//    break;
//  }
//  return mg_l;
//}
//
///**************************************************************************** */
//void Time_NSE3D::mg_solver()
//{
//  System_per_grid& s = this->systems.front();
//  TSquareMatrix3D *sqMat[5];
//  TSquareMatrix **sqmatrices = (TSquareMatrix **)sqMat;
//  TMatrix3D *recMat[4];
//  TMatrix **matrices = (TMatrix **)recMat;
//  MatVecProc *MatVect;
//  DefectProc *Defect;
//  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_TERRIBLY_UNSAFE();
//  switch(TDatabase::ParamDB->NSTYPE)
//  {
//    case 1:
//      ErrThrow("multigrid solver for the nstype 1 is not supported yet");
//      MatVect = MatVect_NSE1;
//      Defect = Defect_NSE1;
//      break;
//    case 2:
//      sqMat[0]  = reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
//      recMat[0] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());
//      recMat[1] = reinterpret_cast<TMatrix3D*>(blocks.at(4).get());
//      recMat[2] = reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
//      recMat[3] = reinterpret_cast<TMatrix3D*>(blocks.at(5).get());
//      MatVect = MatVect_NSE2;
//      Defect = Defect_NSE2;
//      break;
//    case 3:
//      ErrThrow("multigrid solver for the nstype 3 is not supported yet");
//      MatVect = MatVect_NSE3;
//      Defect = Defect_NSE3;
//      break;
//    case 4:
//      sqMat[0]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
//      sqMat[1]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
//      sqMat[2]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(3).get());
//      sqMat[3]=reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
//
//      recMat[0] = reinterpret_cast<TMatrix3D*>(blocks.at(6).get());
//      recMat[1] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
//      recMat[2] = reinterpret_cast<TMatrix3D*>(blocks.at(2).get());
//      recMat[3] = reinterpret_cast<TMatrix3D*>(blocks.at(5).get());
//      MatVect = MatVect_NSE4;
//      Defect = Defect_NSE4;
//      break;
//  }
//
//  int zero_start;
//  int nDof = this->get_size();
//  double *itmethod_rhs, *itmethod_sol;
//  TItMethod *itmethod, *prec;
//  if(TDatabase::ParamDB->SOLVER_TYPE ==1)
//  {
//    switch(TDatabase::ParamDB->SC_SOLVER_SADDLE)
//    {
//      case 11:
//        zero_start = 1;
//        break;
//      case 16:
//        zero_start = 0;
//        break;
//    }
//    switch(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE)
//    {
//      case 5:
//        prec = new TMultiGridIte(MatVect, Defect, nullptr, 0, nDof,
//                                 this->multigrid.get(), zero_start);
//        break;
//      default:
//        ErrThrow("Unknown preconditioner !!!");
//    }
//
//    if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
//    {
//      itmethod_sol = new double[nDof];
//      itmethod_rhs = new double[nDof];
//
//      memcpy(itmethod_sol, s.solution.get_entries(), nDof*SizeOfDouble);
//      memcpy(itmethod_rhs, s.rhs.get_entries(), nDof*SizeOfDouble);
//    }
//    else
//    {
//      itmethod_sol = s.solution.get_entries();
//      itmethod_rhs = s.rhs.get_entries();
//    }
//
//    switch(TDatabase::ParamDB->SC_SOLVER_SADDLE)
//    {
//      case 11:
//        itmethod = new TFixedPointIte(MatVect, Defect, prec, 0, nDof, 0);
//        break;
//      case 16:
//        itmethod = new TFgmresIte(MatVect, Defect, prec, 0, nDof, 0);
//        break;
//      default:
//        ErrThrow("Unknown preconditioner !!!");
//    }
//  }
//
//  switch(TDatabase::ParamDB->SOLVER_TYPE)
//  {
//    case 1:
//      itmethod->Iterate(sqmatrices,matrices,itmethod_sol,itmethod_rhs);
//      break;
//    case 2:
//      break;
//  }
//
//  if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
//  {
//    memcpy(s.solution.get_entries(), itmethod_sol, nDof*SizeOfDouble);
//    memcpy(s.rhs.get_entries(), itmethod_rhs, nDof*SizeOfDouble);
//
//    delete itmethod; delete prec;
//    delete [] itmethod_rhs;
//    delete [] itmethod_sol;
//  }
//}
//
///**************************************************************************** */
//void Time_NSE3D::output(int m, int& image)
//{
//  if(!TDatabase::ParamDB->WRITE_VTK
//    && !TDatabase::ParamDB->MEASURE_ERRORS)
//    return;
//
//  System_per_grid& s = this->systems.front();
//  TFEFunction3D * u1 = s.u.GetComponent(0);
//  TFEFunction3D * u2 = s.u.GetComponent(1);
//
//  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
//       s.p.project_into_L20();
//
//  if(TDatabase::ParamDB->SC_VERBOSE>1)
//  {
//    u1->PrintMinMax();
//    u2->PrintMinMax();
//    s.p.PrintMinMax();
//  }
//
//  if(TDatabase::ParamDB->MEASURE_ERRORS)
//  {
//    double locerr[8];
//    MultiIndex3D allderiv[3]= {D00, D10, D01};
//    const TFESpace3D *v_sp = &this->get_velocity_space();
//    const TFESpace3D *p_sp = &this->get_pressure_space();
//    TAuxParam3D aux;
//    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
//
//    u1->GetErrors(example.get_exact(0), 3, allderiv, 2, L2H1Errors,nullptr,
//                  &aux,1, &v_sp,locerr);
//
//    u2->GetErrors(example.get_exact(1), 3, allderiv, 2, L2H1Errors,nullptr,
//                  &aux,1, &v_sp,locerr+2);
//
//    errors[0] += (locerr[0]*locerr[0]+locerr[2]*locerr[2]
//                  + this->errors[1])*tau*0.5;
//    errors[1] = locerr[0]*locerr[0]+locerr[2]*locerr[2];
//    errors[2] += (locerr[1]*locerr[1]+locerr[3]*locerr[3]
//                  + this->errors[3])*tau*0.5;
//    errors[3] = locerr[1]*locerr[1]+locerr[3]*locerr[3];
//
//    Output::print<1>("L2(u) : ", setprecision(10), sqrt(this->errors[1]));
//    Output::print<1>("H1-semi(u) : ", setprecision(10), sqrt(this->errors[3]));
//
//    Output::print<1>("L2(0,t,L2(u)) : ", sqrt(this->errors[0]));
//    Output::print<1>("L2(0,t,H1-semi(u)) : ", sqrt(this->errors[2]));
//
//    s.p.GetErrors(example.get_exact(2), 3, allderiv, 2, L2H1Errors,
//                  nullptr, &aux, 1, &p_sp, locerr);
//
//    Output::print<1>("L2(p) : ", setprecision(10), locerr[0]);
//    Output::print<1>("H1-semi(p)) : " , setprecision(10), locerr[1] );
//
//    errors[4] += (locerr[0]*locerr[0] + this->errors[5])*tau*0.5;
//    errors[5] = locerr[0]*locerr[0];
//    Output::print<1>("L2(0,t,L2(p)) : ", sqrt(errors[4]) );
//
//    errors[6] += (locerr[1]*locerr[1] + this->errors[7])*tau*0.5;
//    errors[7] = locerr[1]*locerr[1];
//    Output::print<1>("L2(0,t,H1-semi(p)) : ", sqrt(errors[6]) );
//  }
//   delete u1;
//   delete u2;
//
//  if((m==0) || (m/TDatabase::TimeDB->STEPS_PER_IMAGE) )
//  {
//    if(TDatabase::ParamDB->WRITE_VTK)
//    {
//      TOutput3D output(2, 3, 1, 0, NULL);
//      output.AddFEFunction(&s.p);
//      output.AddFEVectFunct(&s.u);
//      std::string filename(TDatabase::ParamDB->OUTPUTDIR);
//      filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
//      if(image<10) filename += ".0000";
//      else if(image<100) filename += ".000";
//      else if(image<1000) filename += ".00";
//      else if(image<10000) filename += ".0";
//      else filename += ".";
//      filename += std::to_string(image) + ".vtk";
//      output.WriteVtk(filename.c_str());
//      image++;
//    }
//  }
//}
///**************************************************************************** */
//std::array< double, int(6) > Time_NSE3D::get_errors()
//{
//  std::array<double, int(6)> error_at_time_points;
//  error_at_time_points[0] = sqrt(errors[1]); // L2 velocity error
//  error_at_time_points[1] = sqrt(errors[3]); // H1 velocity error
//  error_at_time_points[2] = sqrt(errors[5]); // L2 pressure error
//  error_at_time_points[3] = sqrt(errors[7]); // H1 pressure error
//
//  return error_at_time_points;
//}
//
///**************************************************************************** */
