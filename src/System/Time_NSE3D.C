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

//project specific
#include <CoiledPipe.h>


#include <sys/stat.h>

/* *************************************************************************** */
  //TODO  So far of this object only the nonlin it stuff is used - switch entirely!
ParameterDatabase get_default_TNSE3D_parameters()
{
  Output::print<3>("creating a default TNSE3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TNSE3D parameter database");

  //NSE3D requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  return db;
}

/* *************************************************************************** */
Time_NSE3D::System_per_grid::System_per_grid(const Example_TimeNSE3D& example,
                  TCollection& coll, std::pair< int, int > order, 
                  Time_NSE3D::Matrix type
#ifdef _MPI
                  , int maxSubDomainPerDof
#endif
)
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
#ifdef _MPI
     //default construct parallel infrastructure, will be reset in the body
     , parMapperVelocity_(),
     parMapperPressure_(),
     parCommVelocity_(),
     parCommPressure_()
#endif
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
#ifdef _MPI
  velocitySpace_.SetMaxSubDomainPerDof(maxSubDomainPerDof);
  pressureSpace_.SetMaxSubDomainPerDof(maxSubDomainPerDof);

  //Must be reset here, because feSpace needs special treatment
  // This includes copy assignment - all because there is no good
  // way to communicate Maximum number of subdomains per dof to FESpace...
  parMapperVelocity_ = TParFEMapper3D(1, &velocitySpace_);
  parMapperPressure_ = TParFEMapper3D(1, &pressureSpace_);

  parCommVelocity_ = TParFECommunicator3D(&parMapperVelocity_);
  parCommPressure_ = TParFECommunicator3D(&parMapperPressure_);

  //print some information
  parCommVelocity_.print_info();
  parCommPressure_.print_info();

#endif
}

/**************************************************************************** */
Time_NSE3D::Time_NSE3D(const TDomain& domain, const ParameterDatabase& param_db,
                       const Example_TimeNSE3D& ex
#ifdef _MPI
                       , int maxSubDomainPerDof
#endif
)
 : db_(get_default_TNSE3D_parameters()), systems_(), example_(ex),
   solver_(param_db), mg_(nullptr),
   defect_(), old_residual_(), initial_residual_(1e10), errors_(), oldtau_()
{
 // TODO Implement the method "set_parameters" or "Check_parameters". Check
 // if it has to be called here or in the main program (see difference between
 // TNSE2D and NSE3D.
//  this->set_parameters();
  db_.merge(param_db, false);
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

  bool usingMultigrid = solver_.is_using_multigrid();
  if(!usingMultigrid)
  {
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, -4711);

  //CB HERE DO MODIFICATIONS TO COLLECTION DUE TO TWIST!
  if(true) //FIXME! only for twisted pipe example!!
  {
    // ...lie down the cell collection by swapping its vertices x and z coords
    CoiledPipe::swap_x_and_z_coordinates(coll);
    // ...and coil up the pipe by replacing the vertices' coords
    CoiledPipe::coil_pipe_helically(coll);
  }
  //CB END MODIFICATIONS

  #ifdef _MPI
  // create finite element space and function, a matrix, rhs, and solution
  systems_.emplace_back(example_, *coll, velocity_pressure_orders, type,
                      maxSubDomainPerDof);

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

  Output::print("N_Cells     : ", setw(10), coll->GetN_Cells());
  Output::print("h (min,max) : ", setw(10), h_min ," ", setw(12), h_max);
  Output::print("dof Velocity: ", setw(10), 3* n_u);
  Output::print("dof Pressure: ", setw(10), n_p   );
  Output::print("dof all     : ", setw(10), n_dof );
  Output::print("active dof  : ", setw(10), 3*nActive);

  // Initial velocity = interpolation of initial conditions
  TFEFunction3D *u1 = this->systems_.front().u_.GetComponent(0);
  TFEFunction3D *u2 = this->systems_.front().u_.GetComponent(1);
  TFEFunction3D *u3 = this->systems_.front().u_.GetComponent(2);

  u1->Interpolate(example_.get_initial_cond(0));
  u2->Interpolate(example_.get_initial_cond(1));
  u3->Interpolate(example_.get_initial_cond(2));

  #else
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

  Output::print("N_Cells     : ", setw(10), coll->GetN_Cells());
  Output::print("h (min,max) : ", setw(10), h_min ," ", setw(12), h_max);
  Output::print("dof Velocity: ", setw(10), 3* n_u);
  Output::print("dof Pressure: ", setw(10), n_p   );
  Output::print("dof all     : ", setw(10), n_dof );
  Output::print("active dof  : ", setw(10), 3*nActive);

  // Initial velocity = interpolation of initial conditions
  TFEFunction3D *u1 = this->systems_.front().u_.GetComponent(0);
  TFEFunction3D *u2 = this->systems_.front().u_.GetComponent(1);
  TFEFunction3D *u3 = this->systems_.front().u_.GetComponent(2);
  u1->Interpolate(example_.get_initial_cond(0));
  u2->Interpolate(example_.get_initial_cond(1));
  u3->Interpolate(example_.get_initial_cond(2));

#endif
  }
  else
  {
    ErrThrow("No multigrid yet. When implementing, stick e.g. to NSE3D.");
  }
}

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

        // rectangular matrices
        nRectMatrices = 3;
        rectMatrices[0] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get());  // standing B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
        rectMatrices[2] = reinterpret_cast<TMatrix3D*>(blocks.at(11).get());

        nRhs = 3;
        break;
      case 4:
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
      case 14: // TODO: NSType14 has still to be implemented
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
  this->systems_.front().solution_.copy_nonactive(this->systems_.front().rhs_);

  /** After copy_nonactive, the solution vectors needs to be Comm-updated
   * in MPI-case in order to be consistently saved. It is necessary that
   * the vector is consistently saved because it is the only way to
   * ensure that its multiplication with an inconsistently saved matrix
   * (multiplication which appears in the defect and rhs computations)
   * give the correct results.
   * When we call copy_nonactive in MPI-case, we have to remember the following:
   * it can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active.
   * That's why we have to update the values so that the vector becomes consistent again.
   */
  #ifdef _MPI
    double *u1  = this->systems_.front().solution_.block(0);
    double *u2  = this->systems_.front().solution_.block(1);
    double *u3  = this->systems_.front().solution_.block(2);
    double *p   = this->systems_.front().solution_.block(3);
    this->systems_.front().parCommVelocity_.CommUpdate(u1);
    this->systems_.front().parCommVelocity_.CommUpdate(u2);
    this->systems_.front().parCommVelocity_.CommUpdate(u3);
    this->systems_.front().parCommPressure_.CommUpdate(p);
  #endif

  // copy the last right hand side and solution vectors to the old ones
  this->old_rhs_      = this->systems_.front().rhs_;
  this->old_solution_ = this->systems_.front().solution_;
}

/**************************************************************************** */
void Time_NSE3D::assemble_rhs()
{
  // TODO Should it be timesteplength or currenttimesteplength
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
//  const double theta1 = TDatabase::TimeDB->THETA1;
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

  int nRhs = 4;  // number of rhs blocks - TODO for NSType 4 and 14, it is 4
  const TFESpace3D *v_space = &this->get_velocity_space();
  const TFESpace3D *p_space = &this->get_pressure_space();

  // TODO Implement the case NSType 4 and 14 where there's a 4th rhs block
  std::vector<double*> rhsArray(nRhs);
  rhsArray[0] = s.rhs_.block(0);
  rhsArray[1] = s.rhs_.block(1);
  rhsArray[2] = s.rhs_.block(2);
  rhsArray[3] = s.rhs_.block(3);

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

  /** After copy_nonactive, the solution vectors needs to be Comm-updated
     * in MPI-case in order to be consistently saved. It is necessary that
     * the vector is consistently saved because it is the only way to
     * ensure that its multiplication with an inconsistently saved matrix
     * (multiplication which appears in the defect and rhs computations)
     * give the correct results.
     * When we call copy_nonactive in MPI-case, we have to remember the following:
     * it can happen that some slave ACTTIVE DoFs are placed in the block of
     * NON-ACTIVE DoFs (because they are at the interface between processors).
     * Doing copy_nonactive changes then the value of these DOFs,although they are
     * actually active.
     * That's why we have to update the values so that the vector becomes consistent again.
     */
  #ifdef _MPI
    double *u1 = this->systems_.front().solution_.block(0);
    double *u2 = this->systems_.front().solution_.block(1);
    double *u3 = this->systems_.front().solution_.block(2);
    double *p  = this->systems_.front().solution_.block(3);
    this->systems_.front().parCommVelocity_.CommUpdate(u1);
    this->systems_.front().parCommVelocity_.CommUpdate(u2);
    this->systems_.front().parCommVelocity_.CommUpdate(u3);
    this->systems_.front().parCommPressure_.CommUpdate(p);
  #endif

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
  s.massMatrix_.apply_scaled_submatrix(old_solution_, s.rhs_, 3, 3, 1.0);
  // -tau*theta2 * A u^{k-1}
  double factor = -tau*theta2;
  s.matrix_.apply_scaled_submatrix(old_solution_, s.rhs_, 3, 3, factor);

  // scale the BT blocks with time step length
  for(System_per_grid& s : this->systems_)
  {
    if(tau != oldtau_)
    {
      // TODO: change the factor to be THETA1*tau; why??
      factor = /*theta1*/tau;
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

  // Reset old_residual_ for this time step iteration
  // otherwise, ones compares with the old_residual_ from
  // the previous time iteration, which is not correct.
  this->old_residual_ = FixedSizeQueue<10,Residuals>();

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

  Output::print<5>("Assembled the system matrix which will be passed to the ",
                   "solver");
}

/**************************************************************************** */
bool Time_NSE3D::stop_it(unsigned int iteration_counter)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  // compute, update and display defect and residuals
  compute_residuals();

  // stores current norm of the residual. They are normed per default in
  // the class Residuals
  const double normOfResidual        = this->get_full_residual();
//  const double normOfImpulseResidual = this->get_impulse_residual();
//  const double normOfMassResidual    = this->get_mass_residual();
//  const double oldNormOfResidual     = this->old_residual_[1].fullR;
  // hold the residual from up to 10 iterations ago
  const double veryOldNormOfResidual     = this->old_residual_.front().fullResidual;

//  Output::print("nonlinear step  :  " , setw(3), iteration_counter);
//  Output::print("impulse_residual:  " , setw(3), normOfImpulseResidual);
//  Output::print("mass_residual   :  " , setw(3), normOfMassResidual);
//  Output::print("residual        :  " , setw(3), normOfResidual);

  // this is the convergence ratio between actual step and last step
  // TODO : correct oldNormOfResidual to be the residual of last step
  // and print it out
//  if ( iteration_counter > 0 )
//  {
//  Output::print("convergence rate:  " , setw(3), normOfResidual/oldNormOfResidual);
//  }

  if(iteration_counter == 0)
    initial_residual_ = normOfResidual;

  // Parameters for stopping criteria (desired precision epsilon, max number
  // of iteration, convergence rate)
  double epsilon    = db_["nonlinloop_epsilon"];
  size_t max_It     = db_["nonlinloop_maxit"];
  double conv_speed = db_["nonlinloop_slowfactor"];
  bool slow_conv    = false;

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
    this->descale_matrices();
    return true;
   }
   else
    return false;
}

/**************************************************************************** */
void Time_NSE3D::compute_residuals()
{
  System_per_grid& s = this->systems_.front();
  unsigned int number_u_Dof = s.solution_.length(0);
  unsigned int number_p_Dof = s.solution_.length(3);

  // copy rhs to defect and compute defect
  this->defect_ = s.rhs_;
  s.matrix_.apply_scaled_add(s.solution_, defect_,-1.);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    IntoL20Vector3D(&defect_[3*number_u_Dof], number_p_Dof,
                      TDatabase::ParamDB->PRESSURE_SPACE);
  }

  // square norms of the residual components
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  std::vector<double> defect_m(defect_.get_entries_vector()); //copy of defect (entries)
  //Eliminate all non-master rows in defect_m!
  for(int ui = 0; ui < 3; ++ui)
  {//velocity rows
    const int* masters = s.parMapperVelocity_.GetMaster();
    for(size_t i = 0; i<number_u_Dof; ++i)
    {
      if (masters[i]!=my_rank)
      {
        defect_m[number_u_Dof*ui + i]=0;
      }
    }
  }
  {//pressure row
    const int* masters = s.parMapperPressure_.GetMaster();
    for(size_t i = 0; i<number_p_Dof; ++i)
    {
      if (masters[i]!=my_rank)
      {
        defect_m[number_u_Dof*3 + i]=0;
      }
    }
  }
  //TODO write this nicer (std!)
  double impulse_residual = Ddot(3*number_u_Dof, &defect_m.at(0),&defect_m.at(0));
  double mass_residual = Ddot(number_p_Dof, &defect_m.at(3*number_u_Dof),
                              &defect_m.at(3*number_u_Dof));
#else
//  //should not BlockVector be able to do vector*vector?
  double impulse_residual = Ddot(3*number_u_Dof,
                                 &this->defect_[0], &this->defect_[0]);

  double mass_residual    = Ddot(number_p_Dof,
                                 &this->defect_[3*number_u_Dof],
                                 &this->defect_[3*number_u_Dof]);
#endif

  Residuals current_residual(impulse_residual, mass_residual);
  old_residual_.add(current_residual);
}

/**************************************************************************** */
void Time_NSE3D::solve()
{
  System_per_grid& s = systems_.front();

  bool using_multigrid = solver_.is_using_multigrid();

  if(!using_multigrid)
  { // no multigrid
    if(solver_.get_db()["solver_type"].is("direct"))
    {
#ifndef _MPI
      solver_.solve(s.matrix_, s.rhs_, s.solution_);
#endif
#ifdef _MPI
      //two vectors of communicators (const for init, non-const for solving)
      std::vector<const TParFECommunicator3D*> par_comms_init =
      {&s.parCommVelocity_, &s.parCommVelocity_, &s.parCommVelocity_, &s.parCommPressure_};
      std::vector<TParFECommunicator3D*> par_comms_solv =
      {&s.parCommVelocity_, &s.parCommVelocity_, &s.parCommVelocity_, &s.parCommPressure_};

      //set up a MUMPS wrapper
      MumpsWrapper mumps_wrapper(s.matrix_, par_comms_init);

      //kick off the solving process
      mumps_wrapper.solve(s.rhs_, s.solution_, par_comms_solv);
#endif
    }
    else
      solver_.solve(s.matrix_, s.rhs_, s.solution_);
  }
  else  // multigrid preconditioned iterative solver
  {
    ErrThrow("No multigrid yet"); //TODO
  }

  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11, A22 and A33 matrices are
  // reset and assembled again but the non-diagonal blocks are scaled, so
  // for the next iteration we have to descale, see assemble_system()
  this->descale_matrices();

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       s.p_.project_into_L20();

  this->old_solution_ = s.solution_;
}

/**************************************************************************** */
void Time_NSE3D::descale_matrices()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;
  for(System_per_grid& s : this->systems_)
  {
    const FEMatrix& mass_blocks = *s.massMatrix_.get_blocks().at(0).get();
    s.matrix_.add_matrix_actives(mass_blocks, -1.0,
                                 {{0,0}, {1,1}, {2,2}},
                                 {false, false, false});
    const std::vector<std::vector<size_t>>
      cell_positions = {{0,0}, {0,1}, {0,2},
                        {1,0}, {1,1}, {1,2},
                        {2,0}, {2,1}, {2,2}};
    // note: declaring the auxiliary cell_positions is needed by the compiler
    // to sort out the overriding of the function scale_blocks_actives(...,...)
    s.matrix_.scale_blocks_actives(1./factor, cell_positions);
  }
}

/**************************************************************************** */
void Time_NSE3D::output(int m, int &image)
{
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
	bool no_output = !db_["output_write_vtk"] && !db_["output_compute_errors"];
	if(no_output)
		return;

  System_per_grid& s = this->systems_.front();
  TFEFunction3D* u1 = s.u_.GetComponent(0);
  TFEFunction3D* u2 = s.u_.GetComponent(1);
  TFEFunction3D* u3 = s.u_.GetComponent(2);

  if((size_t)db_["verbosity"]> 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    u3->PrintMinMax();
    s.p_.PrintMinMax();
  }

  if((m==0) || (m/TDatabase::TimeDB->STEPS_PER_IMAGE) )
  {
    if(db_["output_write_vtk"])
    {
      // last argument in the following is domain but is never used in this class
      TOutput3D output(5, 5, 2, 1, NULL);
      output.AddFEFunction(&s.p_);
      output.AddFEVectFunct(&s.u_);
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

//    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;

    // Errors in velocity components and pressure
    u1 ->GetErrors(example_.get_exact(0), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &v_space, err_u1);
    u2 ->GetErrors(example_.get_exact(1), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &v_space, err_u2);
    u3 ->GetErrors(example_.get_exact(2), 4, allderiv, 2, L2H1Errors, nullptr,
                  &aux, 1, &v_space, err_u3);
    s.p_.GetErrors(example_.get_exact(3), 4, allderiv, 2, L2H1Errors, nullptr,
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

    errors_[0] = err_u1[0]*err_u1[0] + err_u2[0]*err_u2[0] +
                 err_u3[0]*err_u3[0];  // (L2-norm)^2 for u
    errors_[1] = err_u1[1]*err_u1[1] + err_u2[1]*err_u2[1] +
                 err_u3[1]*err_u3[1];  // (H1-semi)^2 for u
    errors_[2] = err_p[0]*err_p[0];  // (L2-norm)^2 for p
    errors_[3] = err_p[1]*err_p[1];  // (H1-norm)^2 for p

    // TODO : CORRECT THE TIME-CORRECTED NORMS L2 AND H1 AND DISPLAY THEM
//    errors_[4] += (locerr[0]*locerr[0]+locerr[2]*locerr[2]
//                  + this->errors[0])*tau*0.5;
//    errors_[5] += (locerr[1]*locerr[1]+locerr[3]*locerr[3]
//                  + this->errors[1])*tau*0.5;
//    errors[6] += (locerr[0]*locerr[0] + this->errors[2])*tau*0.5;
//    errors[7] += (locerr[1]*locerr[1] + this->errors[3])*tau*0.5;

    // print errors
    if (my_rank == 0 )
    {
      Output::print<1>("L2(u)      : ", setprecision(10), sqrt(this->errors_[0]));
      Output::print<1>("H1-semi(u) : ", setprecision(10), sqrt(this->errors_[1]));
//    Output::print<1>("L2(0,t,L2(u)) : ", sqrt(this->errors[4]));
//    Output::print<1>("L2(0,t,H1-semi(u)) : ", sqrt(this->errors[5]));
      Output::print<1>("L2(p)      : ", setprecision(10), sqrt(this->errors_[2]));
      Output::print<1>("H1-semi(p)): ", setprecision(10), sqrt(this->errors_[3]));
//    Output::print<1>("L2(0,t,L2(p)) : ", sqrt(errors[6]) );
//    Output::print<1>("L2(0,t,H1-semi(p)) : ", sqrt(errors[7]) );
    }
  }
   delete u1;
   delete u2;
   delete u3;

   // do post-processing step depending on what the example implements, if needed
//   example_.do_post_processing(*this);

}

/**************************************************************************** */
const Residuals& Time_NSE3D::get_residuals() const
{
  return old_residual_.back();
}

/**************************************************************************** */
double Time_NSE3D::get_impulse_residual() const
{
  return old_residual_.back().impulsResidual;
}

/**************************************************************************** */
double Time_NSE3D::get_mass_residual() const
{
  return old_residual_.back().massResidual;
}

/**************************************************************************** */
double Time_NSE3D::get_full_residual() const
{
  return old_residual_.back().fullResidual;
}

/**************************************************************************** */
std::array< double, int(6) > Time_NSE3D::get_errors() const
{
  std::array<double, int(6)> error_at_time_points;
  error_at_time_points[0] = sqrt(this->errors_[0]); // L2 velocity error
  error_at_time_points[1] = sqrt(this->errors_[1]); // H1 velocity error
  error_at_time_points[2] = sqrt(this->errors_[2]); // L2 pressure error
  error_at_time_points[3] = sqrt(this->errors_[3]); // H1 pressure error

  return error_at_time_points;
}

/**************************************************************************** */
