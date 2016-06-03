#include <NSE3D.h>
#include <Database.h>
#include <LocalAssembling3D.h>
#include <Assemble3D.h>
#include <MainUtilities.h>

#include <NSE_MGLevel.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <NSE_MGLevel14.h>
#include <LinAlg.h>

#include <ItMethod.h>
#include <FgmresIte.h>
#include <FixedPointIte.h>
#include <MultiGridIte.h>

#include <DirectSolver.h>
#include <Output3D.h>

#include <GridTransfer.h>
#include <Multigrid.h>
#include <Upwind3D.h>

#include <sys/stat.h>

ParameterDatabase get_default_NSE3D_parameters()
{
  Output::print<3>("creating a default NSE3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("NSE3D parameter database");
  
  //NSE3D requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  //stokes case - reduce no nonlin its TODO remove global database dependency
  if (TDatabase::ParamDB->PROBLEM_TYPE == 3)
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

NSE3D::System_per_grid::System_per_grid(const Example_NSE3D& example,
                                    TCollection& coll, std::pair<int, int> order, 
                                    NSE3D::Matrix type
#ifdef _MPI
                                    , int maxSubDomainPerDof
#endif
) :  velocitySpace_(&coll, (char*)"u", (char*)"nse3d velocity", example.get_bc(0), //bd cond at 0 is x velo bc
                    order.first),
     pressureSpace_(&coll, (char*)"p", (char*)"nse3d pressure", example.get_bc(3), //bd condition at 3 is pressure bc
                    order.second),
     matrix_({&velocitySpace_, &velocitySpace_, &velocitySpace_, &pressureSpace_}),
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
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      matrix_ = BlockFEMatrix::NSE3D_Type1(velocitySpace_, pressureSpace_);
      break;                                               
    case 2:                                                
      matrix_ = BlockFEMatrix::NSE3D_Type2(velocitySpace_, pressureSpace_);
      break;                                               
    case 3:                                                
      matrix_ = BlockFEMatrix::NSE3D_Type3(velocitySpace_, pressureSpace_);
      break;                                               
    case 4:                                                
      matrix_ = BlockFEMatrix::NSE3D_Type4(velocitySpace_, pressureSpace_);
      break;                                               
    case 14:                                               
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

void NSE3D::output_problem_size_info() const
{
    const TFESpace3D & velocity_space = this->systems_.front().velocitySpace_;
    const TFESpace3D & pressure_space = this->systems_.front().pressureSpace_;

    size_t nDofu  = velocity_space.GetN_DegreesOfFreedom();
    size_t nDofp  = pressure_space.GetN_DegreesOfFreedom();
    size_t nTotal = 3*nDofu + nDofp;
    size_t nActive= 3*velocity_space.GetActiveBound();

    TCollection* coll = velocity_space.GetCollection();

    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);

    Output::stat("NSE3D", "Mesh data and problem size");
    Output::dash("N_Cells      :  ", setw(10), coll->GetN_Cells());
    Output::dash("h(min, max)  :  ", setw(10), hmin, setw(10), " ", hmax);
    Output::dash("ndof Velocity:  ", setw(10), 3*nDofu );
    Output::dash("ndof Pressure:  ", setw(10), nDofp);
    Output::dash("ndof Total   :  ", setw(10), nTotal );
    Output::dash("nActive      :  ", setw(10), nActive);
}

NSE3D::NSE3D(const TDomain& domain, const ParameterDatabase& param_db,
             const Example_NSE3D& example
#ifdef _MPI
             , int maxSubDomainPerDof
#endif
) : systems_(), example_(example), db(get_default_NSE3D_parameters()),
    solver(param_db), mg_(nullptr),
    defect_(), old_residuals_(),
    initial_residual_(1e10), errors_()
{
  this->db.merge(param_db, false);
  std::pair <int,int> 
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and pressure spaces
  // this function returns a pair which consists of 
  // velocity and pressure order
  this->get_velocity_pressure_orders(velocity_pressure_orders);

  NSE3D::Matrix type;
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:  type = NSE3D::Matrix::Type1;
      break;
    case 2:  type = NSE3D::Matrix::Type2;
      break; 
    case 3:  type = NSE3D::Matrix::Type3; 
      break; 
    case 4:  type = NSE3D::Matrix::Type4;
      break;
    case 14: type = NSE3D::Matrix::Type14;
      break;
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }
  
  bool usingMultigrid = solver.is_using_multigrid();
  if(!usingMultigrid)
  {
    TCollection *coll = domain.GetCollection(It_Finest, 0, -4711);
        
    #ifdef _MPI
    // create finite element space and function, a matrix, rhs, and solution
    systems_.emplace_back(example_, *coll, velocity_pressure_orders, type,
                          maxSubDomainPerDof);
    #else
    // create finite element space and function, a matrix, rhs, and solution
    systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);
    #endif

  }
  else // multigrid
  {
#ifdef _MPI
    ErrThrow("There is no multigrid for NSE3D in MPI yet!");
#endif

    ParameterDatabase database_mg = Multigrid::default_multigrid_database();
    database_mg.merge(param_db, false);

    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    size_t n_levels = database_mg["multigrid_n_levels"];

    std::string mgtype_str = database_mg["multigrid_type"];
    MultigridType mgtype = string_to_multigrid_type(mgtype_str);

    int finest = 0;
    int coarsest = 0;
    if(mgtype == MultigridType::STANDARD)
    {
      finest = domain.get_ref_level();
      coarsest = finest - n_levels + 1;
    }
    else
    if(mgtype == MultigridType::MDML)
    {
      finest = domain.get_ref_level();
      coarsest = finest - n_levels + 2;
      //do the finest algebraic grid in advance
      TCollection *coll = domain.GetCollection(It_EQ, finest, -4711);
#ifdef _MPI
      systems_.emplace_back(example_, *coll, velocity_pressure_orders, type,
                            maxSubDomainPerDof);
#else
      systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);
#endif
      matrices.push_front(&systems_.back().matrix_);
      // set velo/pressure orders to lowest order nonconforming
      velocity_pressure_orders = {-1, 0};
    }
    if(coarsest < 0)
      ErrThrow("More multigrid levels (",n_levels,") than possible due to "
          "refinement and multigrid type requested.");

    for (int grid_no = finest; grid_no >= coarsest; --grid_no)
    {
      TCollection *coll = domain.GetCollection(It_EQ, grid_no, -4711);
#ifdef _MPI
      // create finite element space and function, a matrix, rhs, and solution
      systems_.emplace_back(example_, *coll, velocity_pressure_orders, type,
                            maxSubDomainPerDof);
#else
      // create finite element space and function, a matrix, rhs, and solution
      systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);
#endif
      //prepare input argument for multigrid object
      matrices.push_front(&systems_.back().matrix_);
    }

    // Construct multigrid object
    mg_ = std::make_shared<Multigrid>(database_mg, matrices, mgtype);
  }

  output_problem_size_info();

}

void NSE3D::check_parameters()
{

  // Some implementation/testing constraints on the used discretization.
  if(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE != 0)
  {
    ErrThrow("Every SC_NONLIN_ITE_TYPE_SADDLE except 0 is untested!");
  }
  if(TDatabase::ParamDB->LAPLACETYPE != 0)
  {
    ErrThrow("Every LAPLACETYPE except 0 is untested!");
  }
  if(TDatabase::ParamDB->DISCTYPE != 1)
  {
    ErrThrow("Every DISCTYPE except 1 is untested!");
  }
  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM != 0)
  {
    ErrThrow("Every NSE_NONLINEAR_FORM except 0 is untested!");
  }

}

void NSE3D::get_velocity_pressure_orders(std::pair< int, int >& velocity_pressure_orders)
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
          Output::warn("NSE3D", "The P1/P0 element pair (Q1/Q0 on hexa) is "
              " not stable. Make sure to use stabilization!");
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
}

void NSE3D::assemble_linear_terms()
{
  size_t nFESpace = 2; // spaces used for assembling matrices
  size_t nSqMatrices=10; // maximum no of square matrices (type 14)
  std::vector<TSquareMatrix3D*> sqMatrices(nSqMatrices);
  size_t nReMatrices = 6; // maximum no of rectangular matrices (e.g. type 14)
  std::vector<TMatrix3D*> reMatrices(nReMatrices);
  size_t nRhs=4; // maximum number of right hand sides (e.g. type 4)
  std::vector<double*> rhsArray(nRhs); // right hand side 
  // finite element function used for nonlinear term
  std::vector<TFEFunction3D*> feFunction(4, nullptr);
  
  for(auto &s : this->systems_)
  {

    const TFESpace3D *v_space = &s.velocitySpace_;
    const TFESpace3D *p_space = &s.pressureSpace_;

    // spaces for matrices
    const TFESpace3D *spaces[2] = {v_space, p_space};
    const TFESpace3D *rhsSpaces[4] = {v_space, v_space, v_space, p_space};

    // spaces for right hand side    
    s.rhs_.reset();
    rhsArray[0]=s.rhs_.block(0);
    rhsArray[1]=s.rhs_.block(1);
    rhsArray[2]=s.rhs_.block(2);
    rhsArray[3]=nullptr; //will be reset for type 4 and 14
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_uniquely();
    
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        nSqMatrices = 1;
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        
        nReMatrices = 3;
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[1].get());
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[2].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[3].get());

        nRhs = 3;
        break;
      case 2:
        nSqMatrices = 1;
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        
        nReMatrices = 6;
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[4].get());
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[5].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[6].get());
        reMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks[1].get());
        reMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks[2].get());
        reMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks[3].get());

        nRhs = 3;
        break;
      case 3:
        nSqMatrices = 9;
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
        sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
        sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
        sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
        sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
        sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
        sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
        sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());
        
        nReMatrices = 3;
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[11].get());

        nRhs = 3;
        break;
      case 4:
        nSqMatrices = 9;
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
        sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
        sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
        sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
        sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
        sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
        sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
        sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());
        
        nReMatrices = 6;
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[12].get()); //first the lying B blocks
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[13].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[14].get());
        reMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks[3].get()); //than the standing B blocks
        reMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
        reMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks[11].get());

        //right hand side must be adapted
        nRhs = 4;
        rhsArray[3]=s.rhs_.block(3);
        break;
      case 14:        
        nSqMatrices = 10;
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
        sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
        sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
        sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
        sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
        sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
        sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
        sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());
        sqMatrices[9]=reinterpret_cast<TSquareMatrix3D*>(blocks[15].get());

        nReMatrices = 6;
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[12].get());
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[13].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[14].get());
        reMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
        reMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
        reMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks[11].get());
        
        //right hand side must be adapted
        nRhs = 4;
        rhsArray[3]=s.rhs_.block(3);
        break;
    }// endswitch nstype

    for(unsigned int i=0; i<nSqMatrices; i++)
      sqMatrices[i]->reset();
    for(unsigned int i=0; i<nReMatrices;i++)
      reMatrices[i]->reset();

    // boundary conditions and boundary values
    BoundCondFunct3D * boundContion[4]={
      spaces[0]->getBoundCondition(), spaces[0]->getBoundCondition(),
      spaces[0]->getBoundCondition(), spaces[1]->getBoundCondition() };
      
    std::array<BoundValueFunct3D*, 4> boundValues;
    boundValues[0]=example_.get_bd()[0];
    boundValues[1]=example_.get_bd()[1];
    boundValues[2]=example_.get_bd()[2];
    boundValues[3]=example_.get_bd()[3];

    // finite element functions
    feFunction[0]=s.u_.GetComponent(0);
    feFunction[1]=s.u_.GetComponent(1);
    feFunction[2]=s.u_.GetComponent(2);
    feFunction[3]=&s.p_;

    // local assembling object    
    const LocalAssembling3D la(LocalAssembling3D_type::NSE3D_Linear, 
                         feFunction.data(), example_.get_coeffs());
    
    // assemble now the matrices and right hand side 
    Assemble3D(nFESpace, spaces, 
               nSqMatrices, sqMatrices.data(),
               nReMatrices, reMatrices.data(), 
               nRhs, rhsArray.data(), rhsSpaces,
               boundContion, boundValues.data(), la);
  }// endfor auto grid

  //copy non-actives from rhs to solution on finest grid
  this->systems_.front().solution_.copy_nonactive(systems_.front().rhs_);
}

void NSE3D::assemble_non_linear_term()
{
  size_t nFESpace = 1; // space needed for assembling matrices
  size_t nSqMatrices=3; // no of square matrices (Maximum 3)
  std::vector<TSquareMatrix3D*> sqMatrices(nSqMatrices);
  size_t nReMatrices = 0; // no rectangular matrix in nonlinear iteration
  std::vector<TMatrix3D*> reMatrices{nullptr};
  size_t nRhs=0; // no right hand side to be assembled
  std::vector<double*> rhsArray{nullptr};  
  std::vector<TFEFunction3D*> feFunction(3);
  const TFESpace3D** rhsSpaces{nullptr};
  
  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(systems_.size() > 1)
  {
    for( int block = 0; block < 3 ;++block)
    {
      std::vector<const TFESpace3D*> spaces;
      std::vector<double*> u_entries;
      std::vector<size_t> u_ns_dofs;
      for(auto &s : systems_ )
      {
        spaces.push_back(&s.velocitySpace_);
        u_entries.push_back(s.solution_.block(block));
        u_ns_dofs.push_back(s.solution_.length(block));
      }
      GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
    }
  }

  for(auto &s : this->systems_)
  {
    // spaces for matrices
    const TFESpace3D* spaces[1] = {&s.velocitySpace_};
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_uniquely({{0,0},{1,1},{2,2}});

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        nSqMatrices = 1;
        sqMatrices.resize(nSqMatrices);
        sqMatrices.at(0)=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        break;
      case 3:
      case 4:
      case 14:
        nSqMatrices = 3;
        sqMatrices.at(0)=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        sqMatrices.at(1)=reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
        sqMatrices.at(2)=reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
        break;
    }// endswitch nstype

    for(auto mat : sqMatrices)
    {
      mat->reset();
    }

    // boundary conditions and boundary values
    BoundCondFunct3D * boundCondition[1]={
      spaces[0]->getBoundCondition() };

    std::array<BoundValueFunct3D*, 3> boundValues;
    boundValues[0]=example_.get_bd()[0];
    boundValues[1]=example_.get_bd()[1];
    boundValues[2]=example_.get_bd()[2];

    feFunction[0]=s.u_.GetComponent(0);
    feFunction[1]=s.u_.GetComponent(1);
    feFunction[2]=s.u_.GetComponent(2);

    //decide wether to assemble by upwinding or not
    bool finest_grid = (&s == &systems_.at(0));

    bool upwinding = false;
    if(mg_)
      upwinding = mg_->get_type() == MultigridType::MDML && !finest_grid;

    if(!upwinding)
    {
    // local assembling object    
    const LocalAssembling3D la(LocalAssembling3D_type::NSE3D_NonLinear, 
                         feFunction.data(), example_.get_coeffs());
    
    // assemble now the matrices and right hand side 
    Assemble3D(nFESpace, spaces, 
               nSqMatrices, sqMatrices.data(),
               nReMatrices, reMatrices.data(), 
               nRhs, rhsArray.data(), rhsSpaces,
               boundCondition, boundValues.data(), la);
    }
    else
    {
      double one_over_nu = 1/example_.get_nu(); //the inverse of the example's diffusion coefficient
      for(auto mat : sqMatrices)
      {
        UpwindForNavierStokes3D(
        mat, feFunction[0],feFunction[1],feFunction[2], one_over_nu);
      }
    }

    //TODO: Copying non-actives??

  }// endfor auto grid


}

bool NSE3D::stop_it(unsigned int iteration_counter)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  //compute and update defect and residuals
  compute_residuals();
  
  // the current norm of the residual
  const double normOfResidual = this->get_full_residual();
  // store initial residual, so later we can print the overall reduction
  if(iteration_counter == 0)
    initial_residual_ = normOfResidual;

  // hold the residual from 10 iterations ago
  const double oldNormOfResidual = this->old_residuals_.front().fullResidual;


  size_t max_it = db["nonlinloop_maxit"];
  double conv_speed = db["nonlinloop_slowfactor"];
  bool slow_conv = false;


  if(normOfResidual >= conv_speed*oldNormOfResidual)
    slow_conv = true;

  double limit = db["nonlinloop_epsilon"];
  if (db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= sqrt(this->get_size());
    if(my_rank==0)
     Output::print<1>("stopping tolerance for nonlinear iteration ", limit);
  }

  // check if the iteration has converged, or reached the maximum number of
  // iterations or if convergence is too slow. Then return true otherwise false
  if( (normOfResidual<=limit) || (iteration_counter==max_it) || (slow_conv) )
  {
    if(slow_conv && my_rank==0)
      Output::print<1>(" SLOW !!! ", normOfResidual/oldNormOfResidual);

    // stop iteration
    if(my_rank==0)
    {
      Output::print<1>("\nNonlinear Iterations: ", setw(4), iteration_counter, setprecision(8),
                       " RES : ", normOfResidual, " Reduction : ",
                       normOfResidual/initial_residual_);
    }

    return true;
  }
  else
    return false;
}

/** ************************************************************************ */
void NSE3D::compute_residuals()
{
  System_per_grid& s = this->systems_.front();
  unsigned int n_u_dof = s.solution_.length(0);
  unsigned int n_p_dof = s.solution_.length(3);

  // copy rhs to defect and compute defect
  defect_ = s.rhs_;
  s.matrix_.apply_scaled_add(s.solution_, defect_,-1.);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    IntoL20Vector3D(&defect_[3*n_u_dof], n_p_dof,
                    TDatabase::ParamDB->PRESSURE_SPACE);
  }

  // square of the norms of the residual components
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  std::vector<double> defect_m(defect_.get_entries_vector()); //copy of defect (entries)
  //Eliminate all non-master rows in defect_m!
  for(int ui = 0; ui < 3; ++ui)
  {//velocity rows
    const int* masters = s.parMapperVelocity_.GetMaster();
    for(size_t i = 0; i<n_u_dof; ++i)
    {
      if (masters[i]!=my_rank)
      {
        defect_m[n_u_dof*ui + i]=0;
      }
    }
  }
  {//pressure row
    const int* masters = s.parMapperPressure_.GetMaster();
    for(size_t i = 0; i<n_p_dof; ++i)
    {
      if (masters[i]!=my_rank)
      {
        defect_m[n_u_dof*3 + i]=0;
      }
    }
  }
  //TODO write this nicer (std!)
  double impuls_Residual = Ddot(3*n_u_dof, &defect_m.at(0),&defect_m.at(0));
  double mass_residual = Ddot(n_p_dof, &defect_m.at(3*n_u_dof),
                              &defect_m.at(3*n_u_dof));
#else
  //should not BlockVector be able to do vector*vector?
  double impuls_Residual = Ddot(3*n_u_dof, &this->defect_[0],&this->defect_[0]);
  double mass_residual = Ddot(n_p_dof, &this->defect_[3*n_u_dof],
                              &this->defect_[3*n_u_dof]);
#endif

  Residuals current_residuals(impuls_Residual, mass_residual);
  old_residuals_.add(current_residuals);
}

void NSE3D::solve()
{
  System_per_grid& s = this->systems_.front();
  double damping = this->db["nonlinloop_damping_factor"];
  // store previous solution for damping, it is a pointer so that we can avoid
  // the copy in case of no damping
  std::shared_ptr<BlockVector> old_solution(nullptr);
  if(damping != 1.0)
    old_solution = std::make_shared<BlockVector>(s.solution_);
  
  //determine whether we make use of multigrid
  bool using_multigrid = solver.is_using_multigrid();
  if(!using_multigrid)
  {//no multigrid
    if(this->solver.get_db()["solver_type"].is("direct"))
    {
#ifndef _MPI
      this->solver.solve(s.matrix_, s.rhs_, s.solution_);
#endif
#ifdef _MPI
      if(damping != 1.0)
        Output::warn("NSE3D::solve", "damping in an MPI context is not tested");
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
      this->solver.solve(s.matrix_, s.rhs_, s.solution_);
  }
  else
  {//multigrid preconditioned iterative solver
    solver.solve(s.matrix_, s.rhs_, s.solution_, mg_);
  }

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
   s.p_.project_into_L20();
  }
  
  if(damping != 1.0)
  {
    s.solution_.scale(damping);
    s.solution_.add_scaled(*old_solution, damping);
  }
}

void NSE3D::output(int i)
{
#ifdef _MPI
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
  if(no_output)
    return;
  
  System_per_grid& s=this->systems_.front();
  TFEFunction3D* u1 = s.u_.GetComponent(0);
  TFEFunction3D* u2 = s.u_.GetComponent(1);
  TFEFunction3D* u3 = s.u_.GetComponent(2);
  
  if((size_t)db["verbosity"]> 1)
  {
    u1->PrintMinMax(std::string("u1"));
    u2->PrintMinMax(std::string("u2"));
    u3->PrintMinMax(std::string("u3"));
    s.p_.PrintMinMax(std::string("p"));
  }
  
  // write solution to a vtk file
  if(db["output_write_vtk"])
  {
    // last argument in the following is domain, but is never used in this class
    TOutput3D Output(5, 5, 2, 1, NULL);
    Output.AddFEFunction(&s.p_);
    Output.AddFEVectFunct(&s.u_);
#ifdef _MPI
    char SubID[] = "";
    if(my_rank == 0)
  	  mkdir(db["output_vtk_directory"], 0777);
    std::string dir = db["output_vtk_directory"];
    std::string base = db["output_basename"];
    Output.Write_ParVTK(MPI_COMM_WORLD, 0, SubID, dir, base);
#else
    // Create output directory, if not already existing.
    mkdir(db["output_vtk_directory"], 0777);
    std::string filename = this->db["output_vtk_directory"];
    filename += "/" + this->db["output_basename"].value_as_string();

    if(i >= 0)
      filename += "_" + std::to_string(i);
    filename += ".vtk";
    Output.WriteVtk(filename.c_str());
#endif
  }
  
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    double err_u1[4]; // of these arrays only the two first entries are used,
    double err_u2[4]; // but the evil GetErrors() will corrupt memory if these
    double err_u3[4]; // have not at least size 4
    double err_p[4];


    TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    MultiIndex3D nsAllDerivs[4] = {D000, D100, D010, D001};
    const TFESpace3D *velocity_space = &this->get_velocity_space();
    const TFESpace3D *pressure_space = &this->get_pressure_space();
    
    // errors in first velocity component
    u1->GetErrors(example_.get_exact(0), 4, nsAllDerivs, 2,
                  L2H1Errors, nullptr, &aux, 1, &velocity_space, err_u1);
    // errors in second velocity component
    u2->GetErrors(example_.get_exact(1), 4, nsAllDerivs, 2,
                  L2H1Errors, nullptr, &aux, 1, &velocity_space, err_u2);
    // errors in third velocity component
    u3->GetErrors(example_.get_exact(2), 4, nsAllDerivs, 2,
                  L2H1Errors, nullptr, &aux, 1, &velocity_space, err_u3);
    // errors in pressure
    s.p_.GetErrors(example_.get_exact(3), 4, nsAllDerivs, 2, L2H1Errors,
                   nullptr, &aux, 1, &pressure_space, err_p);
    
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
    for(i=0;i<8;i++)
    {//MPI: sqrt was skipped in GetErrors function - do it here globally!
      err_red[i] = sqrt(err_red[i]);
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
    int my_rank =0;
#endif

    errors_.at(0) = sqrt(err_u1[0]*err_u1[0] + err_u2[0]*err_u2[0] + err_u3[0]*err_u3[0]);//L2
    errors_.at(1) = sqrt(err_u1[1]*err_u1[1] + err_u2[1]*err_u2[1] + err_u3[1]*err_u3[1]);//H1-semi
    errors_.at(2) = err_p[0];
    errors_.at(3) = err_p[1];

    //print errors
    if(my_rank == 0)
    {
      Output::stat("NSE3D", "Measured errors");
      Output::dash("L2(u)     : ", setprecision(10), errors_.at(0));
      Output::dash("H1-semi(u): ", setprecision(10), errors_.at(1));
      Output::dash("L2(p)     : ", setprecision(10), errors_.at(2));
      Output::dash("H1-semi(p): ", setprecision(10), errors_.at(3));
    }
  } // if(this->db["compute_errors"])
  delete u1;
  delete u2;
  delete u3;

  //do postprocessing step depending on what the example implements
  example_.do_post_processing(*this);
}

TFEFunction3D* NSE3D::get_velocity_component(int i)
{
  if(i==0)
    return this->systems_.front().u_.GetComponent(0);
  else if(i==1)
    return this->systems_.front().u_.GetComponent(1);
  else  if(i==2)
    return this->systems_.front().u_.GetComponent(2);
  else
    throw std::runtime_error("There are only three velocity components!");
}

const Residuals& NSE3D::get_residuals() const
{
  return old_residuals_.back();
}

double NSE3D::get_impuls_residual() const
{
  return old_residuals_.back().impulsResidual;
}

double NSE3D::get_mass_residual() const
{
  return old_residuals_.back().massResidual;
}

double NSE3D::get_full_residual() const
{
  return old_residuals_.back().fullResidual;
}

std::array<double, int(4)> NSE3D::get_errors() const
{
  return errors_;
}
