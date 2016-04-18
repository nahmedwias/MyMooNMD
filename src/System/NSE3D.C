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

#endif
}


NSE3D::NSE3D(const TDomain& domain, const Example_NSE3D& example
#ifdef _MPI
             , int maxSubDomainPerDof
#endif
) : systems_(), example_(example), multigrid_(nullptr),
    defect_(), old_residuals_(), initial_residual_(1e10), errors_()
{
  std::pair <int,int> 
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and pressure spaces
  // this function returns a pair which consists of 
  // velocity and pressure order
  this->get_velocity_pressure_orders(velocity_pressure_orders);
  // start with only non-multigrid case

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
  
  bool usingMultigrid = (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5 
                        && TDatabase::ParamDB->SOLVER_TYPE ==1);
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
    
    const TFESpace3D & velocity_space = this->systems_.front().velocitySpace_;
    const TFESpace3D & pressure_space = this->systems_.front().pressureSpace_;
    
    size_t nDofu  = velocity_space.GetN_DegreesOfFreedom();
    size_t nDofp  = pressure_space.GetN_DegreesOfFreedom();
    size_t nTotal = 3*nDofu + nDofp;
    size_t nActive= 3*velocity_space.GetActiveBound();
    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);
    
    Output::print<1>("N_Cells      :  ", setw(10), coll->GetN_Cells());
    Output::print<1>("h(min, max)  :  ",setw(10), hmin, setw(10), " ", hmax);
    Output::print<1>("ndof Velocity:  ", setw(10), 3*nDofu );
    Output::print<1>("ndof Pressure:  ", setw(10), nDofp);
    Output::print<1>("ndof Total   :  ", setw(10), nTotal );
    Output::print<1>("nActive      :  ", setw(10), nActive);

    #endif
  }
  else // multigrid
  {
    size_t n_levels = TDatabase::ParamDB->LEVELS;
    std::vector<TCollection*> collections(n_levels,nullptr);
    for(size_t i =0 ; i< n_levels ; ++i)
    {
      ErrThrow("This loop for multigrid NSE3D is not tested and thus most likely incorrect!");
      collections.at(i)=domain.GetCollection(It_EQ, i, -4711);
    }

    std::vector<double> param(2);
    param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
    param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
    
    this->multigrid_.reset(new TNSE_MultiGrid(1, 2, param.data()));
    this->transposed_B_structures_.resize(n_levels,nullptr);
    
    // constructing systems per grid
    for(auto it : collections)
    {
      #ifdef _MPI
        systems_.emplace_back(example, *it, velocity_pressure_orders, 
                              type, maxSubDomainPerDof);
      #else
        systems_.emplace_back(example, *it, velocity_pressure_orders, 
                              type);
      #endif
    }
    //CLEMENS: this is just for me to be sure about some information
    //You can delete it or modify somehow: b/c the same lines of code 
    // is used above in the case where multigrid is not used 
    {
      TCollection& cellCollection = *collections.front();
      Output::print<1>("N_Cells      :  ", setw(10), cellCollection.GetN_Cells());
      const TFESpace3D & velocity_space = this->systems_.front().velocitySpace_;
      const TFESpace3D & pressure_space = this->systems_.front().pressureSpace_;
      
      size_t nDofu  = velocity_space.GetN_DegreesOfFreedom();
      size_t nDofp  = pressure_space.GetN_DegreesOfFreedom();
      size_t nTotal = 3*nDofu + nDofp;
      size_t nActive= 3*velocity_space.GetActiveBound();
      
      Output::print<1>("ndof Velocity:  ", setw(10), 3*nDofu );
      Output::print<1>("ndof Pressure:  ", setw(10), nDofp);
      Output::print<1>("ndof Total   :  ", setw(10), nTotal );
      Output::print<1>("nActive      :  ", setw(10), nActive);
      double hmin, hmax;
      cellCollection.GetHminHmax(&hmin, &hmax);
      Output::print<1>("h(min, max)  :  ",setw(10), hmin, setw(10), " ", hmax);
    }
    size_t level = 0;
    //Create multigrid-level-objects and add them to the multgrid object.
    // Must be coarsest level first, therefore reverse order iteration.
    for(auto system=systems_.rbegin(); system != systems_.rend(); ++system)
    {
      #ifdef _MPI
        ErrThrow("There is no multigrid in MPI!");
      #else
        multigrid_->AddLevel(this->mg_levels(level, *system));
      #endif
      level++;      
    }
  }
}

void NSE3D::check_parameters()
{
  // this has to do with the relation of UNIFORM_STEPS and LEVELS
  // copied from CD3D, it should actually be unified
  bool usingMultigrid = TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5
                        && TDatabase::ParamDB->SOLVER_TYPE == 1;
  if (!usingMultigrid)
  { //case of direct solve or non-multigrid iterative solve
    if (TDatabase::ParamDB->LEVELS < 1)
    {
      ErrThrow("Parameter LEVELS must be greater or equal 1.");
    }
    TDatabase::ParamDB->UNIFORM_STEPS += TDatabase::ParamDB->LEVELS -1;
    TDatabase::ParamDB->LEVELS = 1;
    Output::print("Non-multigrid solver chosen. Therefore LEVELS -1 was added "
        "to UNIFORM_STEPS and LEVELS set to 1. \n Now: UNIFORM_STEPS = ",
        TDatabase::ParamDB->UNIFORM_STEPS, ".");
  }
  else
  {  // iterative solve with multigrid prec
    if (TDatabase::ParamDB->LEVELS < 2)
    {
      ErrThrow("Parameter LEVELS must be at least 2 for multigrid.");
    }
  }


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
          Output::print<1>("Warning: The P1/P0 element pair (Q1/Q0 on hexa) is "
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
    // local assembling object    
    const LocalAssembling3D la(LocalAssembling3D_type::NSE3D_NonLinear, 
                         feFunction.data(), example_.get_coeffs());
    
    // assemble now the matrices and right hand side 
    Assemble3D(nFESpace, spaces, 
               nSqMatrices, sqMatrices.data(),
               nReMatrices, reMatrices.data(), 
               nRhs, rhsArray.data(), rhsSpaces,
               boundCondition, boundValues.data(), la);

    //TODO: UPWINDING??
    //TODO: Copying non-actives??

  }// endfor auto grid


}

bool NSE3D::stop_it(unsigned int iteration_counter)
{
  //compute and update defect and residuals
  compute_residuals();
  
  // the current norm of the residual
  const double normOfResidual = this->get_full_residual();
  // store initial residual, so later we can print the overall reduction
  if(iteration_counter == 0)
    initial_residual_ = normOfResidual;

  // hold the residual from 10 iterations ago
  const double oldNormOfResidual = this->old_residuals_.front().fullResidual;


  size_t max_it = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  double conv_speed = TDatabase::ParamDB->SC_NONLIN_DIV_FACTOR;
  bool slow_conv = false;


  if(normOfResidual >= conv_speed*oldNormOfResidual)
    slow_conv = true;

  double limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  if (TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE)
  {
    limit *= sqrt(this->get_size());
    Output::print<1>("stopping tolerance for nonlinear iteration ", limit);
  }

  // check if the iteration has converged, or reached the maximum number of
  // iterations or if convergence is too slow. Then return true otherwise false
  if( (normOfResidual<=limit) || (iteration_counter==max_it) || (slow_conv) )
  {
    if(slow_conv)
      Output::print<1>(" SLOW !!! ", normOfResidual/oldNormOfResidual);

    // stop iteration
    Output::print<1>("\nNonlinear Iterations: ", setw(4), iteration_counter, setprecision(8),
                     " RES : ", normOfResidual, " Reduction : ",
                     normOfResidual/initial_residual_);
    // The following line comes from MooNMD and shall be us a reminder to
    // TODO count total number of linear iterations for iterative solvers
    //if(TDatabase::ParamDB->SOLVER_TYPE != 2) // not using direct solver
    //  OutPut(" Linear Iterations Total: " << this->n_linear_iterations);

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

  bool using_multigrid = //determine whether we make use of multigrid
      TDatabase::ParamDB->SOLVER_TYPE == 1 &&
      TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5;

  if(!using_multigrid)
  {//no multigrid
    if(TDatabase::ParamDB->SOLVER_TYPE == 1)
    {
      ErrThrow("You chose a non-multigrid iterative solver. This not available for NSE3D so far.");
    }
    else if(TDatabase::ParamDB->SOLVER_TYPE == 2)
    {
#ifndef _MPI
      // So far only UMFPACK is available as direct solver in sequential case.
      // Actuate it via DirectSolver class.

      /// @todo consider storing an object of DirectSolver in this class
      DirectSolver direct_solver(s.matrix_,
                                 DirectSolver::DirectSolverTypes::umfpack);
      direct_solver.solve(s.rhs_, s.solution_);
#elif _MPI
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
  }
  else
  {//multigrid preconditioned iterative solver is used
    mg_solver();
  }

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
   s.p_.project_into_L20();
  }
}

void NSE3D::output(int i)
{
  if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;
  
  System_per_grid& s=this->systems_.front();
  TFEFunction3D* u1 = s.u_.GetComponent(0);
  TFEFunction3D* u2 = s.u_.GetComponent(1);
  TFEFunction3D* u3 = s.u_.GetComponent(2);
  
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    u3->PrintMinMax();
    s.p_.PrintMinMax();
  }
  
  // write solution to a vtk file
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    // last argument in the following is domain, but is never used in this class
    TOutput3D Output(5, 5, 2, 1, NULL);
    Output.AddFEFunction(&s.p_);
    Output.AddFEVectFunct(&s.u_);
#ifdef _MPI
    char SubID[] = "";
    Output.Write_ParVTK(MPI_COMM_WORLD, 0, SubID);
#else
    std::string filename(TDatabase::ParamDB->OUTPUTDIR);
    filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
    if(i >= 0)
      filename += "_" + std::to_string(i);
    filename += ".vtk";
    Output.WriteVtk(filename.c_str());
#endif
  }
  
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(TDatabase::ParamDB->MEASURE_ERRORS)
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
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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
      Output::print("");
      Output::print<1>("L2(u)     : ", setprecision(10), errors_.at(0));
      Output::print<1>("H1-semi(u): ", setprecision(10), errors_.at(1));
      Output::print<1>("L2(p)     : ", setprecision(10), errors_.at(2));
      Output::print<1>("H1-semi(p): ", setprecision(10), errors_.at(3));
    }
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
  delete u1;
  delete u2;
  delete u3;
}

TNSE_MGLevel* NSE3D::mg_levels(int level, NSE3D::System_per_grid& s)
{
  TNSE_MGLevel * multigridLevel;
  size_t nAuxArray = 2;
  std::vector<double> alpha(2);
  
  if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
        || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE))
     nAuxArray=4;
  else
     nAuxArray=2;
  
  if (level==0)
  {
    alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
    alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
  }
  else
  {
    alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
    alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
  }
  
  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_TERRIBLY_UNSAFE();
  
  TSquareMatrix3D* A11 = reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
  TSquareMatrix3D* A12 = reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
  TSquareMatrix3D* A13 = reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
  TSquareMatrix3D* A21 = reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
  TSquareMatrix3D* A22 = reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
  TSquareMatrix3D* A23 = reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
  TSquareMatrix3D* A31 = reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
  TSquareMatrix3D* A32 = reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
  TSquareMatrix3D* A33 = reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());  
  
  TMatrix3D* B1T = reinterpret_cast<TMatrix3D*>(blocks[3].get());
  TMatrix3D* B2T = reinterpret_cast<TMatrix3D*>(blocks[7].get());
  TMatrix3D* B3T = reinterpret_cast<TMatrix3D*>(blocks[11].get());
  TMatrix3D* B1 = reinterpret_cast<TMatrix3D*>(blocks[12].get());
  TMatrix3D* B2 = reinterpret_cast<TMatrix3D*>(blocks[13].get());
  TMatrix3D* B3 = reinterpret_cast<TMatrix3D*>(blocks[14].get());
  

  std::pair<int, int> velo_pres_code(TDatabase::ParamDB->VELOCITY_SPACE,
                                             TDatabase::ParamDB->PRESSURE_SPACE);
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
    {
      transposed_B_structures_.push_back(B1T->GetStructure().GetTransposed());
      multigridLevel = new TNSE_MGLevel1(level, A11, B1, B2, B3, transposed_B_structures_.back().get(),
                                         s.rhs_.get_entries(),s.solution_.get_entries(), 
                                         nAuxArray, alpha.data(), 
                                         velo_pres_code.first,velo_pres_code.second, 
                                         nullptr, nullptr);
    }
      break;
    case 2:
      multigridLevel = new TNSE_MGLevel2(level, A11, B1, B2, B3, 
                                         B1T, B2T, B3T,
                                         s.rhs_.get_entries(),s.solution_.get_entries(), 
                                         nAuxArray, alpha.data(), 
                                         velo_pres_code.first,velo_pres_code.second, 
                                         nullptr, nullptr);
      break;
    case 3:
    {
      transposed_B_structures_.push_back(B1T->GetStructure().GetTransposed());
      multigridLevel = new TNSE_MGLevel3(level, A11, A12, A13, A21, A22, A23,
                                         A31, A32, A33,
                                         B1, B2, B3, transposed_B_structures_.back().get(),
                                         s.rhs_.get_entries(),s.solution_.get_entries(), 
                                         nAuxArray, alpha.data(), 
                                         velo_pres_code.first,velo_pres_code.second, 
                                         nullptr, nullptr);
    }
      break;
    case 4:
       multigridLevel = new TNSE_MGLevel4(level, 
                                          A11, A12, A13, A21, A22, A23,
                                          A31, A32, A33, B1, B2, B3, 
                                          B1T, B2T, B3T,
                                          s.rhs_.get_entries(),s.solution_.get_entries(), 
                                          nAuxArray, alpha.data(), 
                                          velo_pres_code.first,velo_pres_code.second, 
                                          nullptr, nullptr);
      break;
    case 14:
      TSquareMatrix3D* C   = reinterpret_cast<TSquareMatrix3D*>(blocks[15].get());      
      multigridLevel = new TNSE_MGLevel14(level, 
                                          A11, A12, A13, A21, A22, A23,
                                          A31, A32, A33, C, B1, B2, B3, 
                                          B1T, B2T, B3T,
                                          s.rhs_.get_entries(),s.solution_.get_entries(), 
                                          nAuxArray, alpha.data(), 
                                          velo_pres_code.first,velo_pres_code.second, 
                                          nullptr, nullptr);
      break;
  }
  
  return multigridLevel;
}

void NSE3D::mg_solver()
{
  System_per_grid& s = this->systems_.front();
  std::shared_ptr<TItMethod> itMethod;
  std::shared_ptr<TItMethod> prec;
  int zero_start;  
  TSquareMatrix3D *sqMat[10];
  TSquareMatrix **sqmatrices = (TSquareMatrix **)sqMat;
  TMatrix3D *recMat[6];
  TMatrix **matrices = (TMatrix **)recMat;
  MatVecProc *MatVect;
  DefectProc *Defect;
  
  std::vector<double> itMethodSol;
  std::vector<double> itMethodRhs;

  int nDof = this->get_size();
  
  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_TERRIBLY_UNSAFE();
  
  TSquareMatrix3D *A11=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
  TSquareMatrix3D *A12=reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());
  TSquareMatrix3D *A13=reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());
  TSquareMatrix3D *A21=reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());
  TSquareMatrix3D *A22=reinterpret_cast<TSquareMatrix3D*>(blocks.at(5).get());
  TSquareMatrix3D *A23=reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());
  TSquareMatrix3D *A31=reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());
  TSquareMatrix3D *A32=reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());
  TSquareMatrix3D *A33=reinterpret_cast<TSquareMatrix3D*>(blocks.at(10).get());
  // pressure-pressure block
  TSquareMatrix3D *C = reinterpret_cast<TSquareMatrix3D*>(blocks.at(15).get());
  
  TMatrix3D* B1T=reinterpret_cast<TMatrix3D*>(blocks.at(3).get());
  TMatrix3D* B2T=reinterpret_cast<TMatrix3D*>(blocks.at(7).get());
  TMatrix3D* B3T=reinterpret_cast<TMatrix3D*>(blocks.at(11).get());
  TMatrix3D* B1 =reinterpret_cast<TMatrix3D*>(blocks.at(12).get());
  TMatrix3D* B2 =reinterpret_cast<TMatrix3D*>(blocks.at(13).get());
  TMatrix3D* B3 =reinterpret_cast<TMatrix3D*>(blocks.at(14).get());
  
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      sqMat[0] = A11;
      recMat[0] = B1;
      recMat[1] = B2;
      recMat[2] = B3;
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
      break;
    case 2:
      sqMat[0] = A11;
      
      recMat[0] = B1;
      recMat[1] = B2;
      recMat[2] = B3;
      recMat[3] = B1T;
      recMat[4] = B2T;
      recMat[5] = B3T;
      
      MatVect = MatVect_NSE2;
      Defect = Defect_NSE2;
      break;
    case 3:
      sqMat[0] = A11;
      sqMat[1] = A12;
      sqMat[2] = A13;
      sqMat[3] = A21;
      sqMat[4] = A22;
      sqMat[5] = A23;
      sqMat[6] = A31;
      sqMat[7] = A32;
      sqMat[8] = A33;
      
      recMat[0] = B1;
      recMat[1] = B2;
      recMat[2] = B3;
      MatVect = MatVect_NSE3;
      Defect = Defect_NSE3;
      break;
    case 4:
      sqMat[0] = A11;
      sqMat[1] = A12;
      sqMat[2] = A13;
      sqMat[3] = A21;
      sqMat[4] = A22;
      sqMat[5] = A23;
      sqMat[6] = A31;
      sqMat[7] = A32;
      sqMat[8] = A33;
      
      recMat[0] = B1;
      recMat[1] = B2;
      recMat[2] = B3;
      recMat[3] = B1T;
      recMat[4] = B2T;
      recMat[5] = B3T;

      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
      break;
    case 14:
            sqMat[0] = A11;
      sqMat[1] = A12;
      sqMat[2] = A13;
      sqMat[3] = A21;
      sqMat[4] = A22;
      sqMat[5] = A23;
      sqMat[6] = A31;
      sqMat[7] = A32;
      sqMat[8] = A33;
      sqMat[9] = C;
      
      recMat[0] = B1;
      recMat[1] = B2;
      recMat[2] = B3;
      recMat[3] = B1T;
      recMat[4] = B2T;
      recMat[5] = B3T;

      //FIXME: MatVect_EquOrd_NSE4 is not implemented for 3D case
      // I have to implement that as well
      //MatVect = MatVect_EquOrd_NSE4;
      //Defect = Defect_EquOrd_NSE4;
      break;
  }
  
  
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
        prec=std::make_shared<TMultiGridIte>(MatVect, Defect, nullptr, 0, nDof, 
                                 this->multigrid_.get(), zero_start);
        break;
      default:
        ErrThrow("Unknown preconditioner !!!");
    }
    
    if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
    {
      itMethodSol.resize(nDof);
      itMethodRhs.resize(nDof);
      //FIXME find some other way than memcpy???
      memcpy(itMethodSol.data(), s.solution_.get_entries(), nDof*sizeof(double));
      memcpy(itMethodRhs.data(), s.rhs_.get_entries(), nDof*sizeof(double));
    }
    else
    {
      itMethodSol=s.solution_.get_entries_vector();
      itMethodRhs=s.rhs_.get_entries_vector();
    }
    
    switch(TDatabase::ParamDB->SC_SOLVER_SADDLE)
    {
      case 11:
        itMethod=std::make_shared<TFixedPointIte>(MatVect, Defect, prec.get(), 0, nDof, 0);
        break;
      case 16:
        itMethod=std::make_shared<TFgmresIte>(MatVect, Defect, prec.get(), 0, nDof, 0);
        break;
      default:
        ErrThrow("Unknown preconditioner !!!");
    }
    
    itMethod->Iterate(sqmatrices, matrices, itMethodSol.data(), itMethodRhs.data());
    
    if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE==5)
    {
      s.rhs_ = itMethodRhs.data();
      s.solution_ = itMethodSol.data();
    }
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
