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

NSE3D::SystemPerGrid::SystemPerGrid(const Example_NSE3D& example,
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
     //default construct parallel infrastructrue, will be reset in the body
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

  if(TDatabase::ParamDB->NSTYPE == 1 || TDatabase::ParamDB->NSTYPE == 3)
  {
    ErrThrow("Parallel solve needs explicitely stored BT-blocks, chose NSTYPE 2 or 4.");
  }

  velocitySpace_.SetMaxSubDomainPerDof(maxSubDomainPerDof);
  pressureSpace_.SetMaxSubDomainPerDof(maxSubDomainPerDof);

  //Must be reset here, because feSpace needs special treatment
  // This includes copy assignment - all because there is no good
  // way to communicate Maximum number of subdomains per dof to FESpace...
  parMapperVelocity_ = TParFEMapper3D(3, &velocitySpace_, //magic number
                                      matrix_.get_blocks_TERRIBLY_UNSAFE().at(0)->GetRowPtr(), // The parMapper for velocity expects row and column
                                      matrix_.get_blocks_TERRIBLY_UNSAFE().at(0)->GetKCol());  // non-const pointer from the A matrices
  parMapperPressure_ = TParFEMapper3D(1, &pressureSpace_,
                                      matrix_.get_blocks_TERRIBLY_UNSAFE().at(3)->GetRowPtr(), // The parMapper for velocity expects row and column
                                      matrix_.get_blocks_TERRIBLY_UNSAFE().at(3)->GetKCol());  // non-const pointer from the B-Transposed matrices - I don't know why exactly!
  // TODO The upper thing will lead to extreme trouble if there are no real BT blocks,
  // i.e. NSTYPE 1 or 3. Catch that case in an input-check method AND in the constructor.

  parCommVelocity_ = TParFECommunicator3D(&parMapperVelocity_);
  parCommPressure_ = TParFECommunicator3D(&parMapperPressure_);

#endif
}


NSE3D::NSE3D(std::list<TCollection* > collections, const Example_NSE3D& example
#ifdef _MPI
             , int maxSubDomainPerDof
#endif
) : systems_(), example_(example), multigrid_(nullptr),
    defect(), oldResidual(0), initial_residual(1e10), errors()
{
  std::pair <int,int> 
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and preesure spaces
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
    // Check at least if the collections list contains exactly one Collection.
    if(collections.size() != 1 )
    {
      ErrThrow("Non-multigrid: Expected exactly one collection!");
    }
    // Get the one given collection.
    TCollection& cellCollection = *collections.front();
        
    #ifdef _MPI
    // create finite element space and function, a matrix, rhs, and solution
    systems_.emplace_back(example_, cellCollection, velocity_pressure_orders, type,
                          maxSubDomainPerDof);
    #else
    // create finite element space and function, a matrix, rhs, and solution
    systems_.emplace_back(example_, cellCollection, velocity_pressure_orders, type);
    
    const TFESpace3D & velocity_space = this->systems_.front().velocitySpace_;
    const TFESpace3D & pressure_space = this->systems_.front().pressureSpace_;
    
    size_t nDofu  = velocity_space.GetN_DegreesOfFreedom();
    size_t nDofp  = pressure_space.GetN_DegreesOfFreedom();
    size_t nTotal = 3*nDofu + nDofp;
    size_t nActive= 3*velocity_space.GetActiveBound();
    
    Output::print<1>("N_Cells      :  ", setw(10), cellCollection.GetN_Cells());
    Output::print<1>("ndof Velocity:  ", setw(10), 3*nDofu );
    Output::print<1>("ndof Pressure:  ", setw(10), nDofp);
    Output::print<1>("ndof Total   :  ", setw(10), nTotal );
    Output::print<1>("nActive      :  ", setw(10), nActive);
    double hmin, hmax;
    cellCollection.GetHminHmax(&hmin, &hmax);
    Output::print<1>("h(min, max)  :  ",setw(10), hmin, setw(10), " ", hmax);    
    #endif
  }
  else // multigrid
  {
    size_t nLevels = TDatabase::ParamDB->LEVELS;
    // check
    if(collections.size() != nLevels)
    {
      ErrThrow("Multigrid: Expected ", nLevels, " collections, "
                 , collections.size(), " provided.");
    }
    std::vector<double> param(2);
    param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
    param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
    
    this->multigrid_.reset(new TNSE_MultiGrid(1, 2, param.data()));
    
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
    
    size_t level = 0;
    //Create multigrid-level-objects and add them to the multgrid object.
    // Must be coarsest level first, therefore reverse order iteration.
    for(auto system=systems_.rbegin(); system != systems_.rend(); ++system)
    {
      #ifdef _MPI
        ErrThrow("Clemens has to take care about the MPI implementation");
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

    const TFESpace3D * v_space = &s.velocitySpace_;
    const TFESpace3D * p_space = &s.pressureSpace_;

    // spaces for matrices
    const TFESpace3D* spaces[2] = {v_space, p_space};
    const TFESpace3D* rhsSpaces[4] = {v_space, v_space, v_space, p_space};
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
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_uniquely();
    
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        nSqMatrices = 1;
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());        
        break;
      case 3:
      case 4:
      case 14:
        nSqMatrices = 3;
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());        
        sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
        sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());        
        break;
    }// endswitch nstype
    // boundary conditions and boundary values
    BoundCondFunct3D * boundContion[4]={
      spaces[0]->getBoundCondition(), spaces[0]->getBoundCondition(),
      spaces[0]->getBoundCondition(), spaces[1]->getBoundCondition() };
      
    std::array<BoundValueFunct3D*, 4> boundValues;
    boundValues[0]=example_.get_bd()[0];
    boundValues[1]=example_.get_bd()[1];
    boundValues[2]=example_.get_bd()[2];
    boundValues[3]=example_.get_bd()[3];

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
               boundContion, boundValues.data(), la);
  }// endfor auto grid
}

bool NSE3D::stop_it(unsigned int iteration_counter)
{
  SystemPerGrid& s = this->systems_.front();
  unsigned int nDofU = s.solution_.length(0);
  unsigned int nDofP = s.solution_.length(3);

  defect = s.rhs_;
  s.matrix_.apply_scaled_add(s.solution_, defect, -1);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    IntoL20Vector3D(&defect[3*nDofU], nDofP,
                    TDatabase::ParamDB->PRESSURE_SPACE);
  }
  // compute residuals of mass and momentum
  double residual, impulseResidual, massResidual;
  residual = Ddot(3*nDofU+nDofP, &this->defect[0], &this->defect[0]);
  impulseResidual = Ddot(3*nDofU, &this->defect[0], &this->defect[0]);
  massResidual = Ddot(nDofP,&this->defect[3*nDofU],&this->defect[3*nDofU]);
  
  Output::print("nonlinear step  :  " , setw(3), iteration_counter);
  Output::print("impulse_residual:  " , setw(3), impulseResidual);
  Output::print("mass_residual   :  " , setw(3), massResidual);
  Output::print("residual        :  " , setw(3), sqrt(residual));
  
  if (iteration_counter>0)
  {
    Output::print("rate:           :  " , setw(3), sqrt(residual)/oldResidual);
  }
  oldResidual=sqrt(residual);
  int maxIte = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  double limit=TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  if(TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE)
  {
    limit *=sqrt(this->get_size());
    Output::print("stopping tolerance for nonlinear iteration ", limit);
  }
  if(iteration_counter==0)
    initial_residual = sqrt(residual);
  if ((((sqrt(residual)<=limit)||(iteration_counter==maxIte)))
   && (iteration_counter>=TDatabase::ParamDB->SC_MINIT))
   {
     Output::print("ITE : ", setw(3), iteration_counter, "  RES : ", sqrt(residual), 
                   " Reduction : ",  sqrt(residual)/initial_residual);
     return true;
   }
   else
     return false;
}

void NSE3D::solve()
{

  SystemPerGrid& s = this->systems_.front();

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
#ifdef _SEQ // So far only UMFPACK is available as direct solver in sequential case.
            // Actuate it via DirectSolver class.
      /// @todo consider storing an object of DirectSolver in this class
      DirectSolver direct_solver(s.matrix_,
                                 DirectSolver::DirectSolverTypes::umfpack);
      direct_solver.solve(s.rhs_, s.solution_);
#elif _MPI
    ErrThrow("This is the place to interface the mumps solver!");
#else
    ErrThrow("No OMPONLY solver (PARDISO) interfaced yet.");
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
  
  SystemPerGrid& s=this->systems_.front();
  
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
  {
    s.u_.GetComponent(0)->PrintMinMax();
    s.u_.GetComponent(1)->PrintMinMax();
    s.u_.GetComponent(2)->PrintMinMax();
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
    double err[7];
    TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    MultiIndex3D nsAllDerivs[4] = {D000, D100, D010, D001};
    const TFESpace3D *velocity_space = &this->get_velocity_space();
    const TFESpace3D *pressure_space = &this->get_pressure_space();
    
    // errors in first velocity component
    s.u_.GetComponent(0)->GetErrors(example_.get_exact(0), 4, nsAllDerivs, 2, 
                            L2H1Errors, nullptr, &aux, 1, &velocity_space, err);
    // errors in second velocity component
    s.u_.GetComponent(1)->GetErrors(example_.get_exact(1), 4, nsAllDerivs, 2, 
                            L2H1Errors, nullptr, &aux, 1, &velocity_space, err + 2);
    
    s.u_.GetComponent(2)->GetErrors(example_.get_exact(2), 4, nsAllDerivs, 2, 
                            L2H1Errors, nullptr, &aux, 1, &velocity_space, err + 4);
    
    errors[0] = sqrt(err[0]*err[0] + err[2]*err[2] + err[3]*err[3]);
    errors[1] = sqrt(err[4]*err[4] + err[5]*err[5] + err[6]*err[6]);    
    Output::print<1>("L2(u)     : ", errors[0]);
    Output::print<1>("H1-semi(u): ", errors[1]);
    // errors in pressure
    s.p_.GetErrors(example_.get_exact(3), 4, nsAllDerivs, 2, L2H1Errors, 
                  nullptr, &aux, 1, &pressure_space, err);
    
    errors[2] = err[0];
    errors[3] = err[1];    
    Output::print<1>("L2(p)     : ", errors[2]);
    Output::print<1>("H1-semi(p): ", errors[3]);    
    
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
}

TNSE_MGLevel* NSE3D::mg_levels(int level, NSE3D::SystemPerGrid& s)
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
  
  TStructure structure = B1T->GetStructure();
  std::pair<int, int> velo_pres_code(TDatabase::ParamDB->VELOCITY_SPACE,
                                             TDatabase::ParamDB->PRESSURE_SPACE);
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      multigridLevel = new TNSE_MGLevel1(level, A11, B1, B2, B3, &structure, 
                                         s.rhs_.get_entries(),s.solution_.get_entries(), 
                                         nAuxArray, alpha.data(), 
                                         velo_pres_code.first,velo_pres_code.second, 
                                         nullptr, nullptr);
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
      multigridLevel = new TNSE_MGLevel3(level, A11, A12, A13, A21, A22, A23,
                                         A31, A32, A33,
                                         B1, B2, B3, &structure, 
                                         s.rhs_.get_entries(),s.solution_.get_entries(), 
                                         nAuxArray, alpha.data(), 
                                         velo_pres_code.first,velo_pres_code.second, 
                                         nullptr, nullptr);
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
  SystemPerGrid& s = this->systems_.front();
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
      s.solution_.get_entries_vector() = itMethodSol;
      s.rhs_.get_entries_vector() = itMethodRhs;
    }
}

std::array< double, int(4) > NSE3D::get_errors()
{
  return errors;
}
