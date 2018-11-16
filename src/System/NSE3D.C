#include <NSE3D.h>
#include <Database.h>
#include "LocalAssembling.h"
#include <Assemble3D.h>
#include <MainUtilities.h>
#include <LinAlg.h>
#include <DirectSolver.h>

#include <GridTransfer.h>
#include <Multigrid.h>
#include <Upwind3D.h>
#include <BoundaryAssembling3D.h>
#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!

#include <sys/stat.h>

#ifdef _MPI
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

ParameterDatabase NSE3D::default_NSE_database()
{
  Output::print<5>("creating a default NSE3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("NSE3D parameter database");
  
  //NSE3D requires a nonlinear iteration, set up a nonlinit_database and merge
  db.merge(ParameterDatabase::default_nonlinit_database(), true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);
  
  db.merge(LocalAssembling3D::default_local_assembling_database(), true);
  db.merge(Example3D::default_example_database(), true);
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

NSE3D::System_per_grid::System_per_grid(const Example_NSE3D& example,
                                    TCollection& coll, std::pair<int, int> order, 
                                    NSE3D::Matrix type)
  :  velocitySpace_(new TFESpace3D(&coll, "u", "nse3d velocity", 
                                   example.get_bc(0), //bd cond at 0 is x velo bc
                                   order.first)),
     pressureSpace_(new TFESpace3D(&coll, "p", "nse3d pressure", 
                                   example.get_bc(3), //bd condition at 3 is pressure bc
                                   order.second))
{
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      matrix_ = BlockFEMatrix::NSE3D_Type1(*velocitySpace_, *pressureSpace_);
      break;                                               
    case 2:                                                
      matrix_ = BlockFEMatrix::NSE3D_Type2(*velocitySpace_, *pressureSpace_);
      break;                                               
    case 3:                                                
      matrix_ = BlockFEMatrix::NSE3D_Type3(*velocitySpace_, *pressureSpace_);
      break;                                               
    case 4:                                                
      matrix_ = BlockFEMatrix::NSE3D_Type4(*velocitySpace_, *pressureSpace_);
      break;                                               
    case 14:                                               
      matrix_ = BlockFEMatrix::NSE3D_Type14(*velocitySpace_, *pressureSpace_);
      break;
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }

  rhs_ = BlockVector(matrix_, true);
  solution_ = BlockVector(matrix_, false);

  u_ = TFEVectFunct3D(velocitySpace_.get(), "u", "u", solution_.block(0),
     solution_.length(0), 3);
  p_ = TFEFunction3D(pressureSpace_.get(), "p", "p", solution_.block(3),
     solution_.length(3));

#ifdef _MPI
  //print some information
  velocitySpace_->get_communicator().print_info();
  pressureSpace_->get_communicator().print_info();
#endif
}

void NSE3D::output_problem_size_info() const
{
  int my_rank = 0;
#ifndef _MPI
    auto & velocity_space = *this->systems_.front().velocitySpace_;
    auto & pressure_space = *this->systems_.front().pressureSpace_;

    size_t nDofu  = velocity_space.GetN_DegreesOfFreedom();
    size_t nDofp  = pressure_space.GetN_DegreesOfFreedom();
    size_t nTotal = 3*nDofu + nDofp;
    size_t nActive= 3*velocity_space.GetActiveBound();

    TCollection* coll = velocity_space.GetCollection();

    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);

#else
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    auto velocity_comm = systems_.front().velocitySpace_->get_communicator();
    auto pressure_comm = systems_.front().pressureSpace_->get_communicator();
    int nDofu  = velocity_comm.get_n_global_dof();
    int nDofp  = pressure_comm.get_n_global_dof();
    int nTotal = 3*nDofu + nDofp;

#endif
    if(my_rank ==0)
    {
    Output::stat("NSE3D", "Mesh data and problem size on finest grid");
#ifndef _MPI
    Output::dash("N_Cells      :  ", setw(10), coll->GetN_Cells());
    Output::dash("h(min, max)  :  ", setw(10), hmin, setw(10), " ", hmax);
#endif
    Output::dash("ndof velocity:  ", setw(10), 3*nDofu );
    Output::dash("ndof pressure:  ", setw(10), nDofp);
    Output::dash("ndof total   :  ", setw(10), nTotal );
#ifndef _MPI
    Output::dash("nActive      :  ", setw(10), nActive);
#endif
    }
}

NSE3D::NSE3D(const TDomain& domain, const ParameterDatabase& param_db)
 : NSE3D(domain, param_db, Example_NSE3D(param_db))
{
}

NSE3D::NSE3D(const TDomain& domain, const ParameterDatabase& param_db,
             const Example_NSE3D& example)
  : systems_(), example_(example), db(default_NSE_database()), outputWriter(param_db),
    solver(param_db), defect_(), old_residuals_(), initial_residual_(1e10), 
    errors_()
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
  auto collections = domain.get_grid_collections();
  TCollection *coll = collections.front(); //the finest grid collection
  // create finite element space and function, a matrix, rhs, and solution
  systems_.emplace_back(example_, *coll, velocity_pressure_orders, type);
  
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
    for(int i = n_geo_multigrid_levels; i < n_grids; ++i)
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
    matrices.push_back(&systems_.back().matrix_);

    for(auto coll : collections) // initialize the coarse grid space hierarchy
    {
      systems_.emplace_back(example, *coll, velocity_pressure_orders, type);
      // prepare input argument for multigrid object
      matrices.push_front(&systems_.back().matrix_);
    }
    // initialize the multigrid object with all the matrices on all levels
    mg->initialize(matrices);
  }

  output_problem_size_info();
}

void NSE3D::check_parameters()
{
  if(!db["problem_type"].is(3) && !db["problem_type"].is(5))
  {
    Output::warn<2>("The parameter problem_type doesn't correspond neither to NSE "
        "nor to Stokes. It is now reset to the default value for NSE (=5).");
    db["problem_type"] = 5;
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
}

void NSE3D::get_velocity_pressure_orders(std::pair< int, int >& velocity_pressure_orders)
{
  int velocity_order = velocity_pressure_orders.first;
  int pressure_order = velocity_pressure_orders.second;
  int order = 0;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5: // P_k/Q_k
      order = velocity_order;
      break;
    case 12:
      ErrThrow("Scott-Vogelius P2/P1disc is not stable in 3D even if the final "
               "refinement step is barycentric.");
      break;
    case 13: case 14: case 15:
      // P_k/P_{k-1}^{disc}, k>=3, elements are in general not inf-sup stable on 
      // tetrahedra. If the last refinement step was barycentric refinement, 
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
          // Scott-Vogelius: discontinuous pressure spaces with standard 
          // conforming velocity space. This is not stable on general triangles,
          // be sure to use barycentric refinement.
        case 12: case 13: case 14: case 15:
          pressure_order = -velocity_order+1;
          break;
        case 22: case 23: case 24:
          pressure_order = -(velocity_order-11)*10;
          break;
      }
      break;
    // continuous pressure spaces
    case 1: case 2: case 3: case 4: case 5:
      pressure_order = pressure_order*1;
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

void NSE3D::assemble_linear_terms()
{
  size_t nFESpace = 2; // spaces used for assembling matrices
  size_t nSqMatrices=10; // maximum no of square matrices (type 14)
  std::vector<TSquareMatrix3D*> sqMatrices(nSqMatrices, nullptr);
  size_t nReMatrices = 6; // maximum no of rectangular matrices (e.g. type 14)
  std::vector<TMatrix3D*> reMatrices(nReMatrices, nullptr);
  size_t nRhs=4; // maximum number of right hand sides (e.g. type 4)
  std::vector<double*> rhsArray(nRhs); // right hand side 
  // finite element function used for nonlinear term
  std::vector<TFEFunction3D*> feFunction(4, nullptr);
  
  for(auto &s : this->systems_)
  {

    const TFESpace3D *v_space = s.velocitySpace_.get();
    const TFESpace3D *p_space = s.pressureSpace_.get();

    // spaces for matrices
    const TFESpace3D *spaces[2] = {v_space, p_space};
    const TFESpace3D *rhsSpaces[4] = {v_space, v_space, v_space, p_space};

    // spaces for right hand side    
    s.rhs_.reset();
    rhsArray[0]=s.rhs_.block(0);
    rhsArray[1]=s.rhs_.block(1);
    rhsArray[2]=s.rhs_.block(2);
    rhsArray[3]=s.rhs_.block(3);
    
    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_uniquely();
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[1].get());
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[2].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
        break;
      case 2:
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[4].get());
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[5].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[6].get());
        reMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks[1].get());
        reMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks[2].get());
        reMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
        break;
      case 3:
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
        sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
        sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
        sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
        sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
        sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
        sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
        sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[11].get());
        break;
      case 4:
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
        sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
        sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
        sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
        sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
        sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
        sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[9].get());
        sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[10].get());
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[12].get()); //first the lying B blocks
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[13].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[14].get());
        reMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks[3].get()); //than the standing B blocks
        reMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
        reMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks[11].get());
        break;
      case 14:        
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
        reMatrices[0]=reinterpret_cast<TMatrix3D*>(blocks[12].get());
        reMatrices[1]=reinterpret_cast<TMatrix3D*>(blocks[13].get());
        reMatrices[2]=reinterpret_cast<TMatrix3D*>(blocks[14].get());
        reMatrices[3]=reinterpret_cast<TMatrix3D*>(blocks[3].get());
        reMatrices[4]=reinterpret_cast<TMatrix3D*>(blocks[7].get());
        reMatrices[5]=reinterpret_cast<TMatrix3D*>(blocks[11].get());
        break;
    }// endswitch nstype

    for(unsigned int i=0; i<nSqMatrices; i++)
    {
      if(sqMatrices[i] != nullptr)
        sqMatrices[i]->reset();
    }
    for(unsigned int i=0; i<nReMatrices;i++)
    {
      if(reMatrices[i] != nullptr)
        reMatrices[i]->reset();
    }

    // boundary conditions and boundary values
    BoundCondFunct3D * boundContion[4]={
      spaces[0]->get_boundary_condition(), spaces[0]->get_boundary_condition(),
      spaces[0]->get_boundary_condition(), spaces[1]->get_boundary_condition()};
      
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
    LocalAssembling3D la(this->db, LocalAssembling_type::NSE3D_Linear, 
                         feFunction.data(), example_.get_coeffs());
    
    //HOTFIX: Check the documentation - this ensures that in Galerkin disc,
    // the convective term is not assembled at this point.
    assemble_nse = Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION;

    // assemble now the matrices and right hand side 
    Assemble3D(nFESpace, spaces, 
               nSqMatrices, sqMatrices.data(),
               nReMatrices, reMatrices.data(), 
               nRhs, rhsArray.data(), rhsSpaces,
               boundContion, boundValues.data(), la);
    
    Output::print(" ** START ASSEMBLE PRESSURE BC ON RHS **");

    // get all cells: this is at the moment needed for the boundary assembling
    /// @todo get only the (relevant) boundary cells
    /// e.g., bdCells = coll->get_cells_on_component(i)
    TCollection* coll = v_space->GetCollection();
    std::vector<TBaseCell*> allCells;
    for (int i=0 ; i < coll->GetN_Cells(); i++)
    {
      allCells.push_back(coll->GetCell(i));
    }
	
    BoundaryAssembling3D bi;
    for (int k=0;k<TDatabase::ParamDB->n_neumann_boundary;k++)
    {
      
      double t=TDatabase::TimeDB->CURRENTTIME;
      double PI = acos(-1.0);
      double pressure_of_t = TDatabase::ParamDB->neumann_boundary_value[k]*sin(2*PI*t);

      Output::print(" ** set value ", pressure_of_t,
		    " on boundary ",TDatabase::ParamDB->neumann_boundary_id[k]);
      
      bi.rhs_g_v_n(s.rhs_,v_space,nullptr,
		   allCells,
		   TDatabase::ParamDB->neumann_boundary_id[k],
		   pressure_of_t);
    }
    //delete the temorary feFunctions gained by GetComponent
    for(int i = 0; i<3; ++i)
      delete feFunction[i];
  }// endfor auto grid

  //copy non-actives from rhs to solution on finest grid
  this->systems_.front().solution_.copy_nonactive(systems_.front().rhs_);

/** When we call copy_nonactive in MPI-case, we have to remember the following:
   * it can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active.
   * That's why we have to update the values so that the vector becomes consistent again.
   * This is done here.
   */
#ifdef _MPI
  double *u1 = this->systems_.front().solution_.block(0);
  double *u2 = this->systems_.front().solution_.block(1);
  double *u3 = this->systems_.front().solution_.block(2);
  double *p  = this->systems_.front().solution_.block(3);
  this->systems_.front().velocitySpace_->get_communicator().consistency_update(u1, 3);
  this->systems_.front().velocitySpace_->get_communicator().consistency_update(u2, 3);
  this->systems_.front().velocitySpace_->get_communicator().consistency_update(u3, 3);
  this->systems_.front().pressureSpace_->get_communicator().consistency_update(p, 3);
#endif
}

void NSE3D::assemble_non_linear_term()
{
  size_t nFESpace = 2; // space needed for assembling matrices
  size_t nSqMatrices=10; // no of square matrices (Maximum 3)
  std::vector<TSquareMatrix3D*> sqMatrices(nSqMatrices, nullptr);
  size_t nReMatrices = 6; // no rectangular matrix in nonlinear iteration
  std::vector<TMatrix3D*> reMatrices(nReMatrices, nullptr);
  size_t nRhs=4; // no right hand side to be assembled
  std::vector<double*> rhsArray(nRhs, nullptr);  
  std::vector<TFEFunction3D*> feFunction(3);
  
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
        spaces.push_back(s.velocitySpace_.get());
        u_entries.push_back(s.solution_.block(block));
        u_ns_dofs.push_back(s.solution_.length(block));
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

  for(auto &s : this->systems_)
  {
    //hold the velocity space, we'll need it...
    const TFESpace3D * v_space = s.velocitySpace_.get();
    const TFESpace3D * p_space = s.pressureSpace_.get();
    const TFESpace3D* rhsSpaces[4] = {v_space, v_space, v_space, p_space};
#ifdef _MPI
    //MPI: solution in consistency level 3 (TODO: this might be superfluous here)
    for (size_t bl = 0; bl < s.solution_.n_blocks() ;++bl)
    {
      s.matrix_.get_communicators()[bl]->consistency_update(s.solution_.block(bl), 3);
    }
#endif

    // spaces for matrices
    const TFESpace3D* spaces[2] = {s.velocitySpace_.get(), s.pressureSpace_.get()};
    std::vector<std::vector<size_t>> cells={{0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2}};
    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_uniquely(cells);

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        sqMatrices.at(0)=reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get());
        break;
      case 3: 
      case 4:
      case 14:
        sqMatrices[0]=reinterpret_cast<TSquareMatrix3D*>(blocks[0].get());
        sqMatrices[1]=reinterpret_cast<TSquareMatrix3D*>(blocks[1].get());
        sqMatrices[2]=reinterpret_cast<TSquareMatrix3D*>(blocks[2].get());
        sqMatrices[3]=reinterpret_cast<TSquareMatrix3D*>(blocks[3].get());
        sqMatrices[4]=reinterpret_cast<TSquareMatrix3D*>(blocks[4].get());
        sqMatrices[5]=reinterpret_cast<TSquareMatrix3D*>(blocks[5].get());
        sqMatrices[6]=reinterpret_cast<TSquareMatrix3D*>(blocks[6].get());
        sqMatrices[7]=reinterpret_cast<TSquareMatrix3D*>(blocks[7].get());
        sqMatrices[8]=reinterpret_cast<TSquareMatrix3D*>(blocks[8].get());
        break;
    }// endswitch nstype

    // boundary conditions and boundary values
    BoundCondFunct3D * boundCondition[1]={
      spaces[0]->get_boundary_condition() };

    std::array<BoundValueFunct3D*, 3> boundValues;
    boundValues[0]=example_.get_bd()[0];
    boundValues[1]=example_.get_bd()[1];
    boundValues[2]=example_.get_bd()[2];

    feFunction[0]=s.u_.GetComponent(0);
    feFunction[1]=s.u_.GetComponent(1);
    feFunction[2]=s.u_.GetComponent(2);

    //decide wether to assemble by upwinding or not
    bool finest_grid = (&s == &systems_.at(0));
    bool do_upwinding = (db["space_discretization_type"].is("upwind")
                        || (mdml && !finest_grid))
                        && !is_stokes;
    
    if(!do_upwinding)
    {
      for(auto mat : sqMatrices)
      {
        if(mat != nullptr)
          mat->reset();
      }
      // local assembling object    
      LocalAssembling3D la(this->db, LocalAssembling_type::NSE3D_NonLinear, 
                           feFunction.data(), example_.get_coeffs());
      
      //HOTFIX: Check the documentation!
      assemble_nse = Hotfixglobal_AssembleNSE::WITH_CONVECTION;

      // assemble now the matrices and right hand side 
      Assemble3D(nFESpace, spaces, nSqMatrices, sqMatrices.data(),
                 nReMatrices, reMatrices.data(), nRhs, rhsArray.data(), 
                 rhsSpaces, boundCondition, boundValues.data(), la);
    }
    else
    {
      double one_over_nu = 1/example_.get_nu(); //the inverse of the example's diffusion coefficient
      for(auto mat : sqMatrices)
      {
        UpwindForNavierStokes3D(mat, feFunction[0], feFunction[1],
                                feFunction[2], one_over_nu);
      }
    }

    //delete the temporary feFunctions gained by GetComponent
    for(int i = 0; i<3; ++i)
      delete feFunction[i];

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
  
  // check if minimum number of iterations was performed already
  size_t min_it = db["nonlinloop_minit"];
  if(iteration_counter < min_it)
	  return false;

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
#ifdef _MPI
    //MPI: solution in consistency level 3 (TODO: maybe this is superfluous here
    // (because solution might be in level 3 consistency already)!)
    auto comms = s.matrix_.get_communicators();
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(s.solution_.block(bl), 3);
    }
#endif

  defect_ = s.rhs_;
  s.matrix_.apply_scaled_add(s.solution_, defect_,-1.);

  if(s.matrix_.pressure_projection_enabled())
  {
    TFEFunction3D defect_fctn(s.pressureSpace_.get(),
                              "p_def","pressure defect function",
                              &defect_[3*n_u_dof], n_p_dof);
    defect_fctn.project_into_L20();
  }

  // This is the calculation of the residual, given the defect.
  BlockVector defect_impuls({n_u_dof,n_u_dof,n_u_dof});
  BlockVector defect_mass{n_p_dof};
  //copy the entries (BlockVector offers no functionality to do this more nicely)
  for(size_t i = 0; i<3*n_u_dof ;++i)
    defect_impuls.get_entries()[i] = defect_.get_entries()[i];
  for(size_t i =0 ; i<n_p_dof ; ++i)
    defect_mass.get_entries()[i] = defect_.get_entries()[3*n_u_dof + i];

#ifdef _MPI
  double impuls_residual_square = defect_impuls.norm_global({comms[0],comms[1],comms[2]});
  impuls_residual_square *= impuls_residual_square;
  double mass_residual_square = defect_mass.norm_global({comms[3]});
  mass_residual_square *= mass_residual_square;
#else
  double impuls_residual_square = defect_impuls.norm();
  impuls_residual_square *= impuls_residual_square;
  double mass_residual_square = defect_mass.norm();
  mass_residual_square *= mass_residual_square;
#endif

  Residuals current_residuals(impuls_residual_square, mass_residual_square);
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
  
  // solving:
#ifndef _MPI
  this->solver.solve(s.matrix_, s.rhs_, s.solution_);
#endif
#ifdef _MPI
  if(this->solver.get_db()["solver_type"].is("direct"))
  {
    if(damping != 1.0)
      Output::warn("NSE3D::solve", "damping in an MPI context is not tested");

    //set up a MUMPS wrapper
    MumpsWrapper mumps_wrapper(s.matrix_);

    //kick off the solving process
    mumps_wrapper.solve(s.rhs_, s.solution_);
  }
  else
    this->solver.solve(s.matrix_, s.rhs_, s.solution_); // same as sequential
#endif
  
  if(damping != 1.0)
  {
    s.solution_.scale(damping);
    s.solution_.add_scaled(*old_solution, 1-damping);
  }

  // project pressure if necessary
  if(s.matrix_.pressure_projection_enabled())
    s.p_.project_into_L20();
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
  
  // In multigrid case, print time that was spent on coarse grid
  // (needed for curr project (ParMooN paper, Sep 2016), can be removed after)
  if(solver.is_using_multigrid())
    solver.get_multigrid()->print_coarse_grid_time_total();

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
  outputWriter.add_fe_function(&s.p_);
  outputWriter.add_fe_vector_function(&s.u_);
  outputWriter.write();
  
  
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    double err_u1[4]; // of these arrays only the two first entries are used,
    double err_u2[4]; // but the evil GetErrors() will corrupt memory if these
    double err_u3[4]; // have not at least size 4
    double err_p[4];


    TAuxParam3D aux(1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
    MultiIndex3D nsAllDerivs[4] = {D000, D100, D010, D001};
    const TFESpace3D *velocity_space = &this->get_velocity_space();
    const TFESpace3D *pressure_space = &this->get_pressure_space();
    
    // errors in first velocity component
    u1->GetErrors(example_.get_exact(0), 4, nsAllDerivs, 3,
                  L2H1Errors, nullptr, &aux, 1, &velocity_space, err_u1);
    // errors in second velocity component
    u2->GetErrors(example_.get_exact(1), 4, nsAllDerivs, 3,
                  L2H1Errors, nullptr, &aux, 1, &velocity_space, err_u2);
    // errors in third velocity component
    u3->GetErrors(example_.get_exact(2), 4, nsAllDerivs, 3,
                  L2H1Errors, nullptr, &aux, 1, &velocity_space, err_u3);
    double div_error = s.u_.GetL2NormDivergenceError(example_.get_exact(0),
                                                     example_.get_exact(1), 
                                                     example_.get_exact(2));
    // errors in pressure
    s.p_.GetErrors(example_.get_exact(3), 4, nsAllDerivs, 3, L2H1Errors,
                   nullptr, &aux, 1, &pressure_space, err_p);
    
#ifdef _MPI
    double err_red[9]; //memory for global (across all processes) error
    double err_send[9]; //fill send buffer
    err_send[0]=err_u1[0];
    err_send[1]=err_u1[1];
    err_send[2]=err_u2[0];
    err_send[3]=err_u2[1];
    err_send[4]=err_u3[0];
    err_send[5]=err_u3[1];
    err_send[6]=div_error*div_error;
    err_send[7]=err_p[0];
    err_send[8]=err_p[1];

    MPI_Allreduce(err_send, err_red, 9, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(i=0;i<9;i++)
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
    div_error = err_red[6];
    err_p[0] = err_red[7];
    err_p[1] = err_red[8];
#else
    int my_rank =0;
#endif

    errors_.at(0) = sqrt(err_u1[0]*err_u1[0] + err_u2[0]*err_u2[0] + err_u3[0]*err_u3[0]);//L2
    errors_.at(1) = div_error;
    errors_.at(2) = sqrt(err_u1[1]*err_u1[1] + err_u2[1]*err_u2[1] + err_u3[1]*err_u3[1]);//H1-semi
    errors_.at(3) = err_p[0];
    errors_.at(4) = err_p[1];

    //print errors
    if(my_rank == 0)
    {
      Output::stat("NSE3D", "Measured errors");
      Output::dash("L2(u)     : ", setprecision(14), errors_.at(0));
      Output::dash("L2(div(u)): ", setprecision(14), errors_.at(1));
      Output::dash("H1-semi(u): ", setprecision(14), errors_.at(2));
      Output::dash("L2(p)     : ", setprecision(14), errors_.at(3));
      Output::dash("H1-semi(p): ", setprecision(14), errors_.at(4));
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

std::array<double, int(5)> NSE3D::get_errors() const
{
  return errors_;
}
