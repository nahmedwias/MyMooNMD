#include <Time_CD3D.h>
#include <Database.h>
#include <LocalAssembling3D.h>
#include <Assemble3D.h>
#include <ItMethod.h>
#include <JacobiIte.h>
#include <LinAlg.h>
#include <MultiGridScaIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <Multigrid.h>

#include <MainUtilities.h>

#include <cstring>
#include <sys/stat.h>
//==============================================================================
ParameterDatabase get_default_TCD3D_parameters()
{
  Output::print<3>("creating a default Time_CD3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default Time_CD3D database as well.
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.set_name("Time_CD3d parameter database");
  
  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  parmoon_db.merge(out_db, true);
  
  return parmoon_db;  
}
//==============================================================================
#ifdef _MPI
Time_CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example, TCollection& coll, 
                                             int maxSubDomainPerDof)
: feSpace_(&coll, (char*)"space", (char*)"TCD3D feSpace", example.get_bc(0), 
           TDatabase::ParamDB->ANSATZ_ORDER),
  stiffMatrix_({&feSpace_}), // stiffness matrix (block matrix)
  massMatrix_({&feSpace_}), // mass matrix (block matrix)
  rhs_(stiffMatrix_, true), // rhs hand side vector (filled with zeros)
  solution_(stiffMatrix_, false), // solution vector (filled with zeros)
  old_Au(this->stiffMatrix_, true),
  feFunction_(&feSpace_, (char*)"u", (char*)"u", solution_.get_entries(),solution_.length()),
  parMapper_(),
  parComm_()
{
  //inform the fe space about the maximum number of subdomains per dof
  feSpace_.SetMaxSubDomainPerDof(maxSubDomainPerDof);
  
  stiffMatrix_ = BlockFEMatrix::CD3D(feSpace_);
  massMatrix_ = BlockFEMatrix::CD3D(feSpace_);
  
  parMapper_ = TParFEMapper3D(1, &feSpace_);
  parComm_ = TParFECommunicator3D(&parMapper_);
}
#else /* ***********************************************************************/
Time_CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example, TCollection& coll)
: feSpace_(&coll, (char*)"space", (char*)"TCD3D feSpace", example.get_bc(0), 
           TDatabase::ParamDB->ANSATZ_ORDER),
  stiffMatrix_({&feSpace_}),
  massMatrix_({&feSpace_}),
  rhs_(stiffMatrix_, true),
  solution_(stiffMatrix_, false),
  old_Au(this->stiffMatrix_, true),
  feFunction_(&feSpace_, (char*)"u", (char*)"u", solution_.get_entries(),solution_.length())
{
  stiffMatrix_ = BlockFEMatrix::CD3D(feSpace_);
  massMatrix_ = BlockFEMatrix::CD3D(feSpace_);
}

#endif

//==============================================================================
Time_CD3D::Time_CD3D(std::list<TCollection* >collections, 
			      const ParameterDatabase &param_db,
			      const Example_CD3D& _example
#ifdef _MPI
			      , int maxSubDomainPerDof
#endif
	   )
: systems_(), example_(), multigrid_(nullptr), db(get_default_TCD3D_parameters()),
  solver(param_db), errors_(5,0.0)
{
  this->db.merge(param_db,false); // update this database with given values
  this->checkParameters();
  
  bool usingMultigrid = this->solver.get_db()["solver_type"].is("iterative")
                        && this->solver.get_db()["preconditioner"].is("multigrid");
  if(!usingMultigrid)
  {
    // check at least if the collections list contains exactly one collection
    if(collections.size() != 1)
    {
      ErrThrow("Non-multigrid: expected exactly one collection");
    }
    // the given collection for particular cell
    TCollection& cellCollection = *collections.front();
#ifdef _MPI
    // create finite element spaces, function, matrices, and rhs and solution vectors
    systems_.emplace_back(example_, cellCollection, maxSubDomainPerDof);
#else
    systems_.emplace_back(example_, cellCollection);
#endif    
    // print some useful information
    const TFESpace3D& space = this->systems_.front().feSpace_;
    double hMin, hMax;
    cellCollection.GetHminHmax(&hMin, &hMax);
    Output::print<1>("N_Cells    : ", setw(13), cellCollection.GetN_Cells());
    Output::print<1>("h(min, max): ", setw(13), hMin, " ", setw(13), hMax);
    Output::print<1>("dofs all   : ", setw(13), space.GetN_DegreesOfFreedom());
    Output::print<1>("dof active : ", setw(13), space.GetActiveBound());
    
    // initial concentration
    this->systems_.front().feFunction_.Interpolate(example_.get_initial_cond(0));
  }
  else
  {
    size_t nMgLevels = this->solver.get_db()["n_multigrid_levels"];
    if(collections.size() != nMgLevels)
    {
      ErrThrow("Multigrid: expected ", nMgLevels, " collections ", 
	       collections.size(), "provided.");
    }
    // creating multigrid object
    std::vector<double> param(2);
    param[0] = this->solver.get_db()["damping_factor"];
    param[1] = this->solver.get_db()["damping_factor_finest_grid"];
    multigrid_.reset(new TMultiGrid3D(1,2,param.data()));
    
    // construct all SystemPerGrid and store them
    for(auto it : collections)
    {
#ifdef _MPI
      systems_.emplace_back(example_, *it, maxSubDomainPerDof);
#else
      systems_.emplace_back(example_, *it);
#endif
      systems_.front().feFunction_.Interpolate(example_.get_initial_cond(0));
    }
    // TODO Temporary, one have to change the implementation when
    // the multigrid class is finished (CLEMENS)
    // Coarsest level must be the first, therefore the iteration order is reversed
    size_t level=0;
    for(auto system=systems_.rbegin(); system != systems_.rend(); ++system)
    {
      size_t nAuxArrays = 2;
      if( (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
	|| (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR) )
      {
	nAuxArrays = 4;
      }
      
      TSquareMatrix3D* block 
      = reinterpret_cast<TSquareMatrix3D*>(system->stiffMatrix_.get_blocks_uniquely().at(0).get());
#ifdef _MPI
      TMGLevel3D* multigridLevel = new TMGLevel3D(
        level, block, system->rhs_.get_entries(), 
        system->solution_.get_entries(), 
        &system->parComm_, &system->parMapper_,
	nAuxArrays, NULL);
#else
      TMGLevel3D* multigridLevel = new TMGLevel3D(
	level, block, system->rhs_.get_entries(),
	system->solution_.get_entries(), 
	nAuxArrays, NULL);
#endif
      multigrid_->AddLevel(multigridLevel);
      level++;
    }// end preparing multigrid levels
  }// multigrid case
}

//==============================================================================
void Time_CD3D::SystemPerGrid::descale_stiff_matrix(double tau, double theta1)
{
  if(tau == 0 || theta1 ==0)
  {
    ErrThrow("Do not divide by zero in descaling the stiffMatrix_! ", tau, 
	     " theta1 ", theta1);
  }
  // mass matrix (or matrices M ) are added to stiffMatrix_, subtract first
  const FEMatrix& massBlock = *massMatrix_.get_blocks().at(0).get();
  stiffMatrix_.add_matrix_actives(massBlock, -1, {{0,0}}, {false});
  // dscale stiffMatrix_
  const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
  stiffMatrix_.scale_blocks_actives(1./(tau*theta1), cell_positions);
}

//==============================================================================
void Time_CD3D::SystemPerGrid::update_old_Au()
{
  stiffMatrix_.apply(solution_, old_Au);
}

//==============================================================================
void Time_CD3D::checkParameters()
{
  if(!this->db["problem_type"].is(1))
  {
    this->db["problem_type"] = 1;// set correct problem type
    Output::print("PROBLEM_TYPE set to 1 for (convection-diffusion-reaction" 
                  "for this class Time_CD3D");
  }
  // an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }
#ifdef _MPI
  if(this->solver.get_db()["solver_type"].is("iterative")
     && TDatabase::ParamDB->SC_SOLVER_SCALAR !=11)
  {
    ErrThrow("Only SC_SOLVER_SCALAR: 11 (fixed point iteration) is implemented so far.")
  }
#endif
  // the only preconditioners implemented are Jacobi and multigrid
  if(this->solver.get_db()["solver_type"].is("iterative")
    && !this->solver.get_db()["preconditioner"].is("jacobi")
    && !this->solver.get_db()["preconditioner"].is("multigrid"))
  {
    ErrThrow("Only SC_PRECONDITIONER_SCALAR: 1 (Jacobi) and 5 (multigrid)"
        " are implemented so far.");
  }
}

//==============================================================================
void Time_CD3D::assemble_initial_time()
{
  LocalAssembling3D_type allMatrices = LocalAssembling3D_type::TCD3D;
  for(auto &s : this->systems_)
  {
    TFEFunction3D *feFunction = {&s.feFunction_};
    LocalAssembling3D la(allMatrices, &feFunction, example_.get_coeffs());
    // Assemble stiffness, mass matrices and the rhs. Initially it is independent
    // that which method is used. 
    // 
    call_assembling_routine(s, la, true);
    
    // initialize old_Au
    s.stiffMatrix_.apply(s.solution_, s.old_Au);
  }
  SystemPerGrid& s = this->systems_.front();
  old_rhs = s.rhs_;
}

//==============================================================================
void Time_CD3D::assemble()
{
  // In the case of SUPG: local assemble function itself take care of the 
  // number of matrices. One have to assemble also the weighted mass matrix 
  // which comes from the time discretization of the SUPG method. We have 
  // to assemble the Mass matrix also for each time step due to the convection 
  // field which might also depend on time as well.
  LocalAssembling3D_type stiffMatrixRhs = LocalAssembling3D_type::TCD3DStiffRhs;
  for(auto &s : this->systems_)
  {
    TFEFunction3D *feFunction = {&s.feFunction_};
    LocalAssembling3D la(stiffMatrixRhs, &feFunction, example_.get_coeffs());    
    // call assembling routine 
    switch(TDatabase::ParamDB->DISCTYPE)
    {
      case GALERKIN:
	call_assembling_routine(s, la, false);
	break;
      case SUPG:
	call_assembling_routine(s, la, true);
	break;
    }
  }
  
  if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION == 2)
  {
    ErrThrow("Take care of this later");
  }
  // preparing the right hand side 
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;
  
  SystemPerGrid & s = this->systems_.front();
  
  if(TDatabase::TimeDB->THETA4)
  {
    // scale the right hand side by theta4 and tau
    s.rhs_.scaleActive(tau*theta4);
    // add rhs from previous time step to the rhs scaled with 
    // time step length
    if(theta3 != 0)
      s.rhs_.addScaledActive((this->old_rhs), tau*theta3);
    // save old rhs for the next time steps
    if(theta3)
    {
      this->old_rhs.addScaledActive(s.rhs_, -1./(tau*theta3));
      this->old_rhs.scaleActive(-theta3/theta4);
    }
  }
  else
  {
    if(TDatabase::TimeDB->TIME_DISC == 0)
    {
      ErrThrow("Forward Euler method is not supported. "
          "Choose TDatabase::TimeDB->TIME_DISC as 1 (bw Euler)"
          " or 2 (Crank-Nicoloson)");
    }    
  }
  // For the SUPG method: Mass Matrix = (u, v + delta * bgrad v)
  s.massMatrix_.apply_scaled_add_actives(s.solution_, s.rhs_, 1.0);
  // rhs -= tau*theta2*A_old*uold
  s.rhs_.addScaledActive(s.old_Au, -tau*theta2);
  
  // copy nonactive to the solution vector
  s.solution_.copy_nonactive(s.rhs_);
  
  // preparing the left hand side; stiffness matrix is scaled by
  // tau * theta1 and adding the mass matrix to it
  // descaling is necessary only if the coefficients depends on time
  for(auto &s : this->systems_)
  {
    const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
    s.stiffMatrix_.scale_blocks_actives(tau*theta1, cell_positions);
    
    const FEMatrix& mass_block = *s.massMatrix_.get_blocks().at(0).get();
    s.stiffMatrix_.add_matrix_actives(mass_block, 1.0, {{0,0}}, {false});
  }
}

//==============================================================================
void Time_CD3D::solve()
{
  SystemPerGrid &s = this->systems_.front();
  TSquareMatrix * blocks[1] = {reinterpret_cast<TSquareMatrix*>(
    s.stiffMatrix_.get_blocks_TERRIBLY_UNSAFE().at(0).get()) };
    
#ifdef _MPI
  TParFECommunicator3D* parComm = &s.parComm_;
#endif
  if(this->solver.get_db()["solver_type"].is("iterative"))
  {
    // iterative solver
    int nDof = s.feSpace_.GetN_DegreesOfFreedom();
    TItMethod* preconditioner;
    TItMethod* iterativeSolver;
    
    int zeroStart;
    switch(TDatabase::ParamDB->SC_SOLVER_SCALAR)
    {
      case 11:
	zeroStart = 1;
	break;
      case 16:
	zeroStart = 0;
	break;
    }
    
    switch(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
    {
      case 1: // Jacobi iteration
#ifdef _MPI
	preconditioner = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL, 0, nDof, 
					1, parComm);
#else
	preconditioner = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL, 0, nDof, 1);
#endif
	break;	
      case 5:
	preconditioner = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL, 0, 
					      nDof, multigrid_.get(), zeroStart);
	break;
      default:
	ErrMsg("unknown SC_PRECONDITIONER_SCALAR: " 
	<< TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR);
    }
    
    // preparing solver 
    switch(TDatabase::ParamDB->SC_SOLVER_SCALAR)
    {
      case 11:	
#ifdef _MPI
	iterativeSolver = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, preconditioner,
					     0, nDof, 1, parComm);
#else
	iterativeSolver = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, preconditioner,
					     0, nDof, 1);
#endif
	break; // case 11
      case 16:
#ifdef _MPI
	iterativeSolver = new TFgmresIte(MatVect_Scalar, Defect_Scalar, preconditioner, 
					 0, nDof, 1, parComm);
#else
	iterativeSolver = new TFgmresIte(MatVect_Scalar, Defect_Scalar, preconditioner,
					 0, nDof, 1);
#endif
	break;
      default:
	ErrMsg("unknown solver !!!" << endl);
    }
    
    // solve
    if(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
    {
      std::vector<double> iterativeSol; // = s.solution_.get_entries();
      std::vector<double> iterativeRhs;
      iterativeSol.resize(nDof);
      iterativeRhs.resize(nDof);
      memcpy(iterativeSol.data(), s.solution_.get_entries(), nDof*SizeOfDouble);
      memcpy(iterativeRhs.data(), s.rhs_.get_entries(), nDof*SizeOfDouble);
      
      iterativeSolver->Iterate(blocks, nullptr, iterativeSol.data(), iterativeRhs.data());
      
      memcpy(s.solution_.get_entries(), iterativeSol.data(), nDof*SizeOfDouble);
      memcpy(s.rhs_.get_entries(), iterativeRhs.data(), nDof*SizeOfDouble);
    }
    else
    {
      iterativeSolver->Iterate(blocks, nullptr, s.solution_.get_entries(), s.rhs_.get_entries());
    }
  }
  else if(this->solver.get_db()["solver_type"].is("direct"))
  {
#ifndef _MPI
    this->solver.update_matrix(s.stiffMatrix_);
    this->solver.solve(s.rhs_, s.solution_);
    return;
#elif _MPI
    std::vector<const TParFECommunicator3D*> par_comms_init = {&s.parComm_};
    std::vector<TParFECommunicator3D*> par_comms_solv = {&s.parComm_};

    MumpsWrapper mumps_wrapper(s.stiffMatrix_, par_comms_init);
    mumps_wrapper.solve(s.rhs_, s.solution_, par_comms_solv);
#endif
  }
  else
  {
    ErrThrow("unknown SOLVER_TYPE: Choose either 1 (iterative) or 2 (direct) " );
  }
}

//==============================================================================
void Time_CD3D::output(int m, int& imgage)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  // nothing to do if no error computation or no vtk filess
  // bool noOutput = !db["output_write_vtk"] && !db["output_compute_errors"];
  // if(noOutput)
  //  return;
  SystemPerGrid &s = this->systems_.front();
  s.feFunction_.PrintMinMax();
  
  //write solution for visualization
   if(db["output_write_vtk"])
   {
    //TODO
   }
  
  // compute errors 
//  if(db["output_compute_errors"])
//  {
    MultiIndex3D allDerivatives[4] = { D000, D100, D010, D001 };
    TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    std::vector<double> locError; 
    locError.resize(5);
    const TFESpace3D* space = s.feFunction_.GetFESpace3D();
    
    s.feFunction_.GetErrors(example_.get_exact(0), 4, allDerivatives, 2, L2H1Errors, 
			    example_.get_coeffs(), &aux, 1, &space, locError.data());
#ifdef _MPI
    std::vector<double> errorsReduced;
    errorsReduced.resize(4); // maximum 4; In case of SUPG one need also the SD error
    MPI_Allreduce(locError.data(), errorsReduced.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    errors_.at(0) = locError.at(0);
    errors_.at(1) = locError.at(1);
#else
    int my_rank = 0;
#endif
    
    Output::print<1>("time: ", TDatabase::TimeDB->CURRENTTIME);
    Output::print<1>("  L2: ", locError.at(0));
    Output::print<1>("  H1: ", locError.at(1));
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    
    errors_.at(0) += (locError.at(0) * locError.at(0) 
                    + errors_.at(1) * errors_.at(1))*tau*0.5;
    errors_.at(1) = locError.at(0);    
    
    errors_.at(2) += (locError.at(1) * locError.at(0)
                    + errors_.at(3) * errors_.at(3))*tau*0.5;
    errors_.at(3) = locError.at(1);
    
    if(m==0)
      errors_.at(4) = locError.at(0);
    if(errors_.at(4) < locError.at(0))
      errors_.at(4) = locError.at(0);
    
    if(my_rank == 0)
    {
      Output::print<1>("  L2(0,T;L2)   : ", sqrt(errors_.at(0)));
      Output::print<1>("  L2(0,T;H1)   : ", sqrt(errors_.at(2)));
      Output::print<1>("  Linft(0,T,L2): " , errors_.at(4));
    }
  //}
}

//==============================================================================
void Time_CD3D::call_assembling_routine(Time_CD3D::SystemPerGrid& system, 
					LocalAssembling3D& la, bool assemble_both)
{
  int nFESpaces = 1;
  const TFESpace3D * feSpace = &system.feSpace_;
  
  std::vector<std::shared_ptr<FEMatrix>> stiffBlock = 
                    system.stiffMatrix_.get_blocks_uniquely();
  std::vector<std::shared_ptr<FEMatrix>> massBlock = 
                    system.massMatrix_.get_blocks_uniquely();
  
  int nsqMatrices = 3; // maximum number of matrices
  TSquareMatrix3D *sqMatrices[nsqMatrices];
  if(assemble_both)
  {
    nsqMatrices = 2;    
    sqMatrices[0] = reinterpret_cast<TSquareMatrix3D*>(stiffBlock.at(0).get());
    sqMatrices[0]->reset();
    sqMatrices[1] = reinterpret_cast<TSquareMatrix3D*>(massBlock.at(0).get());
    sqMatrices[1]->reset();
  }
  else
  {
    nsqMatrices = 1; // Mass and Stiffness Matrices at once
    sqMatrices[0] = reinterpret_cast<TSquareMatrix3D*>(stiffBlock.at(0).get());
    sqMatrices[0]->reset();
  }
  
  int nRectMatrices = 0;
  TMatrix3D ** recMatrices{nullptr};
  
  int nRhs = 1;
  double *rhsEntries = system.rhs_.get_entries();
  system.rhs_.reset(); // reset to zeros 
  const TFESpace3D *feSpaceRhs = &system.feSpace_;
  
  BoundCondFunct3D * boundCond = feSpace->getBoundCondition();
  
  BoundValueFunct3D * boundValue[1]{example_.get_bd()[0]};  
  
  Assemble3D(nFESpaces, &feSpace, nsqMatrices, sqMatrices, 
	       nRectMatrices, recMatrices, nRhs, &rhsEntries, &feSpaceRhs, 
	       &boundCond, boundValue, la);
}

//==============================================================================
void Time_CD3D::descale_stiffness()
{
  double theta1 = TDatabase::TimeDB->THETA1;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  // restore the stiffMatrix_ and store the old_Au
  if(!TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION)
  {
    for(auto &s : this->systems_)
    {
      s.descale_stiff_matrix(tau, theta1);
      s.update_old_Au();
    }
  }
  else
  {
    ErrThrow("AFC scheme is not implemented yet" );
  }
}

//==============================================================================
std::array< double, int(3) > Time_CD3D::get_errors() const
{
  std::array<double, int(3)> error_at_time_points;
  error_at_time_points.at(0) = sqrt(errors_.at(0));
  error_at_time_points.at(1) = sqrt(errors_.at(2));
  error_at_time_points.at(2) = sqrt(errors_.at(4));
  return error_at_time_points;
}

//==============================================================================
