#include <Time_CD3D.h>
#include <Database.h>
#include <Assemble3D.h>
#include <LinAlg.h>
#include <Multigrid.h>
#include <AlgebraicFluxCorrection.h>

#include <MainUtilities.h>
#include <TimeDiscretizations.h>

#include <cstring>
#include <sys/stat.h>

#ifdef _MPI
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

//==============================================================================
ParameterDatabase get_default_TCD3D_parameters()
{
  Output::print<5>("creating a default Time_CD3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default Time_CD3D database as well.
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.set_name("Time_CD3d parameter database");
  
  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  parmoon_db.merge(out_db, true);
  
  // a default afc database
  ParameterDatabase afc_db = AlgebraicFluxCorrection::default_afc_database();
  parmoon_db.merge(afc_db, true);
  
  // a default local assembling database
  parmoon_db.merge(LocalAssembling3D::default_local_assembling_database());
  
  parmoon_db.merge(TimeDiscretization::default_TimeDiscretization_database(), true);

  return parmoon_db;  
}
//==============================================================================
#ifdef _MPI
Time_CD3D::SystemPerGrid::SystemPerGrid(const Example_TimeCD3D& example, TCollection& coll, 
                                             int maxSubDomainPerDof)
: feSpace_(&coll, "space", "TCD3D feSpace", example.get_bc(0),
           TDatabase::ParamDB->ANSATZ_ORDER)
{
  //inform the fe space about the maximum number of subdomains per dof
  feSpace_.initialize_parallel(maxSubDomainPerDof);
  
  stiffMatrix_ = BlockFEMatrix::CD3D(feSpace_);
  massMatrix_ = BlockFEMatrix::CD3D(feSpace_);

  rhs_ = BlockVector(stiffMatrix_, true);
  solution_ = BlockVector (stiffMatrix_, false);

  feFunction_ = TFEFunction3D(&feSpace_, "u", "u",
                              solution_.get_entries(),solution_.length());
  solution_m1 = BlockVector(stiffMatrix_, false);
  u_m1 = TFEFunction3D(&feSpace_, "u_m1", "u_m1",
                              solution_m1.get_entries(),solution_m1.length());
  solution_m2 = BlockVector(stiffMatrix_, false);
  u_m2 = TFEFunction3D(&feSpace_, "u_m2", "u_m2",
                              solution_m2.get_entries(),solution_m2.length());
}
#else /* ***********************************************************************/
Time_CD3D::SystemPerGrid::SystemPerGrid(const Example_TimeCD3D& example, 
                                        TCollection& coll)
: feSpace_(&coll, "space", "TCD3D feSpace", example.get_bc(0),
           TDatabase::ParamDB->ANSATZ_ORDER)
{
  stiffMatrix_ = BlockFEMatrix::CD3D(feSpace_);
  massMatrix_ = BlockFEMatrix::CD3D(feSpace_);

  rhs_ = BlockVector(stiffMatrix_, true);
  solution_ = BlockVector (stiffMatrix_, false);

  feFunction_ = TFEFunction3D(&feSpace_, "u", "u",
                              solution_.get_entries(),solution_.length());
  solution_m1 = BlockVector(stiffMatrix_, false);
  u_m1 = TFEFunction3D(&feSpace_, "u_m1", "u_m1",
                              solution_m1.get_entries(),solution_m1.length());
  solution_m2 = BlockVector(stiffMatrix_, false);
  u_m2 = TFEFunction3D(&feSpace_, "u_m2", "u_m2",
                              solution_m2.get_entries(),solution_m2.length());
}

#endif

//==============================================================================
Time_CD3D::Time_CD3D(std::list<TCollection* >collections,
                     const ParameterDatabase &param_db,
                     const Example_TimeCD3D& _example
#ifdef _MPI
                     , int maxSubDomainPerDof
#endif
                     )
:  db(get_default_TCD3D_parameters()),
  solver(param_db), systems_(), example_(_example), errors_({}), 
  outputWriter(param_db), time_stepping_scheme(param_db)
{
  this->db.merge(param_db,false); // update this database with given values
  this->check_and_set_parameters();
  
  bool usingMultigrid = this->solver.is_using_multigrid();
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
    // initial concentration
    this->systems_.front().feFunction_.Interpolate(example_.get_initial_cond(0));
  }
  else
  {
    auto multigrid = this->solver.get_multigrid();
    size_t nMgLevels = multigrid->get_n_geometric_levels();
    if(collections.size() != nMgLevels)
    {
      ErrThrow("Multigrid: expected ", nMgLevels, " collections ", 
	       collections.size(), "provided.");
    }
        std::list<BlockFEMatrix*> matrices;
    // construct all SystemPerGrid and store them
    for(auto it : collections)
    {
#ifdef _MPI
      systems_.emplace_back(example_, *it, maxSubDomainPerDof);
#else
      systems_.emplace_back(example_, *it);
#endif
      systems_.front().feFunction_.Interpolate(example_.get_initial_cond(0));
      matrices.push_front(&systems_.back().stiffMatrix_);
    }
    multigrid->initialize(matrices);
  }// multigrid case
  // print useful information
  this->output_problem_size_info();
  this->rhs_from_time_disc.copy_structure(this->systems_.front().rhs_);
}

//==============================================================================
void Time_CD3D::output_problem_size_info() const
{
  // print some useful information
  const TFESpace3D& space = this->systems_.front().feSpace_;
  double hMin, hMax;
  TCollection *coll = space.GetCollection();
  coll->GetHminHmax(&hMin, &hMax);
#ifndef _MPI
  // sequential
  Output::print<1>("N_Cells    : ", setw(13), coll->GetN_Cells());
  Output::print<1>("h(min, max): ", setw(13), hMin, " ", setw(13), hMax);
  Output::print<1>("dofs all   : ", setw(13), space.GetN_DegreesOfFreedom());
  Output::print<1>("dof active : ", setw(13), space.GetActiveBound());
#else
  // MPI
  int my_rank, size;
  auto comm = TDatabase::ParamDB->Comm;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &size);
  
  double global_h_min = 0, global_h_max = 0;
  MPI_Reduce(&hMin, &global_h_min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&hMax, &global_h_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  
  std::vector<int> n_own_cells(size, 0); // for each process
  std::vector<int> n_halo_cells(size, 0); // for each process
  int local_own_cells = coll->GetN_OwnCells();
  int local_halo_cells = coll->GetN_HaloCells();
  MPI_Gather(&local_own_cells, 1, MPI_INT, // send
             &n_own_cells[0], 1, MPI_INT,  // receive
             0, comm);                     // control
  MPI_Gather(&local_halo_cells, 1, MPI_INT, // send
             &n_halo_cells[0], 1, MPI_INT,  // receive
             0, comm);                      // control
  
  auto par_comm = space.get_communicator();
  int n_local_master_dof = par_comm.GetN_Master();
  std::vector<int> n_masters(size, 0); //for each process
  MPI_Gather(&n_local_master_dof, 1, MPI_INT, // send
             &n_masters[0], 1, MPI_INT,       // receive
             0, comm);                        // control
  if(my_rank == 0)
  {
    Output::stat("Time_CD3D", "information on the fe space");
    size_t sum_cells_total = 0;
    size_t sum_dof_total = 0;
    //FIXME CB: This information output is incorrect!
    for(int i =0; i < size ;++i)
    {
      Output::dash("Process ", i, "\t n_own_cells ", n_own_cells[i], 
                   "\t n_halo_cells ", n_halo_cells[i], 
                   "\t n_master_dof ", n_masters[i]);
      sum_cells_total += n_own_cells[i];
      sum_dof_total += n_masters[i];
    }
    Output::dash("Total number of cells:              ", sum_cells_total);
    Output::dash("Total number of degrees of freedom: ", sum_dof_total);
    Output::dash("Global hMin/hMax: ", hMin, "/", hMax);
  }
#endif
}

//==============================================================================
void Time_CD3D::check_and_set_parameters()
{
  //set problem_type to Time_CD if not yet set
  if(!db["problem_type"].is(2))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 2;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to Time_CD."
          "It is now reset to the correct value for Time_CD (=2).");
      db["problem_type"] = 2;
    }
  }

  // an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }

  //////////////// Algebraic flux correction ////////////
  if(!db["algebraic_flux_correction"].is("none"))
  {//some kind of afc enabled

    if(solver.is_using_multigrid())
      ErrThrow("Multgrid and FEM-FCT are not yet enabled in Time_CD3D"); //TODO

    if(!db["algebraic_flux_correction"].is("fem-fct-cn"))
    {
      db["algebraic_flux_correction"].set("fem-fct-cn");
      Output::print("Parameter 'algebraic_flux_correction' changed to 'fem-fct-cn', ",
                    "which is the only implemented kind of algebraic "
                    "flux correction for TCD problems.");
    }

    if(TDatabase::TimeDB->TIME_DISC !=2)
    {
      ErrThrow("Algebraic flux correction with FEM-FCT is "
          "implemented for Crank-Nicolson time stepping scheme only. "
          "Set TDatabase::TimeDB->TIME_DISC to 2.")
    }

    if(TDatabase::ParamDB->ANSATZ_ORDER != 1)
    {
      ErrThrow("Algebraic flux correction with FEM-FCT does only work for"
          "linear elements. Change ANSATZ_ORDER to 1!");
    }

    //make sure that galerkin discretization is used
    if (!db["space_discretization_type"].is("galerkin"))
    {//some other disctype than galerkin
      db["space_discretization_type"] = "galerkin";
      Output::print("Parameter 'space_discretization_type' changed to 'galerkin' "
          "because Algebraic Flux Correction is enabled.");
    }

    // when using afc, create system matrices as if all dofs were active
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1; //FIXME THIS IS DANGEROUS, does not go well with BlockFEMatrix!
  }

  // Read space discretization type from database and
  // set this->disctype accordingly.
  if (db["space_discretization_type"].is("galerkin"))
    this->disctype = GALERKIN;      // = 1
  else if (db["space_discretization_type"].is("supg"))
    this->disctype = SUPG;          // = 2
  else
    ErrThrow("Unknown space discretization type for TCD3D!");

}

//==============================================================================
void Time_CD3D::assemble_initial_time()
{
  LocalAssembling_type allMatrices = LocalAssembling_type::TCDStiffMassRhs;
  for(auto &s : this->systems_)
  {
    TFEFunction3D *feFunction = {&s.feFunction_};
    LocalAssembling3D la(this->db, allMatrices, &feFunction,
                         example_.get_coeffs(), this->disctype);
    // Assemble stiffness, mass matrices and the rhs. Initially it is independent
    // that which method is used. 
    // 
    call_assembling_routine(s, la, true);
  }
  SystemPerGrid& s = this->systems_.front();
  old_rhs = s.rhs_;
  s.solution_m1 = s.solution_;
  s.solution_m2 = s.solution_m1;  
}

//==============================================================================
void Time_CD3D::assemble()
{
  // In the case of SUPG: local assemble function itself take care of the 
  // number of matrices. One have to assemble also the weighted mass matrix 
  // which comes from the time discretization of the SUPG method. We have 
  // to assemble the Mass matrix also for each time step due to the convection 
  // field which might also depend on time as well.
  LocalAssembling_type la_type = LocalAssembling_type::TCDStiffRhs;
  bool assemble_both = false;
  for(auto &s : this->systems_)
  {
    // call assembling routine 
    if(db["space_discretization_type"].is("supg"))
    {
      la_type = LocalAssembling_type::TCDStiffMassRhs;
    }
    TFEFunction3D *feFunction = {&s.feFunction_};
    LocalAssembling3D la(this->db, la_type, &feFunction,
                         example_.get_coeffs(), this->disctype);
    call_assembling_routine(s, la, assemble_both);
  }
  
  // here the modifications due to time discretization begin
  if (!db["algebraic_flux_correction"].is("none") )
  {
    do_algebraic_flux_correction();
    rhs_from_time_disc = this->systems_.front().rhs_;
    return; // modifications due to time discretization are per-
            // formed inside the afc scheme, so step out here!
  }

  // preparing the right hand side discretized by the used time
  // stepping scheme
  SystemPerGrid& s = this->systems_.front();
  rhs_from_time_disc.reset();
  rhs_from_time_disc = s.rhs_;
  // all matrices are available
  unsigned int n_sols = time_stepping_scheme.n_old_solutions();
  std::vector<BlockVector> old_sols(n_sols);
  old_sols[0] = s.solution_m1;
  if(old_sols.size() == 2)
    old_sols[1] = s.solution_m2;
  std::vector<BlockVector> rhs(2);
  rhs[0] = rhs_from_time_disc;
  rhs[1] = old_rhs;
  // prepare the right hand side from the previous time step
  time_stepping_scheme.prepare_rhs_from_time_disc(s.stiffMatrix_, 
                                                  s.massMatrix_, rhs, old_sols);
  rhs_from_time_disc = rhs[0];
  old_rhs=s.rhs_;
  rhs_from_time_disc.copy_nonactive(s.rhs_);
  
  for(auto &s : this->systems_)
    time_stepping_scheme.prepare_system_matrix(s.stiffMatrix_, s.massMatrix_);
  systems_[0].solution_.copy_nonactive(systems_[0].rhs_);  
}

//==============================================================================
void Time_CD3D::solve()
{
  SystemPerGrid& s = systems_.front();
#ifndef _MPI
  solver.solve(s.stiffMatrix_,rhs_from_time_disc,s.solution_); // sequential
#else // parallel
  if(solver.get_db()["solver_type"].is("direct")) // direct solvers
  {
    MumpsWrapper mumps_wrapper(s.stiffMatrix_);
    mumps_wrapper.solve(s.rhs_, s.solution_);
  }
  else
    solver.solve(s.stiffMatrix_,s.rhs_,s.solution_); // same as sequential
#endif
  // restore stiffness matrix
    if(db["algebraic_flux_correction"].is("none"))
  {
    for(auto &s : this->systems_)
      time_stepping_scheme.reset_linear_matrices(s.stiffMatrix_, s.massMatrix_);
  }
  else
  {
    //in AFC case the stiffness matrix is "ruined" by now -
    Output::print("AFC does not yet reset the stiffness"
        "matrix and old_Au correctly!");
  }
  
  s.solution_m2 = s.solution_m1;
  s.solution_m1 = s.solution_;
}

//==============================================================================
void Time_CD3D::output(int m, int& image)
{
  bool i_am_root = true;
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == 0);
#endif
  // nothing to do if no error computation or no vtk filess
  bool noOutput = !db["output_write_vtk"] && !db["output_compute_errors"];
  if(noOutput)
    return;
  SystemPerGrid &s = this->systems_.front();
  s.feFunction_.PrintMinMax();
  
#ifdef _MPI
  // computing errors as well as writing vtk files requires a minimum 
  // consistency level of 1
  s.feSpace_.get_communicator().consistency_update(s.solution_.get_entries(),1);
#endif // _MPI
  
  //write solution for visualization 
  outputWriter.add_fe_function(&s.feFunction_);
  outputWriter.write(image);
  
  
  // compute errors 
  if(db["output_compute_errors"])
  {
    MultiIndex3D allDerivatives[4] = { D000, D100, D010, D001 };
    TAuxParam3D aux(1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
    std::array<double, 5> locError = {};
    const TFESpace3D* space = s.feFunction_.GetFESpace3D();
    
    s.feFunction_.GetErrors(example_.get_exact(0), 4, allDerivatives, 2, 
                            L2H1Errors, example_.get_coeffs(), &aux, 1, &space,
                            locError.data());
#ifdef _MPI
    /// @todo the GetErrors method in TFEFunction3D should already to the 
    /// communication, it's surprising that in the mpi case the errors are 
    /// squared, while the square root has been taken already in the sequential 
    /// case.
    // global (across all processes) error, L2, H1, Linf, SD-error (for SUPG)
    std::vector<double> errorsReduced(4);
    MPI_Reduce(locError.data(), errorsReduced.data(), 2, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&locError[2], &errorsReduced[2], 1, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
    // correct values only on the root process!
    locError[0] = sqrt(errorsReduced[0]);
    locError[1] = sqrt(errorsReduced[1]);
    locError[2] = errorsReduced[2];
#endif // _MPI
    
    if(i_am_root)
    {
      Output::print<1>("time   : ", time_stepping_scheme.current_time_);
      Output::print<1>("  L2   : ", std::setprecision(14), locError[0]);
      Output::print<1>("  H1   : ", std::setprecision(14), locError[1]);
      Output::print<1>("  L_inf: ", std::setprecision(14), locError.at(2));
    }
    double tau = time_stepping_scheme.get_step_length();
    
    errors_[0] += (locError[0] * locError[0] + errors_[1] * errors_[1])*tau*0.5;
    errors_[1] = locError[0];
    
    errors_[2] += (locError[1] * locError[1] + errors_[3] * errors_[3])*tau*0.5;
    errors_[3] = locError[1];
    
    if(m == 0 || errors_.at(4) < locError[2])
      errors_[4] = locError[2];
    
    if(i_am_root)
    {
      Output::print<1>("  L2(0,T;L2)      : ", sqrt(errors_[0]));
      Output::print<1>("  L2(0,T;H1)      : ", sqrt(errors_[2]));
      Output::print<1>("  L_inf(0,T,L_inf): " , errors_[4]);
    }
  }
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
  
  BoundCondFunct3D * boundCond = feSpace->get_boundary_condition();
  
  BoundValueFunct3D * boundValue[1]{example_.get_bd()[0]};  
  
  Assemble3D(nFESpaces, &feSpace, nsqMatrices, sqMatrices, 
	       nRectMatrices, recMatrices, nRhs, &rhsEntries, &feSpaceRhs, 
	       &boundCond, boundValue, la);
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

// Note: This is almost entirely copy-and-paste from the 2D-case.
// That is a strong indication, that we actually need a better type
// hierarchy in the sense of object-orientation.
void Time_CD3D::do_algebraic_flux_correction()
{
  //determine which kind of afc to use
  if(db["algebraic_flux_correction"].is("default") ||
      db["algebraic_flux_correction"].is("fem-fct-cn"))
  {
    //TODO implement for multigrid!
    SystemPerGrid& s = systems_.front();

    //get references to the relevant objects

    const TFESpace3D& feSpace = s.feSpace_; //space
    const FEMatrix& mass = *s.massMatrix_.get_blocks().at(0).get(); // mass
    FEMatrix& stiff = *s.stiffMatrix_.get_blocks_uniquely().at(0).get(); //stiffness
    //vector entry arrays
#ifdef _MPI
    const TParFECommunicator3D& comm = s.feSpace_.get_communicator();
    // solEntries is needed in consistency 2 for afc.
    // this is done here on the solution blockvector, because
    // solEntries is given as a const and can't be updated
    comm.consistency_update(s.solution_.block(0),2);
#endif
    const std::vector<double>& solEntries = s.solution_.get_entries_vector();
    std::vector<double>& rhsEntries = s.rhs_.get_entries_vector();
    std::vector<double>& oldRhsEntries = old_rhs.get_entries_vector();

    // fill a vector "neumannToDirichlet" with those rows that got
    // internally treated as Neumann although they are Dirichlet
    int firstDiriDof = feSpace.GetActiveBound();
    int nDiri = feSpace.GetN_Dirichlet();
    std::vector<int> neumToDiri(nDiri, 0);
    std::iota(std::begin(neumToDiri), std::end(neumToDiri), firstDiriDof);

    //determine prelimiter from Database
    AlgebraicFluxCorrection::Prelimiter prelim;
    switch((int) db["afc_prelimiter"])
    {
      case 1:
        prelim = AlgebraicFluxCorrection::Prelimiter::MIN_MOD;
        Output::print<2>("FEM-FCT: MinMod prelimiter chosen.");
        break;
      case 2:
        prelim = AlgebraicFluxCorrection::Prelimiter::GRAD_DIRECTION;
        Output::print<2>("FEM-FCT: Gradient direcetion prelimiter chosen.");
        break;
      case 3:
        prelim = AlgebraicFluxCorrection::Prelimiter::BOTH;
        Output::print<2>("FEM-FCT: Double prelimiting (MinMod & Gradient) chosen.");
        break;
      default:
        prelim = AlgebraicFluxCorrection::Prelimiter::NONE;
        Output::print<2>("FEM-FCT: No prelimiting chosen.");
    }

    //get current timestep length from Database
    double delta_t = time_stepping_scheme.get_step_length();

    // apply C-N FEM-FCT
    AlgebraicFluxCorrection::crank_nicolson_fct(
        mass, stiff,
        solEntries,
        rhsEntries, oldRhsEntries,
        delta_t,
        neumToDiri,
        prelim );

    //...and finally correct the entries in the Dirichlet rows
    AlgebraicFluxCorrection::correct_dirichlet_rows(stiff);
    //...and in the right hand side, too
    s.rhs_.copy_nonactive(old_rhs);

#ifdef _MPI
    // this update is necessary here, because copy_nonactive()
    // can disturb the consistency of a vector,
    // see Time_NSE3D::assemble_rhs() for more details
    comm.consistency_update(s.rhs_.block(0),2);
#endif
  }
  else
  {
    ErrThrow("The chosen algebraic flux correction scheme ",
             db["algebraic_flux_correction"]," is unknown to class Time_CD3D.");
  }
}
