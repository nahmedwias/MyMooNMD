#include <CD3D.h>
#include <Example_CD3D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <LocalAssembling3D.h>
#include <Assemble3D.h>
//#include <PostProcessing3D.h>

#include <FixedPointIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>
#include <DirectSolver.h>

#include <MultiGrid3D.h>
#include <MainUtilities.h> // L2H1Errors

#include <Multigrid.h>

#include <sys/stat.h>

#ifdef _MPI
#include <MumpsWrapper.h>
#endif

  bool use_new_solver_debug = false;

ParameterDatabase get_default_CD3D_parameters()
{
  Output::print<3>("creating a default CD3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default CD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("CD3D parameter database");
  
  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  return db;
}

#ifdef _MPI
  CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example,
                                     TCollection& coll, int maxSubDomainPerDof)
   : feSpace_(&coll, (char*)"space", (char*)"cd3d fe_space", example.get_bc(0),
              TDatabase::ParamDB->ANSATZ_ORDER),
     matrix_({&feSpace_}), //system block matrix
     rhs_(matrix_, true), // suitable right hand side vector filled with zeroes
     solution_(matrix_, false), // suitable solution vector filled with zeroes
     feFunction_(&feSpace_, (char*)"c", (char*)"c", solution_.get_entries(),
                 solution_.length()),
     parMapper_(), // will be reset shortly
     parComm_() // will be reset shortly

  {
    //inform the fe space about the maximum number of subdomains per dof
    feSpace_.SetMaxSubDomainPerDof(maxSubDomainPerDof);

    // reset the matrix with named constructor
    matrix_ = BlockFEMatrix::CD3D(feSpace_);

    // Must be reset here, because feSpace needs special treatment
    // This includes copy assignment - all because there is no good
    // way to communicate Maximum number of subdomains per dof to FESpace...
    parMapper_ = TParFEMapper3D(1, &feSpace_);
    parComm_ = TParFECommunicator3D(&parMapper_);
  }
#else
  /* ************************************************************************ */
  CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example,
                                     TCollection& coll)
   : feSpace_(&coll, (char*)"space", (char*)"cd3d fe_space", example.get_bc(0),
              TDatabase::ParamDB->ANSATZ_ORDER),
     matrix_({&feSpace_}), //system block matrix
     rhs_(matrix_, true), // suitable right hand side vector filled with zeroes
     solution_(matrix_, false), // suitable solution vector filled with zeroes
     feFunction_(&feSpace_, (char*)"c", (char*)"c", solution_.get_entries(), solution_.length())
  {
    // reset the matrix with named constructor
    matrix_ = BlockFEMatrix::CD3D(feSpace_);
  }
#endif

  /** ************************************************************************ */
  CD3D::CD3D(std::list<TCollection*> collections,
             const ParameterDatabase& param_db, const Example_CD3D& example
#ifdef _MPI
             ,int maxSubDomainPerDof
#endif
  )
  : systems_(), example_(example), multigrid_(nullptr),
    db(get_default_CD3D_parameters()), solver(param_db), errors_()
  {
    this->db.merge(param_db, false); // update this database with given values
    this->checkParameters();
    // The construction of the members differ, depending on whether
    // a multigrid solver will be used or not.
    bool usingMultigrid = this->solver.is_using_multigrid();

    if (!usingMultigrid)
    {
      //Check at least if the collections list contains exactly one Collection.
      if(collections.size() != 1 )
      {
        ErrThrow("Non-multigrid: Expected exactly one collection!");
      }

      // Get the one given collection.
      TCollection& cellCollection = *collections.front();

#ifdef _MPI
      // create finite element space and function, a matrix, rhs, and solution
      systems_.emplace_back(example_, cellCollection, maxSubDomainPerDof);
#else
      // create finite element space and function, a matrix, rhs, and solution
      systems_.emplace_back(example_, cellCollection);

      // print out some information
      const TFESpace3D & space = this->systems_.front().feSpace_;
      double hMin, hMax;
      cellCollection.GetHminHmax(&hMin, &hMax);
      Output::print<1>("N_Cells    : ", setw(12), cellCollection.GetN_Cells());
      Output::print<1>("h (min,max): ", setw(12), hMin, " ", setw(12), hMax);
      Output::print<1>("dof all    : ", setw(12), space.GetN_DegreesOfFreedom());
      Output::print<1>("dof active : ", setw(12), space.GetN_ActiveDegrees());
#endif
    }
    else
    {// we are using multigrid
      if(use_new_solver_debug)
      {// TEST NEW MULTIGRID
        size_t n_levels = collections.size();
        if(!param_db["multigrid_n_levels"].is(n_levels))
           ErrThrow("Number of collection does not equal number of multigrid levels!");

        ParameterDatabase database_mg = Multigrid::default_multigrid_database();
        database_mg.merge(param_db, false);

        // Construct systems per grid and store them, finest level first
        std::list<BlockFEMatrix*> matrices;
        for (auto coll : collections)
        {
#ifdef _MPI
          systems_.emplace_back(example, *coll, maxSubDomainPerDof);
#else
          systems_.emplace_back(example, *coll);
#endif
          //prepare input argument for multigrid object
          matrices.push_front(&systems_.back().matrix_);
        }

        // Construct multigrid object
        mg_ = std::make_shared<Multigrid>(database_mg, matrices);
      }
      else
      { //THIS IS OLD MULTIGRID
        // number of refinement levels for the multigrid
        size_t nMgLevels = param_db["multigrid_n_levels"];

        // check if we have as many collections as expected multigrid levels
        if(collections.size() != nMgLevels )
        {
          ErrThrow("Multigrid: Expected ", nMgLevels, " collections, "
                   , collections.size(), " provided.");
        }

        // Create multigrid object.
        double *param = new double[2]; // memory leak
        param[0] = this->solver.get_db()["damping_factor"];
        param[1] = this->solver.get_db()["damping_factor_finest_grid"];
        multigrid_.reset(new TMultiGrid3D(1, 2, param));

        // Construct all System(s)PerGrid and store them.
        for(auto it : collections)
        {
#ifdef _MPI
          systems_.emplace_back(example_, *it, maxSubDomainPerDof);
#else
          systems_.emplace_back(example_, *it);
#endif
        }

        // Create multigrid-level-objects and add them to the multgrid object.
        // Must be coarsest level first, therefore reverse order iteration.
        size_t level = 0;
        for(auto system = systems_.rbegin(); system != systems_.rend(); ++system)
        {
          // determine number of auxiliary arrays (cryptic thing needed by Multigrid constructor)
          // (code taken from old ParMooN)
          size_t nAuxArrays = 2;
          if ( (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
              || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR) )
          { // seems this must be 4 when step length control happens
            nAuxArrays= 4;
          }

          //the single block
          TSquareMatrix3D* block = reinterpret_cast<TSquareMatrix3D*>(
              system->matrix_.get_blocks_uniquely().at(0).get());

#ifdef _MPI
          TMGLevel3D* multigridLevel = new TMGLevel3D
              (
                  level, block, system->rhs_.get_entries(),
                  system->solution_.get_entries(),
                  &system->parComm_, &system->parMapper_,
                  nAuxArrays, NULL
              );
#else
          TMGLevel3D* multigridLevel = new TMGLevel3D
              (
                  level, block, system->rhs_.get_entries(),
                  system->solution_.get_entries(),
                  nAuxArrays, NULL
              );
#endif
          multigrid_->AddLevel(multigridLevel);
          level++;
        } //end constructing and adding multigrid levels
      }
    }//end multigrid case
  }

/** ************************************************************************ */
void CD3D::assemble()
{
  //determine the local assembling type to be CD3D
  LocalAssembling3D_type laType = LocalAssembling3D_type::CD3D;
  // this loop has more than one iteration only in case of multigrid
  for(auto & s : systems_)
  {
    TFEFunction3D * feFunctionPtr = &s.feFunction_;

    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling3D laObject(laType, &feFunctionPtr, example_.get_coeffs());

    // assemble the system matrix with given local assembling, solution and rhs
    //s.matrix_.assemble(laObject, s.solution_, s.rhs_);
    call_assembling_routine(s, laObject);
  }

}

/** ************************************************************************ */
void CD3D::solve()
{
  if(use_new_solver_debug)
  {
    SystemPerGrid& s = this->systems_.front();

    //determine whether we make use of multigrid
    bool using_multigrid = this->solver.is_using_multigrid();

    if(!using_multigrid)
    {//no multigrid
      if(this->solver.get_db()["solver_type"].is("direct"))
      {
#ifndef _MPI
        this->solver.solve(s.matrix_, s.rhs_, s.solution_);
#endif
#ifdef _MPI
        //two vectors of communicators (const for init, non-const for solving)
        std::vector<const TParFECommunicator3D*> par_comms_init = {&s.parComm_};
        std::vector<TParFECommunicator3D*> par_comms_solv = {&s.parComm_};

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
    {//multigrid preconditioned iterative solver is used
      // All matrices which mg_'s levels point to must be ready!
      solver.solve(s.matrix_, s.rhs_, s.solution_, mg_);
    }
  }
  else
  {
    // Hold a reference to the system on the finest grid.
    SystemPerGrid& syst = systems_.front();

    // Hold some variable which will be needed in
    // all solver variants.

    //get the block for the solver
    TSquareMatrix* blocks[1] = {reinterpret_cast<TSquareMatrix*>(
        syst.matrix_.get_blocks_TERRIBLY_UNSAFE().at(0).get())};


    // FIXME we could use get_blocks_uniquely() here, but that is
    // intended for assemblers, not solvers - to mark that the basic
    // problem is still the non-const passing of matrices to the solvers,
    // we use the deprecated get_blocks_TERRIBLY_UNSAFE() here on purpose

    double* solutionEntries = syst.solution_.get_entries();
    double* rhsEntries =syst.rhs_.get_entries();
#ifdef _MPI
    TParFECommunicator3D* parComm = &syst.parComm_;
#endif


    if (this->solver.get_db()["solver_type"].is("iterative"))
    { // Iterative solver chosen.

      // Hold/declare some variables which will be needed for all
      // iterative solvers.
      int nDof = syst.feSpace_.GetN_DegreesOfFreedom();
      TItMethod* preconditioner;
      TItMethod* iterativeSolver;

      // Determine and build preconditioner.
      switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
      {
        case 1: //Jacobi Iteration
        {
#ifdef _MPI
          preconditioner = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL, 0, nDof, 1, parComm);
#else
          preconditioner = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL, 0, nDof, 1);
#endif
          break;
        }
        case 5: //TMultiGridScaIte (multgrid iteration)
          preconditioner = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar,
                                                NULL, 0, nDof, multigrid_.get() , 1);
          break;

        default:
          ErrMsg("Unknown SC_PRECONDITIONER_SCALAR: " << TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR);
      } //end building preconditioner

      // Determine and build solver.
      switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
      {
        case 11: //fixed point iteration
        {
#ifdef _MPI
          iterativeSolver = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, preconditioner, 0, nDof, 1, parComm);
#else
          iterativeSolver = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, preconditioner, 0, nDof, 1);
#endif
          break;
        }
        case 16:
        ErrThrow("SC_SOLVER_SCALAR: ", TDatabase::ParamDB->SC_SOLVER_SCALAR,
                 " (FGMRES) has to be implemented.");
          break;
        default:
        ErrThrow("Unknown solver !!!");
      }
      //call the solver
      if(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
      { //using multgrid
        // \todo The multgrid object and the iterative solver work on different copies of the solution and rhs. Is this necessary?

        //make copies on which the iterative solver may work
        double* iterativeSol = new double[nDof];
        double* iterativeRhs = new double[nDof];
        memcpy(iterativeSol, solutionEntries, nDof*SizeOfDouble);
        memcpy(iterativeRhs, rhsEntries, nDof*SizeOfDouble);

        //call the solver
        iterativeSolver->Iterate(blocks, nullptr, iterativeSol, iterativeRhs);

        //copy back & tidy up
        memcpy(solutionEntries, iterativeSol, nDof*SizeOfDouble);
        memcpy(rhsEntries, iterativeRhs, nDof*SizeOfDouble);
        delete[] iterativeSol;
        delete[] iterativeRhs;
      }
      else
      {
        iterativeSolver->Iterate(blocks, nullptr, solutionEntries, rhsEntries);
      }
      delete preconditioner;
      delete iterativeSolver;

    }
    else if (this->solver.get_db()["solver_type"].is("direct"))
    { // Direct solver chosen.
#ifndef _MPI
      this->solver.update_matrix(syst.matrix_);
      this->solver.solve(syst.rhs_, syst.solution_);
      return;
#elif _MPI
      std::vector<const TParFECommunicator3D*> par_comms_init = {&syst.parComm_};
      std::vector<TParFECommunicator3D*> par_comms_solv = {&syst.parComm_};

      MumpsWrapper mumps_wrapper(syst.matrix_, par_comms_init);
      mumps_wrapper.solve(syst.rhs_, syst.solution_, par_comms_solv);
#endif
    }
    else
    {
      ErrThrow("Unknown SOLVER_TYPE. Choose either 1 (iterative) or 2 (direct).");
    }
  }
}

void CD3D::output(int i)
{
#ifdef _MPI
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

	bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
	if(no_output)
		return;

  SystemPerGrid& syst = systems_.front() ;

  // print the value of the largest and smallest entry in the FE vector
  syst.feFunction_.PrintMinMax();

  /*
  // write output
  TOutput3D Output;
  //Output.init();
  Output.addFEFunction(&syst.feFunction_);
  ///@todo parallel output implementation
  #ifdef _MPI
  #else
  Output.write(i,0.0);
  #endif
  */
  
  
  // implementation with the old class TOutput2D
  // write solution to a vtk file
  if(db["output_write_vtk"])
  {
    // last argument in the following is domain, but is never used in this class
    TOutput3D Output(1, 1, 0, 0, NULL);
    Output.AddFEFunction(&syst.feFunction_);
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
    double errors[5];
    TAuxParam3D aux(1, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
    MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
    const TFESpace3D* space = syst.feFunction_.GetFESpace3D();

    syst.feFunction_.GetErrors(example_.get_exact(0), 4, AllDerivatives,
                               2, L2H1Errors, example_.get_coeffs(),
                               &aux, 1, &space, errors);
#ifdef _MPI
    double errorsReduced[4]; //memory for global (across all processes) error

    MPI_Allreduce(errors, errorsReduced, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(i=0;i<2;i++)
      errors[i] = sqrt(errorsReduced[i]);
#else
    int my_rank =0;
#endif

    //store errors
    errors_.at(0) = errors[0]; //L2
    errors_.at(1) = errors[1]; //H1-semi

    //print errors
    if(my_rank == 0)
    {
      Output::print("");
      Output::print( "L2: ", errors_.at(0));
      Output::print( "H1-semi: ", errors_.at(1));
    }
  } // if(this->db["compute_errors"]S)
}

void CD3D::checkParameters()
{
  //check if the correct problem type is set, change eventually
  if (!this->db["problem_type"].is(1))
  {
    this->db["problem_type"] = 1; //set correct problem type
    Output::print("PROBLEM_TYPE set to 1 (convection-diffusion-reaction), "
                  "for this is class CD3D.");
  }

  //an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }

  // the only solving strategy implemented is iterative
#ifdef _MPI
  // among the iterative solvers only 11 (fixed point iteration/
  // Richardson) is working
  if(this->solver.get_db()["solver_type"].is("iterative") 
      && TDatabase::ParamDB->SC_SOLVER_SCALAR != 11)
  {
      ErrThrow("Only SC_SOLVER_SCALAR: 11 (fixed point iteration) is implemented so far.")
  }
#endif // mpi

}

void CD3D::call_assembling_routine(SystemPerGrid& s, LocalAssembling3D& local_assem)
{//FIXME the body of this function was copy and paste

  const TFESpace3D * fe_space = &s.feSpace_;
  BoundCondFunct3D * boundary_conditions = fe_space->getBoundCondition();
  int N_Matrices = 1;
  double * rhs_entries = s.rhs_.get_entries();

  BoundValueFunct3D * non_const_bound_value[1] {example_.get_bd()[0]};

  //fetch stiffness matrix as block
  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix_.get_blocks_uniquely();
  TSquareMatrix3D * block[1]{reinterpret_cast<TSquareMatrix3D*>(blocks.at(0).get())};

  // Do the Assembling!

  // reset right hand side and matrix to zero
  s.rhs_.reset();
  block[0]->reset();
  //and call the method
  Assemble3D(1, &fe_space, N_Matrices, block, 0, NULL, 1, &rhs_entries,
             &fe_space, &boundary_conditions, non_const_bound_value, local_assem);

}

