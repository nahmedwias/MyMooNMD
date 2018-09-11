#include <CD3D.h>
#include <Example_CD3D.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <LinAlg.h>
#include "LocalAssembling.h"
#include <Assemble3D.h>

#include <DirectSolver.h>

#include <MainUtilities.h> // L2H1Errors

#include <Multigrid.h>

#include <sys/stat.h>

#ifdef _MPI
#include <MumpsWrapper.h>
#include "mpi.h"
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

ParameterDatabase get_default_CD3D_parameters()
{
  Output::print<5>("creating a default CD3D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default CD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("CD3D parameter database");
  
  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);
  
  // a default local assembling database
  db.merge(LocalAssembling3D::default_local_assembling_database());

  return db;
}

#ifdef _MPI
  CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example,
                                     TCollection& coll, int maxSubDomainPerDof)
   : feSpace_(new TFESpace3D(&coll, "space", "cd3d fe_space", example.get_bc(0),
              TDatabase::ParamDB->ANSATZ_ORDER))
  {
    //inform the fe space about the maximum number of subdomains per dof
    feSpace_->initialize_parallel(maxSubDomainPerDof);
    feSpace_->get_communicator().print_info();

    // set the matrix with named constructor
    matrix_ = BlockFEMatrix::CD3D(*feSpace_);

    rhs_ = BlockVector(matrix_, true);
    solution_ = BlockVector(matrix_, false);

    feFunction_ = TFEFunction3D(feSpace_.get(), "c", "c",
                                solution_.get_entries(), solution_.length());

  }
#else
  /* ************************************************************************ */
  CD3D::SystemPerGrid::SystemPerGrid(const Example_CD3D& example,
                                     TCollection& coll)
   : feSpace_(new TFESpace3D(&coll, "space", "cd3d fe_space", example.get_bc(0),
              TDatabase::ParamDB->ANSATZ_ORDER))
  {
    // set the matrix with named constructor
    matrix_ = BlockFEMatrix::CD3D(*feSpace_);

    rhs_ = BlockVector(matrix_, true);
    solution_ = BlockVector(matrix_, false);

    feFunction_ = TFEFunction3D(feSpace_.get(), "c", "c", 
                                solution_.get_entries(), solution_.length());
  }
#endif

  /** ************************************************************************ */
  CD3D::CD3D(std::list<TCollection*> collections,
             const ParameterDatabase& param_db, const Example_CD3D& example
#ifdef _MPI
             ,int maxSubDomainPerDof
#endif
  )
  : systems_(), example_(example), db(get_default_CD3D_parameters()),
    outputWriter(param_db), solver(param_db), errors_()
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
#endif
    }
    else
    {// we are using multigrid
        size_t n_levels = collections.size();
        if(!param_db["multigrid_n_levels"].is(n_levels))
           ErrThrow("Number of collection does not equal number of multigrid levels!");

        auto mg = this->solver.get_multigrid();
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
        mg->initialize(matrices);
      }
      // print useful information
      this->output_problem_size_info();
  }

/** ************************************************************************ */
//==============================================================================
void CD3D::output_problem_size_info() const
{
  // print some useful information
  auto& space = *this->systems_.front().feSpace_;
  double hMin, hMax;
  TCollection *coll = space.GetCollection();
  coll->GetHminHmax(&hMin, &hMax);
  Output::print<1>("N_Cells    : ", setw(13), coll->GetN_Cells());
  Output::print<1>("h(min, max): ", setw(13), hMin, " ", setw(13), hMax);
  Output::print<1>("dofs all   : ", setw(13), space.GetN_DegreesOfFreedom());
  Output::print<1>("dof active : ", setw(13), space.GetActiveBound());
}
/** ************************************************************************ */
void CD3D::assemble()
{
  LocalAssembling_type laType = LocalAssembling_type::ConvDiff;
  // this loop has more than one iteration only in case of multigrid
  for(auto & s : systems_)
  {
    TFEFunction3D * feFunctionPtr = &s.feFunction_;

    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling3D laObject(this->db, laType, &feFunctionPtr,
                               example_.get_coeffs());

    // assemble the system matrix with given local assembling, solution and rhs
    //s.matrix_.assemble(laObject, s.solution_, s.rhs_);
    call_assembling_routine(s, laObject);
  }

}

/** ************************************************************************ */
void CD3D::solve()
{
    SystemPerGrid& s = this->systems_.front();
    #ifndef _MPI
    this->solver.solve(s.matrix_, s.rhs_, s.solution_);
    #endif
    #ifdef _MPI
    if(this->solver.get_db()["solver_type"].is("direct"))
    {
      //set up a MUMPS wrapper
      MumpsWrapper mumps_wrapper(s.matrix_);
      //kick off the solving process
      mumps_wrapper.solve(s.rhs_, s.solution_);
    }
    else
      this->solver.solve(s.matrix_, s.rhs_, s.solution_);
    #endif
 }

/** ************************************************************************ */
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

#ifdef _MPI
  // computing errors as well as writing vtk files requires a minimum 
  // consistency level of 1
  syst.feSpace_->get_communicator().consistency_update(
    syst.solution_.get_entries(), 1);
#endif // _MPI
  
  // write solution to a vtk file
  outputWriter.add_fe_function(&syst.feFunction_);
  outputWriter.write();

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    double errors[5];
    TAuxParam3D aux(1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
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
  if(!db["problem_type"].is(1))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 1;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to CD."
          "It is now reset to the correct value for CD (=1).");
      db["problem_type"] = 1;
    }
  }

  //an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }
}

void CD3D::call_assembling_routine(SystemPerGrid& s, LocalAssembling3D& local_assem)
{//FIXME the body of this function was copy and paste

  const TFESpace3D * fe_space = s.feSpace_.get();
  auto * boundary_conditions = fe_space->get_boundary_condition();
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
  Assemble3D(1, &fe_space, N_Matrices, block, 0, nullptr, 1, &rhs_entries,
             &fe_space, &boundary_conditions, non_const_bound_value, local_assem);

}

