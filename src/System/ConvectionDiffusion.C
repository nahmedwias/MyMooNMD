#include "ConvectionDiffusion.h"
#include "ConvDiff.h"
#include "Database.h"
#include "Domain.h"
#include "LocalAssembling.h"
#include "Multigrid.h"
#include "AlgebraicFluxCorrection.h"
#ifdef __2D__
 #include "Assemble2D.h"
 #include "SquareMatrix2D.h"
 #include "AuxParam2D.h"
#else
 #include "Assemble3D.h"
 #include "SquareMatrix3D.h"
 #include "AuxParam3D.h"
#endif
#ifdef _MPI
 #include "ParFECommunicator3D.h"
#endif

/* ************************************************************************* */
template<int d>
ParameterDatabase ConvectionDiffusion<d>::get_default_cd_database()
{
  Output::print<5>("creating a default ConvectionDiffusion parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default CD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("ConvectionDiffusion parameter database");
  
  // a default output database - needed here as long as there's no class handling the output
  db.merge(ParameterDatabase::default_output_database(), true);
  // a default afc database
  db.merge(AlgebraicFluxCorrection::default_afc_database(), true);
  // default local assembling database
  db.merge(LocalAssembling<d>::default_local_assembling_database(), true);

  return db;
}

/* ************************************************************************* */
template<int d>
ConvectionDiffusion<d>::SystemPerGrid::SystemPerGrid(
  const Example_CD& example, TCollection& coll, int ansatz_order)
: fe_space(new FESpace(&coll, "space", "conv diff fe_space", example.get_bc(0),
                       ansatz_order))
{
#ifdef _MPI
  fe_space->get_communicator().print_info();
#endif // _MPI
#ifdef __3D__
  matrix = BlockFEMatrix::CD3D(fe_space);
#else
  matrix = BlockFEMatrix::CD2D(fe_space);
#endif

  rhs = BlockVector(this->matrix, true);
  solution = BlockVector(this->matrix, false);
  fe_function = FEFunction(fe_space, "c", "c", solution.get_entries(),
                           solution.length());
}

/* ************************************************************************* */
template<int d>
ConvectionDiffusion<d>::ConvectionDiffusion(const TDomain& domain,
                                            const ParameterDatabase& param_db)
 : ConvectionDiffusion<d>(domain, param_db, Example_CD(param_db))
{
}

/* ************************************************************************* */
template<int d>
ConvectionDiffusion<d>::ConvectionDiffusion(const TDomain& domain,
                                            const ParameterDatabase& param_db,
                                            const Example_CD& example_cd)
 : systems(), example(example_cd), db(get_default_cd_database()),
   outputWriter(param_db), solver(param_db), errors()
{
  db.merge(param_db, false); // update this database with given values
  set_parameters();
  // The construction of the members differ, depending on whether
  // a multigrid solver will be used or not.
  bool usingMultigrid = solver.is_using_multigrid();
  int ansatz_order = TDatabase::ParamDB->ANSATZ_ORDER;
  
  auto collections = domain.get_grid_collections();
  if(!usingMultigrid)
  {
    // Get the collection on the finest grid
    TCollection& cellCollection = *collections.front();
    // create finite element space and function, a matrix, rhs, and solution
    systems.emplace_back(example, cellCollection, ansatz_order);
  }
  else
  {
    // we are using multigrid
    size_t n_levels = collections.size();
    size_t desirec_multigrid_levels = param_db["multigrid_n_levels"];
    if(desirec_multigrid_levels > n_levels)
      ErrThrow("Not enough collections (", n_levels, ") to use ",
               desirec_multigrid_levels, " multigrid levels!");
    // remove not needed coarser grid from list of collections
    for(size_t i = desirec_multigrid_levels; i < n_levels; ++i)
    {
      collections.pop_back();
    }

    auto mg = this->solver.get_multigrid();
    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    for(auto coll : collections)
    {
      Output::print("creating SystemPerGrid object, n_cells = ", coll->GetN_Cells());
      systems.emplace_back(example, *coll, ansatz_order);
      //prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    mg->initialize(matrices);
  }
  
  outputWriter.add_fe_function(&this->get_function());
  output_problem_size_info();
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::set_parameters()
{
  //set problem_type to CD if not yet set
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
  
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    ErrThrow("Ansatz order 0 is not used in convection diffusion "
             "reaction problems! (Vanishing convection and diffusion term).");
  }
  
  //////////////// Algebraic flux correction ////////////
  if(!db["algebraic_flux_correction"].is("none"))
  {//some kind of afc enabled
    if(!db["algebraic_flux_correction"].is("fem-tvd"))
    {
      db["algebraic_flux_correction"].set("fem-tvd");
      Output::print("Only kind of algebraic flux correction"
          " for CD problems is FEM-TVD (fem-tvd).");
    }
    //make sure that galerkin discretization is used
    if (!db["space_discretization_type"].is("galerkin"))
    {//some other disctype than galerkin
      db["space_discretization_type"] = "galerkin";
      Output::warn<1>("Parameter 'space_discretization_type' changed to 'galerkin' "
          "because Algebraic Flux Correction is enabled.");
    }
    // when using afc, create system matrices as if all dofs were active
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
  }
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::output_problem_size_info() const
{
  // print some useful information
  auto& space = *this->systems.front().fe_space;
  double hMin, hMax;
  auto coll = space.GetCollection();
  coll->GetHminHmax(&hMin, &hMax);
  Output::print<1>("N_Cells    : ", setw(13), coll->GetN_Cells());
  Output::print<1>("h(min, max): ", setw(13), hMin, " ", setw(13), hMax);
  Output::print<1>("dofs all   : ", setw(13), space.GetN_DegreesOfFreedom());
  Output::print<1>("dof active : ", setw(13), space.GetActiveBound());
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::assemble()
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  
  LocalAssembling_type laType = LocalAssembling_type::ConvDiff;
  // this loop has more than one iteration only in case of multigrid
  for(auto & s : systems)
  {
    FEFunction * feFunctionPtr = &s.fe_function;

    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling<d> laObject(this->db, laType, &feFunctionPtr,
                                example.get_coeffs());

    // assemble the system matrix with given local assembling, solution and rhs
    int n_fe_spaces = 1;
    const FESpace * fe_space = s.fe_space.get();
    auto * boundary_conditions = fe_space->get_boundary_condition();
    int n_square_matrices = 1;
    int n_rect_matrices = 0;
    double * rhs_entries = s.rhs.get_entries();
    int n_rhs = 1;
    auto * non_const_bound_value = example.get_bd(0);

    //fetch stiffness matrix as block
    auto blocks = s.matrix.get_blocks_uniquely();
    SquareMatrixD* block = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
    
    // reset right hand side and matrix to zero
    s.rhs.reset();
    block->reset();
    // and call the assmebling method
#ifdef __3D__
    Assemble3D(
#else      
    Assemble2D(
#endif
               n_fe_spaces, &fe_space, n_square_matrices, &block, 
               n_rect_matrices, nullptr, n_rhs, &rhs_entries, &fe_space,
               &boundary_conditions, &non_const_bound_value, laObject);
    
    // copy Dirichlet values from rhs to solution vector (this is not really
    // necessary in case of a direct solver, but also algebraic flux correction
    // assumes the solution has the correct Dirichlet dofs)
    s.solution.copy_nonactive(s.rhs);
  }
  if(!db["algebraic_flux_correction"].is("none"))
  {
    do_algebraic_flux_correction();
  }
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::solve()
{
  double t = GetTime();
  SystemPerGrid& s = this->systems.front();
#ifndef _MPI
  this->solver.solve(s.matrix, s.rhs, s.solution);
#else
  if(this->solver.get_db()["solver_type"].is("direct"))
  {
    MumpsWrapper mumps_wrapper(s.matrix);
    mumps_wrapper.solve(s.rhs, s.solution);
  }
  else
    this->solver.solve(s.matrix, s.rhs, s.solution);
#endif
  t = GetTime() - t;
  Output::print<3>("solving time of a ConvectionDiffusion<", d, "> problem: ",
                   t, " seconds");
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::output(int i)
{
  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  // computing errors as well as writing vtk files requires a minimum 
  // consistency level of 1
  this->systems.front().fe_space->get_communicator().consistency_update(
    this->systems.front().solution.get_entries(), 1);
#endif

  FEFunction & fe_function = this->systems.front().fe_function;
  fe_function.PrintMinMax();

  // write solution
  if(i < 0)
    outputWriter.write();
  else
    outputWriter.write(i);

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    std::array<double, ConvectionDiffusion<d>::n_errors+1> errors;
#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };
#else
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
#endif
    const FESpace* space = this->systems.front().fe_space.get();

    fe_function.GetErrors(example.get_exact(0), d+1, AllDerivatives,
                          ConvectionDiffusion<d>::n_errors,
                          conv_diff_l2_h1_linf_error<d>, example.get_coeffs(),
                          &aux, 1, &space, errors.data());
#ifdef _MPI
    // memory for global (across all processes) error
    std::array<double, ConvectionDiffusion<d>::n_errors> errorsReduced;
    MPI_Reduce(errors.data(), errorsReduced.data(),
               ConvectionDiffusion<d>::n_errors-1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&errors[ConvectionDiffusion<d>::n_errors-1],
               &errorsReduced[ConvectionDiffusion<d>::n_errors-1], 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    for(i = 0; i < ConvectionDiffusion<d>::n_errors-1; i++) // exclude l_inf
      errors[i] = std::sqrt(errorsReduced[i]);
    errors[ConvectionDiffusion<d>::n_errors-1] 
      = errorsReduced[ConvectionDiffusion<d>::n_errors-1];
#endif
    // copy local variable to member variable
    std::copy(errors.begin(), errors.end()-1, this->errors.begin());

    //print errors
    if(my_rank == 0)
    {
      Output::print<1>("L2     : ", setprecision(14), errors[0]);
      Output::print<1>("H1-semi: ", setprecision(14), errors[1]);
      Output::print<1>("SD     : ", setprecision(14), errors[2]);
      Output::print<1>("L_inf  : ", setprecision(14), errors[3]);
    }
  } // if(this->db["compute_errors"])
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::do_algebraic_flux_correction()
{
  for(auto & s : this->systems) // do it on all levels.
  {
    //determine which kind of afc to use
    if(db["algebraic_flux_correction"].is("default") ||
        db["algebraic_flux_correction"].is("fem-tvd"))
    {
      //get pointers/references to the relevant objects
      auto& feSpace = *s.fe_space;
      FEMatrix& one_block = *s.matrix.get_blocks_uniquely().at(0).get();
      const std::vector<double>& solEntries = s.solution.get_entries_vector();
      std::vector<double>& rhsEntries = s.rhs.get_entries_vector();

      // fill a vector "neumannToDirichlet" with those rows that got
      // internally treated as Neumann although they are Dirichlet
      int firstDiriDof = feSpace.GetActiveBound();
      int nDiri = feSpace.GetN_Dirichlet();

      std::vector<int> neumToDiri(nDiri, 0);
      std::iota(std::begin(neumToDiri), std::end(neumToDiri), firstDiriDof);

      // apply FEM-TVD
      AlgebraicFluxCorrection::fem_tvd_algorithm(
          one_block,
          solEntries,rhsEntries,
          neumToDiri);

      //...and finally correct the entries in the Dirchlet rows
      AlgebraicFluxCorrection::correct_dirichlet_rows(one_block);
      //...and in the right hand side, too, assum correct in solution vector
      s.rhs.copy_nonactive(s.solution);
    }
    else
    {
      ErrThrow("The chosen algebraic flux correction scheme is unknown "
               "to class ConvectionDiffusion<", d, ">.");
    }
  }
}

/** ************************************************************************ */
template<int d>
double ConvectionDiffusion<d>::get_L2_error() const
{
  return this->errors[0];
}

/** ************************************************************************ */
template<int d>
double ConvectionDiffusion<d>::get_H1_semi_error() const
{
  return this->errors[1];
}

/** ************************************************************************ */
template<int d>
double ConvectionDiffusion<d>::get_SD_error() const
{
  return this->errors[2];
}

/** ************************************************************************ */
template<int d>
double ConvectionDiffusion<d>::get_L_inf_error() const
{
  return this->errors[3];
}

#ifdef __3D__
template class ConvectionDiffusion<3>;
#else
template class ConvectionDiffusion<2>;
#endif
