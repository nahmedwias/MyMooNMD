#include <CD2D.h>
#include <Database.h>
#include <Multigrid.h>
#include <MainUtilities.h>                        // L2H1Errors
#include <AlgebraicFluxCorrection.h>
#include <PostProcessing2D.h>
#include <LocalAssembling2D.h>
#include <Assemble2D.h>
#include <Upwind.h>
#include <LocalProjection.h>

ParameterDatabase get_default_CD2D_parameters()
{
  Output::print<5>("creating a default CD2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default CD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("CD2D parameter database");

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default afc database
  ParameterDatabase afc_db = AlgebraicFluxCorrection::default_afc_database();
  db.merge(afc_db, true);

  return db;
}


/** ************************************************************************ */
CD2D::System_per_grid::System_per_grid(const Example_CD2D& example,
TCollection& coll, int ansatz_order)
: fe_space(new TFESpace2D(&coll, "space", "cd2d fe_space", example.get_bc(0),
ansatz_order, nullptr))
{
  matrix = BlockFEMatrix::CD2D(*fe_space);

  rhs = BlockVector(this->matrix, true);
  solution = BlockVector(this->matrix, false);

  fe_function = TFEFunction2D(this->fe_space.get(), "c", "c",
    this->solution.get_entries(), this->solution.length());

}


/** ************************************************************************ */
CD2D::CD2D(const TDomain& domain, const ParameterDatabase& param_db,
int reference_id)
: CD2D(domain, param_db, Example_CD2D(param_db), reference_id)
{
}


/** ************************************************************************ */
CD2D::CD2D(const TDomain& domain, const ParameterDatabase& param_db,
const Example_CD2D& example, int reference_id)
: systems(), example(example), db(get_default_CD2D_parameters()),
outputWriter(param_db), solver(param_db), errors()
{
  this->db.merge(param_db, false);                // update this database with given values
  this->set_parameters();
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  // create finite element space and function, a matrix, rhs, and solution
  int ansatz_order = TDatabase::ParamDB->ANSATZ_ORDER;
  this->systems.emplace_back(this->example, *coll, ansatz_order);
  outputWriter.add_fe_function(&this->get_function());

  // print out some information
  auto& space = *this->systems.front().fe_space;
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  Output::print<1>("N_Cells    : ", setw(12), coll->GetN_Cells());
  Output::print<2>("h (min,max): ", setw(12), h_min, " ", setw(12), h_max);
  Output::print<1>("dof all    : ", setw(12), space.GetN_DegreesOfFreedom());
  Output::print<2>("dof active : ", setw(12), space.GetN_ActiveDegrees());
  
  // done with the constructor in case we're not using multigrid
  if(!this->solver.is_using_multigrid())
    return;
  // else multigrid

  auto mg = this->solver.get_multigrid();
  bool mdml = mg->is_using_mdml();
  if(mdml)
  {
    // change the discretization to lowest order
    /// @todo for mdml: is P1/Q1 the correct space on the other grids? Maybe
    /// what we really need is say Q3/P3, Q2/P2, Q1/P1 on the finest grid and
    /// Q1/P1 on all coarser grids.
    if(ansatz_order == -1 || ansatz_order == 1)
    {
      // - using non conforming P1 already, it makes no sense to use another
      //   discretization on the finest grid.
      // - using conforming P1, we don't do another multigrid level with non
      //   conforming P1 elements on the finest grid, because this space is
      //   typically larger that conforming P1.
      // Either way we just do regular multigrid
      mdml = false;
    }
    else
      ansatz_order = -1;
  }
  if(mdml)
    /// @todo mdml for CD2D: We need a special assembling function which
    /// does not assemble the convection term. Instead one then calls an upwind
    /// method.
    ErrThrow("mdml is currently not working.");

  // number of multigrid levels
  size_t n_levels = mg->get_n_geometric_levels();
  // index of finest grid
  int finest = domain.get_ref_level();            // -> there are finest+1 grids
  // index of the coarsest grid used in this multigrid environment
  int coarsest = finest - n_levels + 1;
  if(mdml)
  {
    coarsest++;
  }
  else
  {
    // only for mdml there is another matrix on the finest grid, otherwise
    // the next system to be created is on the next coarser grid
    finest--;
  }
  if(coarsest < 0 )
  {
    ErrThrow("the domain has not been refined often enough to do multigrid "
      "on ", n_levels, " levels. There are only ",
      domain.get_ref_level() + 1, " grid levels.");
  }

  // Construct systems per grid and store them, finest level first
  std::list<BlockFEMatrix*> matrices;
  // matrix on finest grid is already constructed
  matrices.push_back(&systems.back().matrix);
  for (int grid_no = finest; grid_no >= coarsest; --grid_no)
  {
    TCollection *coll = domain.GetCollection(It_EQ, grid_no, reference_id);
    systems.emplace_back(example, *coll, ansatz_order);
    //prepare input argument for multigrid object
    matrices.push_front(&systems.back().matrix);
  }
  mg->initialize(matrices);
}


/** ************************************************************************ */
CD2D::~CD2D()
{
  // delete the collections created during the contructor
  for(auto & s : this->systems)
    delete s.fe_space->GetCollection();
}


/** ************************************************************************ */
void CD2D::set_parameters()
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
  //////////////// Algebraic flux correction ////////////
  
  //How to access the TDatabase from CD2D_AFC?
  if(!db["algebraic_flux_correction"].is("none"))
  {                                               //some kind of afc enabled
    if(!db["algebraic_flux_correction"].is("afc"))
    {
      db["algebraic_flux_correction"].set("afc");
      Output::print("Only kind of algebraic flux correction"
        " for CD problems is AFC (afc).");
    }
    //make sure that galerkin discretization is used
    if (!db["space_discretization_type"].is("galerkin"))
    {                                             //some other disctype than galerkin
      db["space_discretization_type"] = "galerkin";
      Output::warn<1>("Parameter 'space_discretization_type' changed to 'galerkin' "
        "because Algebraic Flux Correction is enabled.");
    }
    // when using afc, create system matrices as if all dofs were active
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
  }
}


/** ************************************************************************ */
void CD2D::assemble()
{
  LocalAssembling2D_type t = LocalAssembling2D_type::ConvDiff;
  bool mdml = this->solver.is_using_multigrid()
    && this->solver.get_multigrid()->is_using_mdml();
  // in case of mdml, we need to change the local assembling, (not yet
  // implemented)

  // this loop has more than one iteration only in case of multigrid
  for(auto & s : this->systems)
  {
    TFEFunction2D * pointer_to_function = &s.fe_function;
    int disc_type_code = 0;

    std::shared_ptr<LocalAssembling2D> la;

    if (db["space_discretization_type"].is("galerkin"))
    {
      disc_type_code = GALERKIN;
    }
    else if (db["space_discretization_type"].is("supg"))
    {
      disc_type_code = SUPG;
    }
    else if (db["space_discretization_type"].is("upwind"))
    {
      disc_type_code = GALERKIN;
    }
    else
      ErrThrow("space_discretization_type ", db["space_discretization_type"].get_name(), " not implemented !");

    {
      // create a local assembling object which is needed to assemble the matrix
      //LocalAssembling2D la(t, &pointer_to_function, example.get_coeffs(),disc_type_code);
      //disc_type_code = (int) db["space_discretization_type"];
      Output::print<4>("assembling discretization ",disc_type_code);
      la = std::make_shared<LocalAssembling2D>(t, &pointer_to_function, example.get_coeffs(),disc_type_code);
    }
    
    // assemble the system matrix with given local assembling, solution and rhs
    const TFESpace2D * fe_space = s.fe_space.get();
    BoundCondFunct2D * boundary_conditions = fe_space->GetBoundCondition();
    int N_Matrices = 1;
    double * rhs_entries = s.rhs.get_entries();

    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();
    TSquareMatrix2D * matrix = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    BoundValueFunct2D * non_const_bound_value[1] {example.get_bd()[0]};

    //Previous Implementation
    {
      s.rhs.reset();
      matrix->reset();
      Output::print<4>("call assemble");
      // assemble
      Assemble2D(1, &fe_space, N_Matrices, &matrix, 0, NULL, 1, &rhs_entries,
                 &fe_space, &boundary_conditions, non_const_bound_value, *la);
    }
    // apply local projection stabilization method
    if(db["space_discretization_type"].is("local_projection")
      && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
    {
      if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
      {
        UltraLocalProjection((void *)&matrix, false);
      }
      else
      {
        ErrThrow("LP_FULL_GRADIENT needs to be one to use LOCAL_PROJECTION");
      }
    }
    bool finest_grid = &systems.front() == &s;
    {
      if ((mdml && !finest_grid) || (db["space_discretization_type"].is("upwind")))
      {
        Output::print<2>("upwind for convection-diffusion equation");
        UpwindForConvDiff(la->GetCoeffFct(), matrix, rhs_entries, fe_space,
          nullptr, nullptr, false);
      }
    }
    // copy Dirichlet values from rhs to solution vector (this is not really
    // necessary in case of a direct solver)
    s.solution.copy_nonactive(s.rhs);
  }
}

/** *********************************************************************** */

void CD2D::solve()
{
  double solving_time = GetTime();
  System_per_grid& s = this->systems.front();
  this->solver.solve(s.matrix, s.rhs, s.solution);
  solving_time = GetTime() - solving_time;
  Output::print("  solving of a CD2D problem done in ", solving_time, " seconds");
}

/** ************************************************************************ */
void CD2D::output(int i)
{
  // print the value of the largest and smallest entry in the finite element
  // vector
  TFEFunction2D & fe_function = this->systems.front().fe_function;
  fe_function.PrintMinMax();
  this->example.do_post_processing(*this);

  // write solution to a vtk file or in case-format
  outputWriter.write(i);

  /*
  // implementation with the old class TOutput2D
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(1, 1, 0, 0, NULL);
    Output.AddFEFunction(&fe_function);

    // Create output directory, if not already existing.
    mkdir(db["output_vtk_directory"], 0777);
    std::string filename = this->db["output_vtk_directory"];
    filename += "/" + this->db["output_basename"].value_as_string();

  if(i >= 0)
  filename += "_" + std::to_string(i);
  filename += ".vtk";
  Output.WriteVtk(filename.c_str());
  }
  */

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(this->db["output_compute_errors"])
  {
    // this should be a little longer than this->errors, because of a bug in
    // FEFunction::GetErrors. Otherwise we could use this->errors directly.
    // Note that we can not write
    // 'constexpr size_t n_errors = errors.max_size();'. The reason is that the
    // method 'max_size' is not marked const in c++11, but it is in c++14. We
    // should switch to that.
    std::array<double, 5> errors;
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D* space = fe_function.GetFESpace2D();

    fe_function.GetErrors(this->example.get_exact(0), 3, AllDerivatives, 4,
      SDFEMErrors, this->example.get_coeffs(), &aux, 1,
      &space, errors.data());

    Output::print<1>("L2     : ", setprecision(14), errors[0]);
    Output::print<1>("H1-semi: ", setprecision(14), errors[1]);
    Output::print<1>("SD     : ", setprecision(14), errors[2]);
    Output::print<1>("L_inf  : ", setprecision(14), errors[3]);
        // copy local variable to member variable
    std::copy(errors.begin(), errors.end()-1, this->errors.begin());
  }
}



/** ************************************************************************ */
double CD2D::get_L2_error() const
{
  return this->errors[0];
}


/** ************************************************************************ */
double CD2D::get_H1_semi_error() const
{
  return this->errors[1];
}


/** ************************************************************************ */
double CD2D::get_SD_error() const
{
  return this->errors[2];
}


/** ************************************************************************ */
double CD2D::get_L_inf_error() const
{
  return this->errors[3];
}
