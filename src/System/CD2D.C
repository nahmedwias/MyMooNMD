#include <CD2D.h>
#include <Database.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <OldSolver.h>
#include <MultiGrid2D.h>
#include <MainUtilities.h> // L2H1Errors
#include <AlgebraicFluxCorrection.h>

#include <LocalAssembling2D.h>
#include <Assemble2D.h>
#include <LocalProjection.h>

#include <numeric>

#include <sys/stat.h>

ParameterDatabase get_default_CD2D_parameters()
{
  Output::print<3>("creating a default CD2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default CD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("CD2D parameter database");
  
  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  return db;
}
/** ************************************************************************ */
CD2D::System_per_grid::System_per_grid(const Example_CD2D& example,
                                       TCollection& coll)
: fe_space(&coll, (char*)"space", (char*)"cd2d fe_space", example.get_bc(0),
           TDatabase::ParamDB->ANSATZ_ORDER, nullptr),
           // TODO CB: Building the matrix here and rebuilding later is due to the
           // highly non-functional class TFEVectFunction2D (and TFEFunction2D,
           // which do neither provide default constructors nor working copy assignments.)
           matrix({&fe_space}),
           rhs(this->matrix, true),
           solution(this->matrix, false),
           fe_function(&this->fe_space, (char*)"c", (char*)"c",
                       this->solution.get_entries(), this->solution.length())
{
  matrix = BlockFEMatrix::CD2D(fe_space);
}
/** ************************************************************************ */
TSquareMatrix2D* CD2D::System_per_grid::get_matrix_pointer()
{
  std::vector<std::shared_ptr<FEMatrix>> blocks =
      matrix.get_blocks_TERRIBLY_UNSAFE();
  return reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
}

/** ************************************************************************ */
CD2D::CD2D(const TDomain& domain, const ParameterDatabase& param_db,
           int reference_id)
 : CD2D(domain, param_db, Example_CD2D(), reference_id)
{
}

/** ************************************************************************ */
CD2D::CD2D(const TDomain& domain, const ParameterDatabase& param_db,
           const Example_CD2D& example, int reference_id)
 : systems(), example(example), multigrid(nullptr), 
   db(get_default_CD2D_parameters()), solver(param_db)
{
  this->db.merge(param_db, false); // update this database with given values
  this->set_parameters();
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  
  // create finite element space and function, a matrix, rhs, and solution
  this->systems.emplace_back(this->example, *coll);
  
  
  // print out some information
  TFESpace2D & space = this->systems.front().fe_space;
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  Output::print<1>("N_Cells    : ", setw(12), coll->GetN_Cells());
  Output::print<2>("h (min,max): ", setw(12), h_min, " ", setw(12), h_max);
  Output::print<1>("dof all    : ", setw(12), space.GetN_DegreesOfFreedom());
  Output::print<2>("dof active : ", setw(12), space.GetN_ActiveDegrees());
  
  
  // done with the conrtuctor in case we're not using multigrid
  if(this->solver.get_db()["solver_type"].is("direct") || 
     !this->solver.get_db()["preconditioner"].is("multigrid"))
    return;
  // else multigrid
  
  // create spaces, functions, matrices on coarser levels
  double *param = new double[2]; // memory leak
  //param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  //param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
  param[0] = this->solver.get_db()["damping_factor"];
  param[1] = this->solver.get_db()["damping_factor_finest_grid"];
  this->multigrid.reset(new TMultiGrid2D(1, 2, param));
  // number of refinement levels for the multigrid
  size_t LEVELS = this->solver.get_db()["n_multigrid_levels"];
  if((int)LEVELS > domain.get_ref_level() + 1)
    LEVELS = domain.get_ref_level() + 1;
  
  // the matrix and rhs side on the finest grid are already constructed 
  // now construct all matrices, rhs, and solutions on coarser grids
  for(int i = LEVELS - 2; i >= 0; i--)
  {
    unsigned int grid = i + domain.get_ref_level() + 1 - LEVELS;
    TCollection *coll = domain.GetCollection(It_EQ, grid, reference_id);
    this->systems.emplace_back(example, *coll);
  }
  
  // create multigrid-level-objects, must be coarsest first
  unsigned int i = 0;
  for(auto it = this->systems.rbegin(); it != this->systems.rend(); ++it)
  {
    //get a non-const pointer to the one block that "matrix" stores
    // TODO must be changed to const as soon as multigrid allows that
    TSquareMatrix2D* one_block = it->get_matrix_pointer();

    TMGLevel2D *multigrid_level = new TMGLevel2D(
      i, one_block, it->rhs.get_entries(),
      it->solution.get_entries(), 2, NULL);
    i++;
    this->multigrid->AddLevel(multigrid_level);
  }
}

/** ************************************************************************ */
CD2D::~CD2D()
{
  // delete the collections created during the contructor
  for(auto & s : this->systems)
    delete s.fe_space.GetCollection();
}

/** ************************************************************************ */
void CD2D::set_parameters()
{
  //////////////// Algebraic flux correction ////////////
  if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION == 1)
  {//some kind of afc enabled
    //make sure that galerkin discretization is used
    if (TDatabase::ParamDB->DISCTYPE !=	1)
    {//some other disctype than galerkin
      TDatabase::ParamDB->DISCTYPE = 1;
      Output::print("DISCTYPE changed to 1 (GALERKIN) because Algebraic Flux ",
                    "Correction is enabled.");
    }
    // when using afc, create system matrices as if all dofs were active
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
  }
  if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION > 1)
  {
    ErrThrow("For CD2D only algebraic flux correction of FEM-TVD is implemented"
        "(ALGEBRAIC_FLUX_CORRECTION: 1).")
  }

}

/** ************************************************************************ */
void CD2D::assemble()
{
  LocalAssembling2D_type t = LocalAssembling2D_type::ConvDiff;

  // this loop has more than one iteration only in case of multigrid
  for(auto & s : this->systems)
  {
    TFEFunction2D * pointer_to_function = &s.fe_function;
    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling2D la(t, &pointer_to_function, example.get_coeffs());

    // assemble the system matrix with given local assembling, solution and rhs
    const TFESpace2D * fe_space = &s.fe_space;
    BoundCondFunct2D * boundary_conditions = fe_space->GetBoundCondition();
    int N_Matrices = 1;
    double * rhs_entries = s.rhs.get_entries();

    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();
    TSquareMatrix2D * matrix = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

    BoundValueFunct2D * non_const_bound_value[1] {example.get_bd()[0]};

      // reset right hand side and matrix to zero (just in case)
      s.rhs.reset();
      matrix->reset();

      // assemble
      Assemble2D(1, &fe_space, N_Matrices, &matrix, 0, NULL, 1, &rhs_entries,
                 &fe_space, &boundary_conditions, non_const_bound_value, la);

      // apply local projection stabilization method
      if(TDatabase::ParamDB->DISCTYPE==LOCAL_PROJECTION
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

      // copy Dirichlet values from rhs to solution vector (this is not really
      // necessary in case of a direct solver)
      s.solution.copy_nonactive(s.rhs);
  }

  // when using afc, do it now
  if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION == 1)
  {
    do_algebraic_flux_correction();
  }

}

/** ************************************************************************ */
void CD2D::solve()
{
  double t = GetTime();
  System_per_grid& s = this->systems.front();
  if(this->solver.get_db()["solver_type"].is("direct") || 
     !this->solver.get_db()["preconditioner"].is("multigrid"))
  {
    this->solver.solve(s.matrix, s.rhs, s.solution);
  }
  else
  {
    // the one matrix stored in this BlockMatrix
    FEMatrix* mat = s.matrix.get_blocks_uniquely().at(0).get();
    // in order for the old method 'Solver' to work we need TSquareMatrix2D
    TSquareMatrix2D *SqMat[1] = { reinterpret_cast<TSquareMatrix2D*>(mat) };
    OldSolver((TSquareMatrix **)SqMat, NULL, s.rhs.get_entries(), 
              s.solution.get_entries(), MatVect_Scalar, Defect_Scalar, 
              this->multigrid.get(), this->get_size(), 0);
  }
  
  t = GetTime() - t;
  Output::print<2>(" solving of a CD2D problem done in ", t, " seconds");
}

/** ************************************************************************ */
void CD2D::output(int i)
{
	bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
	if(no_output)
		return;
  
  // print the value of the largest and smallest entry in the finite element 
  // vector
  TFEFunction2D & fe_function = this->systems.front().fe_function;
  fe_function.PrintMinMax();
  
  // write solution to a vtk file
  if(db["output_write_vtk"])
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
  
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    double errors[5];
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D* space = fe_function.GetFESpace2D();
    
    fe_function.GetErrors(this->example.get_exact(0), 3, AllDerivatives, 4,
                          SDFEMErrors, this->example.get_coeffs(), &aux, 1, 
                          &space, errors);
    
    Output::print<1>("L2     : ", errors[0]);
    Output::print<1>("H1-semi: ", errors[1]);
    Output::print<1>("SD     : ", errors[2]);
    Output::print<1>("L_inf  : ", errors[3]);
  } 
}

/** ************************************************************************ */
void CD2D::do_algebraic_flux_correction()
{
  for(auto & s : this->systems) // do it on all levels - TODO untested for multigrid!
  {
    //determine which kind of afc to use
    switch (TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION)
    {
      case 1: //FEM-TVD
      {
        //get pointers/references to the relevant objects
        TFESpace2D& feSpace = s.fe_space;
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
        break;
      }
      default:
      {
        ErrThrow("The chosen ALGEBRAIC_FLUX_CORRECTION scheme is unknown to class CD2D.");
      }
    }
  }
}
