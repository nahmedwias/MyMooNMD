#include <CD2D.h>
#include <Database.h>
#include <Multigrid.h>
#include <MainUtilities.h> // L2H1Errors
#include <AlgebraicFluxCorrection.h>
#include "LocalAssembling.h"
#include <Assemble2D.h>
#include <Upwind.h>
#include <LocalProjection.h>

#include <Assembler4.h>
#include <BoundaryAssembling2D.h>
#include <AuxParam2D.h>


void Coefficient_Function_CD2D(double *in, double *out)
{
	// coordinates:  x at in[0], y at in[1]
	out[0] = in[2];
	out[1] = in[3];
}

/** ************************************************************************ */
ParameterDatabase CD2D::get_default_CD2D_parameters()
{
  Output::print<5>("creating a default CD2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default CD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("CD2D parameter database");
  
  db.add("coefficient_function_type", (int) 0, "Set the parameter equal to 0 if the"
		  " coefficients are constant, if you want to use a coefficient function that is spatially "
		  "varying and analytically defined in ParMooN_Brinkman2D.C, set this parameter equal to 1 "
		  "and then adjust coeffs via parameters in the example file; if you want to use a coefficient "
		  "function that is spatially varying and defined by a .mesh-file combined with a file describing "
		  "a corresponding FEFunction2D, set this parameter equal to 2; if you want to use two coefficient "
		  "function that are spatially varying and analytically defined, set this parameter equal to 3; if you want to use two coefficient "
		  "functions from files ... .Sol, set this parameter equal to 4 ", (int) 0.,(int) 4.);

  db.add("read_velocity1_function_directory", "." , "This allows to use coefficient_function_type 4 and .... "
  		"The File has to fit with the mesh (refinement level). A default file  is contained in input_files/ .");

  db.add("read_velocity2_function_directory", "." , "This allows to use coefficient_function_type 4 and .... "
  		"The File has to fit with the mesh (refinement level). A default file  is contained in input_files/ .");

  db.add("write_coefficient_function_directory", "." , "");

  db.add("write_velocity1_directory", "." , "This allows to save the computed Brinkmna2d velocity in "
  		"x-direction in a file for later use ( used in coefficient_function_type 3. ");
  db.add("write_velocity2_directory", "." , "This allows to save the computed Brinkmna2d velocity in "
  		"y-direction in a file for later use ( used in coefficient_function_type 3. ");

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default afc database
  ParameterDatabase afc_db = AlgebraicFluxCorrection::default_afc_database();
  db.merge(afc_db, true);
  // default local assembling database
  db.merge(LocalAssembling<2>::default_local_assembling_database(), true);

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
  this->db.merge(param_db, false); // update this database with given values
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
  Output::print<1>("h (min,max): ", setw(12), h_min, " ", setw(12), h_max);
  Output::print<1>("dof all    : ", setw(12), space.GetN_DegreesOfFreedom());
  Output::print<1>("dof active : ", setw(12), space.GetN_ActiveDegrees());

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
  int finest = domain.get_ref_level(); // -> there are finest+1 grids
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

/** ************************************************************************ */
void CD2D::assemble(TFEFunction2D* coefficient_function1, TFEFunction2D* coefficient_function2)
{
	const TFESpace2D *coefficient_function1_space;
	const TFESpace2D *coefficient_function2_space;

 LocalAssembling_type type;
 // use a list of LocalAssembling2D objects
 std::vector< std::shared_ptr <LocalAssembling2D >> la_list;


  type = LocalAssembling_type::ConvDiff;
  int disctype = 1; // Galerkin
  if(db["space_discretization_type"].is("supg") 
    || db["space_discretization_type"].is("sdfem"))
  {
    disctype = 2;
  }
  bool mdml = this->solver.is_using_multigrid()
             && this->solver.get_multigrid()->is_using_mdml();
  // in case of mdml, we need to change the local assembling, (not yet 
  // implemented)

  // this loop has more than one iteration only in case of multigrid
  for(auto & s : this->systems)
  {
	  TFEFunction2D * pointer_to_fe_functions[3] = {&s.fe_function, NULL, NULL};

	  if ((coefficient_function1) && (coefficient_function2))
	  {
	  	coefficient_function1_space = coefficient_function1->GetFESpace2D();
	  	pointer_to_fe_functions[1] = coefficient_function1;
	  	coefficient_function2_space = coefficient_function2->GetFESpace2D();
	  	pointer_to_fe_functions[2] = coefficient_function2;
	  }
	  else if (coefficient_function1)
	  {
	  	coefficient_function1_space = coefficient_function1->GetFESpace2D();
	  	pointer_to_fe_functions[1] = coefficient_function1;
	  }

	  std::shared_ptr <LocalAssembling2D> la(
      new LocalAssembling2D(this->db, type, pointer_to_fe_functions, 
                            this->example.get_coeffs(), disctype));

	  if  ((coefficient_function1) && (coefficient_function2))
	  {
	  	la->setBeginParameter({0});
	  	la->setN_ParamFct(1);
	  	la->setParameterFct({Coefficient_Function_CD2D});
	  	la->setN_FeValues(2);
	  	la->setFeValueFctIndex({1,2});
	  	la->setFeValueMultiIndex({D00, D00});
	  }
	  else if (coefficient_function1)
	  {
	  	// modify la such that it includes the TFEFunction2D coefficient_function
	  	la->setBeginParameter({0});
	  	la->setN_ParamFct(1);
	  	la->setParameterFct({Coefficient_Function_CD2D});
	  	la->setN_FeValues(1);
	  	la->setFeValueFctIndex({1});
	  	la->setFeValueMultiIndex({D00});
	  }

	  la_list.push_back(la);

    // assemble the system matrix with given local assembling, solution and rhs
    std::vector<const TFESpace2D*> fe_spaces;
    fe_spaces.resize(2);
    fe_spaces[0] =  s.fe_space.get();
    const TFESpace2D * fe_space = s.fe_space.get();

    auto * boundary_conditions = fe_spaces[0]->get_boundary_condition();
    int N_Matrices = 1;
    double * rhs_entries = s.rhs.get_entries();

    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();
    TSquareMatrix2D * matrix = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

    BoundValueFunct2D * non_const_bound_value[1] {example.get_bd()[0]};

    // reset right hand side and matrix to zero (just in case)
    s.rhs.reset();
    matrix->reset();

    std::vector<const TFESpace2D*> spaces_for_matrix;
    spaces_for_matrix.resize(1);

    if (coefficient_function1 && coefficient_function2)
    {
    	spaces_for_matrix.resize(3);
    	spaces_for_matrix[1] = coefficient_function1_space;
    	spaces_for_matrix[2] = coefficient_function2_space;
    }
    else if (coefficient_function1)
    {
    	spaces_for_matrix.resize(2);
    	spaces_for_matrix[1] = coefficient_function1_space;
    }

    spaces_for_matrix[0] = s.fe_space.get();

    std::vector<const TFESpace2D*> spaces_for_rhs;
    spaces_for_rhs.resize(1);
    spaces_for_rhs[0] =  s.fe_space.get();

    //Todo: Assembler4 does not produce the same errors as the old Assemble2D(), therefore I added the if else case. Find the reason therefore
    if (coefficient_function1)
    {
    	Assembler4 Ass;
    	Ass.Assemble2D(s.matrix,s.rhs,
    			spaces_for_matrix, spaces_for_rhs,
					example,la_list);
    }
    else
    {
    	// assemble
    	Assemble2D(1, &fe_spaces[0], N_Matrices, &matrix, 0, nullptr, 1, &rhs_entries,
    			&fe_spaces[0], &boundary_conditions, non_const_bound_value, *la);
    }

    for (int k = 0; k < TDatabase::ParamDB->n_neumann_boundary; k++)
    {
    	BoundaryAssembling2D::rhs_g_v_n(s.rhs, s.fe_space.get(),
    			this->example.get_bd(0), //nullptr,                       // g = 1
					TDatabase::ParamDB->neumann_boundary_id[k],               // boundary component
					-1.*TDatabase::ParamDB->neumann_boundary_value[k]);       // mult
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
      if(mdml && !finest_grid)
      {
        UpwindForConvDiff(la->GetCoeffFct(), matrix, rhs_entries, fe_space,
                          nullptr, nullptr, false);
      }

    // copy Dirichlet values from rhs to solution vector (this is not really
    // necessary in case of a direct solver)
    s.solution.copy_nonactive(s.rhs);
      // copy Dirichlet values from rhs to solution vector (this is not really
      // necessary in case of a direct solver)
      s.solution.copy_nonactive(s.rhs);
  }
  // when using afc, do it now
  if(!db["algebraic_flux_correction"].is("none"))
  {
  	do_algebraic_flux_correction();
    do_algebraic_flux_correction();
  }
}

/** ************************************************************************ */
void CD2D::solve()
{
  double t = GetTime();
  System_per_grid& s = this->systems.front();
  this->solver.solve(s.matrix, s.rhs, s.solution);
  
  t = GetTime() - t;
  Output::print<3>(" solving of a CD2D problem done in ", t, " seconds");
}

/** ************************************************************************ */
void CD2D::output(int i)
{
  // print the value of the largest and smallest entry in the finite element 
  // vector
  TFEFunction2D & fe_function = this->systems.front().fe_function;
  fe_function.PrintMinMax();
  
  // write solution to a vtk file or in case-format
  if(i < 0)
    outputWriter.write();
  else
    outputWriter.write(i);

  /*
  // implementation with the old class TOutput2D
  {
    // last argument in the following is domain, but is never used in this class
    TOutput2D Output(1, 1, 0, 0, nullptr);
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
void CD2D::do_algebraic_flux_correction()
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
          "to class CD2D.");
    }
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

