#include <Time_CD2D.h>
#include <Database.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <AlgebraicFluxCorrection.h>
#include <LocalAssembling2D.h>
#include <Assemble2D.h>
#include <LocalProjection.h>
#include <ConvDiff2D.h>
#include <NSE2D_FixPo.h>
#include <FEFunctionInterpolator.h>

/**************************************************************************** */
ParameterDatabase get_default_TCD2D_parameters()
{
  Output::print<5>("creating a default TCD2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default TCD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TCD2D parameter database");

  //TCD2D requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  // a default afc database
  ParameterDatabase afc_db = AlgebraicFluxCorrection::default_afc_database();
  db.merge(afc_db, true);

  return db;
}

/**************************************************************************** */
Time_CD2D::System_per_grid::System_per_grid(const Example_TimeCD2D& example,
                                            TCollection& coll)
: fe_space(&coll, (char*)"space", (char*)"time_cd2d space", example.get_bc(0),
           TDatabase::ParamDB->ANSATZ_ORDER, nullptr),
           // TODO CB: Building the matrix here and rebuilding later is due to the
           // highly non-functional class TFEVectFunction2D (and TFEFunction2D,
           // which do neither provide default constructors nor working copy assignments.)
           stiff_matrix({&fe_space}),
           mass_matrix({&fe_space}),
           rhs(this->stiff_matrix, true),
           solution(this->stiff_matrix, false),
           old_Au(this->stiff_matrix, true),
           fe_function(&this->fe_space, (char*)"c", (char*)"c",
                       this->solution.get_entries(), this->solution.length())
{
  stiff_matrix = BlockFEMatrix::CD2D(fe_space);
  mass_matrix = BlockFEMatrix::CD2D(fe_space);
}

/**************************************************************************** */
void Time_CD2D::System_per_grid::descale_stiff_matrix(double tau, double theta_1)
{
  if (tau==0 || theta_1 == 0)
  {
    ErrThrow("I will not divide by zero in descaling! tau=", tau, "theta_1=", theta_1);
  }
  //subtract the mass matrix (TODO smells like cancellation)...
  const FEMatrix& mass_block = *mass_matrix.get_blocks().at(0).get();
  //subtract the mass matrix...
  stiff_matrix.add_matrix_actives(mass_block, -1.0, {{0,0}}, {false});
  //...and descale the stiffness matrix...
  const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
  stiff_matrix.scale_blocks_actives(1./(tau*theta_1), cell_positions);
}

/**************************************************************************** */
void Time_CD2D::System_per_grid::update_old_Au()
{
  // the stiffness matrix must have been descaled (pure transport operator)
  stiff_matrix.apply(solution,old_Au);
}

/**************************************************************************** */
Time_CD2D::Time_CD2D(const TDomain& domain, const ParameterDatabase& param_db,
		int reference_id)
 : Time_CD2D(domain, param_db, Example_TimeCD2D(param_db), reference_id)
{
  
}

/**************************************************************************** */
Time_CD2D::Time_CD2D(const TDomain& domain, const ParameterDatabase& param_db,
		const Example_TimeCD2D& ex, int reference_id)
 : db(get_default_TCD2D_parameters()), solver(param_db), systems(), example(ex),
   errors(5, 0.0), timeDependentOutput(param_db)
{
  db.merge(param_db, false);
  this->set_parameters();
  // this is the L^inf(L^2(Omega)) error, initialize with some large number
  errors[4] = 1.e10; 
  
  bool usingMultigrid = this->solver.is_using_multigrid();
  if(!usingMultigrid)
  {
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
    // create finite element space and function, a matrix, rhs, and solution
    this->systems.emplace_back(this->example, *coll);
    
    TFESpace2D& space = this->systems.front().fe_space;
    double hmin, hmax;
    coll->GetHminHmax(&hmin, &hmax);
    Output::stat("TCD2D", "Mesh data and problem size");
    Output::dash("N_Cells            :  ", setw(12), coll->GetN_Cells());
    Output::dash("h(min,max)         :  ", setw(12), hmin, " ", setw(12), hmax);
    Output::dash("dof                :  ", setw(12), space.GetN_DegreesOfFreedom());
    Output::dash("active dof         :  ", setw(12), space.GetN_ActiveDegrees());
    
    // old right hand side
    old_rhs.copy_structure(this->systems[0].rhs);
    
    // interpolate the initial data
    TFEFunction2D & fe_function = this->systems.front().fe_function;
    fe_function.Interpolate(example.get_initial_cond(0));
    // add the fe function to the output object. 
    timeDependentOutput.add_fe_function(&fe_function);
  }
  else
  {
    auto multigrid = this->solver.get_multigrid();
    if(multigrid->is_using_mdml())
      ErrThrow("mdml for TCD2D not yet implemented");
    
    size_t nMGLevel = multigrid->get_n_geometric_levels();
    int finest = domain.get_ref_level();
    int coarsest = finest-nMGLevel+1;
    
    std::list<BlockFEMatrix*> matrices;
    for(int grid_no = finest; grid_no >= coarsest; --grid_no)
    {
      TCollection *coll = domain.GetCollection(It_EQ, grid_no, -4711);
      systems.emplace_back(example, *coll);
      
      systems.front().fe_function.Interpolate(example.get_initial_cond(0));
      matrices.push_front(&systems.back().stiff_matrix);
    }
    multigrid->initialize(matrices);
  }
}

/**************************************************************************** */
void Time_CD2D::set_parameters()
{

  if(TDatabase::TimeDB->TIME_DISC == 0)
  {
    ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC 
          << " is not supported");
    throw("TIME_DISC: 0 is not supported. Chose 1 (backward Euler) or 2 (Crank-Nicolson).");
  }

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
    if(!db["algebraic_flux_correction"].is("fem-fct-cn"))
    {
      db["algebraic_flux_correction"].set("fem-fct-cn");
      Output::print("Only kind of algebraic flux correction"
          " for TCD problems is Crank-Nicolson FEM-FCT (fem-fct-cn).");
    }

    if(solver.is_using_multigrid())
      ErrThrow("Multgrid and FEM-FCT are not yet enabled in Time_CD2D"); //TODO

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
    if (TDatabase::ParamDB->DISCTYPE != 1)
    {//some other disctype than galerkin
      TDatabase::ParamDB->DISCTYPE = 1;
      Output::print("DISCTYPE changed to 1 (GALERKIN) because Algebraic Flux ",
                    "Correction is enabled.");
    }

    // when using afc, create system matrices as if all dofs were active
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1; //FIXME THIS IS DANGEROUS, does not go well with BlockFEMatrix!
  }
}

/**************************************************************************** */
void Time_CD2D::assemble_initial_time()
{
  LocalAssembling2D_type mass = LocalAssembling2D_type::TCD2D_Mass;
  LocalAssembling2D_type stiff_rhs = LocalAssembling2D_type::TCD2D;

  for(auto &s : this->systems)
  {
    // assemble mass matrix, stiffness matrix and rhs
    TFEFunction2D* fe_funct[1] = {&s.fe_function}; //wrap up as '**'
    LocalAssembling2D la_mass(mass, fe_funct,
                              this->example.get_coeffs());
    LocalAssembling2D la_a_rhs(stiff_rhs, fe_funct,
                               this->example.get_coeffs());

    call_assembling_routine(s, la_a_rhs, la_mass , true,false);

    // apply local projection stabilization method on stiffness matrix only!
    if(TDatabase::ParamDB->DISCTYPE==LOCAL_PROJECTION
        && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
    {
      if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
      {
        //fetch stiffness matrix as block
        std::vector<std::shared_ptr<FEMatrix>> stiff_blocks =
            s.stiff_matrix.get_blocks_uniquely();
        auto stiff_block = stiff_blocks.at(0).get();
        //call ultra local projection
        UltraLocalProjection((void *)stiff_block, false);
      }
      else
      {
        ErrThrow("LP_FULL_GRADIENT needs to be 1 to use LOCAL_PROJECTION");
      }

    }

    s.stiff_matrix.apply(s.solution,s.old_Au); //put up initial old_Au
  }

  System_per_grid& s = systems.front();
  old_rhs = s.rhs; //put up initial old_rhs
}

/**************************************************************************** */
void Time_CD2D::assemble()
{
  LocalAssembling2D_type stiff_rhs = LocalAssembling2D_type::TCD2D;  
  
  for(auto &s : this->systems)
  {

    TFEFunction2D * pointer_to_function = &s.fe_function;
    // create two local assembling object (second one will only be needed in SUPG case)
    LocalAssembling2D la_a_rhs(stiff_rhs, &pointer_to_function,
                               this->example.get_coeffs());

    if(TDatabase::ParamDB->DISCTYPE == SUPG)
    {
      // In the SUPG case:
      // M = (u,v) + \tau (u,b.grad v)
      LocalAssembling2D_type mass_supg = LocalAssembling2D_type::TCD2D_Mass;;

      LocalAssembling2D la_m_supg(mass_supg, &pointer_to_function,
                               this->example.get_coeffs());

      //call assembling, including mass matrix (SUPG!)
      call_assembling_routine(s, la_a_rhs, la_m_supg , true,false);
    }
    else
    {//call assembling, ignoring mass matrix part (third argument is not relevant)
      call_assembling_routine(s, la_a_rhs, la_a_rhs , false,false);
    }
  }
  
  // here the modifications due to time discretization begin
  if (!db["algebraic_flux_correction"].is("none") )
  {
    do_algebraic_flux_correction();
    return; // modifications due to time discretization are per-
            // formed inside the afc scheme, so step out here!
  }

  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  // preparing rhs 
  System_per_grid& s = this->systems.front();
  
  if(TDatabase::TimeDB->THETA4)
  {
    // scale by time step length and theta4 (only active dofs)
    s.rhs.scaleActive(tau*TDatabase::TimeDB->THETA4);    
    // add old right hand side scaled by time step length and theta3 (only 
    // active dofs)
    if(TDatabase::TimeDB->THETA3 != 0.)
      s.rhs.addScaledActive((this->old_rhs), tau*TDatabase::TimeDB->THETA3);	

    // save old right hand side (only if THETA3 != 0)
    if(TDatabase::TimeDB->THETA3)
    {
      this->old_rhs.addScaledActive(s.rhs, -1./(tau*TDatabase::TimeDB->THETA3));
      this->old_rhs.scaleActive(-TDatabase::TimeDB->THETA3/TDatabase::TimeDB->THETA4);
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
  
  for(auto &s : this->systems)
  {
    // rhs += M*uold
    s.mass_matrix.apply_scaled_add_actives(s.solution, s.rhs, 1.0);
    // rhs -= tau*theta2*A_old*uold
    s.rhs.addScaledActive(s.old_Au, -tau*TDatabase::TimeDB->THETA2);
    
    // preparing the left hand side, i.e., the system matrix
    // stiffness matrix is scaled by tau*THETA1, after solving 
    // the matrix needs to be descaled if the coeffs does not depends
    // on time

    //scale  stiffness matrix...
    const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
    s.stiff_matrix.scale_blocks_actives(tau*TDatabase::TimeDB->THETA1, cell_positions);
    // ...and add the mass matrix
    const FEMatrix& mass_block = *s.mass_matrix.get_blocks().at(0).get();
    s.stiff_matrix.add_matrix_actives(mass_block, 1.0, {{0,0}}, {false});
  }  
  //this->systems[0].rhs.copy_nonactive(this->systems[0].solution);
  systems[0].solution.copy_nonactive(systems[0].rhs);
}

/**************************************************************************** */
void Time_CD2D::descale_stiffness(double tau, double theta_1)
{
    // restore stiffness matrix and store old_Au on all grids
    if(db["algebraic_flux_correction"].is("none"))
    {
      for(auto &s : this->systems)
      {
        s.descale_stiff_matrix(tau, theta_1);
        s.update_old_Au();
      }
    }
    else
    {
      //in AFC case the stiffness matrix is "ruined" by now -
      Output::print("AFC does not yet reset the stiffness"
          "matrix and old_Au correctly!");
    }
}

/**************************************************************************** */
void Time_CD2D::solve()
{
  double t = GetTime();
  System_per_grid& s = this->systems.front();
  
  this->solver.solve(s.stiff_matrix, s.rhs, s.solution);
  
  t = GetTime() - t;
  Output::print<1>("time for solving: ",  t);
  Output::print<2>("solution ", sqrt(Ddot(s.solution.length(),
                                          s.solution.get_entries(),
                                          s.solution.get_entries())) );
}

/**************************************************************************** */
void Time_CD2D::output()
{
  bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
  if(no_output)
    return;
  TFEFunction2D & fe_function = this->systems.front().fe_function;
  fe_function.PrintMinMax();
  
  if(db["output_compute_errors"])
  {
    double loc_e[5];
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D* space = fe_function.GetFESpace2D();
    
    fe_function.GetErrors(this->example.get_exact(0), 3, AllDerivatives, 4,
                          SDFEMErrors, this->example.get_coeffs(), &aux, 1, 
                          &space, loc_e);
    
    Output::print<1>("time: ", TDatabase::TimeDB->CURRENTTIME);
    Output::print<1>("  L2: ", loc_e[0]);
    Output::print<1>("  H1-semi: ", loc_e[1]);
    double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    errors[0] += (loc_e[0]*loc_e[0] + errors[1])*tau*0.5;
    errors[1] = loc_e[0]*loc_e[0];
    Output::print<1>("  L2(0,T;L2) ", sqrt(errors[0]));
    errors[2] += (loc_e[1]*loc_e[1] + errors[3])*tau*0.5;
    errors[3] = loc_e[1]*loc_e[1];
    Output::print<1>("  L2(0,T;H1) ", sqrt(errors[2]));
    
    if(errors[4] < loc_e[0])
      errors[4] = loc_e[0];
    Output::print<1>("  Linfty(0,T;L2) ", errors[4]);
  }

  // write output
  timeDependentOutput.write(TDatabase::TimeDB->CURRENTTIME);

}

/**************************************************************************** */
std::array< double, int(3) > Time_CD2D::get_errors() const
{
  std::array<double, int(3)> error_at_time_points;
  error_at_time_points.at(0) = sqrt(errors.at(0));
  error_at_time_points.at(1) = sqrt(errors.at(2));
  error_at_time_points.at(2) = sqrt(errors.at(4));
  return error_at_time_points;

    
    
}

/**************************************************************************** */

void Time_CD2D::do_algebraic_flux_correction()
{
  //determine which kind of afc to use
  if(db["algebraic_flux_correction"].is("default") ||
      db["algebraic_flux_correction"].is("fem-fct-cn"))
  {
    //TODO implement for multigrid!
    System_per_grid& s = systems.front();

    //get references to the relevant objects

    const TFESpace2D& feSpace = s.fe_space; //space
    const FEMatrix& mass = *s.mass_matrix.get_blocks().at(0).get(); // mass
    FEMatrix& stiff = *s.stiff_matrix.get_blocks_uniquely().at(0).get(); //stiffness
    //vector entry arrays
    const std::vector<double>& solEntries = s.solution.get_entries_vector();
    std::vector<double>& rhsEntries = s.rhs.get_entries_vector();
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
    double delta_t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

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
    s.rhs.copy_nonactive(old_rhs);
  }
  else
  {
    ErrThrow("The chosen algebraic flux correction scheme is unknown "
        "to class Time_CD2D.");
  }
}

void Time_CD2D::call_assembling_routine(
    System_per_grid& s,
    LocalAssembling2D& la_stiff, LocalAssembling2D& la_mass,
    bool assemble_both, bool with_convection_field)
{

  // Assemble mass matrix, stiffness matrix and rhs
  //...variables which are the same for both
  size_t nfespace = 1;
  const TFESpace2D * fe_spaces[3] = {&s.fe_space,nullptr,nullptr};
  // Change the 2 top parameters in case there is a convection field as Parameter
  if (with_convection_field)
  {
    nfespace = 3;
    fe_spaces[1] = la_stiff.get_fe_function(1)->GetFESpace2D();
    fe_spaces[2] = la_stiff.get_fe_function(2)->GetFESpace2D();
  }

  BoundCondFunct2D * boundary_conditions = fe_spaces[0]->GetBoundCondition();
  int N_Matrices = 1;
  double * rhs_entries = s.rhs.get_entries();

  BoundValueFunct2D * non_const_bound_value[1] {example.get_bd()[0]};

  //fetch stiffness matrix as block
  std::vector<std::shared_ptr<FEMatrix>> stiff_blocks = s.stiff_matrix.get_blocks_uniquely();
  TSquareMatrix2D * stiff_block[1]{reinterpret_cast<TSquareMatrix2D*>(stiff_blocks.at(0).get())};

  // Do the Assembling!

  //...stiffness matrix is second
  // reset right hand side and matrix to zero
  s.rhs.reset();
  stiff_block[0]->reset();
  Assemble2D(nfespace, fe_spaces, N_Matrices, stiff_block, 0, NULL, 1, &rhs_entries,
             fe_spaces, &boundary_conditions, non_const_bound_value, la_stiff);

  // If assemble_both is true, we also (re)assemble the mass matrix.
  if (assemble_both)
  {
    std::vector<std::shared_ptr<FEMatrix>> mass_blocks = s.mass_matrix.get_blocks_uniquely();
    TSquareMatrix2D * mass_block[1]{reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get())};

    mass_block[0]->reset();
    Assemble2D(nfespace, fe_spaces, N_Matrices, mass_block, 0, NULL, 0, NULL,
               NULL, &boundary_conditions, non_const_bound_value, la_mass);
  }
}








/* *********** BELOW THIS LINE USER SPECIFIC CODE **************/


// these 4 methods have been written for the case of nonlinear TCD
// they are not used at the moment (14.12.2016)
//void Time_CD2D::assemble_rhs_vector
//(const TFEVectFunct2D* convection_field)
//{
//  LocalAssembling2D_type stiff_rhs = LocalAssembling2D_type::TCD2D;
//
//  for(auto &s : this->systems)
//  {
//
//    TFEFunction2D * pointer_to_function[3] = {&s.fe_function, nullptr, nullptr};
//    // create two local assembling object (second one will only be needed in SUPG case)
//    LocalAssembling2D la_a_rhs(stiff_rhs, pointer_to_function,
//                               this->example.get_coeffs());
//
//
//    if (convection_field)   // this is in case where rhs depends on convection field
//    {// =========================== HERE CODE FOR CONVECTION FIELD
//      cout << "J'AI DETECTE LA PRESENCE D'UN CONVECTION FIELD" << endl;
//      // step 1 : interpolate the given convection_field to our fe space
////      const TFESpace2D& space = s.fe_space;
////      size_t n_dofs = space.GetN_DegreesOfFreedom();
////      std::string name("interpolated velo space");
////      std::string description("interpolated velo space");
////      std::vector<double> interp_funct_values(n_dofs,0.0);
//
////      // set up an interpolator object  (ptr will be shared later)
////      const TFESpace2D* into_space = &s.fe_space;
////      FEFunctionInterpolator interpolator(into_space);
////
////      // length of the values array of the interpolated velo must equal length of the
////      // concentration fe function
////      size_t length_interpolated = s.fe_function.GetLength();
////
////      std::vector<double> entries_velo_x(length_interpolated, 0.0);
////      std::vector<double> entries_velo_y(length_interpolated, 0.0);
//
//      // this awful call is due to the way a TFEVectFunct2D creates new dynamically
//      // allocated TFEFunction2D objects
//      TFEFunction2D* convection_x = convection_field->GetComponent(0);
//      TFEFunction2D* convection_y = convection_field->GetComponent(1);
//
////      TFEFunction2D interpolated_convection_x =
////          interpolator.interpolate(*convection_x, entries_velo_x);
////
////      TFEFunction2D interpolated_convection_y =
////          interpolator.interpolate(*convection_y, entries_velo_y);
////
////      delete convection_x; // call to GetComponent dynamically created fe functs
////      delete convection_y;
//
//      // step 2 - set all the 'parameter'-related values in la_a_rhs accordingly
//
//      // set up the input...
//      std::vector<int> beginParameter = {0};
//
//      //fill up the new fe function array
//      //pointer_to_function[0] = &s.fe_function;
//      pointer_to_function[1] = convection_x;
//      pointer_to_function[2] = convection_y;
//
//      std::vector<int> feValueFctIndex = {1,2}; // to produce first fe value use fe function 1,
//      // for second fe value use function 2
//      std::vector<MultiIndex2D> feValueMultiIndex = {D00, D00}; // to produce first fe value use 0th derivative,
//      // for second fe value as well
//      int N_parameters = 2; // two parameters...
//      int N_feValues = 2;   //..both of which stem from the evaluation of fe fcts
//      int N_paramFct = 1;   // dealing with them is performed by 1 ParamFct
//
//      // chose the parameter function ("in-out function") which shears away
//      // the first to "in" values (x,y) and passes only u_x and u_y
//      std::vector<ParamFct*> parameterFct = {NSParamsVelo};
//
//      // ...and call the corresponding setters
//      la_a_rhs.setBeginParameter(beginParameter);
//      la_a_rhs.setFeFunctions2D(pointer_to_function); //reset - now velo comp included
//      la_a_rhs.setFeValueFctIndex(feValueFctIndex);
//      la_a_rhs.setFeValueMultiIndex(feValueMultiIndex);
//      la_a_rhs.setN_Parameters(N_parameters);
//      la_a_rhs.setN_FeValues(N_feValues);
//      la_a_rhs.setN_ParamFct(N_paramFct);
//      la_a_rhs.setParameterFct(parameterFct);
//      //...this should do the trick
//      //===============================================================END CODE
//    }
//    else
//    {
//      cout << " JE SUIS ICI" << endl;
//      //      TFEFunction2D * pointer_to_function = &s.fe_function;
//    }
//
//
//    if(TDatabase::ParamDB->DISCTYPE == SUPG)
//    {
//      // In the SUPG case:
//      // M = (u,v) + \tau (u,b.grad v)
//      LocalAssembling2D_type mass_supg = LocalAssembling2D_type::TCD2D_Mass;
//
//      LocalAssembling2D la_m_supg(mass_supg, pointer_to_function,
//                                  this->example.get_coeffs());
//
//      //call assembling, including mass matrix (SUPG!)
//      call_assembling_routine(s, la_a_rhs, la_m_supg , true);
//    }
//    else
//    {//call assembling, ignoring mass matrix part (third argument is not relevant)
//      call_assembling_routine(s, la_a_rhs, la_a_rhs , false);
//    }
//  }
//
//
//  // End of the pure assembling operation
//
//
//  // here the modifications due to time discretization begin
//  if ( TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION == 2 )
//  {
//    do_algebraic_flux_correction();
//    return; // modifications due to time discretization are per-
//    // formed inside the afc scheme, so step out here!
//  }
//
//
//  // Start of the "scaling" of the right hand side
//
//
//  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
//  // preparing rhs
//  System_per_grid& s = this->systems.front();
//
//  if(TDatabase::TimeDB->THETA4)
//  {
//    // scale by time step length and theta4 (only active dofs)
//    s.rhs.scaleActive(tau*TDatabase::TimeDB->THETA4);
//    // add old right hand side scaled by time step length and theta3 (only
//    // active dofs)
//    if(TDatabase::TimeDB->THETA3 != 0.)
//      s.rhs.addScaledActive((this->old_rhs), tau*TDatabase::TimeDB->THETA3);
//
//    // save old right hand side (only if THETA3 != 0)
//    if(TDatabase::TimeDB->THETA3)
//    {
//      this->old_rhs.addScaledActive(s.rhs, -1./(tau*TDatabase::TimeDB->THETA3));
//      this->old_rhs.scaleActive(-TDatabase::TimeDB->THETA3/TDatabase::TimeDB->THETA4);
//    }
//  }
//  else
//  {
//    if(TDatabase::TimeDB->TIME_DISC == 0)
//    {
//      ErrThrow("Forward Euler method is not supported. "
//          "Choose TDatabase::TimeDB->TIME_DISC as 1 (bw Euler)"
//          " or 2 (Crank-Nicoloson)");
//    }
//  }
//
//  for(auto &s : this->systems)
//  {
//    // rhs += M*uold
//    s.mass_matrix.apply_scaled_add_actives(s.solution, s.rhs, 1.0);
//    // rhs -= tau*theta2*A_old*uold
//    s.rhs.addScaledActive(s.old_Au, -tau*TDatabase::TimeDB->THETA2);
//  }
//  //this->systems[0].rhs.copy_nonactive(this->systems[0].solution);
//  systems[0].solution.copy_nonactive(systems[0].rhs);
//}
//
//void Time_CD2D::assemble_stiffness_matrix_alone()
//{
//  LocalAssembling2D_type stiff_type = LocalAssembling2D_type::TCD2D;
//
//  for(auto &s : this->systems)
//  {
//    TFEFunction2D * pointer_to_function = &s.fe_function;
//
//
//    // create two local assembling object (second one will only be needed in SUPG case)
//    LocalAssembling2D la(stiff_type, &pointer_to_function,
//                               this->example.get_coeffs());
//
//    la.setAssembleParam(LocalMatrixA_alone);
//
//    // Assemble mass matrix, stiffness matrix and rhs
//    //...variables which are the same for both
//    const TFESpace2D * fe_space = &s.fe_space;
//    BoundCondFunct2D * boundary_conditions = fe_space->GetBoundCondition();
//    int N_Matrices = 1;
//
//    BoundValueFunct2D * non_const_bound_value[1] {example.get_bd()[0]};
//
//    //fetch stiffness matrix as block
//    std::vector<std::shared_ptr<FEMatrix>> stiff_blocks = s.stiff_matrix.get_blocks_uniquely();
//    TSquareMatrix2D * stiff_block[1]{reinterpret_cast<TSquareMatrix2D*>(stiff_blocks.at(0).get())};
//
//    // Do the Assembling!
//
//    // reset matrix to zero
//    stiff_block[0]->reset();
//
//    Assemble2D(1, &fe_space, N_Matrices, stiff_block, 0, NULL, 0, NULL,
//               &fe_space, &boundary_conditions, non_const_bound_value, la);
//  }
//}
//
//void Time_CD2D::assemble_stiffness_matrix_alone_with_convection
//(const TFEVectFunct2D* convection_field)
//{
//  LocalAssembling2D_type stiff_type = LocalAssembling2D_type::TCD2D;
//
//  for(auto &s : this->systems)
//  {
//    TFEFunction2D * pointer_to_function[3] = {&s.fe_function, nullptr, nullptr};
//
//
//    // create two local assembling object (second one will only be needed in SUPG case)
//    LocalAssembling2D la(stiff_type, pointer_to_function,
//                         this->example.get_coeffs());
//
//    la.setAssembleParam(LocalMatrixA_alone);
//
//    if (convection_field)
//    {// =========================== HERE CODE FOR CONVECTION FIELD
//      cout << "J'AI DETECTE LA PRESENCE D'UN CONVECTION FIELD" << endl;
//      // step 1 : interpolate the given convection_field to our fe space
////      const TFESpace2D& space = s.fe_space;
////      size_t n_dofs = space.GetN_DegreesOfFreedom();
////      std::string name("interpolated velo space");
////      std::string description("interpolated velo space");
////      std::vector<double> interp_funct_values(n_dofs,0.0);
//
////      // set up an interpolator object  (ptr will be shared later)
////      const TFESpace2D* into_space = &s.fe_space;
////      FEFunctionInterpolator interpolator(into_space);
////
////      // length of the values array of the interpolated velo must equal length of the
////      // concentration fe function
////      size_t length_interpolated = s.fe_function.GetLength();
////
////      std::vector<double> entries_velo_x(length_interpolated, 0.0);
////      std::vector<double> entries_velo_y(length_interpolated, 0.0);
//
//      // this awful call is due to the way a TFEVectFunct2D creates new dynamically
//      // allocated TFEFunction2D objects
//      TFEFunction2D* convection_x = convection_field->GetComponent(0);
//      TFEFunction2D* convection_y = convection_field->GetComponent(1);
//
////      TFEFunction2D interpolated_convection_x =
////          interpolator.interpolate(*convection_x, entries_velo_x);
////
////      TFEFunction2D interpolated_convection_y =
////          interpolator.interpolate(*convection_y, entries_velo_y);
////
////      delete convection_x; // call to GetComponent dynamically created fe functs
////      delete convection_y;
//
//      // step 2 - set all the 'parameter'-related values in la_a_rhs accordingly
//
//      // set up the input...
//      std::vector<int> beginParameter = {0};
//
//      //fill up the new fe function array
//      //pointer_to_function[0] = &s.fe_function;
//      pointer_to_function[1] = convection_x;
//      pointer_to_function[2] = convection_y;
//
//      std::vector<int> feValueFctIndex = {1,2}; // to produce first fe value use fe function 1,
//      // for second fe value use function 2
//      std::vector<MultiIndex2D> feValueMultiIndex = {D00, D00}; // to produce first fe value use 0th derivative,
//      // for second fe value as well
//      int N_parameters = 2; // two parameters...
//      int N_feValues = 2;   //..both of which stem from the evaluation of fe fcts
//      int N_paramFct = 1;   // dealing with them is performed by 1 ParamFct
//
//      // chose the parameter function ("in-out function") which shears away
//      // the first to "in" values (x,y) and passes only u_x and u_y
//      std::vector<ParamFct*> parameterFct = {NSParamsVelo};
//
//      // ...and call the corresponding setters
//      la.setBeginParameter(beginParameter);
//      la.setFeFunctions2D(pointer_to_function); //reset - now velo comp included
//      la.setFeValueFctIndex(feValueFctIndex);
//      la.setFeValueMultiIndex(feValueMultiIndex);
//      la.setN_Parameters(N_parameters);
//      la.setN_FeValues(N_feValues);
//      la.setN_ParamFct(N_paramFct);
//      la.setParameterFct(parameterFct);
//      //...this should do the trick
//      //===============================================================END CODE
//    }
//    else
//    {
//      cout << " JE SUIS ICI" << endl;
//      //      TFEFunction2D * pointer_to_function = &s.fe_function;
//    }
//
//
//    // Assemble mass matrix, stiffness matrix and rhs
//    //...variables which are the same for both
//    const TFESpace2D * fe_space = &s.fe_space;
//    BoundCondFunct2D * boundary_conditions = fe_space->GetBoundCondition();
//    int N_Matrices = 1;
//
//    BoundValueFunct2D * non_const_bound_value[1] {example.get_bd()[0]};
//
//    //fetch stiffness matrix as block
//    std::vector<std::shared_ptr<FEMatrix>> stiff_blocks = s.stiff_matrix.get_blocks_uniquely();
//    TSquareMatrix2D * stiff_block[1]{reinterpret_cast<TSquareMatrix2D*>(stiff_blocks.at(0).get())};
//
//    // Do the Assembling!
//
//    // reset matrix to zero
//    stiff_block[0]->reset();
//
//    Assemble2D(1, &fe_space, N_Matrices, stiff_block, 0, NULL, 0, NULL,
//               &fe_space, &boundary_conditions, non_const_bound_value, la);
//  }
//}
//
//void Time_CD2D::scale_stiffness_matrix()
//{
//  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
//  for(auto &s : this->systems)
//  {
//    // preparing the left hand side, i.e., the system matrix
//    // stiffness matrix is scaled by tau*THETA1, after solving
//    // the matrix needs to be descaled if the coeffs does not depends
//    // on time
//
//    //scale  stiffness matrix...
//    const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
//    s.stiff_matrix.scale_blocks_actives(tau*TDatabase::TimeDB->THETA1, cell_positions);
//    // ...and add the mass matrix
//    const FEMatrix& mass_block = *s.mass_matrix.get_blocks().at(0).get();
//    s.stiff_matrix.add_matrix_actives(mass_block, 1.0, {{0,0}}, {false});
//  }
//}



// these 2 methods are currently used and up to date
void Time_CD2D::assemble_initial_time_with_convection
(const TFEVectFunct2D* convection_field)
{
  LocalAssembling2D_type mass = LocalAssembling2D_type::TCD2D_Mass;
  LocalAssembling2D_type stiff_rhs = LocalAssembling2D_type::TCD2D;

  for(auto &s : this->systems)
  {
    bool with_convection_field = false;
    // assemble mass matrix, stiffness matrix and rhs
    TFEFunction2D* fe_funct[3] = {&s.fe_function,
                                  nullptr, nullptr}; //wrap up as '**'

    LocalAssembling2D la_mass(mass, fe_funct,
                              this->example.get_coeffs());
    LocalAssembling2D la_a_rhs(stiff_rhs, fe_funct,
                               this->example.get_coeffs());



    if (convection_field)    // assembles initial RHS and Stiffness with convection field
    {
      with_convection_field = true;
      // HERE IS THE CODE TO SET UP THE CONVECTION FIELD FOR LOCAL ASSEMBLING OBJECT

      // this awful call is due to the way a TFEVectFunct2D creates new dynamically
      // allocated TFEFunction2D objects
      TFEFunction2D* convection_x = convection_field->GetComponent(0);
      TFEFunction2D* convection_y = convection_field->GetComponent(1);

      // we assume convection_x and convection_y are not too exotic,
      // and have the same space...
      int Ndof_convection     = convection_x->GetFESpace2D()->GetN_DegreesOfFreedom();
      int Ndof_thisfefunction = s.fe_space.GetN_DegreesOfFreedom();


      // step 1: check if the spaces "this->systems.fe_space" and
      // "convection_field->GetFESpace2D()" are the same. If yes, no interpolation needed
      // otherwise, do the interpolation
      if (Ndof_convection == Ndof_thisfefunction) //this is a simple check condition...
      {
        Output::info<1>("Time_CD2D", "The spaces of the convection field and the "
            "scalar field are the same ==> There will be no interpolation.");
//        Output::print<3>("Degres of freedoms of scalar field = " , Ndof_thisfefunction);
//        Output::print<3>("Degres of freedoms of velo x compo = " , Ndof_convection);

        //fill up the new fe function array
        fe_funct[1] = convection_x;
        fe_funct[2] = convection_y;

      }
      else  // do the interpolation
      {
        Output::warn<1>("Time_CD2D", "The spaces of the convection field and the "
            "scalar field are not the same ==> adapting assemble method");
        fe_funct[1] = convection_x;
        fe_funct[2] = convection_y;
//        // set up an interpolator object  (ptr will be shared later)
//        const TFESpace2D* into_space = &s.fe_space;
//        FEFunctionInterpolator interpolator(into_space);
//
//        // length of the values array of the interpolated velo must equal length of the
//        // concentration fe function
//        size_t length_interpolated = s.fe_function.GetLength();
//
//        std::vector<double> temporary_x(length_interpolated, 0.0);
//        std::vector<double> temporary_y(length_interpolated, 0.0);
//
//        this->entries_velo_x = temporary_x;
//        this->entries_velo_y = temporary_y;
//
//        TFEFunction2D interpolated_convection_x =
//            interpolator.interpolate(*convection_x, this->entries_velo_x);
//
//        TFEFunction2D interpolated_convection_y =
//            interpolator.interpolate(*convection_y, this->entries_velo_y);
//
//        //fill up the new fe function array
//        fe_funct[1] = &interpolated_convection_x;
//        fe_funct[2] = &interpolated_convection_y;
//      delete convection_x; // call to GetComponent dynamically created fe functs
//      delete convection_y;
      }

      // step 2 - set all the 'parameter'-related values in la_a_rhs accordingly
      // set up the input...
      std::vector<int> beginParameter = {0};
      std::vector<int> feValueFctIndex = {1,2}; // to produce first fe value use fe function 1,
      // for second fe value use function 2
      std::vector<MultiIndex2D> feValueMultiIndex = {D00, D00}; // to produce first fe value use 0th derivative,
      // for second fe value as well
      int N_parameters = 2; // two parameters...
      int N_feValues = 2;   //..both of which stem from the evaluation of fe fcts
      int N_paramFct = 1;   // dealing with them is performed by 1 ParamFct

      // chose the parameter function ("in-out function") which shears away
      // the first to "in" values (x,y) and passes only u_x and u_y
      std::vector<ParamFct*> parameterFct = {NSParamsVelo};

      // ...and call the corresponding setters
      la_a_rhs.setBeginParameter(beginParameter);
      la_a_rhs.setFeFunctions2D(fe_funct); //reset - now velo comp included
      la_a_rhs.setFeValueFctIndex(feValueFctIndex);
      la_a_rhs.setFeValueMultiIndex(feValueMultiIndex);
      la_a_rhs.setN_Parameters(N_parameters);
      la_a_rhs.setN_FeValues(N_feValues);
      la_a_rhs.setN_ParamFct(N_paramFct);
      la_a_rhs.setParameterFct(parameterFct);
      //...this should do the trick
      //===============================================================END CODE
    }


    // do the actual assembling
    call_assembling_routine(s, la_a_rhs, la_mass , true, with_convection_field);

    // apply local projection stabilization method on stiffness matrix only!
    if(TDatabase::ParamDB->DISCTYPE==LOCAL_PROJECTION
        && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
    {
      if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
      {
        //fetch stiffness matrix as block
        std::vector<std::shared_ptr<FEMatrix>> stiff_blocks =
            s.stiff_matrix.get_blocks_uniquely();
        auto stiff_block = stiff_blocks.at(0).get();
        //call ultra local projection
        UltraLocalProjection((void *)stiff_block, false);
      }
      else
      {
        ErrThrow("LP_FULL_GRADIENT needs to be 1 to use LOCAL_PROJECTION");
      }

    }

    s.stiff_matrix.apply(s.solution,s.old_Au); //put up initial old_Au
  }

  System_per_grid& s = systems.front();
  old_rhs = s.rhs; //put up initial old_rhs
}

void Time_CD2D::assemble_with_convection
(const TFEVectFunct2D* convection_field)
{
  LocalAssembling2D_type stiff_rhs = LocalAssembling2D_type::TCD2D;

  for(auto &s : this->systems)
  {
    bool with_convection_field = false;
    TFEFunction2D* fe_funct[3] = {&s.fe_function,
                                  nullptr, nullptr};

    // create two local assembling object (second one will only be needed in SUPG case)
    LocalAssembling2D la_a_rhs(stiff_rhs, fe_funct,
                               this->example.get_coeffs());

    if (convection_field)    // assembles initial RHS and Stiffness with convection field
    {
      with_convection_field = true;
      // HERE IS THE CODE TO SET UP THE CONVECTION FIELD FOR LOCAL ASSEMBLING OBJECT

      // this awful call is due to the way a TFEVectFunct2D creates new dynamically
      // allocated TFEFunction2D objects
      TFEFunction2D* convection_x = convection_field->GetComponent(0);
      TFEFunction2D* convection_y = convection_field->GetComponent(1);

      // we assume convection_x and convection_y are not too exotic,
      // and have the same space...
      int Ndof_convection     = convection_x->GetFESpace2D()->GetN_DegreesOfFreedom();
      int Ndof_thisfefunction = s.fe_space.GetN_DegreesOfFreedom();


      // step 1: check if the spaces "this->systems.fe_space" and
      // "convection_field->GetFESpace2D()" are the same. If yes, no interpolation needed
      // otherwise, do the interpolation
      if (Ndof_convection == Ndof_thisfefunction) //this is a simple check condition...
      {
        Output::info<1>("Time_CD2D", "The spaces of the convection field and the "
            "scalar field are the same ==> There will be no interpolation.");
//        Output::print<3>("Degres of freedoms of scalar field = " , Ndof_thisfefunction);
//        Output::print<3>("Degres of freedoms of velo x compo = " , Ndof_convection);

        //fill up the new fe function array
        fe_funct[1] = convection_x;
        fe_funct[2] = convection_y;

      }
      else  // do the interpolation
      {
        Output::warn<1>("Time_CD2D", "The spaces of the convection field and the "
            "scalar field are not the same ==> adapting assemble method");

        fe_funct[1] = convection_x;
        fe_funct[2] = convection_y;
//        // set up an interpolator object  (ptr will be shared later)
//        const TFESpace2D* into_space = &s.fe_space;
//        FEFunctionInterpolator interpolator(into_space);
//
//        // length of the values array of the interpolated velo must equal length of the
//        // concentration fe function
//        size_t length_interpolated = s.fe_function.GetLength();
//
//        std::vector<double> temporary_x(length_interpolated, 0.0);
//        std::vector<double> temporary_y(length_interpolated, 0.0);
//
//        this->entries_velo_x = temporary_x;
//        this->entries_velo_y = temporary_y;
//
//        TFEFunction2D interpolated_convection_x =
//            interpolator.interpolate(*convection_x, this->entries_velo_x);
//
//        TFEFunction2D interpolated_convection_y =
//            interpolator.interpolate(*convection_y, this->entries_velo_y);
//
//        //fill up the new fe function array
//        fe_funct[1] = &interpolated_convection_x;
//        fe_funct[2] = &interpolated_convection_y;
//
////      delete convection_x; // call to GetComponent dynamically created fe functs
////      delete convection_y;
      }

      // step 2 - set all the 'parameter'-related values in la_a_rhs accordingly
      // set up the input...
      std::vector<int> beginParameter = {0};
      std::vector<int> feValueFctIndex = {1,2}; // to produce first fe value use fe function 1,
      // for second fe value use function 2
      std::vector<MultiIndex2D> feValueMultiIndex = {D00, D00}; // to produce first fe value use 0th derivative,
      // for second fe value as well
      int N_parameters = 2; // two parameters...
      int N_feValues = 2;   //..both of which stem from the evaluation of fe fcts
      int N_paramFct = 1;   // dealing with them is performed by 1 ParamFct

      // chose the parameter function ("in-out function") which shears away
      // the first to "in" values (x,y) and passes only u_x and u_y
      std::vector<ParamFct*> parameterFct = {NSParamsVelo};

      // ...and call the corresponding setters
      la_a_rhs.setBeginParameter(beginParameter);
      la_a_rhs.setFeFunctions2D(fe_funct); //reset - now velo comp included
      la_a_rhs.setFeValueFctIndex(feValueFctIndex);
      la_a_rhs.setFeValueMultiIndex(feValueMultiIndex);
      la_a_rhs.setN_Parameters(N_parameters);
      la_a_rhs.setN_FeValues(N_feValues);
      la_a_rhs.setN_ParamFct(N_paramFct);
      la_a_rhs.setParameterFct(parameterFct);
      //...this should do the trick
      //===============================================================END CODE
    }

    if(TDatabase::ParamDB->DISCTYPE == SUPG)
    {
      // In the SUPG case:
      // M = (u,v) + \tau (u,b.grad v)
      LocalAssembling2D_type mass_supg = LocalAssembling2D_type::TCD2D_Mass;;

      LocalAssembling2D la_m_supg(mass_supg, fe_funct,
                               this->example.get_coeffs());

      //call assembling, including mass matrix (SUPG!)
      call_assembling_routine(s, la_a_rhs, la_m_supg , true,with_convection_field);
    }
    else
    {//call assembling, ignoring mass matrix part (third argument is not relevant)
      call_assembling_routine(s, la_a_rhs, la_a_rhs , false,with_convection_field);
    }
  }

  // here the modifications due to time discretization begin
  if ( !db["algebraic_flux_correction"].is("none")  )
  {
    do_algebraic_flux_correction();
    return; // modifications due to time discretization are per-
            // formed inside the afc scheme, so step out here!
  }

  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  // preparing rhs
  System_per_grid& s = this->systems.front();

  if(TDatabase::TimeDB->THETA4)
  {
    // scale by time step length and theta4 (only active dofs)
    s.rhs.scaleActive(tau*TDatabase::TimeDB->THETA4);
    // add old right hand side scaled by time step length and theta3 (only
    // active dofs)
    if(TDatabase::TimeDB->THETA3 != 0.)
      s.rhs.addScaledActive((this->old_rhs), tau*TDatabase::TimeDB->THETA3);

    // save old right hand side (only if THETA3 != 0)
    if(TDatabase::TimeDB->THETA3)
    {
      this->old_rhs.addScaledActive(s.rhs, -1./(tau*TDatabase::TimeDB->THETA3));
      this->old_rhs.scaleActive(-TDatabase::TimeDB->THETA3/TDatabase::TimeDB->THETA4);
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

  for(auto &s : this->systems)
  {
    // rhs += M*uold
    s.mass_matrix.apply_scaled_add_actives(s.solution, s.rhs, 1.0);
    // rhs -= tau*theta2*A_old*uold
    s.rhs.addScaledActive(s.old_Au, -tau*TDatabase::TimeDB->THETA2);

    // preparing the left hand side, i.e., the system matrix
    // stiffness matrix is scaled by tau*THETA1, after solving
    // the matrix needs to be descaled if the coeffs does not depends
    // on time

    //scale  stiffness matrix...
    const std::vector<std::vector<size_t>> cell_positions = {{0,0}};
    s.stiff_matrix.scale_blocks_actives(tau*TDatabase::TimeDB->THETA1, cell_positions);
    // ...and add the mass matrix
    const FEMatrix& mass_block = *s.mass_matrix.get_blocks().at(0).get();
    s.stiff_matrix.add_matrix_actives(mass_block, 1.0, {{0,0}}, {false});
  }
  //this->systems[0].rhs.copy_nonactive(this->systems[0].solution);
  systems[0].solution.copy_nonactive(systems[0].rhs);
}

