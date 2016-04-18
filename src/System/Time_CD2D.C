#include <Time_CD2D.h>
#include <Database.h>
#include <MultiGrid2D.h>
#include <Output2D.h>
#include <OldSolver.h>
#include <DirectSolver.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <AlgebraicFluxCorrection.h>

#include <LocalAssembling2D.h>
#include <Assemble2D.h>
#include <LocalProjection.h>

#include <numeric>


/**************************************************************************** */
Time_CD2D::System_per_grid::System_per_grid(const Example_CD2D& example,
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
  stiff_matrix.scale_blocks_actives(1./(tau*theta_1), {{0,0}});
}

/**************************************************************************** */

void Time_CD2D::System_per_grid::update_old_Au()
{
  // the stiffness matrix must have been descaled (pure transport operator)
  stiff_matrix.apply(solution,old_Au);
}

/**************************************************************************** */
TSquareMatrix2D* Time_CD2D::System_per_grid::get_stiff_matrix_pointer()
{
  std::vector<std::shared_ptr<FEMatrix>> blocks =
      stiff_matrix.get_blocks_TERRIBLY_UNSAFE();
  return reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
}

/**************************************************************************** */
Time_CD2D::Time_CD2D(const TDomain& domain, int reference_id)
 : Time_CD2D(domain, Example_CD2D(), reference_id)
{
  
}

/**************************************************************************** */
Time_CD2D::Time_CD2D(const TDomain& domain, const Example_CD2D& ex,
                     int reference_id)
 : systems(), example(ex), multigrid(nullptr), errors(5, 0.0)
{
  this->set_parameters();
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  // create finite element space and function, a matrix, rhs, and solution
  this->systems.emplace_back(this->example, *coll);
  
  TFESpace2D& space = this->systems.front().fe_space;
  double hmin, hmax;
  coll->GetHminHmax(&hmin, &hmax);
  Output::print<1>("N_Cells    : ", setw(12), coll->GetN_Cells());
  Output::print<1>("h(min,max) : ", setw(12), hmin, " ", setw(12), hmax);
  Output::print<1>("dof        : ", setw(12), space.GetN_DegreesOfFreedom());
  Output::print<1>("active dof : ", setw(12), space.GetN_ActiveDegrees());
  
  // old right hand side
  old_rhs.copy_structure(this->systems[0].rhs);
  
  // interpolate the initial data
  this->systems.front().fe_function.Interpolate(example.get_initial_cond(0));
  
  // done with the constructor in case we're not using multigrid
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR != 5 
    || TDatabase::ParamDB->SOLVER_TYPE != 1)
    return;
  // else multigrid
  
  double *param = new double[2];
  param[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  param[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
  this->multigrid.reset(new TMultiGrid2D(1, 2, param));
  // number of refinement levels for the multigrid
  int LEVELS = TDatabase::ParamDB->LEVELS;
  if(LEVELS > domain.get_ref_level() + 1)
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
    TSquareMatrix2D* stiff_block = it->get_stiff_matrix_pointer();

    TMGLevel2D *multigrid_level = new TMGLevel2D(
      i, stiff_block, it->rhs.get_entries(),
      it->solution.get_entries(), 2, NULL);
    i++;
    this->multigrid->AddLevel(multigrid_level);
  }
}

/**************************************************************************** */
void Time_CD2D::set_parameters()
{
  if(TDatabase::ParamDB->EXAMPLE < 101)
  {
    ErrMsg("Example " << TDatabase::ParamDB->EXAMPLE 
           << "is not supported for time dependent problem");
    exit(1);
  }
  
  if(TDatabase::TimeDB->TIME_DISC == 0)
  {
    ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC 
          << " is not supported");
    throw("TIME_DISC: 0 is not supported. Chose 1 (backward Euler) or 2 (Crank-Nicolson).");
  }

  //////////////// Algebraic flux correction ////////////
  if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION != 0)
  {//some kind of afc enabled
    if(TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION != 2 )
    {
      TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION = 2;
      Output::print("Only kind of algebraic flux correction"
          " for TCD problems is Crank-Nicolson FEM-FCT (2).");
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

    // so far only direct solver
    if (TDatabase::ParamDB->SOLVER_TYPE != 2)
    {
      ErrThrow("With AFC only direct solver is working so far!");
    }

    //make sure that galerkin discretization is used
    if (TDatabase::ParamDB->DISCTYPE != 1)
    {//some other disctype than galerkin
      TDatabase::ParamDB->DISCTYPE = 1;
      Output::print("DISCTYPE changed to 1 (GALERKIN) because Algebraic Flux ",
                    "Correction is enabled.");
    }

    // when using afc, create system matrices as if all dofs were active
    TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
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

    call_assembling_routine(s, la_a_rhs, la_mass , true);

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
      call_assembling_routine(s, la_a_rhs, la_m_supg , true);
    }
    else
    {//call assembling, ignoring mass matrix part (third argument is not relevant)
      call_assembling_routine(s, la_a_rhs, la_a_rhs , false);
    }
  }
  
  // here the modifications due to time discretization begin
  if ( TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION == 2 )
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
    s.stiff_matrix.scale_blocks_actives(tau*TDatabase::TimeDB->THETA1, {{0,0}});
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
    if(!TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION)
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
  if(TDatabase::ParamDB->SOLVER_TYPE == 2) // use direct solver
  {
    /// @todo consider storing an object of DirectSolver in this class
    DirectSolver direct_solver(s.stiff_matrix, 
                               DirectSolver::DirectSolverTypes::umfpack);
    direct_solver.solve(s.rhs, s.solution);
  }
  else
  {
    // the one matrix stored in this BlockMatrix
    FEMatrix* mat = s.stiff_matrix.get_blocks_uniquely().at(0).get();
    // in order for the old method 'Solver' to work we need TSquareMatrix2D
    TSquareMatrix2D *sqMat[1] = { reinterpret_cast<TSquareMatrix2D*>(mat) };
    OldSolver((TSquareMatrix **)sqMat, NULL, s.rhs.get_entries(), 
              s.solution.get_entries(), MatVect_Scalar, Defect_Scalar, 
              this->multigrid.get(), s.solution.length(), 0);
  }
  
  t = GetTime() - t;
  Output::print<1>("time for solving: ",  t);
  Output::print<2>("solution ", sqrt(Ddot(s.solution.length(),
                                          s.solution.get_entries(),
                                          s.solution.get_entries())) );

}

/**************************************************************************** */
void Time_CD2D::output(int m, int& image)
{
  if(!TDatabase::ParamDB->WRITE_VTK && !TDatabase::ParamDB->MEASURE_ERRORS)
    return;
  
  TFEFunction2D & fe_function = this->systems.front().fe_function;
  fe_function.PrintMinMax();
  
  if(TDatabase::ParamDB->MEASURE_ERRORS)
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
    
    if(m==0)
      errors[4]= loc_e[0];
    
    if(errors[4] < loc_e[0])
      errors[4] = loc_e[0];
    Output::print<1>("  Linfty(0,T;L2) ", errors[4]);
  }

  if((m==1) || (m%TDatabase::TimeDB->STEPS_PER_IMAGE == 0))
  {
    if(TDatabase::ParamDB->WRITE_VTK)
    {
      TOutput2D Output(1, 1, 0, 0, NULL);
      Output.AddFEFunction(&fe_function);
      std::string filename(TDatabase::ParamDB->OUTPUTDIR);
      filename += "/" + std::string(TDatabase::ParamDB->BASENAME);
      if(image<10) filename += ".0000";
      else if(image<100) filename += ".000";
      else if(image<1000) filename += ".00";
      else if(image<10000) filename += ".0";
      else filename += ".";
      filename += std::to_string(image) + ".vtk";
      Output.WriteVtk(filename.c_str());
      image++;
    }
  }
}

/**************************************************************************** */

void Time_CD2D::do_algebraic_flux_correction()
{
  //determine which kind of afc to use
  switch (TDatabase::ParamDB->ALGEBRAIC_FLUX_CORRECTION)
  {
    case 2: //Crank Nicolson FEM-FCT
    { //TODO implement for multigrid!
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
      switch(TDatabase::ParamDB->FEM_FCT_PRELIMITING)
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
      break;
    }
    default:
    {
      ErrThrow("The chosen ALGEBRAIC_FLUX_CORRECTION scheme is unknown to class Time_CD2D.");
    }
  }
}

void Time_CD2D::call_assembling_routine(
    System_per_grid& s,
    LocalAssembling2D& la_stiff, LocalAssembling2D& la_mass,
    bool assemble_both)
{

  // Assemble mass matrix, stiffness matrix and rhs
  //...variables which are the same for both
  const TFESpace2D * fe_space = &s.fe_space;
  BoundCondFunct2D * boundary_conditions = fe_space->GetBoundCondition();
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
  Assemble2D(1, &fe_space, N_Matrices, stiff_block, 0, NULL, 1, &rhs_entries,
             &fe_space, &boundary_conditions, non_const_bound_value, la_stiff);

  // If assemble_both is true, we also (re)assemble the mass matrix.
  if (assemble_both)
  {
    std::vector<std::shared_ptr<FEMatrix>> mass_blocks = s.mass_matrix.get_blocks_uniquely();
    TSquareMatrix2D * mass_block[1]{reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get())};

    mass_block[0]->reset();
    Assemble2D(1, &fe_space, N_Matrices, mass_block, 0, NULL, 0, NULL,
               NULL, &boundary_conditions, non_const_bound_value, la_mass);
  }
}
