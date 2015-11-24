#include <Time_CD2D.h>
#include <Database.h>
#include <MultiGrid2D.h>
#include <Output2D.h>
#include <Solver.h>
#include <LinAlg.h>
#include <MainUtilities.h>



/**************************************************************************** */
Time_CD2D::System_per_grid::System_per_grid(const Example_CD2D& example,
                                            TCollection& coll)
 : fe_space(&coll, (char*)"space", (char*)"time_cd2d space", example.get_bc(0),
            TDatabase::ParamDB->ANSATZ_ORDER, nullptr),
   Stiff_matrix(this->fe_space, example.get_bd(0)),
   Mass_Matrix(this->fe_space, example.get_bd(0), true),
   rhs(this->Stiff_matrix, true),
   solution(this->Stiff_matrix, false),
   fe_function(&this->fe_space, (char*)"u", (char*)"u", 
               this->solution.get_entries(), this->solution.length())
{
}

/**************************************************************************** */
Time_CD2D::Time_CD2D(const TDomain& domain, int reference_id)
 : Time_CD2D(domain, *(new Example_CD2D()), reference_id)
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
  
  // done with the conrtuctor in case we're not using multigrid
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
    TMGLevel2D *multigrid_level = new TMGLevel2D(
      i, it->Stiff_matrix.get_matrix(), it->rhs.get_entries(), 
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
           << "does not supported for time dependent problem");
    exit(1);
  }
  
  if(TDatabase::TimeDB->TIME_DISC == 0)
  {
    ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC 
          << " does not supported");
    throw("TIME_DISC: 0 does not supported");
  }  
}

/**************************************************************************** */
void Time_CD2D::assemble_initial_time()
{
  LocalAssembling2D_type m_rhs = TCD2D_Mass_Rhs_Galerkin;
  LocalAssembling2D_type stiff_rhs = TCD2D_Stiff_Rhs_Galerkin;
  
  if(TDatabase::ParamDB->DISCTYPE==SUPG)
  {
    m_rhs = TCD2D_Mass_Rhs_SUPG;
    stiff_rhs = TCD2D_Stiff_Rhs_SUPG;
    ErrMsg("DISCTYPE " <<TDatabase::ParamDB->DISCTYPE << " does not supported yet"
          << "One have to implement differently in the local assemble function");
    exit(1);
  }
  
  for(auto &s : this->systems)
  {
    TFEFunction2D * pointer_to_function = &s.fe_function;
    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling2D la_m_rhs(m_rhs, &pointer_to_function,
                               this->example.get_coeffs());
    LocalAssembling2D la_a_rhs(stiff_rhs, &pointer_to_function,
                               this->example.get_coeffs());
    // assemble mass matrix, stiffness matrix and rhs
    s.Mass_Matrix.Assemble(la_m_rhs, s.solution, s.rhs);
    s.Stiff_matrix.Assemble(la_a_rhs, s.solution, s.rhs);
  } 
  // copy the current right hand side vector to the old_rhs 
  this->old_rhs = this->systems.front().rhs;
}

/**************************************************************************** */
void Time_CD2D::assemble()
{
  LocalAssembling2D_type stiff_rhs = TCD2D_Stiff_Rhs_Galerkin;
  if(TDatabase::ParamDB->DISCTYPE==SUPG)
    stiff_rhs = TCD2D_Stiff_Rhs_SUPG;
  
  for(auto &s : this->systems)
  {
    TFEFunction2D * pointer_to_function = &s.fe_function;
    // create a local assembling object
    LocalAssembling2D la_a_rhs(stiff_rhs, &pointer_to_function,
                               this->example.get_coeffs());
    s.Stiff_matrix.Assemble(la_a_rhs,s.solution,s.rhs);
  }
  
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  // preparing rhs 
  System_per_grid& s = this->systems.front();
  if(TDatabase::TimeDB->THETA4)
  {
    // scale by time step length and theta4 (only active dofs)
    s.rhs.scale(tau*TDatabase::TimeDB->THETA4, 0);
    // add old right hand side scaled by time step length and theta3 (only 
    // active dofs)
    if(TDatabase::TimeDB->THETA3 != 0.)
      s.rhs.add_scaled((this->old_rhs), tau*TDatabase::TimeDB->THETA3);	

    // save old right hand side (only if THETA3 != 0)
    if(TDatabase::TimeDB->THETA3)
    {
      this->old_rhs.add_scaled(s.rhs, -1./(tau*TDatabase::TimeDB->THETA3));
      this->old_rhs.scale(-TDatabase::TimeDB->THETA3/TDatabase::TimeDB->THETA4);
    }    
  }
  else
  {
    if(TDatabase::TimeDB->TIME_DISC == 0)
    {
      ErrMsg("Bacward Euler method does not implemented");
      throw("Bacward Euler method does not implemented");
    }
  }
  
  for(auto &s : this->systems)
  {
    // rhs += M*uold
    s.Mass_Matrix.apply_scaled_add(s.solution.get_entries(),
                                   s.rhs.get_entries(), 1.0);
    // rhs += tau*theta3*A*uold
    s.Stiff_matrix.apply_scaled_add(s.solution.get_entries(),
                                    s.rhs.get_entries(),
                                    -tau*TDatabase::TimeDB->THETA2);
    
    // preparing the left hand side, i.e., the system matrix
    // stiffness matrix is scaled by tau*THETA1, after solving 
    // the matrix needs to be descaled if the coeffs does not depends
    // on time
    s.Stiff_matrix.scale_active(tau*TDatabase::TimeDB->THETA1);
    s.Stiff_matrix.add_scaled_active(s.Mass_Matrix, 1.0);
  }
  
  this->systems[0].rhs.copy_nonactive(this->systems[0].solution);  
}

/**************************************************************************** */
void Time_CD2D::solve()
{
  double t = GetTime();
  System_per_grid& s = this->systems.front();
  TSquareMatrix2D *sqMat[1] = {s.Stiff_matrix.get_matrix()};
  Solver((TSquareMatrix **)sqMat, NULL, s.rhs.get_entries(), 
         s.solution.get_entries(), MatVect_Scalar, Defect_Scalar, 
         this->multigrid.get(), s.solution.length(), 0);
  
  t = GetTime() - t;
  Output::print<1>("time for solving: ",  t);
  Output::print<2>("solution ", sqrt(Ddot(s.solution.length(),
                                          s.solution.get_entries(),
                                          s.solution.get_entries())) );

  // descale the stiffness matrix  
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  for(auto &s : this->systems)
  {
    s.Stiff_matrix.add_scaled_active(s.Mass_Matrix, -1.0);
    s.Stiff_matrix.scale_active(1./(tau*TDatabase::TimeDB->THETA1));
  }
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
    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
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
