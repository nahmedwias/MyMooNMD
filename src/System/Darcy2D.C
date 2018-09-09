/** ************************************************************************ 
* @brief     source file for Darcy2D
* @author    Ulrich Wilbrandt,
* @date      15.03.15
 ************************************************************************  */
#include <Database.h>
#include <Darcy2D.h>
#include <Assemble2D.h>
#include <MainUtilities.h>
#include <AuxParam2D.h>


ParameterDatabase get_default_Darcy2D_parameters()
{
  Output::print<5>("creating a default Darcy2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default Darcy2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Darcy2D parameter database");
  
  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  return db;
}


/** ************************************************************************ */
Darcy2D::System_per_grid::System_per_grid(const Example_Darcy2D& example,
                                          TCollection& coll)
 : velocity_space(new TFESpace2D(&coll, "u", "Darcy velocity",
                                 example.get_bc(0), 
                                 TDatabase::ParamDB->VELOCITY_SPACE, nullptr)),
   pressure_space(new TFESpace2D(&coll, "p", "Darcy pressure",
                                 example.get_bc(1),
                                 TDatabase::ParamDB->PRESSURE_SPACE, nullptr))

{
  /// @note witout this call, we could switch save the FESpaces as const, ie, 
  /// we could use std::shared_ptr<const TFESpace2D>.
  velocity_space->SetAsDGSpace();
  pressure_space->SetAsDGSpace();
  
  matrix = BlockFEMatrix::Darcy2D(*this->velocity_space, *this->pressure_space);

  rhs = BlockVector(this->matrix, true);
  solution = BlockVector(this->matrix, false);

  u = TFEFunction2D(this->velocity_space.get(), "u", "u",
                    this->solution.block(0), this->solution.length(0));
  p = TFEFunction2D(this->pressure_space.get(), "p", "p",
                    this->solution.block(1), this->solution.length(1));
}

/** ************************************************************************ */
Darcy2D::Darcy2D(const TDomain& domain, const ParameterDatabase& param_db, 
                 int reference_id)
 : Darcy2D(domain, param_db, Example_Darcy2D(param_db), reference_id)
{
  // note that the way we construct the example above will produce a memory 
  // leak, but that class is small.
}

/** ************************************************************************ */
Darcy2D::Darcy2D(const TDomain& domain, const ParameterDatabase& param_db,
                 const Example_Darcy2D ex, int reference_id)
 : systems(), example(ex), db(get_default_Darcy2D_parameters()),
   outputWriter(param_db), solver(param_db), errors()
{
  // get the parameters to control the behavior of this class
  this->db.merge(param_db, false);
  // make sure all parameters in the database are set consistently
  this->set_parameters();
  
  // a collection is basically only an array of cells, which is needed to create
  // a finite element space
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  
  // create finite element spaces and functions, a matrix, rhs, and solution
  this->systems.emplace_back(example, *coll);
  
  outputWriter.add_fe_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  
  // print out some information on the finite element space
  auto& v_space = *this->systems.front().velocity_space;
  auto& p_space = *this->systems.front().pressure_space;
  int n_u = v_space.GetN_DegreesOfFreedom();
  int n_u_active = v_space.GetActiveBound();
  int n_p = p_space.GetN_DegreesOfFreedom();
  int n_dof = n_u + n_p;
  
  Output::print<1>(" dof velocity (vector-valued) : ", setw(5), n_u);
  Output::print<1>(" active dof velocity          : ", setw(5), n_u_active);
  Output::print<1>(" dof pressure                 : ", setw(5), n_p);
  Output::print<1>(" dof all                      : ", setw(5), n_dof);
}

/** ************************************************************************ */
Darcy2D::~Darcy2D()
{
  // delete the collections created during the contructor
  for(auto & s : this->systems)
    delete s.velocity_space->GetCollection();
}

/** ************************************************************************ */
void Darcy2D::set_parameters()
{
  if(TDatabase::ParamDB->PRESSURE_SPACE == -4711)
  {
    switch(TDatabase::ParamDB->VELOCITY_SPACE)
    {
      case 1000: // Raviart-Thomas, order 0
        TDatabase::ParamDB->PRESSURE_SPACE = 0;
        break;
      case 1001: // Raviart-Thomas, order 1
        TDatabase::ParamDB->PRESSURE_SPACE = -11;
        break;
      case 1002: // Raviart-Thomas, order 2
        TDatabase::ParamDB->PRESSURE_SPACE = -12;
        break;
      case 1003: // Raviart-Thomas, order 3
        TDatabase::ParamDB->PRESSURE_SPACE = -13;
        break;
      case 1011: // Brezzi-Douglas-Marini, order 1
        TDatabase::ParamDB->PRESSURE_SPACE = 0;
        break;
      case 1012: // Brezzi-Douglas-Marini, order 2
        TDatabase::ParamDB->PRESSURE_SPACE = -110;
        break;
      case 1013: // Brezzi-Douglas-Marini, order 3
        TDatabase::ParamDB->PRESSURE_SPACE = -120;
        break;
      default:
        ErrMsg("unknown velocity space for Darcy2D");
        throw(std::runtime_error("unknown velocity space for Darcy2D"));
        break;
    }
  }
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("multigrid not yet implemented for Darcy2D");
  }
}

/** ************************************************************************ */
void Darcy2D::assemble()
{
  for(System_per_grid& s : this->systems)
  {
    // the class LocalAssembling2D which we will need next, requires an array of
    // pointers to finite element functions, i.e. TFEFunction2D **.
    TFEFunction2D *fe_functions[2] = { &s.u, &s.p };
    // create a local assembling object which is needed to assemble the matrices
    LocalAssembling2D la(this->db, LocalAssembling_type::Darcy,
                         fe_functions, this->example.get_coeffs());
    
    // everything which follows within this for loop only has the goal to call
    // Assemble2D_VectFE at the end. 
    const TFESpace2D * v_space = s.velocity_space.get();
    const TFESpace2D * p_space = s.pressure_space.get();

    const size_t n_fe_spaces = 2;

    const TFESpace2D *fespmat[2] = {v_space, p_space};
    const size_t n_sq_mat = 2;
    TSquareMatrix2D *sq_matrices[n_sq_mat]{};

    const size_t n_rect_mat = 2;
    TMatrix2D *rect_matrices[n_rect_mat]{};

    const size_t n_rhs = 2;
    double *RHSs[n_rhs] = {s.rhs.block(0), s.rhs.block(1)};
    
    BoundCondFunct2D * boundary_conditions[n_fe_spaces] = {
      v_space->GetBoundCondition(), p_space->GetBoundCondition() };
    std::array<BoundValueFunct2D*, n_fe_spaces> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    
    std::vector<std::shared_ptr<FEMatrix>> blocks 
      = s.matrix.get_blocks_uniquely();
    sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
    rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get());
    rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
    
    Assemble2D_VectFE(n_fe_spaces, fespmat, n_sq_mat, sq_matrices, n_rect_mat, 
                      rect_matrices, n_rhs, RHSs, fespmat, la, 
                      boundary_conditions, non_const_bound_values.data());
  }
  // copy Dirichlet values from rhs to solution vector (this is not really 
  // necessary in case of a direct solver)
  this->systems.front().solution.copy_nonactive(this->systems.front().rhs);
} // void Darcy2D::Assemble

/** ************************************************************************ */
void Darcy2D::solve()
{
  double t = GetTime();
  System_per_grid& s = this->systems.front();
  
  this->solver.solve(s.matrix, s.rhs, s.solution);
  
  if(s.matrix.pressure_projection_enabled())
    s.p.project_into_L20();
  
  t = GetTime() - t;
  Output::print<2>(" solving of a Darcy2D problem done in ", t, " seconds");
}

/** ************************************************************************ */
void Darcy2D::output(int i)
{
  System_per_grid & s = this->systems.front();
  if((size_t)db["verbosity"]> 1)
  {
    s.u.PrintMinMax();
    s.p.PrintMinMax();
  }
  
  if(i < 0)
    outputWriter.write();
  else
    outputWriter.write(i);
  if(db["output_compute_errors"])
  {
    DoubleFunct2D *const *Exact = &(example.get_exact())[0];
    ErrorMethod2D *L2DivH1 = L2DivH1Errors;
    // unfortunatly this needs to be longer than one would expect. This is 
    // because TFEFunction2D::GetErrors needs this array to be of size
    // N_Errors+1 (here 2+1).
    std::array<double, 6> errors;
    s.u.GetErrorsForVectorValuedFunction(Exact, L2DivH1, errors.data());
    
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
    const TFESpace2D * pointer_to_p_space = s.pressure_space.get();
    s.p.GetErrors(example.get_exact(2), 3, AllDerivatives, 2, L2H1Errors,
                  nullptr, &aux, 1, &pointer_to_p_space, errors.data() + 3);

    Output::print<1>(" L2(u):      ", setprecision(14), errors[0]);
    Output::print<1>(" L2(div(u)): ", setprecision(14), errors[1]);
    Output::print<1>(" H1-semi(u): ", setprecision(14), errors[2]);
    Output::print<1>(" L2(p):      ", setprecision(14), errors[3]);
    Output::print<1>(" H1-semi(p): ", setprecision(14), errors[4]);
    
    // copy 
    std::copy(errors.begin(), errors.end()-1, this->errors.begin());
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
}

/** ************************************************************************ */
double Darcy2D::getL2VelocityError() const
{
  return this->errors[0];
}

/** ************************************************************************ */
double Darcy2D::getL2DivergenceError() const
{
  return this->errors[1];
}

/** ************************************************************************ */
double Darcy2D::getH1SemiVelocityError() const
{
  return this->errors[2];
}

/** ************************************************************************ */
double Darcy2D::getL2PressureError() const
{
  return this->errors[3];
}

/** ************************************************************************ */
double Darcy2D::getH1SemiPressureError() const
{
  return this->errors[4];
}




