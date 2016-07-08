#include <Brinkman2D.h>
#include <MainUtilities.h> // GetVelocityAndPressureSpace
#include <Database.h>
#include <LinAlg.h> // DDot
#include <DirectSolver.h>
#include <Upwind.h>

#include<Assemble2D.h>

ParameterDatabase get_default_Brinkman2D_parameters()
{
  Output::print<3>("creating a default Brinkman2D parameter database");

  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Brinkman2D parameter database");

  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);


  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);

  return db;
}

/** ************************************************************************ */

Brinkman2D::System_per_grid::System_per_grid (const Example_Brinkman2D& example,
               TCollection& coll, std::pair<int,int> velocity_pressure_orders,
               Brinkman2D::Matrix type)
 : velocity_space(&coll, (char*)"u", (char*)"Brinkman velocity", example.get_bc(0),
                  velocity_pressure_orders.first, nullptr),
   pressure_space(&coll, (char*)"p", (char*)"Brinkman pressure", example.get_bc(2),
                  velocity_pressure_orders.second, nullptr)
// TODO CB: Building the matrix here and rebuilding later is due to the
// highly non-functional class TFEVectFunction2D (and TFEFunction2D,
// which do neither provide default constructors nor working copy assignments.)
   ,matrix({&velocity_space, &velocity_space, &pressure_space}),
   rhs(matrix, true),
   solution(matrix, false),
   u(&velocity_space, (char*)"u", (char*)"u", solution.block(0),
     solution.length(0), 2),
   p(&pressure_space, (char*)"p", (char*)"p", solution.block(2),
     solution.length(2))
{
    

    
// rebuild the matrix due to NSE type. We must be sure, that the rhs and solution
// vector which we built above do fit the new matrix, too!
  switch (type)
  {
    case Brinkman2D::Matrix::Type1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
    break;
    case Brinkman2D::Matrix::Type2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
    break;
    case Brinkman2D::Matrix::Type3:
     matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
    break;
    case Brinkman2D::Matrix::Type4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
    break;
    case Brinkman2D::Matrix::Type14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
     break;
    default:
      ErrThrow("Unknown NSE type given to constructor of NSE2D::System_per_grid.");
  }

}

/** ************************************************************************ */
Brinkman2D::Brinkman2D(const TDomain& domain, const ParameterDatabase& param_db,
                       int reference_id)
 : Brinkman2D(domain, param_db, Example_Brinkman2D(param_db["example"],param_db),
              reference_id)
{
  // note that the way we construct the example above will produce a memory 
  // leak, but that class is small.
  // FIXME Find a workaround - we do not like memory leaks at all,
  // because they pollute our valgrind tests.
}

/** ************************************************************************ */

Brinkman2D::Brinkman2D(const TDomain & domain, const ParameterDatabase& param_db,
                       const Example_Brinkman2D & e,
             unsigned int reference_id)
    : db(get_default_Brinkman2D_parameters()), outputWriter(param_db),
      systems(), example(e), solver(param_db), defect(), oldResiduals(),
      initial_residual(1e10), errors()
{
  db.merge(param_db,false);
  
  std::pair <int,int>
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and preesure spaces
  // this function returns a pair which consists of 
  // velocity and pressure order
  this->get_velocity_pressure_orders(velocity_pressure_orders);
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  
 
  // we use always Matrix Type 14
    this->systems.emplace_back(example, *coll, velocity_pressure_orders, Brinkman2D::Matrix::Type14);
  
  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);

  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());

  // print out some information
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_u_active = this->get_velocity_space().GetN_ActiveDegrees();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 2 * n_u + n_p; // total number of degrees of freedom
 
    
  double h_min, h_max;
  coll->GetHminHmax(&h_min, &h_max);
  Output::print<1>("N_Cells            : ", setw(10), coll->GetN_Cells());
  Output::print<1>("h (min,max)        : ", setw(10), h_min, " ", setw(12),
                   h_max);
  Output::print<1>("dof velocity       : ", setw(10), 2* n_u);
  Output::print<1>("dof velocity active: ", setw(10), 2* n_u_active);
  Output::print<1>("dof pressure       : ", setw(10), n_p);
  Output::print<1>("dof all            : ", setw(10), n_dof);
  

  // done with the constructor in case we're not using multigrid
  if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE != 5
    || TDatabase::ParamDB->SOLVER_TYPE != 1)
    return;
  // else multigrid

}

/** ************************************************************************ */
Brinkman2D::~Brinkman2D()
{
}

/** ************************************************************************ */


void Brinkman2D::get_velocity_pressure_orders(std::pair <int,int>
                 &velocity_pressure_orders)
{
//  int velocity_order = velocity_pressure_orders.first;
//  int pressure_order = velocity_pressure_orders.second;
   Output::print("velocity space", setw(10), TDatabase::ParamDB->VELOCITY_SPACE);
  Output::print("pressure space", setw(10), TDatabase::ParamDB->PRESSURE_SPACE);
}

void Brinkman2D::set_parameters()
{
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
  {
    //Assemble2DSlipBC does not work, and is not implemented yet
    ErrThrow("Set INTERNAL_SLIP_WITH_FRICTION to 0, this feature is not yet included.");
  }
}

/** ************************************************************************ */


void Brinkman2D::assemble()
{
    
   
  for(System_per_grid& s : this->systems)
  {
    s.rhs.reset(); //right hand side reset (TODO: is that necessary?)

      
    const TFESpace2D * v_space = &s.velocity_space;
    const TFESpace2D * p_space = &s.pressure_space;

    // declare the variables which Assemble2D needs and each brinkmantype has to fill
    size_t N_FESpaces = 2;

    const TFESpace2D *fespmat[2] = {v_space, p_space};
    size_t n_sq_mat;
    TSquareMatrix2D *sq_matrices[5]{nullptr};//it's five pointers maximum (Type14)

    size_t n_rect_mat;
    TMatrix2D *rect_matrices[4]{nullptr};//it's four pointers maximum (Types 2, 4, 14)

    size_t N_Rhs = 2; //is 3 if NSE type is 4 or 14
    double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), nullptr}; //third place gets only filled
    const TFESpace2D *fesprhs[3] = {v_space, v_space, nullptr};  // if NSE type is 4 or 14
    BoundCondFunct2D * boundary_conditions[3] = {
      v_space->GetBoundCondition(), v_space->GetBoundCondition(),
      p_space->GetBoundCondition() };
    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    //same for all: the local asembling object
    TFEFunction2D *fe_functions[3] =
      { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    LocalAssembling2D la(Brinkman2D_Galerkin1, fe_functions,
                     this->example.get_coeffs());

    std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely();
  
// Note: We use only Type 14 for Brinkman (for now)
      //CB DEBUG
        if(blocks.size() != 9)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        //END DEBUG
        n_sq_mat = 5;
        sq_matrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sq_matrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sq_matrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sq_matrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        sq_matrices[4] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());

        n_rect_mat = 4;
        rect_matrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rect_matrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rect_matrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rect_matrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

        RHSs[2] = s.rhs.block(2); // NSE type 14 includes pressure rhs
        fesprhs[2]  = p_space;
        N_Rhs = 3;

      
 Output::print("assemble");
 
      
    // call the assemble method with the information that has been patched together
    Assemble2D(N_FESpaces, fespmat, n_sq_mat, sq_matrices,
               n_rect_mat, rect_matrices, N_Rhs, RHSs, fesprhs,
               boundary_conditions, non_const_bound_values.data(), la);
      Output::print("assemble");
      
    // do upwinding TODO remove dependency of global values
    if((TDatabase::ParamDB->DISCTYPE == UPWIND)
       && !(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 3))
    {
      switch(TDatabase::ParamDB->BrinkmanTYPE)
      {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                                la.get_fe_function(0), la.get_fe_function(1));
          Output::print<3>("UPWINDING DONE : level ");
          break;
        case 3:
        case 4:
        case 14:
          // do upwinding with two matrices
          Output::print<3>("UPWINDING DONE : level ");
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                                la.get_fe_function(0), la.get_fe_function(1));
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[3],
                                la.get_fe_function(0), la.get_fe_function(1));
          break;
      } // endswitch
    } // endif

    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);

    // TODO Maybe we have to explicitely set non-actives in non-diagonal blocks
    // to zero here, that was done in former code, but maybe we can move it to the solver part

    //tidy up
    delete fe_functions[0];
    delete fe_functions[1];
      

  }
}



/** ************************************************************************ */
bool Brinkman2D::stopIt(unsigned int iteration_counter)
{
  // compute the residuals with the current matrix and solution
  this->computeNormsOfResiduals();
  // the current norm of the residual
  const double normOfResidual = this->getFullResidual();
  // store initial residual, so later we can print the overall reduction
  if(iteration_counter == 0)
    initial_residual = normOfResidual;
  // the residual from 10 iterations ago
  const double oldNormOfResidual = this->oldResiduals.front().fullResidual;
  
  const unsigned int Max_It = db["nonlinloop_maxit"];
  const double convergence_speed = db["nonlinloop_slowfactor"];
  bool slow_conv = false;
  
  
  if(normOfResidual >= convergence_speed*oldNormOfResidual)
    slow_conv = true;
  
  double limit = db["nonlinloop_epsilon"];
  if ( db["nonlinloop_scale_epsilon_with_size"] )
  {
    limit *= sqrt(this->get_size());
    Output::print<1>("stopping tolerance for nonlinear iteration ", limit);
  }

  // check if the iteration has converged, or reached the maximum number of
  // iterations or if convergence is too slow. Then return true otherwise false
  if( (normOfResidual<=limit) || (iteration_counter==Max_It) || (slow_conv) )
  {
    if(slow_conv)
      Output::print<1>(" SLOW !!! ", normOfResidual/oldNormOfResidual);
    // stop iteration
    Output::print<1>(" ITE : ", setw(4), iteration_counter, setprecision(8),
                     " RES : ", normOfResidual, " Reduction : ",
                     normOfResidual/initial_residual);
    return true;
  }
  else
    return false;
}


/** ************************************************************************ */
void Brinkman2D::computeNormsOfResiduals()
{
  System_per_grid& s = this->systems.front();
  unsigned int n_u_dof = s.solution.length(0);
  unsigned int n_p_dof = s.solution.length(2);
  
  // copy rhs to defect
  this->defect = s.rhs;
  s.matrix.apply_scaled_add(s.solution, defect,-1.);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    IntoL20FEFunction(&defect[2*n_u_dof], n_p_dof, &this->get_pressure_space(),
                      TDatabase::ParamDB->VELOCITY_SPACE, 
                      TDatabase::ParamDB->PRESSURE_SPACE);
  }
  
  // square of the norms of the residual components
  double impuls_Residual = Ddot(2*n_u_dof, &this->defect[0],&this->defect[0]);
  double mass_residual = Ddot(n_p_dof, &this->defect[2*n_u_dof],
                              &this->defect[2*n_u_dof]);
  
  Residuals currentResiduals(impuls_Residual, mass_residual);
  oldResiduals.add(currentResiduals);
}

/** ************************************************************************ */
void Brinkman2D::solve()
{
  System_per_grid& s = systems.front();
  double damping = db["nonlinloop_damping_factor"];
  // store previous solution for damping, it is a pointer so that we can avoid
  // the copy in case of no damping
  std::shared_ptr<BlockVector> old_solution(nullptr);
  if(damping != 1.0)
    old_solution = std::make_shared<BlockVector>(s.solution);

  solver.solve(s.matrix,s.rhs, s.solution);

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    s.p.project_into_L20();

}

/** ************************************************************************ */
void Brinkman2D::output(int i)
{

	bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
	if(no_output)
		return;

  System_per_grid& s = this->systems.front();
  TFEFunction2D* u1 = s.u.GetComponent(0);
  TFEFunction2D* u2 = s.u.GetComponent(1);
  
  // print the value of the largest and smallest entry in the finite element 
  // vector
  if((size_t)db["verbosity"]> 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    s.p.PrintMinMax();
  }
  
  // write solution to a vtk file
  if(db["output_write_vtk"])
  {
    outputWriter.write(i);
  }
  
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    double err[4];
    TAuxParam2D NSEaux_error;
    MultiIndex2D NSAllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D *velocity_space = &this->get_velocity_space();
    const TFESpace2D *pressure_space = &this->get_pressure_space();
    
    // errors in first velocity component
    u1->GetErrors(example.get_exact(0), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &velocity_space, err);
    // errors in second velocity component
    u2->GetErrors(example.get_exact(1), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &velocity_space, err + 2);
    
    errors.at(0) = sqrt(err[0]*err[0] + err[2]*err[2]);
    errors.at(1) = sqrt(err[1]*err[1] + err[3]*err[3]);    
    Output::print<1>("L2(u)     : ", errors[0]);
    Output::print<1>("H1-semi(u): ", errors[1]);
    // errors in pressure
    s.p.GetErrors(example.get_exact(2), 3, NSAllDerivatives, 2, L2H1Errors, 
                  nullptr, &NSEaux_error, 1, &pressure_space, err);
    
    errors.at(2) = err[0];
    errors.at(3) = err[1];    
    Output::print<1>("L2(p)     : ", errors[2]);
    Output::print<1>("H1-semi(p): ", errors[3]);    
    
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
  delete u1;
  delete u2;
}

/** ************************************************************************ */
std::array< double, int(6) > Brinkman2D::get_errors()
{
  return errors;
}

/** ************************************************************************ */
const Brinkman2D::Residuals& Brinkman2D::getResiduals() const
{
  return this->oldResiduals.back();
}

/** ************************************************************************ */
double Brinkman2D::getImpulsResidual() const
{
  return this->oldResiduals.back().impulsResidual;
}

/** ************************************************************************ */
double Brinkman2D::getMassResidual() const
{
  return this->oldResiduals.back().massResidual;
}

/** ************************************************************************ */
double Brinkman2D::getFullResidual() const
{
  return this->oldResiduals.back().fullResidual;
}

/** ************************************************************************ */
Brinkman2D::Residuals::Residuals()
 : impulsResidual(1e10), massResidual(1e10), fullResidual(1e10)
{}

/** ************************************************************************ */
Brinkman2D::Residuals::Residuals(double imR, double maR)
 : impulsResidual(sqrt(imR)), massResidual(sqrt(maR)),
   fullResidual(sqrt(imR + maR))
{}

/** ************************************************************************ */
std::ostream& operator<<(std::ostream& s, const Brinkman2D::Residuals& n)
{
  s << setw(14) << n.impulsResidual << "\t" << setw(14)
    << n.massResidual << "\t" << setw(14) << n.fullResidual;
  return s;
}

