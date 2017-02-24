#include <Time_NSE2D_BDF.h>
#include <Database.h>
#include <Assemble2D.h>
#include <GridTransfer.h>
#include <LinAlg.h>

/* *************************************************************************** */
  //TODO  So far of this object only the nonlin it stuff is used - switch entirely!
ParameterDatabase get_default_TNSE2D_parameters()
{
  Output::print<5>("creating a default TNSE2D parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default TNSE2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TNSE2D parameter database");

  //Time_NSE2D requires a nonlinear iteration, set up a nonlinit_database and merge
  ParameterDatabase nl_db = ParameterDatabase::default_nonlinit_database();
  db.merge(nl_db,true);

  // a default output database - needed here as long as there's no class handling the output
  ParameterDatabase out_db = ParameterDatabase::default_output_database();
  db.merge(out_db, true);
  
  ParameterDatabase space_disc_db = ParameterDatabase::default_space_disc_database();
  db.merge(out_db, true);
  
  ParameterDatabase turbulence_db = ParameterDatabase::default_turbulence_model_database();
  db.merge(turbulence_db, true);


  ParameterDatabase time_disc_db = ParameterDatabase::default_time_database();
  db.merge(time_disc_db, true);
  
  return db;
}

Time_NSE2D_BDF::System_per_grid::System_per_grid(const Example_TimeNSE2D& example, 
   TCollection& coll, std::pair< int, int > order, Time_NSE2D_BDF::Matrix type)
: velocity_space(&coll, (char*)"u", (char*)"velocity space",  example.get_bc(0),
                  order.first, nullptr),
   pressure_space(&coll, (char*)"p", (char*)"pressure space", example.get_bc(2),
                  order.second, nullptr),
   matrix({&velocity_space, &velocity_space, &pressure_space}),
   Mass_Matrix({&velocity_space, &velocity_space}),
   rhs(matrix, true),
   solution(matrix, false),
   u(&velocity_space, (char*)"u", (char*)"u", solution.block(0), 
     solution.length(0), 2),
   p(&pressure_space, (char*)"p", (char*)"p", this->solution.block(2),
     solution.length(2)), 
   MatrixK({&velocity_space, &velocity_space, &pressure_space}),
   solution_m1(matrix, false),
   u_m1(&velocity_space, (char*)"u", (char*)"u", solution_m1.block(0), 
        solution_m1.length(0), 2),
   solution_m2(matrix, false), 
   u_m2(&velocity_space, (char*)"u", (char*)"u", solution_m2.block(0), 
        solution_m2.length(0), 2),
   p_old(&pressure_space, (char*)"p", (char*)"p", this->solution_m1.block(2),
     solution_m1.length(2)),
   combined_old_sols(matrix, false),
   comb_old_u(&velocity_space, (char*)"u", (char*)"u", combined_old_sols.block(0), 
        combined_old_sols.length(0), 2), 
   extrapolate_sol(matrix, false),
   extrapolate_u(&velocity_space, (char*)"u", (char*)"u", extrapolate_sol.block(0), 
        extrapolate_sol.length(0), 2)
{
  // Mass Matrix
  // Output::increaseVerbosity(5);
  //solution_m2 = solution_m1;
  Mass_Matrix = BlockFEMatrix::Mass_NSE2D(velocity_space);
      
  switch(type)
  {
    case Time_NSE2D_BDF::Matrix::Type1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
      break;
    case Time_NSE2D_BDF::Matrix::Type2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
      break;
    case Time_NSE2D_BDF::Matrix::Type3:
      matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
      break;
    case Time_NSE2D_BDF::Matrix::Type4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);

      MatrixK = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      break;
    case Time_NSE2D_BDF::Matrix::Type14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);

      MatrixK = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      break;
  }
}

Time_NSE2D_BDF::Time_NSE2D_BDF(const TDomain& domain, const ParameterDatabase& param_db,
                       const Example_TimeNSE2D& ex, int reference_id)
: db(get_default_TNSE2D_parameters()), outputWriter(param_db), systems(),
   example(ex), solver(param_db), defect(), oldResiduals(), 
   initial_residual(1e10), errors(10,0.), oldtau(0.0), 
   pre_step_time_db(ParameterDatabase::default_time_database())
{
  db.merge(param_db);
  this->set_parameters();
  
  std::pair <int,int> velo_pres_order(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  this->get_velocity_pressure_orders(velo_pres_order);
  
  Time_NSE2D_BDF::Matrix type;
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case  1: type = Matrix::Type1;  break;
    case  2: type = Matrix::Type2;  break;
    case  3: type = Matrix::Type3;  break;
    case  4: type = Matrix::Type4;  break;
    case 14: type = Matrix::Type14; break;
    default:
      ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
               " That NSE Block Matrix Type is unknown to class NSE2D.");
  }
  
  bool usingMultigrid = this->solver.is_using_multigrid();
  // create the collection of cells from the domain (finest grid)
  TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
  // create finite element space, functions, matrices, rhs and solution
  // at the finest grid
  this->systems.emplace_back(example, *coll, velo_pres_order, type);
  
  TFEFunction2D * u1 = this->systems.front().u.GetComponent(0);
  TFEFunction2D * u2 = this->systems.front().u.GetComponent(1);
  
  u1->Interpolate(example.get_initial_cond(0));
  u2->Interpolate(example.get_initial_cond(1));

  if(usingMultigrid)
  {
    auto mg = solver.get_multigrid();
    bool mdml = mg->is_using_mdml();
    if(mdml)
    {
      // change the discretization on the coarse grids to lowest order 
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velo_pres_order = {-1, -4711};
      this->get_velocity_pressure_orders(velo_pres_order);
    }
    //determine multigrid levels
    int second_grid;
    int coarsest_grid = domain.get_ref_level() -
                         mg->get_n_geometric_levels() + 1;
    if(mdml)
      second_grid = domain.get_ref_level();
    else 
      second_grid = domain.get_ref_level() - 1;
    
    if(coarsest_grid < 0 )
    {
      ErrThrow("The domain has not been refined often enough to do multigrid "
          "on ", mg->get_n_geometric_levels(), " geometric levels (",
          mg->get_n_algebraic_levels()," algebraic levels). There are"
          " only ", domain.get_ref_level() + 1, " geometric grid levels.");
    }
    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    //prepare input argument for multigrid object
    matrices.push_back(&systems.back().matrix);
    for(int grid_no = second_grid; grid_no >= coarsest_grid; --grid_no)
    {
      TCollection *coll = domain.GetCollection(It_EQ, grid_no, reference_id);
      systems.emplace_back(example, *coll, velo_pres_order,
                            type);
      matrices.push_front(&systems.back().matrix);
    }
    mg->initialize(matrices);
  }
  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);
  
  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  
  // print out the information (cells, dofs, etc)
  this->output_problem_size_info();
}

/**************************************************************************** */
void Time_NSE2D_BDF::set_parameters()
{
  if(!db["problem_type"].is(6))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 6;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to Time_NSE."
          "It is now reset to the correct value for Time_NSE (=6).");
      db["problem_type"] = 6;
    }
  }
  if(TDatabase::TimeDB->TIME_DISC == 0)
  {
    ErrMsg("TIME_DISC: " << TDatabase::TimeDB->TIME_DISC 
          << " does not supported");
    throw("TIME_DISC: 0 is not supported");
  }
  if(TDatabase::ParamDB->DISCTYPE == SUPG)
  {
    if(TDatabase::TimeDB->TIME_DISC !=5 && TDatabase::TimeDB->TIME_DISC !=1)
    {
      ErrThrow("TIME_DISC: " , TDatabase::TimeDB->TIME_DISC, 
               " does not supported for SUPG method" , 
               " Only (BDF2) TIME_DISC: 5 can be used");
    }
  }
  
  // set the parameters in the global database
  if(db["disctype"].is("galerkin"))
    TDatabase::ParamDB->DISCTYPE = 1;
  else if (db["disctype"].is("supg"))
  {
    if(db["ansatz_test_extrapolate"].is("only_velocity_test"))
      TDatabase::ParamDB->DISCTYPE = -2;
    else
      TDatabase::ParamDB->DISCTYPE = 2;
  }
  else if (db["disctype"].is("residual_based_vms"))
  {
    if(db["ansatz_test_extrapolate"].is("only_velocity_test"))
      TDatabase::ParamDB->DISCTYPE = -101;
    else
      TDatabase::ParamDB->DISCTYPE = 101;
  }
  else
  {
    ErrThrow("disctype ", db["disctype"], " is not supported " );
  }
  // set time disc parameters
  if(db["time_discretization"].is("backward_euler"))
  {
    // db["theta1"] = 1.0;
    // db["theta2"] = 0.0;
    // db["theta3"] = 0.0;
    // db["theta4"] = 1.0;
    
    TDatabase::TimeDB->TIME_DISC = 1;    
  }
  else if (db["time_discretization"].is("crank_nicolson"))
  {
    // db["theta1"] = 0.5;
    // db["theta2"] = 0.5;
    // db["theta3"] = 0.5;
    // db["theta4"] = 0.5;
    TDatabase::TimeDB->TIME_DISC = 2;
  }
  else if (db["time_discretization"].is("fractional_step"))
  {
    ErrThrow("Not yet supported ");
    // db["theta1"] = 0.0;
    // db["theta2"] = 0.0;
    // db["theta3"] = 0.0;
    // db["theta4"] = 0.0;
    TDatabase::TimeDB->TIME_DISC = 3;
  }  
  else if (db["time_discretization"].is("bdf_two"))
  {
//     // at first step one uses the b
//     if(current_step_ == 1)
//       TDatabase::TimeDB->TIME_DISC = 1;
//     else
//     {
//       TDatabase::TimeDB->THETA1 = 2./3.;
//       TDatabase::TimeDB->THETA2 = 0.0;
//       TDatabase::TimeDB->THETA3 = 0.0;
//       TDatabase::TimeDB->THETA4 = 2./3.;
//     }
    TDatabase::TimeDB->TIME_DISC == 5;
  }
  
  TDatabase::TimeDB->TIMESTEPLENGTH = db["time_step_length"];
  TDatabase::TimeDB->ENDTIME = db["time_end"];
  TDatabase::TimeDB->STARTTIME = db["time_start"];
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = db["current_time_step_length"];
}

/**************************************************************************** */
void Time_NSE2D_BDF::get_velocity_pressure_orders(std::pair< int, int > &velo_pres_order)
{
  int velocity_order = velo_pres_order.first;
  int pressure_order = velo_pres_order.second;
  int order = 0;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5:
    case 12: case 13: case 14: case 15:
      if(velocity_order > 10)
        order = velocity_order-10;
      else
        order = velocity_order;
      break;
    case -1: case -2: case -3: case -4: case -5:
    case -101:
      order = velocity_order;
      break;
    // conforming fe spaces with bubbles on triangles 
    case 22: case 23: case 24:
      order = velocity_order;
      break;
      // discontinuous spaces 
    case -11: case -12: case -13:
      order = velocity_order*10;
      break;
  }
  TDatabase::ParamDB->VELOCITY_SPACE = order;
  velo_pres_order.first = order;
  switch(pressure_order)
  {
    case -4711:
      switch(velocity_order)
      {
        case -1:
        case -2:
        case -3:
        case -4:
          // nonconforming pw (bi)linear velo/ pw constant pressure
          // conforming pw (bi)linear velo/ pw constant pressure (not stable !!!)
          pressure_order = -velocity_order-1;
          break; 
        case 1: // discontinuous space 
          pressure_order = 0;
          break;
        case 2: case 3: case 4: case 5:
        // standard conforming velo and continuous pressure
          pressure_order = velocity_order-1;
          break;
          // discontinuous pressure spaces 
          // standard conforming velo and discontinuous pressure
          // this is not stable on triangles !!!
        case 12: case 13: case 14: case 15:
          pressure_order = -(velocity_order-1)*10;
          break;
        case 22: case 23: case 24:
          pressure_order = -(velocity_order-11)*10;
          break;
      }
      break;
    case 1:case 2: case 3: case 4: case 5:
      // pressure order is chosen correctly
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
    default:
      ErrThrow("pressure space is not chosen properly ", pressure_order);
  }
  TDatabase::ParamDB->PRESSURE_SPACE  = pressure_order;
  velo_pres_order.second = pressure_order;
  
  Output::print("velocity space", setw(10), velo_pres_order.first);
  Output::print("pressure space", setw(10), velo_pres_order.second);
}


void Time_NSE2D_BDF::assemble_initial_time()
{
  for(auto &s : this->systems)
  {
    call_assembling_routine(s, LocalAssembling2D_type::TNSE2D);
    // copy the solution 
    s.solution_m1 = s.solution;
    s.solution_m2 = s.solution;
  }
  old_solution = systems.front().solution;
}

void Time_NSE2D_BDF::assemble_rhs(bool rhs_assemble)
{
  if(db["time_discretization"].is("bdf_two") && current_step_ == 2)
  {
    TDatabase::TimeDB->TIME_DISC = 5;
    if(current_step_ == 2)
    {
      TDatabase::TimeDB->THETA1 = 2./3.;
      TDatabase::TimeDB->THETA2 = 0.0;
      TDatabase::TimeDB->THETA3 = 0.0;
      TDatabase::TimeDB->THETA4 = 2./3.;
      
      Output::print<1>("New parameters for the BDF's schemes");
      Output::print<1>("Theta1: ", TDatabase::TimeDB->THETA1);
      Output::print<1>("Theta2: ", TDatabase::TimeDB->THETA2);
      Output::print<1>("Theta3: ", TDatabase::TimeDB->THETA3);
      Output::print<1>("Theta4: ", TDatabase::TimeDB->THETA4);
    }
  }
  System_per_grid& s = this->systems.front();
  if(rhs_assemble)
    call_assembling_routine(s, TNSE2D_Rhs);
  s.solution.copy_nonactive(s.rhs);
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double t4 = TDatabase::TimeDB->THETA4;
  // scale the current right hand side with the factor tau*theta4
  s.rhs.scaleActive(tau*t4);
  s.matrix.get_blocks().at(0)->Print("M");
  if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
    this->modify_slip_bc(true);
  exit(0);
  // mass matrix contribution to the right hand side
  if(db["time_discretization"].is("bdf_two") && current_step_ > 1)
  {
    if(db["disctype"].is("residual_based_vms"))
    {
      s.MatrixK.apply_scaled_submatrix(s.solution_m1, s.rhs, 2, 2, 4./3.);
      s.MatrixK.apply_scaled_submatrix(s.solution_m2, s.rhs, 2, 2, -1./3.);
    }
    else
    {
      s.Mass_Matrix.apply_scaled_submatrix(s.solution_m1, s.rhs, 2, 2, 4./3.);
      s.Mass_Matrix.apply_scaled_submatrix(s.solution_m2, s.rhs, 2, 2, -1./3.);
    }
  }
  else
  {
    if(db["disctype"].is("residual_based_vms"))
      s.MatrixK.apply_scaled_submatrix(old_solution, s.rhs, 2, 2, 1.0);
    else
      s.Mass_Matrix.apply_scaled_submatrix(old_solution, s.rhs, 2, 2, 1.0);
  }
  
  for(auto &s: this->systems)
  {
    if(db["ansatz_test_extrapolate"].is("only_velocity_test"))
    {
      // in this case one have to re-assemble the mass matrix and 
      // b-blocks, so the scaling with the factor times step length
      // is important
      double t1 = TDatabase::TimeDB->THETA1;
      s.matrix.scale_blocks_actives(tau*t1, {{0,2}, {1,2}, {2,0}, {2,1}});
      // TODO: C-block depends on solution because of the stabilization parameter
      // check and re-assemble again
      if(TDatabase::ParamDB->NSTYPE == 14)/// && current_step_ ==1)
        s.matrix.scale_blocks_actives(tau*t1, {{2,2}});
    }
  }
  // copy back the dirichlet dof's from solution to rhs
  s.rhs.copy_nonactive(s.solution);  
}

void Time_NSE2D_BDF::assemble_nonlinear_term()
{
  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(systems.size() > 1)
  {
    for( int b = 0; b < 2 ;++b)
    {
      std::vector<const TFESpace2D*> spaces;
      std::vector<double*> u_entries;
      std::vector<size_t> u_ns_dofs;
      for(auto &s : systems )
      {
        spaces.push_back(&s.velocity_space);
        u_entries.push_back(s.solution.block(b));
        u_ns_dofs.push_back(s.solution.length(b));
      }
      GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
    }
  }
  
  for(System_per_grid &s : this->systems)
    this->call_assembling_routine(s, LocalAssembling2D_type::TNSE2D_NL);
  
  if( (db["disctype"].is("supg") || db["disctype"].is("residual_based_vms"))
     && (db["ansatz_test_extrapolate"].is("no_extrapolation")) )
  {
    // rhs and the nonlinear matrices are assemble during the function call 
    // "assemble_nonlinear_term()". The right-hand side for the solver needs 
    // to be re-assemble
    this->assemble_rhs(0);
  }
  
  if( (db["disctype"].is("supg") || db["disctype"].is("residual_based_vms"))
     && (db["ansatz_test_extrapolate"].is("only_velocity_test")) 
     && (TDatabase::ParamDB->NSTYPE == 14))
  {
    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    double t1 = TDatabase::TimeDB->THETA1;
    for(System_per_grid& s : this->systems)
      s.matrix.scale_blocks_actives(tau*t1, {{2,0}, {2,1}});
  }
}

void Time_NSE2D_BDF::assemble_system()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double t1 = TDatabase::TimeDB->THETA1;
  
  for(System_per_grid &s : systems)
  {
    const std::vector<std::vector<size_t>> pos = {{0,0}, {0,1}, {1,0}, {1,1}};
    s.matrix.scale_blocks_actives(t1*tau, pos);
    if(db["disctype"].is("residual_based_vms"))
    {
      const FEMatrix& M00 = *s.MatrixK.get_blocks().at(0).get();
      s.matrix.add_matrix_actives(M00, 1.0, {{0,0}}, {false});
      
      const FEMatrix& M01 = *s.MatrixK.get_blocks().at(1).get();
      s.matrix.add_matrix_actives(M01, 1.0, {{0,1}}, {false});
      
      const FEMatrix& M10 = *s.MatrixK.get_blocks().at(3).get();
      s.matrix.add_matrix_actives(M10, 1.0, {{1,0}}, {false});
      
      const FEMatrix& M11 = *s.MatrixK.get_blocks().at(4).get();
      s.matrix.add_matrix_actives(M11, 1.0, {{1,1}}, {false});
    }
    else
    {
      // add mass matrix to the scaled stiff matrix
      const FEMatrix& M = *s.Mass_Matrix.get_blocks().at(0).get();
      
      s.matrix.add_matrix_actives(M, 1.0, {{0,0}, {1,1}}, {false, false});
    }
  }
}

bool Time_NSE2D_BDF::stopIte(unsigned int it_counter)
{
  System_per_grid& s = this->systems.front();
  unsigned int nuDof = s.solution.length(0);
  unsigned int npDof = s.solution.length(2);
  unsigned int sc_minit = db["nonlinloop_minit"];
  
  this->defect = s.rhs; 
  s.matrix.apply_scaled_add(s.solution, defect,-1.);
  // 
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20FEFunction(&defect[2*nuDof], npDof, &this->get_pressure_space(),
                      TDatabase::ParamDB->VELOCITY_SPACE, 
                      TDatabase::ParamDB->PRESSURE_SPACE);
  double residual =  Ddot(2*nuDof+npDof, &this->defect[0], &this->defect[0]);
  double impulse_residual = Ddot(2*nuDof, &this->defect[0],
         &this->defect[0]);
  double mass_residual    = Ddot(npDof,&this->defect[2*nuDof],
         &this->defect[2*nuDof]);
  Output::print<1>("nonlinear step : ", it_counter);
  if(it_counter > 0)
  {
    Output::print("Residuals  :  " , setw(15), impulse_residual, 
                setw(15), mass_residual, setw(15), sqrt(residual),
                setw(15), sqrt(residual)/oldResidual);
  }
  else
  {
    Output::print("Residuals  :  " , setw(15), impulse_residual, 
                setw(15), mass_residual, setw(15), sqrt(residual));
  }
  
  oldResidual = sqrt(residual);
  if(it_counter == 0)
    initial_residual = sqrt(residual);
  
  size_t Max_It = db["nonlinloop_maxit"];
  double limit = db["nonlinloop_epsilon"];
  if (db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= sqrt(this->get_size());
    Output::print("stopping tolerance for nonlinear iteration ", limit);
  }
  
  if ( (((sqrt(residual)<=limit)||(it_counter==Max_It))) )
   {
     for(System_per_grid& s: this->systems)
     {
       s.solution_m2 = s.solution_m1;
       s.solution_m1 = s.solution;
     }
     this->old_solution = s.solution;
     
     Output::print("ITE : ", setw(3), it_counter, "  RES : ", sqrt(residual), 
                   " Reduction : ",  sqrt(residual)/initial_residual);
     
     if(imex_scheme(0) && it_counter >0)
     {
       return true;
     }
     else
     {
       // descale the matrices, since only the diagonal A block will 
       // be reassembled in the next time step
       this->deScaleMatrices();
       return true;
     }     
   }
   else
     return false;  
}

void Time_NSE2D_BDF::solve()
{
  System_per_grid& s = this->systems.front();

  solver.solve(s.matrix, s.rhs, s.solution);
  
  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11 and A22 matrices are 
  // reset and assembled again but the A12 and A21 are scaled, so
  // for the next iteration we have to descale, see assemble_system()
  this->deScaleMatrices();

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       s.p.project_into_L20();
}

void Time_NSE2D_BDF::deScaleMatrices()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double t1 = TDatabase::TimeDB->THETA1;  

  for(System_per_grid& s : this->systems)
  {
    if(db["disctype"].is("residual_based_vms"))
    {
      const FEMatrix& M00 = *s.MatrixK.get_blocks().at(0).get();
      s.matrix.add_matrix_actives(M00, -1.0, {{0,0}}, {false});
      
      const FEMatrix& M01 = *s.MatrixK.get_blocks().at(1).get();
      s.matrix.add_matrix_actives(M01, -1.0, {{0,1}}, {false});
      
      const FEMatrix& M10 = *s.MatrixK.get_blocks().at(3).get();
      s.matrix.add_matrix_actives(M10, -1.0, {{1,0}}, {false});
      
      const FEMatrix& M11 = *s.MatrixK.get_blocks().at(4).get();
      s.matrix.add_matrix_actives(M11, -1.0, {{1,1}}, {false});
    }
    else 
    {
      const FEMatrix& M = *s.Mass_Matrix.get_blocks().at(0).get();
      s.matrix.add_matrix_actives(M, -1.0, {{0,0}, {1,1}}, {false, false});
    }
    
    const std::vector<std::vector<size_t>>pos = {{0,0}, {0,1}, {1, 0}, {1, 1}};
    s.matrix.scale_blocks_actives(1./(t1*tau), pos);    
  }
}


void Time_NSE2D_BDF::call_assembling_routine(Time_NSE2D_BDF::System_per_grid& s, 
 LocalAssembling2D_type type)
{
    // set arrays of spaces for matrices and rhs 
  std::vector<const TFESpace2D*> spaces_mat;
  std::vector<const TFESpace2D*> spaces_rhs;
  std::vector<TFEFunction2D*> fefunctios;
  // call to routine to set arrays 
  set_arrays(s, spaces_mat, spaces_rhs, fefunctios);
  
  // prepare matrices and rhs for assembling 
  std::vector<TSquareMatrix2D*> sqMatrices;
  std::vector<TMatrix2D*> rectMatrices;
  std::vector<double*> rhs_array;
  // call the routine to prepare the matrices
  set_matrices_rhs(s, type, sqMatrices, rectMatrices, rhs_array);
  // boundary conditions and boundary values array
  // boundary conditions:
  std::vector<const BoundCondFunct2D*> bc(3);
  bc[0]=s.velocity_space.GetBoundCondition();
  bc[1]=bc[0];
  bc[2]=s.pressure_space.GetBoundCondition();
  
  // boundary values:
  std::vector<BoundValueFunct2D*>bv(3);
  bv[0]=example.get_bd(0);
  bv[1]=example.get_bd(1);
  bv[2]=example.get_bd(2);
  
    // local assembling settings
  LocalAssembling2D la(type, fefunctios.data(),
                           this->example.get_coeffs());

  // assemble all the matrices and right hand side 
  Assemble2D(spaces_mat.size(), spaces_mat.data(), 
               sqMatrices.size(), sqMatrices.data(), 
               rectMatrices.size(), rectMatrices.data(), 
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(), 
               bc.data(), bv.data(), la);
}

void Time_NSE2D_BDF::set_arrays(Time_NSE2D_BDF::System_per_grid& s, std::vector<const TFESpace2D*> &spaces, 
std::vector< const TFESpace2D* >& spaces_rhs, std::vector< TFEFunction2D*> &functions)
{
  spaces.resize(2);
  spaces_rhs.resize(2);
  
  spaces[0] = &s.velocity_space;
  spaces[1] = &s.pressure_space;
  
  spaces_rhs[0] = &s.velocity_space;
  spaces_rhs[1] = &s.velocity_space;
  
  if(TDatabase::ParamDB->NSTYPE == 14)
  {
    spaces_rhs.resize(3);
    spaces_rhs[2] = &s.pressure_space;
  }
  
  bool is_imex = imex_scheme(0);
  functions.resize(2);
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  if(!is_imex)
  {
    functions[0] = s.u.GetComponent(0);
    functions[1] = s.u.GetComponent(1);
    // extrapolation is used in the weighted test functions occurs
    //  due the supg/rbVMS method 
    if( (db["disctype"].is("supg") || db["disctype"].is("residual_based_vms") )
       && (db["ansatz_test_extrapolate"].is("only_velocity_test")) ) 
    { 
      // The extrapolation is only used in the veloicyt/pressure test functions
      // Depending on the time order the extrapolation is used: i.e.,
      // BDF1: constant extrapolation
      // BDF2: linear extrapolation
      if(db["time_discretization"].is("backward_euler") || current_step_ <=1)
      {
        if(db["disctype"].is("supg"))
        {
          functions.resize(6);
          // these will be used in the assembling process: 
          // for example in the assembly of right hand side for the 
          // pressure (u^n-1/dt, \nabla q): for the BDF2 scheme this 
          // will be 2 u^n-1 - 1/2 u^n-2 ::
          // Also this old solution will be used to compute the residual
          // in the method residual_based_vms
          functions[2] = s.u_m1.GetComponent(0);
          functions[3] = s.u_m1.GetComponent(1);
          
          s.combined_old_sols.reset();
          s.combined_old_sols = s.solution_m1;
          
          functions[4] = s.comb_old_u.GetComponent(0);
          functions[5] = s.comb_old_u.GetComponent(1);
        }
        if(db["disctype"].is("residual_based_vms"))
        {
          // for the residual in the nonlinear term
          // constant or linear extrapolation of the solution
          // is used
          s.extrapolate_sol.reset();
          s.extrapolate_sol = s.solution_m1;
          functions.resize(9);
          functions[2] = s.extrapolate_u.GetComponent(0);
          functions[3] = s.extrapolate_u.GetComponent(1);
          functions[4] = &s.p_old;
          
          //For the residual based vms scheme fill this
          //with the time derivative which is used in the residual 
          //computation 
          s.combined_old_sols.reset();
          s.combined_old_sols = s.solution;
          s.combined_old_sols.add_scaled(s.solution_m1, -1.);
          s.combined_old_sols.scale(1./tau);
          
          functions[5] = s.comb_old_u.GetComponent(0);
          functions[6] = s.comb_old_u.GetComponent(1);
          
          // also the solution from previous time steps
          // is used for the computation of right-hand
          // side fro pressure term
          s.combined_old_sols.reset();
          s.combined_old_sols = s.solution_m1;
          functions[7] = s.comb_old_u.GetComponent(0);
          functions[8] = s.comb_old_u.GetComponent(1);
        }
      }
      else if(db["time_discretization"].is("bdf_two"))// BDF2 case 
      {
        functions.resize(6);
        if(db["disctype"].is("supg"))
        {
          // 2 u^{n-1} - u^{n-2}
          s.extrapolate_sol.reset();
          s.extrapolate_sol = s.solution_m1;
          s.extrapolate_sol.scale(2.);
          s.extrapolate_sol.add_scaled(s.solution_m2, -1.);
          functions[2] = s.extrapolate_u.GetComponent(0);
          functions[3] = s.extrapolate_u.GetComponent(1);
          
          // now in the case of equal-order: 
          // linear combination of previous solution is 
          // necessary for the pressure part in rhs
          // it's devided by the time step within the 
          // local assemble routine
          s.combined_old_sols.reset();
          s.combined_old_sols = s.solution_m1;
          s.combined_old_sols.scale(2.);
          s.combined_old_sols.add_scaled(s.solution_m2, -0.5);
          functions[4] = s.comb_old_u.GetComponent(0);
          functions[5] = s.comb_old_u.GetComponent(1);
        }
        else if(db["disctype"].is("residual_based_vms"))
        {
          // for the residual in the nonlinear term
          // constant or linear extrapolation of the solution
          // is used
          functions.resize(9);
          
          s.extrapolate_sol.reset();
          s.extrapolate_sol = s.solution_m1;
          s.extrapolate_sol.scale(2.);
          s.extrapolate_sol.add_scaled(s.solution_m2, -0.1);
          
          functions[2] = s.extrapolate_u.GetComponent(0);
          functions[3] = s.extrapolate_u.GetComponent(1);
          functions[4] = &s.p_old;
          
          //For the residual based vms scheme fill this
          //with the time derivative which is used in the residual 
          //computation 
          s.combined_old_sols.reset();
          s.combined_old_sols = s.solution_m1;
          s.combined_old_sols.scale(1./tau);
          s.combined_old_sols.add_scaled(s.solution_m2, -.5/tau);
          
          functions[5] = s.comb_old_u.GetComponent(0);
          functions[6] = s.comb_old_u.GetComponent(1);
          
          // also the solution from previous time steps
          // is used for the computation of right-hand
          // side fro pressure term
          s.combined_old_sols.reset();
          s.combined_old_sols = s.solution_m1;
          s.combined_old_sols.scale(2.);
          s.combined_old_sols.add_scaled(s.solution_m2, -0.5);
          functions[7] = s.comb_old_u.GetComponent(0);
          functions[8] = s.comb_old_u.GetComponent(1);
        }        
      }
    }
    else //
    {
      ErrThrow("Implement the linear extrapolation");
    }
  }
  else
  {
    ErrThrow("not tested so far");
  }
}

void Time_NSE2D_BDF::set_matrices_rhs(Time_NSE2D_BDF::System_per_grid& s, 
 LocalAssembling2D_type type, std::vector< TSquareMatrix2D* >& sqMat, 
 std::vector< TMatrix2D* >& reMat, std::vector< double* >& rhs_array)
{
  // common in all type
  std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix.get_blocks_uniquely();
  std::vector<std::shared_ptr<FEMatrix>> Mass_SUPG
        = s.Mass_Matrix.get_blocks_uniquely();
  std::vector<std::shared_ptr<FEMatrix>> Mass_RBVMS
        = s.MatrixK.get_blocks_uniquely();
  int nstype = TDatabase::ParamDB->NSTYPE;
  switch(type){
    case TNSE2D:
      switch(nstype){
        case 4:
          if(db["disctype"].is("galerkin"))
          {

          }
          else if(db["disctype"].is("supg"))
          {
            // rhs-array fill
            rhs_array.resize(2);
            rhs_array[0] = s.rhs.block(0);
            rhs_array[1] = s.rhs.block(1);
            s.rhs.reset();
            // 4 A blocks and a Mass Matrix
            sqMat.resize(5);
            sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());  
            sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
            
            sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(Mass_SUPG.at(0).get());
            
            reMat.resize(4);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
            reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());   
          }
          else if(db["disctype"].is("residual_based_vms"))
          {
            // rhs-array fill
            rhs_array.resize(2);
            rhs_array[0] = s.rhs.block(0);
            rhs_array[1] = s.rhs.block(1);
            s.rhs.reset();
            // 4 A blocks and a Mass Matrix
            sqMat.resize(8);
            sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());  
            sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
            
            sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(0).get());
            sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(1).get());
            sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(3).get());  
            sqMat[7] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(4).get());
            
            reMat.resize(4);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
            reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());  
          }
          else
          {
            ErrThrow("disctype: ", db["disctype"], " for nstype ", 
                     nstype, " is not supported");
          }
          break;
        case 14:
          if(db["disctype"].is("galerkin"))
          {
          }
          else if(db["disctype"].is("supg"))
          {
            // rhs-array fill
            rhs_array.resize(3);
            rhs_array[0] = s.rhs.block(0);
            rhs_array[1] = s.rhs.block(1);
            rhs_array[2] = s.rhs.block(2);
            s.rhs.reset();
            // 4 A blocks and a Mass Matrix
            sqMat.resize(6);
            sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());  
            sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
            
            sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(Mass_SUPG.at(0).get());
            sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
            
            reMat.resize(4);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
            reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());   
          }
          else if(db["disctype"].is("residual_based_vms"))
          {
            // rhs-array fill
            rhs_array.resize(3);
            rhs_array[0] = s.rhs.block(0);
            rhs_array[1] = s.rhs.block(1);
            rhs_array[2] = s.rhs.block(2);
            s.rhs.reset();
            // 4 A blocks and a Mass Matrix
            sqMat.resize(9);
            sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());  
            sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
            
            sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(0).get());
            sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(1).get());
            sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(3).get());  
            sqMat[7] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(4).get());
            
            sqMat[8] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
            
            reMat.resize(4);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
            reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());  
          }
          else
          {
            ErrThrow("disctype: ", db["disctype"], " for nstype ", 
                     nstype, " is not supported");
          }
          break; // case 14 
      } // switch nstype
      break; // case TNSE2D
//-----------------
    case TNSE2D_NL: 
      switch(nstype){
        case 4:
          if(db["disctype"].is("galerkin"))
          {
          }
          else if(db["disctype"].is("supg"))
          {
            if(db["ansatz_test_extrapolate"].is("only_velocity_test"))
            {
              sqMat.resize(2);
              sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
              sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
              rhs_array.resize(0);
            }
            else
            {
              ErrThrow("fully version is already implemented: only needs to be merged");
            }
          }
          else if(db["disctype"].is("residual_based_vms"))
          {
            if(db["ansatz_test_extrapolate"].is("only_velocity_test"))
            {
              sqMat.resize(4);
              sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
              sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
              sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
              sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
              rhs_array.resize(0);
            }
            else
            {
              ErrThrow("not supported yet");
            }
          }
          else
          {
            ErrThrow("disctype: ", db["disctype"], " for nstype ", 
                     nstype, " is not supported");
          }
          break;
        case 14:
          if(db["disctype"].is("galerkin"))
          {
          }
          else if(db["disctype"].is("supg"))
          {
            if(db["ansatz_test_extrapolate"].is("only_velocity_test"))
            {
              sqMat.resize(2);
              sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
              sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
              reMat.resize(2);
              reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); 
              reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            
              rhs_array.resize(0);
            }
            else
            {
              ErrThrow("fully version is already implemented: only needs to be merged");
            }
          }
          else if(db["disctype"].is("residual_based_vms"))
          {
            if(db["ansatz_test_extrapolate"].is("only_velocity_test"))
            {
              sqMat.resize(4);
              sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
              sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
              sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
              sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
              
              reMat.resize(2);
              reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get());
              reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
              rhs_array.resize(0);
            }
            else
            {
              ErrThrow("not supported yet");
            }
          }
          else
          {
            ErrThrow("disctype: ", db["disctype"], " for nstype ", 
                     nstype, " is not supported");
          }
          break; // case 14 
      } // switch nstype
      break; // case TNSE2D_NL
//-----------------
    case TNSE2D_Rhs: 
      switch(nstype){
        case 4:
          rhs_array.resize(2);
          rhs_array[0] = s.rhs.block(0);
          rhs_array[1] = s.rhs.block(1);
          s.rhs.reset();
          
          if(db["disctype"].is("galerkin"))
          {
          }
          else if(db["disctype"].is("supg"))
          {
            // set extra matrices to be assembled together with
            // right-hand side: 
            sqMat.resize(1);
            sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(Mass_SUPG.at(0).get());
            
            reMat.resize(4);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
            reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());   
          }
          else if(db["disctype"].is("residual_based_vms"))
          {
            // rhs-array fill
            rhs_array.resize(2);
            rhs_array[0] = s.rhs.block(0);
            rhs_array[1] = s.rhs.block(1);
            s.rhs.reset();
            // 4 A blocks and a Mass Matrix
            sqMat.resize(4);
            sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(3).get());  
            sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(4).get());
            
            reMat.resize(4);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); 
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); 
            reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());  
          }
          else
          {
            ErrThrow("disctype: ", db["disctype"], " for nstype ", 
                     nstype, " is not supported");
          }
          break;
        case 14:
          rhs_array.resize(3);
          rhs_array[0] = s.rhs.block(0);
          rhs_array[1] = s.rhs.block(1);
          rhs_array[2] = s.rhs.block(2);
          s.rhs.reset();
          
          if(db["disctype"].is("galerkin"))
          {
          }
          else if(db["disctype"].is("supg"))
          {
            // set extra matrices to be assembled together with
            // right-hand side: 
            sqMat.resize(2);
            sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(Mass_SUPG.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
            
            reMat.resize(4);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
            reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());   
          }
          else if(db["disctype"].is("residual_based_vms"))
          {
            // set extra matrices to be assembled together with
            // right-hand side: 
            sqMat.resize(5);
            sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(0).get());
            sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(1).get());
            sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(3).get());  
            sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(Mass_RBVMS.at(4).get());
            
            sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
            
            reMat.resize(4);
            reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
            reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
            reMat[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
            reMat[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());   
          }
          else
          {
            ErrThrow("disctype: ", db["disctype"], " for nstype ", 
                     nstype, " is not supported");
          }
          break; // case 14 
      } // switch nstype
      break; // case TNSE2D_Rhs    
  }
  for(auto sm : sqMat)
    sm->reset();
  for(auto rm : reMat)
    rm->reset();
  
}

bool Time_NSE2D_BDF::imex_scheme(bool print_info)
{
  return false;
}

void Time_NSE2D_BDF::output(int m)
{
  bool no_output = !db["output_write_vtk"] && !db["output_compute_errors"];
  if(no_output)
    return;

  System_per_grid& s = this->systems.front();
  TFEFunction2D * u1 = s.u.GetComponent(0);
  TFEFunction2D * u2 = s.u.GetComponent(1);

  if((size_t)db["verbosity"]> 1)
  {
    u1->PrintMinMax();
    u2->PrintMinMax();
    s.p.PrintMinMax();
  }

  if(db["output_compute_errors"])
  {
    double locerr[8];
    MultiIndex2D allderiv[3]= {D00, D10, D01};
    const TFESpace2D *v_sp = &this->get_velocity_space();
    const TFESpace2D *p_sp = &this->get_pressure_space();
    TAuxParam2D aux;
    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
    
    u1->GetErrors(example.get_exact(0), 3, allderiv, 2, L2H1Errors,nullptr,
                  &aux,1, &v_sp,locerr);
    
    u2->GetErrors(example.get_exact(1), 3, allderiv, 2, L2H1Errors,nullptr,
                  &aux,1, &v_sp,locerr+2);

    errors[0] += (locerr[0]*locerr[0]+locerr[2]*locerr[2] 
                  + this->errors[1])*tau*0.5;
    errors[1] = locerr[0]*locerr[0]+locerr[2]*locerr[2];
    errors[2] += (locerr[1]*locerr[1]+locerr[3]*locerr[3] 
                  + this->errors[3])*tau*0.5;
    errors[3] = locerr[1]*locerr[1]+locerr[3]*locerr[3];  

    Output::print<1>("L2(u) : ", setprecision(10), sqrt(this->errors[1]));
    Output::print<1>("H1-semi(u) : ", setprecision(10), sqrt(this->errors[3]));

    Output::print<1>("L2(0,t,L2(u)) : ", sqrt(this->errors[0]));
    Output::print<1>("L2(0,t,H1-semi(u)) : ", sqrt(this->errors[2]));

    s.p.GetErrors(example.get_exact(2), 3, allderiv, 2, L2H1Errors, 
                  nullptr, &aux, 1, &p_sp, locerr);

    Output::print<1>("L2(p) : ", setprecision(10), locerr[0]);
    Output::print<1>("H1-semi(p)) : " , setprecision(10), locerr[1] );

    errors[4] += (locerr[0]*locerr[0] + this->errors[5])*tau*0.5;
    errors[5] = locerr[0]*locerr[0];
    Output::print<1>("L2(0,t,L2(p)) : ", sqrt(errors[4]) );
    
    errors[6] += (locerr[1]*locerr[1] + this->errors[7])*tau*0.5;
    errors[7] = locerr[1]*locerr[1];
    Output::print<1>("L2(0,t,H1-semi(p)) : ", sqrt(errors[6]) );
  }
  delete u1;
  delete u2;

  // example.do_post_processing(*this);
  
  if((m==0) || (m/TDatabase::TimeDB->STEPS_PER_IMAGE) )
  {
    if(db["output_write_vtk"])
    {
      outputWriter.write(TDatabase::TimeDB->CURRENTTIME);
    }
  }
}

void Time_NSE2D_BDF::output_problem_size_info() const
{
  int n_u = this->get_velocity_space().GetN_DegreesOfFreedom();
  int n_u_active = this->get_velocity_space().GetN_ActiveDegrees();
  int n_p = this->get_pressure_space().GetN_DegreesOfFreedom();
  int n_dof = 2 * n_u + n_p; // total number of degrees of freedom
  
  double h_min, h_max;
  TCollection * coll = this->get_velocity_space().GetCollection();
  coll->GetHminHmax(&h_min, &h_max);
  Output::stat("NSE2D", "Mesh data and problem size");
  Output::dash("cells              :  ", setw(10), coll->GetN_Cells());
  Output::dash("h (min, max)       :  ", setw(10), h_min, setw(10), " ", h_max);
  Output::dash("dof velocity       :  ", setw(10), 2*n_u );
  Output::dash("dof velocity active:  ", setw(10), 2*n_u_active);
  Output::dash("dof pressure       :  ", setw(10), n_p);
  Output::dash("dof all            :  ", setw(10), n_dof);
}

void Time_NSE2D_BDF::modify_slip_bc(bool BT_Mass)
{
  // modification of the matrices due to the 
  // slip type boundary conditions: If the mass matrices,
  // the off-diagonal A-blocks , and the BT's block,
  // are unchanged during the time iteration, then this modification
  // is done only once in the time loop. However, in the SUPG
  // and residual based VMS method these matrices are also 
  // updated during the time steps, so modification of all 
  // of them including the right-hand side is necessary.The 
  // modification of the diagonal A-blocks are necessary 
  // in any case.
  if(TDatabase::ParamDB->NSTYPE < 4)
  {
    ErrThrow("Slip with friction b.c. is only implemented for NSTYPE 4");
  }
  std::vector<const TFESpace2D*> spaces_mat(1);
  std::vector<double*> rhs_array(2);
  std::vector<const TFESpace2D*> rhs_space(2);
  
  for(System_per_grid& s: this->systems)
  {
    spaces_mat[0] = &s.velocity_space;
    rhs_space[0] = spaces_mat[0];
    rhs_space[1] = spaces_mat[0];
    
    rhs_array[0] = s.rhs.block(0);
    rhs_array[1] = s.rhs.block(1);
    
    std::vector<std::shared_ptr<FEMatrix>> blocks
         = s.matrix.get_blocks_uniquely();
    
    std::vector<std::shared_ptr<FEMatrix>> mass_blocks;
    
    mass_blocks = s.Mass_Matrix.get_blocks_uniquely();
    if(db["disctype"].is("residual_based_vms"))
    {
      mass_blocks = s.MatrixK.get_blocks_uniquely();
    }
    
    std::vector<const BoundCondFunct2D*> bc(3);
    bc[0]=s.velocity_space.GetBoundCondition();
    bc[1]=bc[0];
    bc[2]=s.pressure_space.GetBoundCondition();
    // boundary values:
    std::vector<BoundValueFunct2D*>bv(3);
    bv[0]=example.get_bd(0);
    bv[1]=example.get_bd(1);
    bv[2]=example.get_bd(2);
    
    std::vector<TSquareMatrix2D*> sqMat;
    std::vector<TMatrix2D*> reMat;
    sqMat.resize(5);
    sqMat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
    sqMat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
    sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());  
    sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());

     if(db["disctype"].is("residual_based_vms"))
     {
       sqMat.resize(8);
       sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
       sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(1).get());
       sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(3).get());  
       sqMat[7] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());
     }
     if(BT_Mass)
       sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());

     reMat.resize(0);
     if(BT_Mass)
     {
       reMat.resize(2);    
       reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
       reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
     }
     // update the matrices and right hand side
     Assemble2DSlipBC(spaces_mat.size(), spaces_mat.data(), 
                   sqMat.size(), sqMat.data(), reMat.size(), reMat.data(), 
                   rhs_array.size(), rhs_array.data(), rhs_space.data(), 
                   bc.data(), bv.data(), s.u.GetComponent(0), s.u.GetComponent(1));
  }

}
