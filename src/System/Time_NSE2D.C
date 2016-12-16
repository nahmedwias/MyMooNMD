#include <Time_NSE2D.h>
#include <Database.h>
#include <Assemble2D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <GridTransfer.h>
#include <LocalAssembling2D.h>
#include <FEFunctionInterpolator.h>

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

  return db;
}
/* *************************************************************************** */

/**************************************************************************** */
Time_NSE2D::System_per_grid::System_per_grid(const Example_TimeNSE2D& example, 
                  TCollection& coll, std::pair< int, int > order, 
                  Time_NSE2D::Matrix type)
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
     solution.length(2))
{
  // Mass Matrix
  // Output::increaseVerbosity(5);

  Mass_Matrix = BlockFEMatrix::Mass_NSE2D(velocity_space);
      
  switch(type)
  {
    case Time_NSE2D::Matrix::Type1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
      break;
    case Time_NSE2D::Matrix::Type2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
      break;
    case Time_NSE2D::Matrix::Type3:
      matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
      break;
    case Time_NSE2D::Matrix::Type4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      break;
    case Time_NSE2D::Matrix::Type14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
      break;
  }
}

/**************************************************************************** */
Time_NSE2D::Time_NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
                       int reference_id)
  : Time_NSE2D(domain, param_db, Example_TimeNSE2D(param_db), reference_id)
{
  
}

/**************************************************************************** */
Time_NSE2D::Time_NSE2D(const TDomain& domain, const ParameterDatabase& param_db,
                       const Example_TimeNSE2D& ex, int reference_id)
 : db(get_default_TNSE2D_parameters()), outputWriter(param_db), systems(),
   example(ex), solver(param_db), defect(), oldResidual(0), 
   initial_residual(1e10), errors(10,0.), oldtau(0.0)
{
  db.merge(param_db);
  this->set_parameters();
  
  std::pair <int,int> velo_pres_order(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  this->get_velocity_pressure_orders(velo_pres_order);
  
  Time_NSE2D::Matrix type;
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
  
  if(!usingMultigrid)
  {
    // create the collection of cells from the domain (finest grid)
    TCollection *coll = domain.GetCollection(It_Finest, 0, reference_id);
    this->systems.emplace_back(example, *coll, velo_pres_order, type);
    
    TFEFunction2D * u1 = this->systems.front().u.GetComponent(0);
    TFEFunction2D * u2 = this->systems.front().u.GetComponent(1);
    
    u1->Interpolate(example.get_initial_cond(0));
    u2->Interpolate(example.get_initial_cond(1));
  }
  else
  {
    auto multigrid = this->solver.get_multigrid();
    
    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    size_t n_levels = multigrid->get_n_geometric_levels();
    int finest = domain.get_ref_level();
    int coarsest = finest - n_levels + 1;
    for (int grid_no = finest; grid_no >= coarsest; --grid_no)
    {
      TCollection *coll = domain.GetCollection(It_EQ, grid_no, reference_id);
      systems.emplace_back(example, *coll, velo_pres_order,
                            type);
      //prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    multigrid->initialize(matrices);
  }
  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);
  
  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  
  // print out the information (cells, dofs, etc)
  this->output_problem_size_info();
}

/**************************************************************************** */
void Time_NSE2D::set_parameters()
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
}

/**************************************************************************** */
void Time_NSE2D::get_velocity_pressure_orders(std::pair< int, int > &velo_pres_order)
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
    // discontinuous spaces
    case 1:case 2: case 3: case 4: case 5:
      pressure_order = -(velocity_order-1)*10;
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
  }
  TDatabase::ParamDB->PRESSURE_SPACE  = pressure_order;
  velo_pres_order.second = pressure_order;
  
  Output::print("velocity space", setw(10), velo_pres_order.first);
  Output::print("pressure space", setw(10), velo_pres_order.second);
}

/**************************************************************************** */
void Time_NSE2D::assemble_initial_time()
{
  for(System_per_grid& s : this->systems)
  {
    s.rhs.reset();
    const TFESpace2D * velo_space = &s.velocity_space;
    const TFESpace2D * pres_space = &s.pressure_space;
    // variables which are same for all nstypes
    size_t n_fe_spaces = 2;
    const TFESpace2D *fespmat[2] = {velo_space, pres_space};
    
    size_t n_square_matrices = 6; // 
    TSquareMatrix2D *sqMatrices[6]{nullptr}; // maximum number of square matrices
    
    size_t n_rect_matrices = 4; // maximum number of rectangular matrices
    TMatrix2D *rectMatrices[4]{nullptr}; // maximum number of pointers
    
    size_t nRhs = 2; //is 3 if NSE type is 4 or 14
    double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), nullptr}; //third place gets only filled
    const TFESpace2D *fe_rhs[3] = {velo_space, velo_space, nullptr};  // if NSE type is 4 or 14
    
    BoundCondFunct2D * boundary_conditions[3] = {
      velo_space->GetBoundCondition(), velo_space->GetBoundCondition(),
      pres_space->GetBoundCondition() };
      
    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    TFEFunction2D *fe_functions[3] =
      { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };

    LocalAssembling2D la(TNSE2D, fe_functions, 
                         this->example.get_coeffs());
    std::vector<std::shared_ptr<FEMatrix>> blocks 
         = s.matrix.get_blocks_uniquely();
    std::vector<std::shared_ptr<FEMatrix>> mass_blocks
         = s.Mass_Matrix.get_blocks_uniquely();
    
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        if(blocks.size() != 3)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        n_square_matrices = 2;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        // mass matrix 
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        // rectangular matrices
        n_rect_matrices = 2;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get());
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
        break;
      case 2:
        if(blocks.size() != 5)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        n_square_matrices = 2;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        // mass matrix
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        // rectangular matrices
        n_rect_matrices = 4;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(3).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(4).get());
        rectMatrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get()); //than the standing B blocks
        rectMatrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
        break;
      case 3:
        if(blocks.size() != 6)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        n_square_matrices = 5;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());        
        sqMatrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());        
        sqMatrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());        
        // mass matrices
        sqMatrices[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());        
        // ErrThrow("not tested yet!!!, kindly remove one mass matrix from the LocalAssembling2D routine");
        // rectangular matrices
        n_rect_matrices = 2;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
        break;
      case 4:
        if(blocks.size() != 8)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        n_square_matrices = 5;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());        
        sqMatrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());        
        sqMatrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());        
        // mass matrices
        sqMatrices[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        // rectangular matrices
        n_rect_matrices = 4;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rectMatrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rectMatrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

        RHSs[2] = s.rhs.block(2); // NSE type 4 includes pressure rhs
        fe_rhs[2] = pres_space;
        nRhs = 3;

        break;
      case 14:
        if(blocks.size() != 9)
        {
          ErrThrow("Wrong blocks.size() ", blocks.size());
        }
        n_square_matrices = 6;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sqMatrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sqMatrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        // mass matrices
        sqMatrices[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        // C block pressure pressure
        sqMatrices[5] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
        // rectangular matrices  
        n_rect_matrices = 4;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rectMatrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rectMatrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

        RHSs[2] = s.rhs.block(2); // NSE type 14 includes pressure rhs
        fe_rhs[2]  = pres_space;
        nRhs = 3;

        break;
      default:
        ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
               " That NSE Block Matrix Type is unknown to class Time_NSE2D.");
    }
    // assemble all the matrices and right hand side 
    Assemble2D(n_fe_spaces, fespmat, n_square_matrices, sqMatrices, 
               n_rect_matrices, rectMatrices, nRhs, RHSs, fe_rhs, 
               boundary_conditions, non_const_bound_values.data(), la);
    // copy nonactives
    s.solution.copy_nonactive(s.rhs);
  }
  
// this piece of code is just to check A matrix
//  BlockVector testing1 = this->systems.front().solution;
//  testing1 = 1;
//  const BlockVector testing2= testing1;
//  this->systems.front().matrix.apply(testing2,testing1);
//  testing1.write("vecteur_test_initial");


  // copy the current right hand side vector to the old_rhs 
  this->old_rhs = this->systems.front().rhs; 
  this->old_solution = this->systems.front().solution;
}

/**************************************************************************** */
void Time_NSE2D::assemble_rhs()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  const double theta2 = TDatabase::TimeDB->THETA2;
  const double theta3 = TDatabase::TimeDB->THETA3;
  const double theta4 = TDatabase::TimeDB->THETA4;
  
  System_per_grid& s = this->systems.front();
  // reset the right hand side
  s.rhs.reset();
  // assembling of the right hand side 
  TFEFunction2D *fe_functions[3] = 
  { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
  
  LocalAssembling2D la(TNSE2D_Rhs, fe_functions,
                           this->example.get_coeffs());
  
  int N_Rhs = 3;
  const TFESpace2D * v_space = &this->get_velocity_space();
  const TFESpace2D * p_space = &this->get_pressure_space();
  
  double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), s.rhs.block(2)};
  
  const TFESpace2D *fespmat[2] = {v_space, p_space};
  const TFESpace2D *fesprhs[3] = {v_space, v_space, p_space};
  
  BoundCondFunct2D * boundary_conditions[3] = {
             v_space->GetBoundCondition(), v_space->GetBoundCondition(), 
              p_space->GetBoundCondition() };
  
   std::array<BoundValueFunct2D*, 3> non_const_bound_values;
   non_const_bound_values[0] = this->example.get_bd(0);
   non_const_bound_values[1] = this->example.get_bd(1);
   non_const_bound_values[2] = this->example.get_bd(2);
  
   Assemble2D(1, fespmat, 0, nullptr,
              0, nullptr, N_Rhs, RHSs, fesprhs,
              boundary_conditions, non_const_bound_values.data(), la);
   // copy the non active to the solution vector
   // since the rhs vector will be passed to the solver
   // and is modified with matrix vector multiplication
   // which also uses the non-actives
   s.solution.copy_nonactive(s.rhs);
   
  // now it is this->systems[i].rhs = f^k
  // scale by time step length and theta4 (only active dofs)  
  s.rhs.scaleActive(tau*theta4);
  // add rhs from previous time step 
  if(theta3 != 0)
  {    
    s.rhs.addScaledActive((this->old_rhs), tau*theta3);
    
    // now it is this->systems[i].rhs = tau*theta3*f^{k-1} + tau*theta4*f^k
    // next we want to set old_rhs to f^k (to be used in the next time step)
    this->old_rhs.addScaledActive(s.rhs, -1./(tau*theta3));
    this->old_rhs.scaleActive(-theta3/theta4);
    this->old_rhs.copy_nonactive(s.rhs);
  }  
  // FIXME FInd other solution than this submatrix method.
  // M u^{k-1}
  s.Mass_Matrix.apply_scaled_submatrix(old_solution, s.rhs, 2, 2, 1.0);
  // -tau*theta2 * A u^{k-1}
  double factor = -tau*theta2;
  s.matrix.apply_scaled_submatrix(old_solution, s.rhs, 2, 2, factor);
  
  // scale the BT blocks with time step length
  for(System_per_grid& s : this->systems)
  {
    if(tau != oldtau)
    {
      // TODO: change the factor to be THETA1*tau;
      factor = /*TDatabase::TimeDB->THETA1**/tau;
      if(this->oldtau != 0.0)
      {
        factor /= this->oldtau;
        Output::print<1>("change in tau", this->oldtau, "->", tau);
      }
      // scale the BT transposed blocks with the current time step
      const std::vector<std::vector<size_t>> cell_positions = {{0,2}, {1,2}};
	s.matrix.scale_blocks(factor, cell_positions);      
      if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
      {
        const std::vector<std::vector<size_t>> cell_positions_t = {{2,0}, {2,1}};
	s.matrix.scale_blocks(factor, cell_positions_t);
      }
    }
  }
  this->oldtau = tau;
  // copy non active from solution into rhs vector
  s.rhs.copy_nonactive(s.solution);  

  Output::print<5>("assembled the system right hand side ");  
}

/**************************************************************************** */
void Time_NSE2D::assemble_system()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;
  
  for(System_per_grid& s : this->systems)
  {
    const std::vector<std::vector<size_t>>
      cell_positions = {{0,0}, {0,1}, {1, 0}, {1, 1}};
    // note: declaring the auxiliary cell_positions is needed by the compiler
    // to sort out the overriding of the function scale_blocks_actives(...,...)
    s.matrix.scale_blocks_actives(factor, cell_positions);
    const FEMatrix& mass_bloks = *s.Mass_Matrix.get_blocks().at(0).get();
    s.matrix.add_matrix_actives(mass_bloks, 1.0, {{0,0}, {1,1}}, {false, false});
  }
  Output::print<5>("Assembled the system matrix which will be passed to the ", 
                   "solver");
}

/**************************************************************************** */
void Time_NSE2D::assemble_nonlinear_term()
{
  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(systems.size() > 1)
  {
    for( int block = 0; block < 2 ;++block)
    {
      std::vector<const TFESpace2D*> spaces;
      std::vector<double*> u_entries;
      std::vector<size_t> u_ns_dofs;
      for(auto &s : systems )
      {
        spaces.push_back(&s.velocity_space);
        u_entries.push_back(s.solution.block(block));
        u_ns_dofs.push_back(s.solution.length(block));
      }
      GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
    }
  }
  
  for(System_per_grid& s : this->systems)
  {
    const TFESpace2D *velocity_space = &s.velocity_space;
    size_t n_fe_spaces = 1;
    const TFESpace2D *fespmat[1]={velocity_space};
    
    size_t n_square_matrices;
    TSquareMatrix2D* sqMatrices[2]{nullptr};
    
    size_t n_rect_matrices = 0;
    TMatrix2D** rectMatrices=nullptr;
    
    BoundCondFunct2D * boundary_conditions[1] 
       = {velocity_space->GetBoundCondition() };
    
     std::array<BoundValueFunct2D*, 3> non_const_bound_values;
     non_const_bound_values[0] = this->example.get_bd(0);
     non_const_bound_values[1] = this->example.get_bd(1);
     non_const_bound_values[2] = this->example.get_bd(2);
     
     
    TFEFunction2D *fe_functions[3] = 
      { s.u.GetComponent(0), s.u.GetComponent(1), &s.p };
    LocalAssembling2D la_nonlinear(TNSE2D_NL, fe_functions,
                                   this->example.get_coeffs());
    
    
    std::vector<std::shared_ptr<FEMatrix>> blocks 
         = s.matrix.get_blocks_uniquely({{0,0},{1,1}});
    
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        n_square_matrices = 1;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        break;
      case 3:
      case 4:
      case 14:
        n_square_matrices = 2;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        break;
      default:
        ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
               " That NSE Block Matrix Type is unknown to class Time_NSE2D.");
    }
    // reset matrices to zero
    for(size_t m=0; m<n_square_matrices; m++)
      sqMatrices[m]->reset();
    
    Assemble2D(n_fe_spaces, fespmat, n_square_matrices, sqMatrices,
               n_rect_matrices, rectMatrices, 0, nullptr, nullptr, 
               boundary_conditions, non_const_bound_values.data(),
               la_nonlinear);
  }
  Output::print<5>("Assembled the nonlinear matrix only ");
}

/**************************************************************************** */
bool Time_NSE2D::stopIte(unsigned int it_counter)
{//TODO This has no "slow convergence criterion yet!"
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
  
  Output::print("nonlinear step  :  " , setw(3), it_counter);
  Output::print("impulse_residual:  " , setw(3), impulse_residual);
  Output::print("mass_residual   :  " , setw(3), mass_residual);
  Output::print("residual        :  " , setw(3), sqrt(residual));
  
  if (it_counter>0)
  {
  Output::print("rate:           :  " , setw(3), sqrt(residual)/oldResidual);
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
  
  if ((((sqrt(residual)<=limit)||(it_counter==Max_It)))
   && (it_counter>=sc_minit))
   {
     Output::print("ITE : ", setw(3), it_counter, "  RES : ", sqrt(residual), 
                   " Reduction : ",  sqrt(residual)/initial_residual);
     // descale the matrices, since only the diagonal A block will 
     // be reassembled in the next time step
     this->deScaleMatrices();     
     return true;
   }
   else
     return false;
}

/**************************************************************************** */
void Time_NSE2D::solve()
{
  System_per_grid& s = this->systems.front();
  
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("multigrid solver is not tested yet")
  }
  solver.solve(s.matrix,s.rhs, s.solution);
  
  // Important: We have to descale the matrices, since they are scaled
  // before the solving process. Only A11 and A22 matrices are 
  // reset and assembled again but the A12 and A21 are scaled, so
  // for the next iteration we have to descale, see assemble_system()
  this->deScaleMatrices();

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       s.p.project_into_L20();

  this->old_solution = s.solution;
}

/**************************************************************************** */
void Time_NSE2D::deScaleMatrices()
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  double factor = tau*TDatabase::TimeDB->THETA1;  
  for(System_per_grid& s : this->systems)
  {
    const FEMatrix& mass_bloks = *s.Mass_Matrix.get_blocks().at(0).get();
    s.matrix.add_matrix_actives(mass_bloks, -1.0, {{0,0}, {1,1}}, {false, false});
    const std::vector<std::vector<size_t>>
      cell_positions = {{0,0}, {0,1}, {1, 0}, {1, 1}};
    // note: declaring the auxiliary cell_positions is needed by the compiler
    // to sort out the overriding of the function scale_blocks_actives(...,...)
    s.matrix.scale_blocks_actives(1./factor, cell_positions);
  }  
}

/**************************************************************************** */
void Time_NSE2D::output(int m)
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
  
  if((m==0) || (m/TDatabase::TimeDB->STEPS_PER_IMAGE) )
  {
    if(db["output_write_vtk"])
    {
      outputWriter.write(TDatabase::TimeDB->CURRENTTIME);
    }
  }
}

/**************************************************************************** */
void Time_NSE2D::output_problem_size_info() const
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

/**************************************************************************** */
std::array< double, int(6) > Time_NSE2D::get_errors()
{
  std::array<double, int(6)> error_at_time_points;
  error_at_time_points[0] = sqrt(errors[1]); // L2 velocity error
  error_at_time_points[1] = sqrt(errors[3]); // H1 velocity error
  error_at_time_points[2] = sqrt(errors[5]); // L2 pressure error
  error_at_time_points[3] = sqrt(errors[7]); // H1 pressure error
  
  return error_at_time_points;
}

/**************************************************************************** */










/* *********** BELOW THIS LINE USER SPECIFIC CODE **************/
/** ************************************************************************ */
void Time_NSE2D::assemble_initial_time_withfields(TFEFunction2D* rho_field,
                     TFEFunction2D* mu_field)
{
  for(System_per_grid& s : this->systems)
  {
    s.rhs.reset();
    const TFESpace2D * velo_space = &s.velocity_space;
    const TFESpace2D * pres_space = &s.pressure_space;
    // variables which are same for all nstypes
    size_t n_fe_spaces = 2;
    const TFESpace2D *fespmat[2] = {velo_space, pres_space};

    size_t n_square_matrices = 6; //
    TSquareMatrix2D *sqMatrices[6]{nullptr}; // maximum number of square matrices

    size_t n_rect_matrices = 4; // maximum number of rectangular matrices
    TMatrix2D *rectMatrices[4]{nullptr}; // maximum number of pointers

    size_t nRhs = 2; //is 3 if NSE type is 4 or 14
    double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), nullptr}; //third place gets only filled
    const TFESpace2D *fe_rhs[3] = {velo_space, velo_space, nullptr};  // if NSE type is 4 or 14

    BoundCondFunct2D * boundary_conditions[3] = {velo_space->GetBoundCondition(),
                                                 velo_space->GetBoundCondition(),
                                                 pres_space->GetBoundCondition()};

    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    non_const_bound_values[2] = example.get_bd()[2];

    TFEFunction2D *fe_functions[5] =
    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p, nullptr, nullptr };

    LocalAssembling2D la(TNSE2D, fe_functions,
                         this->example.get_coeffs());

    // the following should add fluid property fluids
    // to obtain a dimensional formulation of NSE
    if (rho_field != nullptr && mu_field != nullptr)
    {
      // HERE IS THE CODE TO SET UP THE RHO AND MU FIELDS FOR LOCAL ASSEMBLING OBJECT

      // we assume rho and mu fields are not too exotic,
      // and have the same space...
      int Ndof_rho     = rho_field->GetFESpace2D()->GetN_DegreesOfFreedom();
      int Ndof_velocity = velo_space->GetN_DegreesOfFreedom();

      // step 1: check if the velocity space and "rho_field->GetSFESpace2d()"
      // are the same. If yes, no interpolation needed
      // otherwise, do the interpolation
      if (Ndof_rho == Ndof_velocity) //this is a simple check condition...
      {
        Output::info<1>("Time_NSE2D", "The spaces of the velocity field and the "
                   "scalar field (rho) are the same ==> There will be no interpolation.");
//        Output::print<3>("Degres of freedoms of scalar field rho= " , Ndof_rho);
//        Output::print<3>("Degres of freedoms of velocity   = " , Ndof_velocity);

        //fill up the new fe function array
        fe_functions[3] = rho_field;
        fe_functions[4] = mu_field;
      }
      else  // do the interpolation
      {
        Output::warn<1>("Time_NSE2D", "The spaces of the velocity field and the "
                        "scalar field are not the same ==> Starting interpolation... this may take some time");

        // set up an interpolator object  (ptr will be shared later)
        // we interpolate into velo_space
        FEFunctionInterpolator interpolator(velo_space);

        // length of the values array of the interpolated rho
        // velo must equal length of the velocity components
        size_t length_interpolated = s.u.GetComponent(0)->GetLength();

        std::vector<double> temporary_rho(length_interpolated, 0.0);
        std::vector<double> temporary_mu(length_interpolated, 0.0);

        this->entries_rho_scalar_field = temporary_rho;
        this->entries_mu_scalar_field = temporary_mu;

        TFEFunction2D interpolated_rho_field =
            interpolator.interpolate(*rho_field,
                                     this->entries_rho_scalar_field);

        TFEFunction2D interpolated_mu_field =
            interpolator.interpolate(*mu_field,
                                     this->entries_mu_scalar_field);

        //fill up the new fe function array
        fe_functions[3] = &interpolated_rho_field;
        fe_functions[4] = &interpolated_mu_field;

      }

      // step 2 - set all the 'parameter'-related values in la accordingly
      // set up the input...
      la.setBeginParameter({0});
      la.setFeFunctions2D(fe_functions); //reset - now velo comp included
      la.setFeValueFctIndex({0,1,3,4});
      la.setFeValueMultiIndex({D00,D00,D00,D00});
      la.setN_Parameters(4);
      la.setN_FeValues(4);
      la.setN_ParamFct(1);
      la.setParameterFct_string("TimeNSParamsVelo_dimensional");

      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          la.setAssembleParam_string("TimeNSType1Galerkin_dimensional");  //this is for dimensional NSE
          break;
        case 3:
          switch(TDatabase::ParamDB->LAPLACETYPE)
          {
            case 1:
              // Assembling routine for NSType 3 with DD
              la.setAssembleParam_string("TimeNSType3GalerkinDD_dimensional");
              break;
            default:
              ErrThrow("NSTYPE 3 works only with LAPLACETYPE 1, please correct input parameter!");
          }
          break;
        default:
          ErrThrow("Time_NSE2D: NSType 2 and 4 are not implemented yet for the Dimensional NSE!");
      }
      //...this should do the trick
    }


    std::vector<std::shared_ptr<FEMatrix>> blocks
    = s.matrix.get_blocks_uniquely();
    std::vector<std::shared_ptr<FEMatrix>> mass_blocks
    = s.Mass_Matrix.get_blocks_uniquely();

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        if(blocks.size() != 3)
        {ErrThrow("Wrong blocks.size() ", blocks.size());}
        n_square_matrices = 2;
        sqMatrices[0]   = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sqMatrices[1]   = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        n_rect_matrices = 2;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get());
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
        break;
      case 2:
        if(blocks.size() != 5)
        {ErrThrow("Wrong blocks.size() ", blocks.size());}
        n_square_matrices = 2;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        n_rect_matrices = 4;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(3).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(4).get());
        rectMatrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(1).get()); //than the standing B blocks
        rectMatrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get());
        break;
      case 3:
        if(blocks.size() != 6)
        {ErrThrow("Wrong blocks.size() ", blocks.size());}
        n_square_matrices = 5;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sqMatrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sqMatrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        sqMatrices[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        // ErrThrow("not tested yet!!!, kindly remove one mass matrix from the LocalAssembling2D routine");
        n_rect_matrices = 2;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
        break;
      case 4:
        if(blocks.size() != 8)
        {ErrThrow("Wrong blocks.size() ", blocks.size());}
        n_square_matrices = 5;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sqMatrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sqMatrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        sqMatrices[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        n_rect_matrices = 4;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rectMatrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rectMatrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

        RHSs[2] = s.rhs.block(2); // NSE type 4 includes pressure rhs
        fe_rhs[2] = pres_space;
        nRhs = 3;

        break;
      case 14:
        if(blocks.size() != 9)
        {ErrThrow("Wrong blocks.size() ", blocks.size());}
        n_square_matrices = 6;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        sqMatrices[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
        sqMatrices[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(4).get());
        sqMatrices[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
        // C block pressure pressure
        sqMatrices[5] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(8).get());
        n_rect_matrices = 4;
        rectMatrices[0] = reinterpret_cast<TMatrix2D*>(blocks.at(6).get()); //first the lying B blocks
        rectMatrices[1] = reinterpret_cast<TMatrix2D*>(blocks.at(7).get());
        rectMatrices[2] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //than the standing B blocks
        rectMatrices[3] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());

        RHSs[2] = s.rhs.block(2); // NSE type 14 includes pressure rhs
        fe_rhs[2]  = pres_space;
        nRhs = 3;

        break;
      default:
        ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
                 " That NSE Block Matrix Type is unknown to class Time_NSE2D.");
    }
    // assemble all the matrices and right hand side
    Assemble2D(n_fe_spaces, fespmat, n_square_matrices, sqMatrices,
               n_rect_matrices, rectMatrices, nRhs, RHSs, fe_rhs,
               boundary_conditions, non_const_bound_values.data(), la);
    // copy nonactives
    s.solution.copy_nonactive(s.rhs);
  }

// this piece of code is just to check A matrix
//  BlockVector testing1 = this->systems.front().solution;
//  testing1 = 1;
//  const BlockVector testing2= testing1;
//  this->systems.front().matrix.apply(testing2,testing1);
//  testing1.write("vecteur_test_DIMENSIONAL");

  // copy the current right hand side vector to the old_rhs
  this->old_rhs = this->systems.front().rhs;
  this->old_solution = this->systems.front().solution;
}



void Time_NSE2D::assemble_nonlinear_term_withfields(TFEFunction2D* rho_field,
                                                    TFEFunction2D* mu_field)
{
  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(systems.size() > 1)
  {
    for( int block = 0; block < 2 ;++block)
    {
      std::vector<const TFESpace2D*> spaces;
      std::vector<double*> u_entries;
      std::vector<size_t> u_ns_dofs;
      for(auto &s : systems )
      {
        spaces.push_back(&s.velocity_space);
        u_entries.push_back(s.solution.block(block));
        u_ns_dofs.push_back(s.solution.length(block));
      }
      GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
    }
  }

  for(System_per_grid& s : this->systems)
  {
    const TFESpace2D *velocity_space = &s.velocity_space;
    size_t n_fe_spaces = 1;
    const TFESpace2D *fespmat[1]={velocity_space};

    size_t n_square_matrices;
    TSquareMatrix2D* sqMatrices[2]{nullptr};

    size_t n_rect_matrices = 0;
    TMatrix2D** rectMatrices=nullptr;

    BoundCondFunct2D * boundary_conditions[1]
                                           = {velocity_space->GetBoundCondition() };

    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = this->example.get_bd(0);
    non_const_bound_values[1] = this->example.get_bd(1);
    non_const_bound_values[2] = this->example.get_bd(2);


    TFEFunction2D *fe_functions[5] =
    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p, nullptr, nullptr };
    LocalAssembling2D la_nonlinear(TNSE2D_NL, fe_functions,
                                   this->example.get_coeffs());


    // the following should add fluid property fluids
    // to obtain a dimensional formulation of NSE
    if (rho_field != nullptr && mu_field != nullptr)
    {
      // HERE IS THE CODE TO SET UP THE RHO AND MU FIELDS FOR LOCAL ASSEMBLING OBJECT

      // we assume rho and mu fields are not too exotic,
      // and have the same space...
      int Ndof_rho     = rho_field->GetFESpace2D()->GetN_DegreesOfFreedom();
      int Ndof_velocity = velocity_space->GetN_DegreesOfFreedom();

      // step 1: check if the velocity space and "rho_field->GetSFESpace2d()"
      // are the same. If yes, no interpolation needed
      // otherwise, do the interpolation
      if (Ndof_rho == Ndof_velocity) //this is a simple check condition...
      {
        Output::info<1>("Time_NSE2D", "The spaces of the velocity field and the "
                        "scalar field (rho) are the same ==> There will be no interpolation.");
//        Output::print<3>("Degres of freedoms of scalar field rho= " , Ndof_rho);
//        Output::print<3>("Degres of freedoms of velocity   = " , Ndof_velocity);

        //fill up the new fe function array
        fe_functions[3] = rho_field;
        fe_functions[4] = mu_field;
      }
      else  // do the interpolation
      {
        Output::warn<1>("Time_NSE2D", "The spaces of the velocity field and the "
                        "scalar field are not the same ==> Starting interpolation... this may take some time");

        // set up an interpolator object  (ptr will be shared later)
        // we interpolate into velo_space
        FEFunctionInterpolator interpolator(velocity_space);

        // length of the values array of the interpolated rho
        // velo must equal length of the velocity components
        size_t length_interpolated = s.u.GetComponent(0)->GetLength();

        std::vector<double> temporary_rho(length_interpolated, 0.0);
        std::vector<double> temporary_mu(length_interpolated, 0.0);

        this->entries_rho_scalar_field = temporary_rho;
        this->entries_mu_scalar_field = temporary_mu;

        TFEFunction2D interpolated_rho_field =
            interpolator.interpolate(*rho_field,
                                     this->entries_rho_scalar_field);

        TFEFunction2D interpolated_mu_field =
            interpolator.interpolate(*mu_field,
                                     this->entries_mu_scalar_field);

        //fill up the new fe function array
        fe_functions[3] = &interpolated_rho_field;
        fe_functions[4] = &interpolated_mu_field;

      }

      // step 2 - set all the 'parameter'-related values in la accordingly
      // set up the input...
      la_nonlinear.setBeginParameter({0});
      la_nonlinear.setFeFunctions2D(fe_functions); //reset - now velo comp included
      la_nonlinear.setFeValueFctIndex({0,1,3,4});
      la_nonlinear.setFeValueMultiIndex({D00,D00,D00,D00});
      la_nonlinear.setN_Parameters(4);
      la_nonlinear.setN_FeValues(4);
      la_nonlinear.setN_ParamFct(1);
      la_nonlinear.setParameterFct_string("TimeNSParamsVelo_dimensional");

      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          la_nonlinear.setAssembleParam_string("TimeNSType1_2NLGalerkin_dimensional");  //this is for dimensional NSE
          break;
        case 3:
          switch(TDatabase::ParamDB->LAPLACETYPE)
          {
            case 1:
              // Assembling routine for NSType 3 with DD
              la_nonlinear.setAssembleParam_string("TimeNSType3_4NLGalerkinDD_dimensional");
              break;
            default:
              ErrThrow("NSTYPE 3 works only with LAPLACETYPE 1, please correct input parameter!");
          }
          break;
            default:
              ErrThrow("Time_NSE2D: NSType 2 and 4 are not implemented yet for the Dimensional NSE!");
      }
      //...this should do the trick
    }


    std::vector<std::shared_ptr<FEMatrix>> blocks
    = s.matrix.get_blocks_uniquely({{0,0},{1,1}});

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        n_square_matrices = 1;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
        break;
      case 3:
      case 4:
      case 14:
        n_square_matrices = 2;
        sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());

        sqMatrices[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
        break;
      default:
        ErrThrow("TDatabase::ParamDB->NSTYPE = ", TDatabase::ParamDB->NSTYPE ,
                 " That NSE Block Matrix Type is unknown to class Time_NSE2D.");
    }
    // reset matrices to zero
    for(size_t m=0; m<n_square_matrices; m++)
      sqMatrices[m]->reset();

    Assemble2D(n_fe_spaces, fespmat, n_square_matrices, sqMatrices,
               n_rect_matrices, rectMatrices, 0, nullptr, nullptr,
               boundary_conditions, non_const_bound_values.data(),
               la_nonlinear);
  }
  Output::print<5>("Assembled the nonlinear matrix only ");
}



void Time_NSE2D::assemble_rhs_withfields(TFEFunction2D* rho_field,
                                         TFEFunction2D* mu_field)
{
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  const double theta2 = TDatabase::TimeDB->THETA2;
  const double theta3 = TDatabase::TimeDB->THETA3;
  const double theta4 = TDatabase::TimeDB->THETA4;

  System_per_grid& s = this->systems.front();
  const TFESpace2D * velo_space = &s.velocity_space;
  // reset the right hand side
  s.rhs.reset();
  // assembling of the right hand side
  TFEFunction2D *fe_functions[5] =
  { s.u.GetComponent(0), s.u.GetComponent(1), &s.p, nullptr, nullptr };

  LocalAssembling2D la(TNSE2D_Rhs, fe_functions,
                       this->example.get_coeffs());


  // the following should add fluid property fluids
  // to obtain a dimensional formulation of NSE
  if (rho_field != nullptr && mu_field != nullptr)
  {
    // HERE IS THE CODE TO SET UP THE RHO AND MU FIELDS FOR LOCAL ASSEMBLING OBJECT

    // we assume rho and mu fields are not too exotic,
    // and have the same space...
    int Ndof_rho     = rho_field->GetFESpace2D()->GetN_DegreesOfFreedom();
    int Ndof_velocity = s.velocity_space.GetN_DegreesOfFreedom();

    // step 1: check if the velocity space and "rho_field->GetSFESpace2d()"
    // are the same. If yes, no interpolation needed
    // otherwise, do the interpolation
    if (Ndof_rho == Ndof_velocity) //this is a simple check condition...
    {
      Output::info<1>("Time_NSE2D", "The spaces of the velocity field and the "
                 "scalar field (rho) are the same ==> There will be no interpolation.");
//        Output::print<3>("Degres of freedoms of scalar field rho= " , Ndof_rho);
//        Output::print<3>("Degres of freedoms of velocity   = " , Ndof_velocity);

      //fill up the new fe function array
      fe_functions[3] = rho_field;
      fe_functions[4] = mu_field;
    }
    else  // do the interpolation
    {
      Output::warn<1>("Time_NSE2D", "The spaces of the velocity field and the "
                      "scalar field are not the same ==> Starting interpolation... this may take some time");

      // set up an interpolator object  (ptr will be shared later)
      // we interpolate into velo_space
      FEFunctionInterpolator interpolator(velo_space);

      // length of the values array of the interpolated rho
      // velo must equal length of the velocity components
      size_t length_interpolated = s.u.GetComponent(0)->GetLength();

      std::vector<double> temporary_rho(length_interpolated, 0.0);
      std::vector<double> temporary_mu(length_interpolated, 0.0);

      this->entries_rho_scalar_field = temporary_rho;
      this->entries_mu_scalar_field = temporary_mu;

      TFEFunction2D interpolated_rho_field =
          interpolator.interpolate(*rho_field,
                                   this->entries_rho_scalar_field);

      TFEFunction2D interpolated_mu_field =
          interpolator.interpolate(*mu_field,
                                   this->entries_mu_scalar_field);

      //fill up the new fe function array
      fe_functions[3] = &interpolated_rho_field;
      fe_functions[4] = &interpolated_mu_field;

    }

    // step 2 - set all the 'parameter'-related values in la accordingly
    // set up the input...
    la.setBeginParameter({0});
    la.setFeFunctions2D(fe_functions); //reset - now velo comp included
    la.setFeValueFctIndex({0,1,3,4});
    la.setFeValueMultiIndex({D00,D00,D00,D00});
    la.setN_Parameters(4);
    la.setN_FeValues(4);
    la.setN_ParamFct(1);
    la.setParameterFct_string("TimeNSParamsVelo_dimensional");

    la.setAssembleParam_string("TimeNSRHS_dimensional");  //this is for dimensional NSE
    //...this should do the trick
  }



  int N_Rhs = 3;
  const TFESpace2D * v_space = &this->get_velocity_space();
  const TFESpace2D * p_space = &this->get_pressure_space();

  double *RHSs[3] = {s.rhs.block(0), s.rhs.block(1), s.rhs.block(2)};

  const TFESpace2D *fespmat[2] = {v_space, p_space};
  const TFESpace2D *fesprhs[3] = {v_space, v_space, p_space};

  BoundCondFunct2D * boundary_conditions[3] = {
                                               v_space->GetBoundCondition(), v_space->GetBoundCondition(),
                                               p_space->GetBoundCondition() };

  std::array<BoundValueFunct2D*, 3> non_const_bound_values;
  non_const_bound_values[0] = this->example.get_bd(0);
  non_const_bound_values[1] = this->example.get_bd(1);
  non_const_bound_values[2] = this->example.get_bd(2);

  Assemble2D(1, fespmat, 0, nullptr,
             0, nullptr, N_Rhs, RHSs, fesprhs,
             boundary_conditions, non_const_bound_values.data(), la);
  // copy the non active to the solution vector
  // since the rhs vector will be passed to the solver
  // and is modified with matrix vector multiplication
  // which also uses the non-actives
  s.solution.copy_nonactive(s.rhs);

  // now it is this->systems[i].rhs = f^k
  // scale by time step length and theta4 (only active dofs)
  s.rhs.scaleActive(tau*theta4);
  // add rhs from previous time step
  if(theta3 != 0)
  {
    s.rhs.addScaledActive((this->old_rhs), tau*theta3);

    // now it is this->systems[i].rhs = tau*theta3*f^{k-1} + tau*theta4*f^k
    // next we want to set old_rhs to f^k (to be used in the next time step)
    this->old_rhs.addScaledActive(s.rhs, -1./(tau*theta3));
    this->old_rhs.scaleActive(-theta3/theta4);
    this->old_rhs.copy_nonactive(s.rhs);
  }
  // FIXME FInd other solution than this submatrix method.
  // M u^{k-1}
  s.Mass_Matrix.apply_scaled_submatrix(old_solution, s.rhs, 2, 2, 1.0);
  // -tau*theta2 * A u^{k-1}
  double factor = -tau*theta2;
  s.matrix.apply_scaled_submatrix(old_solution, s.rhs, 2, 2, factor);

  // scale the BT blocks with time step length
  for(System_per_grid& s : this->systems)
  {
    if(tau != oldtau)
    {
      // TODO: change the factor to be THETA1*tau;
      factor = /*TDatabase::TimeDB->THETA1**/tau;
      if(this->oldtau != 0.0)
      {
        factor /= this->oldtau;
        Output::print<1>("change in tau", this->oldtau, "->", tau);
      }
      // scale the BT transposed blocks with the current time step
      const std::vector<std::vector<size_t>> cell_positions = {{0,2}, {1,2}};
      s.matrix.scale_blocks(factor, cell_positions);
      if(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
      {
        const std::vector<std::vector<size_t>> cell_positions_t = {{2,0}, {2,1}};
        s.matrix.scale_blocks(factor, cell_positions_t);
      }
    }
  }
  this->oldtau = tau;
  // copy non active from solution into rhs vector
  s.rhs.copy_nonactive(s.solution);

  Output::print<5>("assembled the system right hand side ");

}


void Time_NSE2D::assemble_massmatrix_withfields(TFEFunction2D* rho_field)
{
  for(System_per_grid& s : this->systems)
  {
    const TFESpace2D *velocity_space = &s.velocity_space;

    BoundCondFunct2D * boundary_conditions[1]
                                           = {velocity_space->GetBoundCondition() };

    std::array<BoundValueFunct2D*, 3> non_const_bound_values;
    non_const_bound_values[0] = this->example.get_bd(0);
    non_const_bound_values[1] = this->example.get_bd(1);
    non_const_bound_values[2] = this->example.get_bd(2);

    TFEFunction2D *fe_functions[4] =
    { s.u.GetComponent(0), s.u.GetComponent(1), &s.p, nullptr };

    LocalAssembling2D la_mass(TNSE2D_Mass, fe_functions,
                                   this->example.get_coeffs());

    // the following should add fluid property fluids
    // to obtain a dimensional formulation of NSE
    if (rho_field != nullptr)
    {
      // HERE IS THE CODE TO SET UP THE RHO AND MU FIELDS FOR LOCAL ASSEMBLING OBJECT

      // we assume rho and mu fields are not too exotic,
      // and have the same space...
      int Ndof_rho     = rho_field->GetFESpace2D()->GetN_DegreesOfFreedom();
      int Ndof_velocity = velocity_space->GetN_DegreesOfFreedom();

      // step 1: check if the velocity space and "rho_field->GetSFESpace2d()"
      // are the same. If yes, no interpolation needed
      // otherwise, do the interpolation
      if (Ndof_rho == Ndof_velocity) //this is a simple check condition...
      {
        Output::info<1>("Time_NSE2D", "The spaces of the velocity field and the "
                        "scalar field (rho) are the same ==> There will be no interpolation.");
        //        Output::print<3>("Degres of freedoms of scalar field rho= " , Ndof_rho);
        //        Output::print<3>("Degres of freedoms of velocity   = " , Ndof_velocity);

        //fill up the new fe function array
        fe_functions[3] = rho_field;
      }
      else  // do the interpolation
      {
        Output::warn<1>("Time_NSE2D", "The spaces of the velocity field and the "
                        "scalar field are not the same ==> Starting interpolation... this may take some time");

        // set up an interpolator object  (ptr will be shared later)
        // we interpolate into velo_space
        FEFunctionInterpolator interpolator(velocity_space);

        // length of the values array of the interpolated rho
        // velo must equal length of the velocity components
        size_t length_interpolated = s.u.GetComponent(0)->GetLength();

        std::vector<double> temporary_rho(length_interpolated, 0.0);

        this->entries_rho_scalar_field = temporary_rho;

        TFEFunction2D interpolated_rho_field =
            interpolator.interpolate(*rho_field,
                                     this->entries_rho_scalar_field);

        //fill up the new fe function array
        fe_functions[3] = &interpolated_rho_field;
      }

      // step 2 - set all the 'parameter'-related values in la_mass accordingly
      // set up the input...
      la_mass.setBeginParameter({0});
      la_mass.setFeFunctions2D(fe_functions); //reset - now velo comp included
      la_mass.setFeValueFctIndex({0,1,3});
      la_mass.setFeValueMultiIndex({D00,D00,D00});
      la_mass.setN_Parameters(3);
      la_mass.setN_FeValues(3);
      la_mass.setN_ParamFct(1);
      la_mass.setParameterFct_string("TimeNSParamsVelo_dimensional");

      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          // the following line is normally done in the above la_mass constructor
          //la_mass.setAssembleParam_string("TimeNSType1GalerkinMass_dimensional"); //this is for dimensional NSE
          break;
        case 3:
          switch(TDatabase::ParamDB->LAPLACETYPE)
          {
            case 1:
              // Assembling routine for NSType 3 with DD
              la_mass.setAssembleParam_string("TimeNSType3GalerkinDDMass_dimensional");
              break;
            default:
              ErrThrow("NSTYPE 3 works only with LAPLACETYPE 1, please correct input parameter!");
          }
          break;
            default:
              ErrThrow("Time_NSE2D: NSType 2 and 4 are not implemented yet for the Dimensional NSE!");
      }
      //...this should do the trick
    }
    else
    {
      Output::warn<1>("TNSE2D", "Using this method without a valid rho_field "
                      "won't give correct results...");
      ErrThrow("Impossible to use this method, see warning message...");
    }

    size_t n_fe_spaces = 1;
    const TFESpace2D *fespmat[1]={velocity_space};

    size_t n_square_matrices = 1;
    TSquareMatrix2D* sqMatrices[1];

    std::vector<std::shared_ptr<FEMatrix>> mass_blocks
         = s.Mass_Matrix.get_blocks_uniquely();

    sqMatrices[0] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());
    sqMatrices[0]->reset();

    Assemble2D(n_fe_spaces, fespmat, n_square_matrices, sqMatrices,
               0, nullptr, 0, nullptr, nullptr,
               boundary_conditions, non_const_bound_values.data(),
               la_mass);
  }
}

