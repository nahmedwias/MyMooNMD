#include <Saddle_point_preconditioner.h>
#ifdef __2D__
#include <LocalAssembling2D.h>
#include <Assemble2D.h>
#else
#include <LocalAssembling3D.h>
#include <Assemble3D.h>
#endif
#include <MainUtilities.h>
#include <MooNMD_Io.h>
#include <algorithm>

#include <Joint.h>
#include <BoundEdge.h>
#include <BoundFace.h>

/* ************************************************************************** */
Saddle_point_preconditioner::Saddle_point_preconditioner(
  const BlockFEMatrix& m, type t,  bool direct_velocity_solve)
  : spp_type(t), lsc_strategy(direct_velocity_solve ? 0 : 1), M(&m),
    velocity_block(), gradient_block(nullptr), divergence_block(nullptr),
    velocity_solver(nullptr), inverse_diagonal(), 
    velocity_space(&m.get_row_space(0)), pressure_space(nullptr),
    damping_factor(1.0), Poisson_solver_matrix(nullptr),
    Poisson_solver(nullptr), up_star(m), bdryCorrectionMatrix_(),
    poissonMatrixBdry_(nullptr), poissonSolverBdry_(nullptr)
{
  Output::print<3>("constructing a Saddle_point_preconditioner");
  if(lsc_strategy > 0)
  {
    Output::print("WARNING: solving systems within LSC using some iterative "
                  "routine requires a flexible solver.");
  }
  bool pressure_correction_in_matrix = this->M->pressure_correction_enabled();
  this->M->disable_pressure_correction();
  
  // number of block columns and rows
  unsigned int n_rows = this->M->get_n_cell_rows();
  unsigned int n_cols = this->M->get_n_cell_columns();
  pressure_space = &m.get_row_space(n_rows - 1);
  if(n_rows < 2 || n_cols != n_rows)
    ErrThrow("can not create a Saddle_point_preconditioner with this matrix");
  if(n_rows == 2 && this->spp_type == Saddle_point_preconditioner::type::bd_lsc)
    ErrThrow("boundary corrected LSC only available for (Navier-) Stokes type "
             "problems. This seems to be a Darcy type problem.\nIt is unclear "
             "how the boundary correction has to be implemented for H(div) "
             "elements");
  
  //velocity block K
  // take all blocks except from last row and last column
  this->velocity_block = M->get_sub_blockfematrix(0, n_rows-2);
  this->fill_inverse_diagonal();
  
  {
    // velocity solver database
    ParameterDatabase vs_db = Solver<>::default_solver_database();
    if(lsc_strategy == 0)
    {
      // use umfpack direct solver
      vs_db["solver_type"] = "direct";
      vs_db["direct_solver_type"] = "umfpack";
    }
    else
    {
      vs_db["solver_type"] = "iterative";
      vs_db["iterative_solver_type"] = "fgmres";
      vs_db["preconditioner"] = "jacobi";
      vs_db["max_n_iterations"] = 100;
      vs_db["residual_tolerance"] = 1.0e-15; // hardly ever reached
      vs_db["residual_reduction"] = 0.01;    // the actual stopping criterion
      vs_db["gmres_restart"] = 10;
      vs_db["damping_factor"] = 1.0; // no damping
    }
    this->velocity_solver.reset(new Solver<BlockFEMatrix, BlockVector>(vs_db));
    this->velocity_solver->update_matrix(this->velocity_block);
  }
  
  //gradient block B^T
  // take last column without last block (which would be pressure-pressure)
  this->gradient_block = this->M->get_combined_submatrix({0,n_cols-1}, 
                                                         {n_rows-2,n_cols-1});
  //divergence block B
  this->divergence_block = this->M->get_combined_submatrix({n_rows-1,0}, 
                                                           {n_rows-1,n_cols-2});
  
  if(this->spp_type == Saddle_point_preconditioner::type::simple)
  {
    // scale the gradient block B^T with the approximation of the mass-matrix
    this->gradient_block->scale(&inverse_diagonal[0], true);
  }
  // construct an approximation to the Schur complement matrix
  this->Poisson_solver_matrix = this->compute_Poisson_solver_matrix();
  this->Poisson_solver.reset(new DirectSolver(*this->Poisson_solver_matrix,
                             DirectSolver::DirectSolverTypes::umfpack));
  u_star = BlockVector(velocity_block);
  p_star = BlockVector(*this->Poisson_solver_matrix);
  u_tmp = BlockVector(u_star);
  p_tmp = BlockVector(p_star);
  
  if(this->spp_type == Saddle_point_preconditioner::type::bd_lsc)
  {
    //boundary corrected LSC
    //construct correctionMatrixBdry (anew...) and fill with ones
    bdryCorrectionMatrix_ = std::vector<double>(this->inverse_diagonal.size(),
                                                1.);
    //set up the bdryCorrectionMatrix_
    computeBdryCorrectionMatrix(m);
    //set up the poissonMatrixBdry_
    computePoissonMatrixBdry();
    //and wrap the latter into a solver
    //poissonMatrixBdry_.PrintFull("pMB");
    poissonSolverBdry_.reset(
      new DirectSolver(*poissonMatrixBdry_,
                       DirectSolver::DirectSolverTypes::umfpack));
  }
  
  if(pressure_correction_in_matrix)
    this->M->enable_pressure_correction();
}

/* ************************************************************************** */
Saddle_point_preconditioner::Saddle_point_preconditioner(
  const BlockMatrix& m, type t,  bool direct_velocity_solve)
{
  ErrThrow("Creating a Saddle_point_preconditioner with a BlockMatrix is not "
           "possible, you need a BlockFEMatrix.");
  // otherwise the get_combined_matrix methods are possibly not giving the 
  // correct behavior.
}

/* ************************************************************************** */
void Saddle_point_preconditioner::update(const BlockFEMatrix & m)
{
  if(this->M != &m)
    ErrThrow("cannot update Saddle_point_preconditioner with different matrix");
  
  // we assume the velocity block has changed but not the other blocks!
  // number of block columns and rows
  unsigned int n_rows = this->M->get_n_cell_rows();
  //unsigned int n_cols = this->M->get_n_cell_columns();
  
  // velocity block K, delete it first because it might have changed, then
  // create a new one
  // take all blocks except from last row and last column
  this->velocity_block = M->get_sub_blockfematrix(0, n_rows-2);
  
  this->velocity_solver->update_matrix(this->velocity_block);
  
  // we assume the blocks involving pressure did not change
  if(this->spp_type == Saddle_point_preconditioner::type::simple)
  { // SIMPLE
    // this is not well implemented, we need the method 
    // TMatrix* TMatrix::multiply(TMatrix*, double)
    // to take the product structure, this way we could speed things up here.
    // Instead all members are deleted and recreated
    
    // descale the gradient_block:
    unsigned int n_diagonal_entries = inverse_diagonal.size();
    for(unsigned int d = 0; d < n_diagonal_entries; ++d)
      this->inverse_diagonal[d] = 1.0 / this->inverse_diagonal[d];
    this->gradient_block->scale(&inverse_diagonal[0], true);
    
    // recompute the inverse_diagonal (of velocity_block)
    this->fill_inverse_diagonal();
    
    //scale the gradient block B^T with the inverse_diagonal
    this->gradient_block->scale(&inverse_diagonal[0], true);
    
    //recompute the Poisson_solver_matrix and then create new Poisson_solver
    this->Poisson_solver_matrix = this->compute_Poisson_solver_matrix();
    this->Poisson_solver.reset(
      new DirectSolver(*this->Poisson_solver_matrix,
                       DirectSolver::DirectSolverTypes::umfpack));
  }
  else if(this->spp_type == Saddle_point_preconditioner::type::lsc 
          || this->spp_type == Saddle_point_preconditioner::type::bd_lsc)
  {
    // nothing more needs to be done here
    
    // members which do not change and therefore do not need to be recomputed:
    // - gradient and divergence block
    // - inverse_diagonal 
    // - Poisson_solver_matrix
    // - Poisson_solver
    // - velocity_space, pressure_space
  }
}

/* ************************************************************************** */
void Saddle_point_preconditioner::update(const BlockMatrix& m)
{
  ErrThrow("Updating a Saddle_point_preconditioner with a BlockMatrix is not "
           "possible, you need a BlockFEMatrix.");
  // otherwise the get_combined_matrix methods are possibly not giving the 
  // correct behavior.
}

/* ************************************************************************** */
void solve_velocity(std::shared_ptr<Solver<BlockFEMatrix>> velocity_solver, 
                    const BlockVector& rhs, BlockVector& sol)
{
  unsigned int verbosity = Output::getVerbosity();
  Output::setVerbosity(1);
  velocity_solver->solve(rhs, sol);
  Output::setVerbosity(verbosity);
}

/* ************************************************************************** */
void Saddle_point_preconditioner::apply(const BlockVector &z,
                                        BlockVector &r) const
{
  r = z;
  Output::print<5>("Saddle_point_preconditioner::solve");
  if(r.n_blocks() == 0)
  {
    // only the standard constructor has been called for r so far
    r.copy_structure(z); // copy the structure of z, all entries are zero
  }
  else
  {
    r.reset(); // set all entries to zero
  }
  up_star = 0.;
  
  switch(this->spp_type)
  {
    case Saddle_point_preconditioner::type::simple:
    {
      //Output::print<5>("SIMPLE Saddle_point_preconditioner::solve");
      // up_star includes du* and dp*
      // ----------------------------------------------------------------------
      // solve A du* = r_u
      u_tmp = z.get_entries();
      solve_velocity(this->velocity_solver, u_tmp, u_star);
      
      // ----------------------------------------------------------------------
      // solve ^S dp* = r_p - B du*
      // use BlockVector r to store the right hand side for the system 
      // involving the Schur approximation
      unsigned int n_blocks = z.n_blocks();
      p_tmp = z.block(n_blocks - 1);
      // do r_p += -B * du*
      
      this->divergence_block->multiply(u_star.get_entries(), 
                                       p_tmp.get_entries(), -1.);
      this->Poisson_solver->solve(p_tmp, p_star);
      
      // ----------------------------------------------------------------------
      // dp = dp*
      r.copy(p_star.get_entries(), n_blocks - 1);
      // du = du* - D^-1 B^T dp*
      // there are zeros in the u-block of r currently
      // store the result of -B^T dp* in u-block of r
      gradient_block->multiply(p_star.get_entries(), r.get_entries(), -1.);
      // scale with D^-1, loop over all velocity entries
      for(unsigned int i = 0, n_v = gradient_block->GetN_Rows(); i < n_v; ++i)
      {
        r[i] *= this->inverse_diagonal[i];
      }
      // add du*
      for(unsigned int i = 0; i < n_blocks - 1; ++i)
        r.add(u_star.block(i), i);
      
      // ----------------------------------------------------------------------
      // update r = z + omega (dp, du)
      r *= this->damping_factor;
      r += z;
      break;
    }
    case Saddle_point_preconditioner::type::lsc:
    case Saddle_point_preconditioner::type::bd_lsc:
    {
      // a short name to tell if the correction at the boundary should be done
      bool correct_boundary = 
        this->spp_type == Saddle_point_preconditioner::type::bd_lsc;
      //if(correct_boundary) 
      //  Output::print<5>("bd_LSC Saddle_point_preconditioner::solve");
      //else
      //  Output::print<5>("LSC Saddle_point_preconditioner::solve");
      // -----------------------------------------------------------------------
      // The cases lsc and bd_lsc only differ in Step 1 and substep 2.2!
      unsigned int n_blocks = z.n_blocks();
      
      // Step 1: Poisson solver: solve S_p p_star = p_tmp = z_p, where 
      // S_p = B Q^(-1) B^T
      p_tmp = z.block(n_blocks - 1);
      
      if(correct_boundary)
        this->poissonSolverBdry_->solve(p_tmp, p_star);
      else
        this->Poisson_solver->solve(p_tmp, p_star);
      
      // Step 2. Update p_tmp = -B Q^(-1) K Q^(-1) B^T p_star
      
      // Substep 2.1. Compute B^T p_star and store it in u_tmp
      u_tmp = 0.;
      gradient_block->multiply(p_star.get_entries(), u_tmp.get_entries(), -1.0);
      
      // Substep 2.2. Compute Q^(-1) B^T p_star and store it in u_tmp
      if(correct_boundary)
      {
        for(unsigned int i = 0, n_v = gradient_block->GetN_Rows(); i < n_v; ++i)
        {
          u_tmp[i] *= this->bdryCorrectionMatrix_[i];
        }
      }
      else
      {
        for(unsigned int i = 0, n_v = gradient_block->GetN_Rows(); i < n_v; ++i)
        {
          u_tmp[i] *= this->inverse_diagonal[i];
        }
      }
      
      
      // Substep 2.3. Compute  K Q^(-1) B^T p_star and store it in u_star
      velocity_block.apply(u_tmp, u_star);
      
      
      //Substep 2.4. Compute Q^(-1) K Q^(-1) B^T p_star and store it in u_star
      for(unsigned int i = 0, n_v = gradient_block->GetN_Rows(); i < n_v; ++i)
      {
        u_star[i] *= this->inverse_diagonal[i];
      }
      
      // Subset 2.5. Compute B Q^(-1) K Q^(-1) B^T p_star and store it in p_tmp
      p_tmp.reset();
      divergence_block->multiply(u_star.get_entries(), p_tmp.get_entries(), 1.);
      
      //Step 3: 2. solver solve S_p p_star = p_tmp and store it in p_star
      Poisson_solver->solve(p_tmp, p_star);
      
      //Step 4. Update u_star = z_u - B^T p_star
      u_star = z.get_entries();
      gradient_block->multiply(p_star.get_entries(), u_star.get_entries(), -1.);
      
      //Step 5. Solve K u_tmp = u_star
      solve_velocity(this->velocity_solver, u_star, u_tmp);
      
      r.reset();
      // copy u_tmp and p_star back to r
      for(unsigned int i = 0; i < n_blocks - 1; ++i)
        r.add(u_tmp.block(i), i);
      r.copy(p_star.block(0), n_blocks-1);
      /// TODO damping? this->damping_factor
       
      //IntoL20FEFunction(r.block(n_blocks-1), r.length(n_blocks-1),
      //                  pressure_space,
      //                  TDatabase::ParamDB->VELOCITY_SPACE,
      //                  TDatabase::ParamDB->PRESSURE_SPACE);
      break;
    }
    default:
      ErrThrow("unknown saddle point preconditioner type");
      break;
  }
}

/* ************************************************************************** */
/* LSC === */
std::shared_ptr<BlockMatrix> 
Saddle_point_preconditioner::compute_Poisson_solver_matrix() const
{
  std::shared_ptr<TMatrix> ret(nullptr);
  // get the matrix, put it into a shared_ptr, essentially this is the matrix
  // which should be returned
  switch(this->spp_type)
  {
    case Saddle_point_preconditioner::type::simple:
      ret.reset(this->divergence_block->multiply(this->gradient_block.get(),
                                                 1.));
      break;
    case Saddle_point_preconditioner::type::lsc: // original LSC
    case Saddle_point_preconditioner::type::bd_lsc: // bdry corrected LSC
      ret.reset(
        divergence_block->multiply_with_transpose_from_right(inverse_diagonal));
      break;
    default:
      ErrThrow("unknown preconditioner for saddle point problems ");
      break;
  }
  // put the shared_ptr into a vector
  std::vector<std::shared_ptr<TMatrix>> ret_as_vector{ret};
  // create a BlockMatrix and return it
  return std::make_shared<BlockMatrix>(1, 1, ret_as_vector);
}
/* === LSC */

/* ************************************************************************** */
// local assembling routine to assemble a mass velocity matrix
// This needs to go into some assembling class!
bool has_vector_valued_basis_functions;
void local_assembling_velocity_mass(double Mult, double *coeff, double *param,
                                    double hK, double **OrigValues, 
                                    int *N_BaseFuncts, double ***LocMatrices, 
                                    double **LocRhs)
{
  // assemble (u,v) at a specific quadrature point in a specific cell
  double ** MatrixM = LocMatrices[0];
  const double * u_values = OrigValues[0]; // the values of u (and v)
  
  const int N_U = N_BaseFuncts[0];
  for(int i = 0; i < N_U; i++)
  {
    double *MatrixM_row  = MatrixM[i];
    double test00 = u_values[i];
    
    for(int j = 0; j < N_U; j++)
    {
      double ansatz00 = u_values[j];
      MatrixM_row[j] += Mult*(ansatz00*test00);
    }                            // endfor j
  }                              // endfor i
  if(has_vector_valued_basis_functions)
  {
    for(int i = 0; i < N_U; i++)
    {
      double *MatrixM_row  = MatrixM[i];
      double test00_y = u_values[i + N_U];
      
      for(int j = 0; j < N_U; j++)
      {
        double ansatz00_y = u_values[j+N_U];
        MatrixM_row[j] += Mult*(ansatz00_y*test00_y);
      }                            // endfor j
    }                              // endfor i
  }
}

void Saddle_point_preconditioner::fill_inverse_diagonal()
{
  // number of entries on diagonal in velocity_block
  size_t n_diagonal_entries = this->velocity_block.get_n_total_rows();
  this->inverse_diagonal.resize(n_diagonal_entries, 0.);
  
  //mass-matrix Q and its approximation
  if(this->spp_type == Saddle_point_preconditioner::type::simple)
  {
    // SIMPLE
    for(unsigned int d = 0; d < n_diagonal_entries; ++d)
    {
      this->inverse_diagonal[d] = 1.0 / this->velocity_block.get(d, d);
    }
  }
  else
  {
    // either lsc or bd_lsc
    // first: assemble a velocity mass matrix
    
#ifdef __2D__
    typedef MultiIndex2D MultiIndex;
    MultiIndex v = D00; // value, no derivatives
    typedef CoeffFct2D CoeffFct;
    typedef AssembleFctParam2D AssembleFctParam;
    typedef ManipulateFct2D ManipulateFct;
    typedef TFEFunction2D FEFunction;
    typedef LocalAssembling2D LocalAssembling;
    typedef TSquareMatrix2D SqMat;
    typedef TMatrix2D RectMat;
    typedef TFESpace2D FESpace;
    typedef BoundCondFunct2D BoundCondFunct;
    typedef BoundValueFunct2D BoundValueFunct;
#else
    typedef MultiIndex3D MultiIndex;
    MultiIndex v = D000; // value, no derivatives
    typedef CoeffFct3D CoeffFct;
    typedef AssembleFctParam3D AssembleFctParam;
    typedef ManipulateFct3D ManipulateFct;
    typedef TFEFunction3D FEFunction;
    typedef LocalAssembling3D LocalAssembling;
    typedef TSquareMatrix3D SqMat;
    typedef TMatrix3D RectMat;
    typedef TFESpace3D FESpace;
    typedef BoundCondFunct3D BoundCondFunct;
    typedef BoundValueFunct3D BoundValueFunct;
#endif
    
    // create an approriate LocalAssembling2D object:
    // We will assemble only one velocity block and then create a BlockMatrix
    // where this one block appears twice (three times in 3D)
    int n_terms = 1; // only 1 term (u,v)
    std::vector<MultiIndex> derivatives(1, v); // that one term is the value
    std::vector<int> fe_space_numbers(1, 0); // only 1 space with index 0
    std::vector<int> row_space(1, 0); // 1 matrix 
    std::vector<int> column_space(1, 0); // 1 matrix
    std::vector<int> rhs_space(1, 0); // not needed
    CoeffFct * coeff = nullptr; // no coefficients
    AssembleFctParam* local_assembling_function = 
      local_assembling_velocity_mass;
    ManipulateFct * manipulate = nullptr; // not needed
    int n_matrices = 1; // assemble only one matrix
    int n_rhs = 0; // no right hand side
    int n_parameter = 0; // not needed
    std::vector<ParamFct*> parameter_functions; // no parameter functions
    std::vector<int> begin_parameter; // not needed
    int n_parameters = 0; // not needed
    FEFunction ** fe_functions = nullptr;
    int n_fe_values = 0; // no other fe values are needed
    std::vector<int> fe_value_function_index; // not needed
    std::vector<MultiIndex> fe_value_multi_index; // not needed
    LocalAssembling la(n_terms, derivatives, fe_space_numbers, row_space, 
                       column_space, rhs_space, coeff, 
                       local_assembling_function, manipulate, n_matrices,
                       n_rhs, n_parameter, parameter_functions, 
                       begin_parameter, n_parameters, fe_functions, 
                       n_fe_values, fe_value_function_index, 
                       fe_value_multi_index);
    
    // prepare a call to Assemble2D
    int n_fe_spaces = 1;
    SqMat *sq_matrices[1];
    
    //FESpace* neumann_velocity_space = 
    //  this->velocity_space->copy_with_all_Neumann_boundary();
    //const FESpace* v_space = neumann_velocity_space;
    const FESpace* v_space = this->velocity_space;
    
    SqMat one_block_of_mass_matrix(v_space);
    
    sq_matrices[0] =  &one_block_of_mass_matrix;
    int n_rect_mat = 0; // only square matrices
    RectMat *rect_matrices[0];
    double **rhs = nullptr; // right hand sides
    const FESpace ** fe_spaces_rhs = nullptr;
    // which boundary conditions are correct here? We could as well use the ones
    // from the velocity_space. We use three even in 2D, where two would be ok.
    BoundCondFunct * boundary_conditions[3] = {
      BoundConditionNoBoundCondition, BoundConditionNoBoundCondition, 
      BoundConditionNoBoundCondition };
    BoundValueFunct* non_const_bound_values[3] = {
      BoundaryValueHomogenous, BoundaryValueHomogenous, 
      BoundaryValueHomogenous };
   
#ifdef __2D__
    if(this->velocity_space->get_fe(0).GetBaseFunct2D()->GetBaseVectDim() != 1)
      has_vector_valued_basis_functions = true;
    else
      has_vector_valued_basis_functions = false;
    
    Assemble2D(n_fe_spaces, &v_space, n_matrices, sq_matrices,
               n_rect_mat, rect_matrices, n_rhs, rhs, fe_spaces_rhs,
               boundary_conditions, non_const_bound_values, la);
    BlockFEMatrix mass_matrix({v_space, v_space});
    mass_matrix.replace_blocks(one_block_of_mass_matrix,
                               {{0,0}, {1,1}}, {false, false});
#else
    if(this->velocity_space->get_fe(0).GetBaseFunct3D()->GetBaseVectDim() != 1)
      has_vector_valued_basis_functions = true;
    else
      has_vector_valued_basis_functions = false;
    
    Assemble3D(n_fe_spaces, &v_space, n_matrices, sq_matrices,
               n_rect_mat, rect_matrices, n_rhs, rhs, fe_spaces_rhs,
               boundary_conditions, non_const_bound_values, la);
    BlockFEMatrix mass_matrix({v_space, v_space, v_space});
    mass_matrix.replace_blocks(one_block_of_mass_matrix,
                               {{0,0}, {1,1}, {2,2}}, {false, false, false});
#endif
    
    //number of entries on diagonal in mass matrix
    for(unsigned int d = 0; d < n_diagonal_entries; ++d)
    {
      this->inverse_diagonal[d] = 1.0 / mass_matrix.get(d, d);
    }
    //delete neumann_velocity_space;
  }
}


/* ************************************************************************** */
/* Extra methods for boundary corrected LSC */

/** Sets up bdryCorrectionMatrix_. */
#ifdef __2D__
void Saddle_point_preconditioner::computeBdryCorrectionMatrix(
    const BlockFEMatrix& m)
{
  //the problem parameter epsilon - chosen 0.1 as proposed by the authors
  double epsilon = 0.1;
  
  // just a simple check
  if(this->velocity_space != &m.get_row_space(0))
    ErrThrow("bd_lsc preconditioner: matrix has wrong fe space");
  
  //Hold a reference to the TCollection underlying the problem
  const TCollection& collection = *this->velocity_space->GetCollection();
  
  //loop over cells in the collection
  for(int iCells = 0; iCells < collection.GetN_Cells(); ++iCells)
  {
    //hold a reference to the current cell
    const TBaseCell& cell = *collection.GetCell(iCells);
    //loop over cell's joints
    for(int iJoints = 0; iJoints < cell.GetN_Joints(); ++iJoints)
    {
      //find out if this cell has a boundary joint
      if(!cell.GetJoint(iJoints)->InnerJoint())
      {
        //this is a joint on the boundary, so a TBoundEdge (in 2D)
        // but is it on a Dirichlet boundary part?
        const TBoundEdge& edge = *dynamic_cast<const TBoundEdge*>(cell.GetJoint(
            iJoints));
        //get ID of the boundary component -- looks like we can work with this 
        // component!
        int componentID = edge.GetBoundComp()->GetID();
        BoundCond componentType;
        // the parameter on the boundary component which refers to the center
        // of 'edge'
        double t_c = (edge.GetStartParameter() + edge.GetEndParameter()) / 2.;
        // the 0 for velo space, componentID for global ID of the boundary 
        // component,
        this->velocity_space->GetBoundCondition()(componentID, t_c, 
                                                  componentType);
        if(componentType == DIRICHLET)
        {
          //OutPut("Found a cell with a Dirchlet boundary edge." << endl);
          // we found a cell with an edge on the dirichlet bondary. time to 
          // start working:
          
          // fetch and store coordinates of beginning and ending of the edge
          //(ASSUMING THE EDGE IS A LINE (class TBdLine))
          double startX, startY, endX, endY;
          edge.GetBoundComp()->GetXYofT(edge.GetStartParameter(), startX,
                                        startY);
          edge.GetBoundComp()->GetXYofT(edge.GetEndParameter(), endX, endY);
          // calculate unit normal on the edge (orientation is of no interest 
          // here)
          double normalX = endY - startY;
          double normalY = startX - endX;
          double norm = sqrt(normalX * normalX + normalY * normalY);
          normalX /= norm;
          normalY /= norm;
          
          //those are the values to be written into the matrix
          double horizontalDofsCorrectionValue = std::max(fabs(normalX),
                                                          epsilon);
          double verticalDofsCorrectionValue = std::max(fabs(normalY),
                                                        epsilon);
          
          //now write the relevant matrix entries
          //loop over velocity dofs belonging to cell
          for(int index = velocity_space->GetBeginIndex()[iCells];
              index < velocity_space->GetBeginIndex()[iCells + 1]; ++index)
          {
            int iVeloDof = velocity_space->GetGlobalNumbers()[index];
            
            //here we rely on the correct ordering of the velo dofs!
            //this is a horizontal (x-direction) velo dof
            bdryCorrectionMatrix_[iVeloDof] = horizontalDofsCorrectionValue;
            //OutPut("Horizontal velo dof " << iVeloDof << " in Cell " << iCells
            //       << " put correction to " << horizontalDofsCorrectionValue
            //       << endl);
            //this is a vertical (y-direction) velo dof
            bdryCorrectionMatrix_.at(velocity_space->GetN_DegreesOfFreedom()
                + iVeloDof) = verticalDofsCorrectionValue;
          } //end loop over velocity dofs
        } // end if dirichlet boundary
      } //end if boundary joint
    } //end loop over joints
  } //end loop over cells
  
  //to store the correct matrix "H^{-1}= D*D_Q^{-1}" we have to scale with the
  // diagonal mass inverse
  std::transform(bdryCorrectionMatrix_.begin(), bdryCorrectionMatrix_.end(),
                 inverse_diagonal.begin(), bdryCorrectionMatrix_.begin(),
                 std::multiplies<double>());
}

#endif // __2D__
#ifdef __3D__
void Saddle_point_preconditioner::computeBdryCorrectionMatrix(
    const BlockFEMatrix& m)
{
  //the problem parameter epsilon - chosen 0.1 as proposed by the authors
  double epsilon = 0.1;
  
  // just a simple check
  if(this->velocity_space != &m.get_row_space(0))
    ErrThrow("bd_lsc preconditioner: matrix has wrong fe space");
  
  //Hold a reference to the TCollection underlying the problem
  const TCollection& collection = *this->velocity_space->GetCollection();
  
  //loop over cells in the collection
  for(int iCells = 0; iCells < collection.GetN_Cells(); ++iCells)
  {
    //hold a reference to the current cell
    TBaseCell& cell = *collection.GetCell(iCells);
    //loop over cell's joints
    for(int iJoints = 0; iJoints < cell.GetN_Joints(); ++iJoints)
    {
      
      //find out if this cell has a boundary joint
      if(!cell.GetJoint(iJoints)->InnerJoint())
      {
        // this joint is on the boundary, check if it is on a Dirichlet boundary
        BoundCond componentType;
        const int * face_vertex, *face_vertex_length;
        int m_l; // max_length
        cell.GetShapeDesc()->GetFaceVertex(face_vertex, face_vertex_length,m_l);
        // compute center of face in order to then evaluate the boundary 
        // condition there
        double x = 0, y = 0, z = 0;
        for (int v = 0; v < face_vertex_length[iJoints]; ++v)
        {
          TVertex* vertex = cell.GetVertex(face_vertex[iJoints*m_l+v]);
          x += vertex->GetX();
          y += vertex->GetY();
          z += vertex->GetZ();
        }
        x /= face_vertex_length[iJoints];
        y /= face_vertex_length[iJoints];
        z /= face_vertex_length[iJoints];
        // the 0 for velo space, 1 would be pressure
        this->velocity_space->getBoundCondition()(x, y, z, componentType);
        if(componentType == DIRICHLET)
        {
          //OutPut("Found a cell with a Dirchlet boundary edge." << endl);
          // we found a cell with an edge on  dirichlet bondary. time to start 
          // working:
          
          // find the normal of this face
          // Assuming this face is a PLANE or WALL (class TBdPlane or TBdWall)
          double normalX, normalY, normalZ;
          {
            TBoundFace& face =
                *dynamic_cast<TBoundFace*>(cell.GetJoint(iJoints));
            BoundTypes bt = face.GetBoundComp()->GetType();
            if(bt == Cylinder || bt == Sphere)
              ErrThrow("getting normal on cylinder or sphere possibly not "
                       "working");
            
            if(face_vertex_length[iJoints] < 3) 
              ErrThrow("a joint in 3D has less than three vertices???");
            double ux, uy, uz, vx, vy, vz;
            // helping doubles to evaluate coordinates for some vertices
            double X = 0, Y = 0, Z = 0;
            cell.GetVertex(face_vertex[iJoints * m_l])->GetCoords(X, Y, Z);
            ux = vx = -X;
            uy = vy = -Y;
            uz = vz = -Z;
            cell.GetVertex(face_vertex[iJoints * m_l + 1])->GetCoords(X, Y, Z);
            ux += X;
            uy += Y;
            uz += Z;
            cell.GetVertex(face_vertex[iJoints * m_l + 2])->GetCoords(X, Y, Z);
            vx += X;
            vy += Y;
            vz += Z;
            normalX = uy*vz - uz*vy;
            normalY = uz*vx - ux*vz;
            normalZ = ux*vy - uy*vx;
            double norm = sqrt(
                normalX * normalX + normalY * normalY + normalZ * normalZ);
            normalX /= norm;
            normalY /= norm;
            normalZ /= norm;
          }
          
          //those are the values to be written into the matrix
          double x_DofsCorrectionValue = std::max(fabs(normalX), epsilon);
          double y_DofsCorrectionValue = std::max(fabs(normalY), epsilon);
          double z_DofsCorrectionValue = std::max(fabs(normalZ), epsilon);
          
          //now write the relevant matrix entries
          //loop over velocity dofs belonging to cell
          for(int index = velocity_space->GetBeginIndex()[iCells];
              index < velocity_space->GetBeginIndex()[iCells + 1]; ++index)
          {
            int iVeloDof = velocity_space->GetGlobalNumbers()[index];
            
            //here we rely on the correct ordering of the velo dofs!
            //this is a horizontal (x-direction) velo dof
            bdryCorrectionMatrix_[iVeloDof] = x_DofsCorrectionValue;
            //Output::print("Horizontal velo dof ", iVeloDof, " in Cell ",
            //              iCells, " put correction to ",
            //              horizontalDofsCorrectionValue);
            //this is a vertical (y-direction) velo dof
            bdryCorrectionMatrix_[velocity_space->GetN_DegreesOfFreedom()
                + iVeloDof] = y_DofsCorrectionValue;
            bdryCorrectionMatrix_[2*velocity_space->GetN_DegreesOfFreedom()
                            + iVeloDof] = z_DofsCorrectionValue;
          } //end loop over velocity dofs
        } // end if dirichlet boundary
      } //end if boundary joint
    } //end loop over joints
  } //end loop over cells
  
  //to store the correct matrix "H^{-1}= D*D_Q^{-1}" we have to scale with the 
  // diagonal mass inverse
  std::transform(bdryCorrectionMatrix_.begin(), bdryCorrectionMatrix_.end(),
                 inverse_diagonal.begin(), bdryCorrectionMatrix_.begin(),
                 std::multiplies<double>());
}

#endif // 3D

void Saddle_point_preconditioner::computePoissonMatrixBdry()
{
  if(this->spp_type == Saddle_point_preconditioner::type::bd_lsc)
  {
    bool transpose;
    TMatrix* ret = divergence_block->multiply_with_transpose_from_right(
      bdryCorrectionMatrix_,
      Poisson_solver_matrix->get_block(0, 0, transpose)->GetStructure());
    poissonMatrixBdry_.reset(ret);
  }
  else
    ErrThrow("Call of Saddle_point_preconditioner::setUpPoissonSolverBdry() "
             "for wrong preconditioner.");
}

/* End extra methods for boundary corrected LSC */

