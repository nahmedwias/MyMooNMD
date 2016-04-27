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
Saddle_point_preconditioner::Saddle_point_preconditioner(const BlockFEMatrix& m,
                                                         type t)
  : spp_type(t),
    lsc_strategy(0), M(&m), velocity_block(), gradient_block(nullptr),
    divergence_block(nullptr), Schur_complement(nullptr),
    velocity_solver(nullptr), Schur_solver(nullptr), inverse_diagonal(),
    velocity_space(&m.get_row_space(0)), pressure_space(&m.get_row_space(2)),
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
#ifdef __2D__
  unsigned int dim = 2;
#else
  unsigned int dim = 3;
#endif
  bool pressure_correction_in_matrix = this->M->pressure_correction_enabled();
  this->M->disable_pressure_correction();
  
  // number of block columns and rows
  unsigned int n_rows = this->M->get_n_cell_rows();
  unsigned int n_cols = this->M->get_n_cell_columns();
  if(n_rows != (dim+1) || n_cols != (dim+1))
    ErrThrow("currently we only support a BlockMatrix for Stokes/Navier-Stokes "
             "problems.");
  
  //velocity block K
  // take all blocks except from last row and last column
  this->velocity_block = M->get_sub_blockfematrix(0, n_rows-2);
  this->fill_inverse_diagonal();
  
  if(lsc_strategy == 0)
  {
    this->velocity_solver.reset(
      new DirectSolver(this->velocity_block,
                       DirectSolver::DirectSolverTypes::umfpack));
  }
  
  //gradient block B^T
  // take last column without last block (which would be pressure-pressure)
  this->gradient_block = this->M->get_combined_submatrix({0,n_cols-1}, 
                                                         {n_rows-2,n_cols-1});
  //divergence block B
  this->divergence_block = this->M->get_combined_submatrix({n_rows-1,0}, 
                                                           {n_rows-1,n_cols-2});
  
  if(this->spp_type == Saddle_point_preconditioner::type::simple)
  { // SIMPLE
    //scale the gradient block B^T with the approximation of the mass-matrix
    this->gradient_block->scale(&inverse_diagonal[0], true);
    // construct an approximation to the Schur complement matrix
    this->Schur_complement = this->compute_Schur_complement_approximation();
    this->Schur_solver.reset(
      new DirectSolver(*this->Schur_complement,
                       DirectSolver::DirectSolverTypes::umfpack));
  }
  else
  { // original LSC (bdry corrected LSC treated the same here)
    //construct an approximation to the Poisson solve matrix
    this->Poisson_solver_matrix = this->compute_Poisson_solver_matrix();
    this->Poisson_solver.reset(
      new DirectSolver(*this->Poisson_solver_matrix,
                       DirectSolver::DirectSolverTypes::umfpack));
    u_star = BlockVector(velocity_block);
    p_star = BlockVector(*this->Poisson_solver_matrix);
    u_tmp = BlockVector(u_star);
    p_tmp = BlockVector(p_star);
  }
  
  if(this->spp_type == Saddle_point_preconditioner::type::bd_lsc)
  {
    //boundary corrected LSC
    //construct correctionMatrixBdry (anew...) and fill with ones
    bdryCorrectionMatrix_ = std::vector<double>(
        dim * velocity_space->GetN_DegreesOfFreedom(), 1.);
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
Saddle_point_preconditioner::Saddle_point_preconditioner(const BlockMatrix& m,
                                                         type t)
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
  
  if(lsc_strategy == 0)
  {
    try
    {
      this->velocity_solver.reset(
        new DirectSolver(this->velocity_block, 
                         DirectSolver::DirectSolverTypes::umfpack));
    }
    catch(...)
    {
      ErrThrow("something is wrong with the velocity matrix");
    }
  }
  
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
    
    for(unsigned int d = 0; d < n_diagonal_entries; ++d)
    {
      this->inverse_diagonal[d] = 1.0 / this->velocity_block.get(d, d);
    }
    
    // delete old objects which are no longer valid 
    this->Schur_complement.reset();
    this->Schur_solver.reset();
    
    //scale the gradient block B^T with the approximation of the mass-matrix
    this->gradient_block->scale(&inverse_diagonal[0], true);
    // construct an approximation to the Schur complement matrix
    this->Schur_complement = this->compute_Schur_complement_approximation();
    this->Schur_solver.reset(
      new DirectSolver(*this->Schur_complement,
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
  ErrThrow("Updateing a Saddle_point_preconditioner with a BlockMatrix is not "
           "possible, you need a BlockFEMatrix.");
  // otherwise the get_combined_matrix methods are possibly not giving the 
  // correct behavior.
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
      this->solve_velocity_block(z.get_entries(), up_star.block(0));
      //u_tmp = z.get_entries();
      //this->velocity_solver->solve(u_tmp, u_star);
      
      // ----------------------------------------------------------------------
      // solve ^S dp* = r_p - B du*
      // use BlockVector r to store the right hand side for the system 
      // involving the Schur approximation
      unsigned int n_blocks = z.n_blocks();
      r.copy(z.block(n_blocks - 1), n_blocks - 1); // copy last block, i.e. r_p
      //p_tmp = z.block(n_blocks - 1);
      // do r_p += -B * du*
      this->divergence_block->multiply(up_star.block(0), r.block(n_blocks - 1),
                                       -1.);
      this->Schur_solver->solve(r.block(n_blocks - 1),
                               up_star.block(n_blocks - 1));
      //this->divergence_block->multiply(u_star.get_entries(), 
      //                                 p_tmp.get_entries(), -1.);
      //this->Schur_solver->solve(p_tmp, p_star);
      
      
      // ----------------------------------------------------------------------
      // dp = dp*
      r.copy(up_star.block(n_blocks - 1), n_blocks - 1);
      //r.copy(p_star.get_entries(), n_blocks - 1);
      // du = du* - D^-1 B^T dp*
      // there are zeros in the u-block of r currently
      // store the result of -B^T dp* in u-block of r
      gradient_block->multiply(up_star.block(n_blocks - 1), r.block(0), -1.);
      //gradient_block->multiply(p_star.get_entries(), r.get_entries(), -1.);
      // scale with D^-1, loop over all velocity entries
      for(unsigned int i = 0, n_v = gradient_block->GetN_Rows(); i < n_v; ++i)
      {
        r[i] *= this->inverse_diagonal[i];
      }
      // add du*
      r.add(up_star.block(0), 0);
      r.add(up_star.block(1), 1);
      //r.add(u_star.block(0), 0);
      //r.add(u_star.block(1), 1);
#ifdef __3D__
      //r.add(u_star.block(2), 2);
#endif
      
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
      this->velocity_solver->solve(u_star, u_tmp);
      
      r.reset();
      // copy u_tmp and p_star back to r
      r.copy(u_tmp.block(0), 0);
      r.copy(u_tmp.block(1), 1);
#ifdef __3D__
      r.copy(u_tmp.block(2), 2);
#endif
      r.copy(p_star.block(0), n_blocks-1);
       
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
std::shared_ptr<TMatrix>
Saddle_point_preconditioner::compute_Schur_complement_approximation() const
{
  switch(this->spp_type)
  {
    case Saddle_point_preconditioner::type::simple:
      return std::shared_ptr<TMatrix>(
        this->divergence_block->multiply(this->gradient_block.get(), 1.));
    default:
      ErrThrow("unknown preconditioner for saddle point problems ");
  }
  return nullptr;
}

/* ************************************************************************** */
/* LSC === */
std::shared_ptr<BlockMatrix> 
Saddle_point_preconditioner::compute_Poisson_solver_matrix() const
{
  switch(this->spp_type)
  {
    case Saddle_point_preconditioner::type::lsc: // original LSC
    case Saddle_point_preconditioner::type::bd_lsc: // bdry corrected LSC
    {
      // construct matrix with pointer to structure of divergence_block but own
      // "entries" array. this is basically transposing, but with already known
      // structure
      TMatrix divBlockModified(*divergence_block);
      for(int gradRow = 0; gradRow < gradient_block->GetN_Rows(); ++gradRow)
      {
        int segmentStart = gradient_block->GetRowPtr()[gradRow];
        int segmentEnd = gradient_block->GetRowPtr()[gradRow + 1];
        
        for(int index = segmentStart; index < segmentEnd; ++index)
        {
          //set according entry in divBlockModified
          divBlockModified.set(gradient_block->GetKCol()[index], gradRow,
                               gradient_block->GetEntries()[index]);
        }
      }
      
      // compute ret as B*D*B^T for B: divergence_block and D: inverse_diagonal
      //TMatrix* ret = divBlockModified.multiply_with_transpose_from_right(
      //    inverse_diagonal);
      // get the matrix, put it into a shared_ptr, essentially this is the 
      // matrix which should be returned
      std::shared_ptr<TMatrix> 
        ret(divergence_block->multiply_with_transpose_from_right(
          inverse_diagonal));
      // put the shared_ptr into a vector
      std::vector<std::shared_ptr<TMatrix>> ret_as_vector{ret};
      // create a BlockMatrix and return it
      return std::make_shared<BlockMatrix>(1, 1, ret_as_vector);
      break;
    }
    default:
      ErrThrow("unknown preconditioner for saddle point problems ");
      break;
  }
  return nullptr;
}
/* === LSC */

/* ************************************************************************** */
// local assembling routine to assemble a mass velocity matrix
// This needs to go into some assembling class!
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
      MatrixM_row[j] += Mult*(ansatz00*test00);;
    }                            // endfor j
  }                              // endfor i
}

void Saddle_point_preconditioner::fill_inverse_diagonal()
{
  // number of entries on diagonal in velocity_block
  size_t n_diagonal_entries = this->velocity_block.get_n_total_rows();
  this->inverse_diagonal.reserve(n_diagonal_entries);
  
  //mass-matrix Q and its approximation
  if(this->spp_type == Saddle_point_preconditioner::type::simple)
  { // SIMPLE
    for(unsigned int d = 0; d < n_diagonal_entries; ++d)
    {
      this->inverse_diagonal.push_back(1.0 / this->velocity_block.get(d, d));
    }
  }
  else
  {
    Output::decreaseVerbosity(1);
    // either lsc or bd_lsc
    // first: assemble a velocity mass matrix
    
#ifdef __2D__
    // create an approriate LocalAssembling2D object:
    // We will assemble only one velocity block and then create a BlockMatrix
    // where this one block appears twice (three times in 3D)
    int n_terms = 1; // only 1 term (u,v)
    std::vector<MultiIndex2D> derivatives(1, D00); // that one term is the value
    std::vector<int> fe_space_numbers(1, 0); // only 1 space with index 0
    std::vector<int> row_space(1, 0); // 1 matrix 
    std::vector<int> column_space(1, 0); // 1 matrix
    std::vector<int> rhs_space(1, 0); // not needed
    CoeffFct2D * coeff = nullptr; // no coefficients
    AssembleFctParam2D* local_assembling_function = 
      local_assembling_velocity_mass;
    ManipulateFct2D * manipulate = nullptr; // not needed
    int n_matrices = 1; // assemble only one matrix
    int n_rhs = 0; // no right hand side
    int n_parameter = 0; // not needed
    std::vector<ParamFct*> parameter_functions; // no parameter functions
    std::vector<int> begin_parameter; // not needed
    int n_parameters = 0; // not needed
    TFEFunction2D ** fe_functions = nullptr;
    int n_fe_values = 0; // no other fe values are needed
    std::vector<int> fe_value_function_index; // not needed
    std::vector<MultiIndex2D> fe_value_multi_index; // not needed
    LocalAssembling2D la(n_terms, derivatives, fe_space_numbers, row_space, 
                         column_space, rhs_space, coeff, 
                         local_assembling_function, manipulate, n_matrices,
                         n_rhs, n_parameter, parameter_functions, 
                         begin_parameter, n_parameters, fe_functions, 
                         n_fe_values, fe_value_function_index, 
                         fe_value_multi_index);
    
    // prepare a call to Assemble2D
    int n_fe_spaces = 1;
    TSquareMatrix2D *sq_matrices[1];
    
    //TFESpace2D* neumann_velocity_space = 
    //  this->velocity_space->copy_with_all_Neumann_boundary();
    //const TFESpace2D* v_space = neumann_velocity_space;
    const TFESpace2D* v_space = this->velocity_space;
    
    TSquareMatrix2D one_block_of_mass_matrix(v_space);
    
    sq_matrices[0] =  &one_block_of_mass_matrix;
    int n_rect_mat = 0; // only square matrices
    TMatrix2D *rect_matrices[0];
    double **rhs = nullptr; // right hand sides
    const TFESpace2D ** fe_spaces_rhs = nullptr;
    // which boundary conditions are correct here? We could as well use the ones
    // from the velocity_space
    BoundCondFunct2D * boundary_conditions[2] = {
      BoundConditionNoBoundCondition, BoundConditionNoBoundCondition};
    BoundValueFunct2D* non_const_bound_values[2] = {
      BoundaryValueHomogenous, BoundaryValueHomogenous };
    
    Assemble2D(n_fe_spaces, &v_space, n_matrices, sq_matrices,
               n_rect_mat, rect_matrices, n_rhs, rhs, fe_spaces_rhs,
               boundary_conditions, non_const_bound_values, la);
    
    BlockFEMatrix mass_matrix({v_space, v_space});
    mass_matrix.replace_blocks(one_block_of_mass_matrix,
                               {{0,0}, {1,1}}, {false, false});
    
    //number of entries on diagonal in mass matrix
    for(unsigned int d = 0; d < n_diagonal_entries; ++d)
    {
      this->inverse_diagonal.push_back(1.0 / mass_matrix.get(d, d));
      //Output::print("inverse_diagonal[", d, "] = ", inverse_diagonal[d]);
    }
    //delete neumann_velocity_space;
#else
    ErrThrow("not yet implemented in 3D");
#endif
    
    Output::increaseVerbosity(1);
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
          double horizontalDofsCorrectionValue = std::max(std::abs(normalX),
                                                          epsilon);
          double verticalDofsCorrectionValue = std::max(std::abs(normalY),
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
            bdryCorrectionMatrix_[velocity_space->GetN_DegreesOfFreedom()
                + iVeloDof] = verticalDofsCorrectionValue;
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
          double x_DofsCorrectionValue = std::max(std::abs(normalX), epsilon);
          double y_DofsCorrectionValue = std::max(std::abs(normalY), epsilon);
          double z_DofsCorrectionValue = std::max(std::abs(normalZ), epsilon);
          
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
    // construct matrix with pointer to structure of divergence_block but own
    // "entries" array
    // this is basically transposing, but with already known structure
    TMatrix divBlockModified(*divergence_block);
    for(int gradRow = 0; gradRow < gradient_block->GetN_Rows(); ++gradRow)
    {
      
      int segmentStart = gradient_block->GetRowPtr()[gradRow];
      int segmentEnd = gradient_block->GetRowPtr()[gradRow + 1];
      
      for(int index = segmentStart; index < segmentEnd; ++index)
      {
        //set according entry in divBlockModified
        divBlockModified.set(gradient_block->GetKCol()[index], gradRow,
                             gradient_block->GetEntries()[index]);
      }
    }
    
    //a pointer to a new, heap constructed matrix B*H^{-1}*B^T
    //TMatrix* ret = divBlockModified.multiply_with_transpose_from_right(
    //    bdryCorrectionMatrix_, *Poisson_solver_matrix->GetStructure());
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

/* ************************************************************************** */
void Saddle_point_preconditioner::solve_velocity_block(const double* rhs,
                                                       double * solution) const
{
  switch(this->lsc_strategy)
  {
    case 0: // umfpack
      this->velocity_solver->solve(rhs, solution);
      break;
    case 1:
    case 2:
//     {
//       TSquareMatrix * mat = new TSquareMatrix(
//           new TSquareStructure(this->velocity_block->GetN_Rows(),
//                                this->velocity_block->GetN_Entries(),
//                                this->velocity_block->GetKCol(),
//                                this->velocity_block->GetRowPtr()),
//           this->velocity_block->GetEntries());
//       
//       if(0)
//       {
//         // call AMG library to be able to use different solvers and 
//         // preconditioners
//         // this produces a lot of output
//         TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 3;//AMG_SSOR, reset later
//         TDatabase::ParamDB->SC_SMOOTHER_SCALAR = 3; // AMG_SSOR, not reset
//         Solver(mat, rhs, solution);
//         TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 20;
//       }
//       else // a possible way to do the same without the amg library
//       {
//         bool suppress_output = true;//TDatabase::ParamDB->SC_VERBOSE < 3;
//         std::ofstream * lStream;
//         std::streambuf * lBufferOld;
//         if(suppress_output)
//         {
//           lStream = new std::ofstream("/dev/null");
//           lBufferOld = std::cout.rdbuf();
//           std::cout.rdbuf(lStream->rdbuf());
//         }
//         
//         OutPut("velocity solver\n");
//         TItMethod * prec;
//         // solve velocity equation using FGMRES preconditioned with SSOR or BCGS
//         if(this->lsc_strategy == 1)
//           prec = new TSSORIte(MatVect_Scalar, Defect_Scalar, NULL, 0,
//                               this->velocity_block->GetN_Rows(), 1);
//         else
//           prec = new TBcgs(MatVect_Scalar, Defect_Scalar, NULL, 0,
//                            this->velocity_block->GetN_Rows(), 1);
//         
//         TFgmresIte vel_solver(MatVect_Scalar, Defect_Scalar, prec, 0,
//                               this->velocity_block->GetN_Rows(), 1);
//         //vel_solver.SetRedFactor(1e-2);
//         vel_solver.SetMaxit(100);
//         int n_it = vel_solver.Iterate(&mat, NULL, solution, rhs);
//         delete prec;
//         if(suppress_output)
//         {
//           std::cout.rdbuf( lBufferOld );
//           delete lStream;
//         }
//         OutPut("  number of iterations for velocity solve " << n_it << endl);
//       }
//       break;
//     }
    default:
      ErrThrow("unknown value for lsc_strategy ", this->lsc_strategy);
      break;
  }
}
