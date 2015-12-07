#include <algorithm>
#include <StokesDarcy2D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <FEDatabase2D.h>
#include <preconditioner.h>
#include <Output2D.h>
#include <DirectSolver.h>
#include <solvers/solver.h>


/** ************************************************************************ */
StokesDarcy2D::StokesDarcy2D(
    std::map< InterfaceCondition, StokesProblem* > ns_problems,
    std::map< InterfaceCondition, DarcyPrimal* > d_problems)
 : big_matrix(nullptr), big_rhs(nullptr), big_solution(nullptr)
{
  //Output::print<1>("StokesDarcy2D constructor");
  int solution_strategy = TDatabase::ParamDB->StoDa_solutionStrategy;
  error = 1e10; // set to large value at the beginning
  initialError = 0; // needs to be set correctly
  bigResidual = 1e10; // set to large value at the beginning
  // if no big system is used, set bigResidual to 0 (i.e., no restriction)
  if(solution_strategy <= 0)
    bigResidual = 0.0;
  initialBigResidual = 0; // needs to be set correctly
  //s = NULL; // solver object
  f_prec = NULL; // preconditioner object
  p_prec = NULL; // preconditioner object
  StokesFirst = (TDatabase::ParamDB->StoDa_StokesFirst == 1) ? true : false;
  eta_hom = NULL;
  C = NULL;
  CT = NULL;
  
  switch(TDatabase::ParamDB->StoDa_problemType)
  {
    case 0: // Neumann-Neumann
      coupled_stokes = ns_problems.at(Neumann);
      coupled_darcy = d_problems.at(Neumann);
      if(solution_strategy == 1 || solution_strategy == -1
         || solution_strategy == 2 || solution_strategy == -2)
      {
        iteration_stokes = coupled_stokes;
        iteration_darcy = coupled_darcy;
        if(TDatabase::ParamDB->StoDa_updatingStrategy == 3) // C-RR
        {
          iteration_stokes = ns_problems.at(Robin);
          iteration_darcy = d_problems.at(Robin);
        }
      }
      else if(solution_strategy == 3)
      { // Stecklov-Poincare
        if(StokesFirst)
        {
          iteration_stokes = ns_problems.at(Neumann);
          iteration_darcy = d_problems.at(Dirichlet);
          p_prec = d_problems.at(Neumann);
        }
        else
        {
          iteration_stokes = ns_problems.at(Dirichlet);
          iteration_darcy = d_problems.at(Neumann);
          f_prec = ns_problems.at(Neumann);
        }
      }
      else if (solution_strategy == 0)
      { // only direct solution
        iteration_stokes = NULL;
        iteration_darcy = NULL;
      }
      else
        ErrThrow("unkonwn solution strategy, set 'StoDa_solutionStrategy' to ",
                 "either 0, 1, -1, 2  or -2.");
      break;
    case 1: // Robin-Robin
      coupled_stokes = ns_problems.at(Robin);
      coupled_darcy = d_problems.at(Robin);
      iteration_stokes = coupled_stokes;
      iteration_darcy = coupled_darcy;
      break;
    case 2: // weak Robin-Robin
      coupled_stokes = ns_problems.at(weakRobin);
      coupled_darcy = d_problems.at(weakRobin);
      iteration_stokes = coupled_stokes;
      iteration_darcy = coupled_darcy;
      break;
    case 4: // Dirichlet-Dirichlet
      if(solution_strategy != 0)
        ErrThrow("Dirichlet-Dirichlet does not work iteratively");
      coupled_stokes = ns_problems.at(Dirichlet);
      coupled_darcy = d_problems.at(Dirichlet);
      break;
    case 3:
    case 5:
    case 6:
      ErrThrow("these problem Types are not yet implemented");
      break;
    default:
      ErrThrow("unknown problem type");
      break;
  }
}

/** ************************************************************************ */
StokesDarcy2D::~StokesDarcy2D()
{
  Output::print<3>("StokesDarcy2D destructor");
  //delete s;
  delete eta_hom;
  // the StokesProblem and DarcyPrimal are destroyed in main, which is where
  // they were created
}

/** ************************************************************************ */
template<class problem>
void check_linearity(problem& p, InterfaceFunction& eta)
{
  // check p(-eta) = -p(eta) + 2 p(0)
  // set to something other than zero
  eta = 5.0 / ((eta.length() - 1.0) / 2.0);
  InterfaceFunction eta2(eta);
  eta2.reset();
  
  eta *= -1.0;
  p->solve(eta);
  // eta2 = p(-eta)
  p->map_solution_to_interface(eta2, 1.0);
  
  eta *= -1.0;
  p->solve(eta);
  // eta2 = p(-eta) + p(eta)
  p->map_solution_to_interface(eta2, 1.0);
  
  eta.reset();
  p->solve(eta);
  // eta2 = p(-eta) + p(eta) - 2 p(0)
  p->map_solution_to_interface(eta2, -2.0);
  
  if(eta2.norm() > 1e-12 * sqrt(eta2.length()))
    Error("norm which should be zero " << eta2.norm() << endl);
  
  // check p(eta + eta2) = p(eta) + p(eta2) - p(0)
  // set to something other than zero
  eta = 5.0 / ((eta.length() - 1.0) / 2.0);
  eta2 = 1.0;
  InterfaceFunction eta3(eta);
  eta3.reset();
  
  p->solve(eta);
  // eta3 = p(eta)
  p->map_solution_to_interface(eta3, 1.0);
  
  p->solve(eta2);
  // eta3 = p(eta) + p(eta2)
  p->map_solution_to_interface(eta3, 1.0);
  
  eta += eta2;
  p->solve(eta);
  // eta3 = p(eta) + p(eta2) - p(eta + eta2)
  p->map_solution_to_interface(eta3, -1.0);
  
  eta.reset();
  p->solve(eta);
  // eta3 = p(eta) + p(eta2) - p(eta + eta2) - p(0)
  p->map_solution_to_interface(eta3, -1.0);
  
  if(eta3.norm() > 1e-12 * sqrt(eta3.length()))
    Error("norm which should be zero " << eta3.norm() << endl);
}

/** ************************************************************************ */
void StokesDarcy2D::solveDirect()
{
  Output::print<1>("\nSolving the coupled system directly");
  // full system matrix is 
  // bigMatrix = ( S  CT )
  //             ( C  D  )
  GetCouplingMatrix(); // assemble C and CT
  c_stokes()->setCouplingMatrix(CT);
  c_darcy()->setCouplingMatrix(C);
  SolveOneSystem();
}

/** ************************************************************************ */
void StokesDarcy2D::solve_fixed_point(InterfaceFunction& eta) const
{
  const double theta =
      StokesFirst ?
          TDatabase::ParamDB->StoDa_theta_f : TDatabase::ParamDB->StoDa_theta_p;
  InterfaceFunction eta_damp; // store old eta for damping
  if(theta != 1.0) // damping
  {
    eta_damp = eta;
    eta_damp *= 1-theta;
  }
  
  switch(TDatabase::ParamDB->StoDa_updatingStrategy)
  {
    case 1: // fixed point iteration (Neumann-Neumann or Robin-Robin)
    case 4: // fixed point iteration D-RR
      if(StokesFirst)
      {
        // Interface integrals & solving for the Stokes part
        stokes()->solve(eta); 
        stokes()->update(eta);
        // Interface integrals & solving for the Darcy part
        darcy()->solve(eta);
        darcy()->update(eta); // includes damping
      }
      else
      {
        // Interface integrals & solving for the Darcy part
        darcy()->solve(eta);
        darcy()->update(eta);
        // Interface integrals & solving for the Stokes part
        stokes()->solve(eta); 
        stokes()->update(eta); // includes damping
      }
      break;
    case 2: // Newton
      ErrThrow("Newton method for fixed point equation not yet implemented");
      break;
    default:
      ErrThrow("solving fixed point equation. Set the value of ", 
               "'StoDa_updatingStrategy' to either 1, 2, or 4");
      break;
  }
  
  if(theta != 1.0) // damping
  {
    eta *= theta;
    eta += eta_damp;
  }
}

/** ************************************************************************ */
void StokesDarcy2D::solve(const InterfaceFunction& eta_r,
                          InterfaceFunction& eta_z) const
{
  // solve "M eta_z = eta_r" if "M" is this preconditioner 
  Output::print<1>("apply preconditioner");
  // no preconditioner would mean eta_z = eta_r
  if(eta_z.length() == 0)
    eta_z = eta_r; // this done here just in case eta_z is not yet initialized
  eta_z.reset();
  
  InterfaceFunction * eta_hom_prec = NULL;
  
  if(StokesFirst)
  {
    eta_hom_prec = p_prec->get_homogeneous_solution();
    p_prec->solve(eta_r);
    p_prec->map_solution_to_interface(eta_z, -1.0);
  }
  else
  {
    // eta_hom is solution to  H_f,N( 0 )
    eta_hom_prec = f_prec->get_homogeneous_solution();
    f_prec->solve(eta_r);
    f_prec->map_solution_to_interface(eta_z, -1.0);
  }
  
  // homogeneous solution should have been initialized in 
  // 'StokesDarcy2D::solve_Stecklov_Poincare(eta)'
  if(eta_hom_prec == NULL)
  {
    ErrThrow("homogeneous solution of preconditioner not available");
  }
  // substract homogeneous solution to make this preconditioner linear (not 
  // only affine)
  eta_z -= *eta_hom_prec;
  
  if(StokesFirst)
    eta_z *= TDatabase::ParamDB->StoDa_theta_f;
  else
    eta_z *= TDatabase::ParamDB->StoDa_theta_p;
}

/** ************************************************************************ */
void StokesDarcy2D::solve_coupled_system(InterfaceFunction& eta_f,
                                         InterfaceFunction& eta_p) const 
{
  bool GaussSeidel = (TDatabase::ParamDB->StoDa_algorithm == 1) ? true : false;
  if(GaussSeidel)
  {
    if(StokesFirst)
    {
      // Interface integrals & solving for the Stokes part
      stokes()->solve(eta_f); 
      stokes()->update(eta_p, &eta_f); // includes damping
    }
    {
      // Interface integrals & solving for the Darcy part
      darcy()->solve(eta_p);
      darcy()->update(eta_f, &eta_p); // includes damping
    }
    if(!StokesFirst)
    {
      // Interface integrals & solving for the Stokes part
      stokes()->solve(eta_f); 
      stokes()->update(eta_p, &eta_f); // includes damping
    }
  }
  else
  { // Jacobi
    stokes()->solve(eta_f);
    darcy()->solve(eta_p);
    
    darcy()->update(eta_f, &eta_p); // includes damping
    stokes()->update(eta_p, &eta_f); // includes damping
  }
}

/** ************************************************************************ */
bool StokesDarcy2D::stopIteration(int it)
{
  // diff will be the (weighted) sum of the relative differences of the 
  // solution vectors. Initially it is set to a value such that a convergence 
  // criterion is not fulfilled.
  double diff = TDatabase::ParamDB->StoDa_relDiff_solution;
  // true if convergence is slow, i.e. if succesive iterates are identical up  
  // to machine precision
  bool slow = false;
  //if(it) // compute reldiff
  {
    int n_DOF;
    double diffDarcyP, diffStokesU, diffStokesP;
    const double a1 = TDatabase::ParamDB->StoDa_relDiff_factor1; // StokesU
    const double a2 = TDatabase::ParamDB->StoDa_relDiff_factor2; // StokesP
    const double a3 = TDatabase::ParamDB->StoDa_relDiff_factor3; // DarcyP
    n_DOF = darcy()->getP().GetFESpace2D()->GetN_DegreesOfFreedom();
    diffDarcyP = relabsDiff(n_DOF, darcy()->getSolOld(), darcy()->getSol());
    
    n_DOF = 2 * stokes()->get_velocity_space().GetN_DegreesOfFreedom();
    diffStokesU = relabsDiff(n_DOF, stokes()->get_u_old().GetValues(),
                             stokes()->get_velocity().GetValues());
    
    n_DOF = stokes()->get_pressure_space().GetN_DegreesOfFreedom();
    diffStokesP = relabsDiff(n_DOF, stokes()->get_p_old().GetValues(),
                             stokes()->get_pressure().GetValues());
    
    if(diffDarcyP + diffStokesU + diffStokesP < 1e-14)
      slow = true;
    //diff = sqrt(a1*POW(diffStokesU,2)
    //          + a2*POW(diffStokesP,2)
    //          + a3*POW(diffDarcyP,2));
    diff = a1 * diffStokesU + a2 * diffStokesP + a3 * diffDarcyP;
    Output::print<1>(" DIFFERENCE in Stokes velocity vector ", diffStokesU);
    Output::print<1>(" DIFFERENCE in Stokes pressure vector ", diffStokesP);
    Output::print<1>(" DIFFERENCE in Darcy pressure vector ", diffDarcyP);
    Output::print<1>(" DIFFERENCE in overall solution ", diff);
  }
  
  //===========================================================================
  // error on interface ( e^k )
  const double error_old = error;
  error = ErrorOnInterface(*stokes(), *darcy(), it);
  const double relError = 1 - error / error_old;
  Output::print<1>(" Error on interface ", error);
  if(it)
    Output::print<1>(" relative difference of error on interface ", relError);
  if(it == 0)
    initialError = error;
  
  //if(relError < 1e-14)
  //  slow = true;
  
  //===========================================================================
  // compute defect of large system
  if(TDatabase::ParamDB->StoDa_solutionStrategy >= 0)
  { // compute defect of big system:
    // r = ( S  CT ) ( u,p ) _ ( f_f )
    //     ( C  D  ) ( phi )   ( f_p )
    BlockVector r(*big_matrix);
    r.copy(stokes()->get_solution().get_entries(), 0);
    r.copy(darcy()->getSol(), 1);
    r *= *big_matrix;
    r -= *big_rhs;
    bigResidual = r.norm();
    Output::print<1>(" Norm of big system defect ", bigResidual);
  }
  
  //===========================================================================
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
    //compute residuals and print results to console
    ComputeResiduals(*stokes(), *darcy());
  //===========================================================================
  // write output and measure errors
  WriteVtk_and_measureErrors(*stokes(), *darcy(), NULL, it);
  
  //===========================================================================
  // check convergence criteria
  bool stopIteration = false;
  if( //abs(relError) < TDatabase::ParamDB->StoDa_relDiff_interfaceError && 
  diff < TDatabase::ParamDB->StoDa_relDiff_solution
  && bigResidual < TDatabase::ParamDB->StoDa_bigResidual)
  {
    Output::print<1>("Converged after ", it+1, " iterations");
    stopIteration = true;
  }
  // check divergence criterion
  if(error > 1e10)
  {
    Output::print<1>("DIVERGED after ", it+1, " iterations !!!!!");
    stopIteration = true;
  }
  // check if error on interface is zero up to machine precision or so
  if(0)
    if(error < 1e-14) // (almost) never going to happen
    {
      Output::print<1>("Converged after ", it+1,
                       " iterations with interface error ", error);
      stopIteration = true;
    }
  // check if error on interface increased
  if(0)if(relError < -1e-1)
  {
    Output::print<1>("Error increased after ", it+1, " iterations.");
    stopIteration = true;
  }
  // check if maximum number of iterations is reached
  if(it == TDatabase::ParamDB->StoDa_nIterations - 1)
  {
    Output::print<1>("Reached maximum number of ",
                     TDatabase::ParamDB->StoDa_nIterations, " iterations!!");
    stopIteration = true;
  }
  // check if convergence is slow
  if(slow && !stopIteration)
  {
    Output::print<1>("Iteration is stuck after ", it+1, " iterations!!");
    stopIteration = true;
  }
  
  if(stopIteration)
  {
    Output::print<1>("\n");
    if(TDatabase::ParamDB->StoDa_solutionStrategy >= 1)
    {
      // big system has been solved
      Output::print<1>("big residual changed from ", initialBigResidual, " to ",
                       bigResidual, " quotient ",
                       bigResidual/initialBigResidual, 
                       "\n  average reduction factor per iteration ",
                       exp(log(bigResidual/initialBigResidual)/it));
    }
    Output::print<1>("interface error changed from ", initialError, " to ",
                     error, " quotient ", error/initialError,
                     "\n  average reduction factor per iteration ",
                     exp(log(error/initialError)/it));
  }
  return stopIteration;
}

/** ************************************************************************ */
void StokesDarcy2D::SolveOneSystem()
{
  // full system matrix is 
  // bigMatrix = ( S  CT )
  //             ( C  D  )
  // combine all matrices into one big matrix
  CombineBigMatrix(); // create 'big_matrix'
  big_solution.reset(new BlockVector(*big_matrix));
  big_rhs.reset(new BlockVector(*big_matrix));
  big_rhs->copy(c_stokes()->get_rhs().get_entries(), 0);
  big_rhs->copy(c_darcy()->get_rhs().get_entries(), 1);
  
  // write out some information on the big matrix system
  Output::print<1>(" Big system: Size: ", big_matrix->n_total_rows(), " x ",
                   big_matrix->n_total_cols(), "\tnumber of entries: ",
                   big_matrix->n_total_entries());
  
  
  //big_matrix->get_combined_matrix()->PrintFull("M");
  //big_matrix->get_combined_matrix()->info(3);
  //c_stokes()->getComposedMatForBigSystem()->PrintFull("S");
  //c_stokes()->get_matrix().get_A_block(0)->PrintFull("A11");
  //this->CT->PrintFull("CT");
  //c_darcy()->getComposedMatForBigSystem()->PrintFull("D");
  //this->C->PrintFull("C");
  
  typedef preconditioner <BlockMatrix, BlockVector> block_prec;
  typedef solver <BlockMatrix, BlockVector, block_prec> block_solver;
  
  /// create preconditioner and solver object for direct solution of big system.
  block_prec big_prec(big_matrix.get());
  block_solver big_s(big_matrix.get(), big_solution.get(), big_rhs.get(),
                     &big_prec);
  big_s.solve();
  
  c_stokes()->setDirectSol(big_solution->block(0));
  c_darcy()->setDirectSol(big_solution->block(1));
  if((TDatabase::ParamDB->StoDa_updatingStrategy == 3
      && TDatabase::ParamDB->StoDa_solutionStrategy == 1) // C-RR
     || TDatabase::ParamDB->StoDa_solutionStrategy == 3) // Stecklov Poincare 
  {
    stokes()->setDirectSol(big_solution->block(0));
    darcy()->setDirectSol(big_solution->block(1));
  }
  // else c_stokes() == stokes()  (same with darcy) or only direct solution
  
  if(TDatabase::ParamDB->SC_VERBOSE > 0)
  {
    BlockVector r(*big_rhs); // copy structure
    r = *big_solution; // copy values
    r *= *big_matrix;
    r -= *big_rhs;
    Output::print<1>("  residual of big system: ", norm(r));
  }
  // compute initial big residual
  {
    BlockVector r(*big_matrix);
    r.copy(c_stokes()->get_solution().get_entries(), 0);
    r.copy(c_darcy()->get_solution().get_entries(), 1);
    r *= *big_matrix;
    r -= *big_rhs;
    initialBigResidual = norm(r);
  }
}

/** ************************************************************************ */
void StokesDarcy2D::GetCouplingMatrix()
{
  // only the structure of C is build here, at the end we transpose this 
  // structure to get the structure of CT
  // Note: this ignores Dirichlet rows. i.e. there are entries in these rows.
  //       However during assembling nothing is written in these rows. The 
  //       advantage of doing it this way is: we can simply call transpose.
  
  // velocity and pressure spaces for Stokes and Darcy
  const TFESpace2D* vSpaceStokes = &c_stokes()->get_velocity_space();
  const TFESpace2D* pSpaceStokes = &c_stokes()->get_pressure_space();
  const TFESpace2D* pSpaceDarcy = c_darcy()->getP().GetFESpace2D();
  // number of Stokes velocity degrees of freedom for each component
  const int N_SveloDOF = vSpaceStokes->GetN_DegreesOfFreedom();
  
  // everything needed to call the constructor of TStructure
  const int n_rows = c_darcy()->get_rhs().length();
  const int n_cols = c_stokes()->get_size();
  int N_entries = 0;
  int *cols, *rows = new int[n_rows + 1];
  
  // for each row we have one vector containing all (Stokes) dofs which couple 
  // with this dof. For all dofs which do not belong to a cell adjacent to the
  // interface, the corresponding vector will be of length 0
  std::vector<std::vector<int> > coupledDOF(n_rows);
  
  // array containing the number of basis functions for all finite elements
  int const * const N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  
  // Id of the Stokes space
  const int S_id = vSpaceStokes->GetCollection()->GetCell(0)->GetReference_ID();
  
  // loop over the interface edges
  for(unsigned int iEdge = 0; iEdge < c_stokes()->getInterface().size(); 
      ++iEdge)
  {
    const TInnerInterfaceJoint * thisEdge = c_stokes()->getInterface()[iEdge];
    TBaseCell *d_cell, *s_cell = thisEdge->GetNeighbour(0);
    if(s_cell->GetReference_ID() != S_id)
    {
      d_cell = s_cell;
      s_cell = thisEdge->GetNeighbour(1);
    }
    else
      d_cell = thisEdge->GetNeighbour(1);
    
    const FE2D pDarcyFEId = pSpaceDarcy->GetFE2D(0, d_cell);
    const int N_pDarcyBaseFunct = N_BaseFunct[pDarcyFEId];
    int const* const pDarcyDOF = pSpaceDarcy->GetGlobalDOF(
        d_cell->GetCellIndex());
    
    const FE2D uStokesFEId = vSpaceStokes->GetFE2D(0, s_cell);
    const int N_uStokesBaseFunct = N_BaseFunct[uStokesFEId];
    int const* const uStokesDOF = vSpaceStokes->GetGlobalDOF(
        s_cell->GetCellIndex());
    
    const FE2D pStokesFEId = pSpaceStokes->GetFE2D(0, s_cell);
    const int N_pStokesBaseFunct = N_BaseFunct[pStokesFEId];
    int const* const pStokesDOF = pSpaceStokes->GetGlobalDOF(
        s_cell->GetCellIndex());
    
    for(int jD = 0; jD < N_pDarcyBaseFunct; jD++)
    {
      const int testDOF = pDarcyDOF[jD];
      // loop over all degrees of freedom of this Sokes cell
      for(int jS = 0; jS < N_uStokesBaseFunct; jS++)
      {
        const int ansatz_vDOF = uStokesDOF[jS];
        // first Stokes velocity component
        coupledDOF[testDOF].push_back(ansatz_vDOF);
        // second Stokes velocity component
        coupledDOF[testDOF].push_back(ansatz_vDOF + N_SveloDOF);
      }
      for(int jS = 0; jS < N_pStokesBaseFunct; jS++)
      {
        const int ansatz_pDOF = pStokesDOF[jS];
        // Stokes pressure component
        coupledDOF[testDOF].push_back(ansatz_pDOF + 2 * N_SveloDOF);
      }
    }
  }
  
  // fill the array 'rows' and compute 'N_entries'
  rows[0] = 0;
  for(int i = 0; i < n_rows; i++)
  {
    // sort the Stokes dofs for the i-th Darcy dof
    std::sort(coupledDOF[i].begin(), coupledDOF[i].end());
    // remove all duplicates 
    std::vector<int>::iterator it = std::unique(coupledDOF[i].begin(),
                                                coupledDOF[i].end());
    // resize to the size without duplicates
    coupledDOF[i].resize(std::distance(coupledDOF[i].begin(), it));
    N_entries += coupledDOF[i].size();
    rows[i + 1] = N_entries;
  }
  
  cols = new int[N_entries];
  // fill the array 'cols'
  for(int i = 0; i < n_rows; i++)
  {
    for(int j = 0; j < rows[i + 1] - rows[i]; j++)
    {
      cols[rows[i] + j] = coupledDOF[i].at(j);
    }
  }
  // generate sparse matrix
  std::shared_ptr<TStructure> Cstructure(new TStructure(n_rows, n_cols, 
                                                        N_entries, cols, rows));
  C.reset(new TMatrix(Cstructure)); // empty matrix
  CT.reset(new TMatrix(Cstructure->GetTransposed())); // empty matrix
  // Assemble the coupling matrix
  AssembleCouplingMatrices();
  
  reorder_rows(CT);
  if(TDatabase::ParamDB->StoDa_periodicBoundary)
  {
    c_darcy()->makePeriodicBoundary(C);
    c_stokes()->makePeriodicBoundary(CT);
  }  
  // remove all zeros, which are explicitly saved in the sparsity pattern
  CT->remove_zeros(-1);//1e-14);
  C->remove_zeros(-1);//1e-14);
}

/** ************************************************************************ */
void StokesDarcy2D::AssembleCouplingMatrices() const
{
  Output::print<1>("Assemble interface integrals into coupling matrices");
  // velocity and pressure spaces for Stokes and Darcy
  const TFESpace2D & vSpaceStokes = c_stokes()->get_velocity_space();
  const TFESpace2D & pSpaceStokes = c_stokes()->get_pressure_space();
  const TFESpace2D & pSpaceDarcy = c_darcy()->get_space();
   
  
  // number of Stokes velocity degrees of freedom for each component
  const int N_SveloDOF = vSpaceStokes.GetN_DegreesOfFreedom();
  const int N_Sactive = vSpaceStokes.GetActiveBound();
  
  // number of Darcy velocity degrees of freedom for each component
  const int N_Dactive =  pSpaceDarcy.GetActiveBound();
  
  // Id of the Stokes space
  const int S_id = vSpaceStokes.GetCollection()->GetCell(0)->GetReference_ID();
  
  // everything needed for the local assembling
  local_edge_assembling l;
  
  // loop over the interface edges
  for(unsigned int iEdge = 0; iEdge < c_stokes()->getInterface().size(); 
      ++iEdge)
  {
    const TInnerInterfaceJoint * thisEdge = c_stokes()->getInterface()[iEdge];
    TBaseCell *d_cell = thisEdge->GetNeighbour(0);
    TBaseCell *s_cell = thisEdge->GetNeighbour(1);
    if(s_cell->GetReference_ID() != S_id)
    { // swap
      TBaseCell *temp_cell = d_cell;
      d_cell = s_cell;
      s_cell = temp_cell;
    }
    
    const FE2D pDarcyFEId = pSpaceDarcy.GetFE2D(0, d_cell);
    TFE2D *pD = TFEDatabase2D::GetFE2D(pDarcyFEId);
    const int N_pDBf = pD->GetN_DOF(); // number of Darcy pressure basis functs
    int const* const pDarcyDOF = pSpaceDarcy.GetGlobalDOF(
        d_cell->GetCellIndex());
    // Basis functions for Darcy pressure
    TBaseFunct2D * pBfD = pD->GetBaseFunct2D();
    RefTrans2D refTransD = pD->GetRefTransID();
    
    const FE2D uStokesFEId = vSpaceStokes.GetFE2D(0, s_cell);
    TFE2D *uS = TFEDatabase2D::GetFE2D(uStokesFEId);
    // number of Stokes velocity basis functions
    const int N_uSBf = uS->GetN_DOF(); 
    int *uStokesDOF = vSpaceStokes.GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * uBfS = uS->GetBaseFunct2D();
    RefTrans2D refTransS = uS->GetRefTransID();
    
    const FE2D pStokesFEId = pSpaceStokes.GetFE2D(0, s_cell);
    TFE2D *pS = TFEDatabase2D::GetFE2D(pStokesFEId);
    // number of Stokes pressure basis functions
    const int N_pSBf = pS->GetN_DOF(); 
    int *pStokesDOF = pSpaceStokes.GetGlobalDOF(s_cell->GetCellIndex());
    // Basis functions for Stokes velocity
    TBaseFunct2D * pBfS = pS->GetBaseFunct2D();
    
    // polynomial degree of finite element, needed for choosing an 
    // appropriate quadrature formula
    int fe_degree = pBfD->GetPolynomialDegree();
    fe_degree = MAX(fe_degree, uBfS->GetPolynomialDegree());
    // get the quadrature formula
    QuadFormula1D QFId; // quadrature formula id
    int N_LinePoints; // number of qudrature points
    double *LineWeights, *zeta; // quadrature weights and points (on [-1,1])
    {
      QFId = TFEDatabase2D::GetQFLineFromDegree(2 * fe_degree);
      TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(QFId);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      // qf1 no longer needed, only local scope
    }
    
    // make sure all functions & derivatives are available for this quadrature
    pBfD->MakeRefElementData(QFId);
    uBfS->MakeRefElementData(QFId);
    pBfS->MakeRefElementData(QFId);
    
    // compute length of the edge
    const double hE = thisEdge->GetLength();
    
    // index of this edge in the adjacent cells in Stokes and Darcy subdomain
    const int eID = thisEdge->GetIndexInNeighbor(d_cell); // edge Index Darcy
    const int eIS = thisEdge->GetIndexInNeighbor(s_cell); // edge Index Stokes
    // normal (pointing out of the Stokes subdomain) 
    double nx, ny;
    getNormal(s_cell, thisEdge, nx, ny);
    // tangent: we take the normal and rotate it counterclockwise
    double tx = -ny;
    double ty = nx;
    
    l.nt.set(nx, ny, tx, ty);
    l.hE = hE;
    
    
    // values of all functions and derivatives on reference element 
    double **pDref = TFEDatabase2D::GetJointValues2D(pBfD->GetID(), QFId, eID);
    double **pDxiref = TFEDatabase2D::GetJointDerivatives2D(pBfD->GetID(), QFId,
                                                            eID, D10);
    double **pDetaref = TFEDatabase2D::GetJointDerivatives2D(pBfD->GetID(),
                                                             QFId, eID, D01);
    
    double **uSref = TFEDatabase2D::GetJointValues2D(uBfS->GetID(), QFId, eIS);
    double **uSxiref = TFEDatabase2D::GetJointDerivatives2D(uBfS->GetID(), QFId,
                                                            eIS, D10);
    double **uSetaref = TFEDatabase2D::GetJointDerivatives2D(uBfS->GetID(),
                                                             QFId, eIS, D01);
    double **pSref = TFEDatabase2D::GetJointValues2D(pBfS->GetID(), QFId, eIS);
    
    // values and derivatives of all functions at one quadrature point
    double *pDorig = new double[N_pDBf];
    double *pDxorig = new double[N_pDBf];
    double *pDyorig = new double[N_pDBf];
    double *uSorig = new double[N_uSBf];
    double *uSxorig = new double[N_uSBf];
    double *uSyorig = new double[N_uSBf];
    double *pSorig = new double[N_pSBf];
    
    // loop over all quadrature points
    for(int k = 0; k < N_LinePoints; k++)
    {
      TFEDatabase2D::SetCellForRefTrans(s_cell, refTransS);
      TFEDatabase2D::GetOrigValues(refTransS, zeta[k], uBfS, eIS, uSref[k],
                                   uSxiref[k], uSetaref[k], uSorig, uSxorig,
                                   uSyorig);
      TFEDatabase2D::GetOrigValues(refTransS, zeta[k], pBfS, eIS, pSref[k],
                                   NULL, NULL, pSorig, NULL, NULL);
      TFEDatabase2D::SetCellForRefTrans(d_cell, refTransD);
      // DIRTY HACK! 
      // usually the quadrature points go through an edge counterclockwise
      // however here this is only true for one subdomain. From the other we
      // go clockwise. Since the quadrature points are symmetric in [-1,1]
      // (assumption!), we simply reverse their order for the Darcy subdomain
      // by reversing the sign.
      TFEDatabase2D::GetOrigValues(refTransD, -zeta[k], pBfD, eID,
                                   pDref[N_LinePoints - k - 1],
                                   pDxiref[N_LinePoints - k - 1],
                                   pDetaref[N_LinePoints - k - 1], pDorig,
                                   pDxorig, pDyorig);
      
      // quadrature weight and determinant of tranformation
      l.qw = LineWeights[k] * (hE / 2);
      
      // loop over all Darcy pressure basis functions
      for(int row = 0; row < N_pDBf; row++)
      {
        l.at.setTest(pDorig[row], pDxorig[row], pDyorig[row]);
        l.testDOF = pDarcyDOF[row];
        
        // Stokes velocity ansatz functions
        for(int col = 0; col < N_uSBf; col++)
        {
          l.at.setAnsatz(uSorig[col], uSxorig[col], uSyorig[col]);
          l.ansatzDOF = uStokesDOF[col];
          
          localAssemble_pressure_velocity(l);
          
          const int AnsatzDOF = uStokesDOF[col];
          if(AnsatzDOF >= N_Sactive)
          {
            l.m.m[1][AnsatzDOF].clear();
            l.m.m[1][AnsatzDOF + N_SveloDOF].clear();
          }
        }
        // Stokes pressure ansatz functions
        for(int col = 0; col < N_pSBf; col++)
        {
          // x- and y-derivatives not needed
          l.at.setAnsatz(pSorig[col], 1e10, 1e10);
          l.ansatzDOF = pStokesDOF[col];
          
          localAssemble_pressure_pressure(l);
        }
        
        const int TestDOF = pDarcyDOF[row];
        if(TestDOF >= N_Dactive)
        {
          l.m.m[0][TestDOF].clear();
        }
      }
    }
    delete[] pDorig;
    delete[] pDxorig;
    delete[] pDyorig;
    delete[] uSorig;
    delete[] uSxorig;
    delete[] uSyorig;
    delete[] pSorig;
  }
  
  C->add(l.m.m[0]);
  CT->add(l.m.m[1]);
}

/** ************************************************************************ */
void StokesDarcy2D::reorder_rows(std::shared_ptr<TMatrix> m) const
{
  c_stokes()->reorder_rows(m);
}

/** ************************************************************************ */
void StokesDarcy2D::localAssemble_pressure_velocity(
    local_edge_assembling &l) const
{
  const double nu = 1.0 / TDatabase::ParamDB->RE_NR;
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const int N_SveloDOF=c_stokes()->get_velocity_space().GetN_DegreesOfFreedom();
  const double nx = l.nt.nx, ny = l.nt.ny;
  const double q = l.at.v, qx = l.at.vx, qy = l.at.vy; // Darcy pressure
  const double u = l.at.u, ux = l.at.ux, uy = l.at.uy; // Stokes velocity
  const int tDOF = l.testDOF;
  const int aDOF = l.ansatzDOF;
  const double qw = l.qw; // quadrature weight
  switch(c_stokes()->getTypeOf_bci())
  {
    case Neumann:
    {
      //Output::print<1>("localAssemble_pressure_velocity ", tDOF, " ", aDOF,
      //                 "  ---  ", q, "  ", u, "  ", qw, "\t\t\t",
      //                 l.m.m[0][3][5 + N_SveloDOF]);
      l.m.m[0][tDOF][aDOF] += -q * u * nx * qw; // first component
      l.m.m[0][tDOF][aDOF + N_SveloDOF] += -q * u * ny * qw; // second comp
      l.m.m[1][aDOF][tDOF] += q * u * nx * qw; // first component
      l.m.m[1][aDOF + N_SveloDOF][tDOF] += q * u * ny * qw; // second comp
      break;
    }
    case Robin:
    {
      const double gamma_f = TDatabase::ParamDB->StoDa_gamma_f;
      const double gamma_p = TDatabase::ParamDB->StoDa_gamma_p;
      // Darcy pressure test function, Stokes velocity ansatz function
      double val = q * (-u + 2 * nu * (ux * nx + uy * ny) / gamma_p) * qw;
      l.m.m[0][tDOF][aDOF] += val * nx;
      l.m.m[0][tDOF][aDOF + N_SveloDOF] += val * ny;
      // Stokes velocity test function, Darcy pressure ansatz function
      val = (q + K * (qx * nx + qy * ny) * gamma_f) * u * qw;
      l.m.m[1][aDOF][tDOF] += val * nx;
      l.m.m[1][aDOF + N_SveloDOF][tDOF] += val * ny;
      break;
    }
    case weakRobin:
    {
      // could be a different gamma for Stokes and Darcy
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      const double gamma_p = TDatabase::ParamDB->StoDa_gamma_p;
      const double gamma_f = TDatabase::ParamDB->StoDa_gamma_f;
      // Darcy pressure test function, Stokes velocity ansatz function
      double val = (-gamma_p * u + 2 * nu * (ux * nx + uy * ny))
                   * (q - gamma * l.hE * K * (qx * nx + qy * ny)) * qw;
      val /= (gamma_p + gamma * l.hE);
      l.m.m[0][tDOF][aDOF] += val * nx; // first component
      l.m.m[0][tDOF][aDOF + N_SveloDOF] += val * ny; // second component
      // Stokes velocity test function, Darcy pressure ansatz function
      val = (q + gamma_f * K * (qx * nx + qy * ny))
            * (u - gamma * l.hE * 2 * nu * (ux * nx + uy * ny)) * qw;
      val /= (gamma_f + gamma * l.hE);
      l.m.m[1][aDOF][tDOF] += val * nx; // first component
      l.m.m[1][aDOF + N_SveloDOF][tDOF] += val * ny; // second component
      break;
    }
    case DirichletSTAB:
    {
      const double gamma_p = TDatabase::ParamDB->StoDa_gamma_p;
      double val = q
                   * (-u * l.hE / (gamma_p + l.hE)
                      + 2 * nu * (ux * nx + uy * ny) * K * gamma_p / l.hE)
                   * qw;
      l.m.m[0][tDOF][aDOF] += val * nx; // first component
      l.m.m[0][tDOF][aDOF + N_SveloDOF] += val * ny; // second component
                                           
      const double gamma_f = TDatabase::ParamDB->StoDa_gamma_f;
      val = (K * (qx * nx + qy * ny) * nu * gamma_f / l.hE
             + q * l.hE / (gamma_f + l.hE))
            * u * qw;
      l.m.m[1][aDOF][tDOF] += val * nx; // first component
      l.m.m[1][aDOF + N_SveloDOF][tDOF] += val * ny; // second component
      break;
    }
    case Dirichlet:
    {
      // could be a different gamma for Stokes and Darcy
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      double val = gamma * q * nu / l.hE;
      //        gamma ^    psi     2 nu n.DD(u).n / h
      // first component
      l.m.m[0][tDOF][aDOF] += val * (2 * ux * nx + uy * ny) * nx * qw; 
      // second component
      l.m.m[0][tDOF][aDOF + N_SveloDOF] += val * (2 * uy * ny + ux * nx) * ny
                                           * qw;
      val = -gamma * (-K * (qx * nx + qy * ny)) * u * qw / l.hE;
      l.m.m[1][aDOF][tDOF] += val * nx; // first component
      l.m.m[1][aDOF + N_SveloDOF][tDOF] += val * ny; // second component
      break;
    }
    default:
      ErrThrow("unsupported type of Problem\n");
      break;
  }
}

/** ************************************************************************ */
void StokesDarcy2D::localAssemble_pressure_pressure(local_edge_assembling &l)
 const
{
  const double K = TDatabase::ParamDB->SIGMA_PERM;
  const int N_SveloDOF=c_stokes()->get_velocity_space().GetN_DegreesOfFreedom();
  const double nx = l.nt.nx, ny = l.nt.ny;
  const double q = l.at.v, qx = l.at.vx, qy = l.at.vy; // Darcy pressure
  const double p = l.at.u; //, px = l.at.ux, py = l.at.uy; // Stokes pressure
  const int tDOF = l.testDOF;
  const int aDOF = l.ansatzDOF + 2 * N_SveloDOF;
  const double qw = l.qw; // quadrature weight
  
  switch(c_stokes()->getTypeOf_bci())
  {
    case Neumann:
    {
      break;
    }
    case Robin:
    {
      const double gamma_p = TDatabase::ParamDB->StoDa_gamma_p;
      const double val = qw * q * (-p) / gamma_p;
      l.m.m[0][tDOF][aDOF] += val;
      break;
    }
    case weakRobin:
    {
      // could be a different gamma for Stokes and Darcy
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      const double gamma_p = TDatabase::ParamDB->StoDa_gamma_p;
      const double gamma_f = TDatabase::ParamDB->StoDa_gamma_f;
      // Darcy pressure test function, Stokes pressure ansatz function
      double val = (q - gamma * l.hE * K * (qx * nx + qy * ny)) * (-p) * qw;
      val /= (gamma_p + gamma * l.hE);
      l.m.m[0][tDOF][aDOF] += val;
      // Stokes pressure test function, Darcy pressure ansatz function
      val = (q + gamma_f * K * (qx * nx + qy * ny)) * (gamma * l.hE * p) * qw;
      val /= (gamma_f + gamma * l.hE);
      l.m.m[1][aDOF][tDOF] += val;
      break;
    }
    case Dirichlet:
    {
      const double gamma = abs(TDatabase::ParamDB->StoDa_weakGamma);
      const double val = qw * q * (-p) * gamma / l.hE;
      l.m.m[0][tDOF][aDOF] += val;
      break;
    }
    case DirichletSTAB:
    {
      const double val = qw * q * (-p) * K * TDatabase::ParamDB->StoDa_gamma_p
                         / l.hE;
      l.m.m[0][tDOF][aDOF] += val;
      break;
    }
    default:
      ErrThrow("unsupported type of Problem");
      break;
  }
}

/** ************************************************************************ */
void StokesDarcy2D::CombineBigMatrix()
{
  /* combine the matrices into one big matrix
   * ( S  CT )
   * ( C  D  )
   * */
   // Stokes Matrix
  std::shared_ptr<TMatrix> S = c_stokes()->getComposedMatForBigSystem();
  // Darcy  Matrix
  std::shared_ptr<TMatrix> D = c_darcy()->getComposedMatForBigSystem(); 
  
  if(S == NULL || D == NULL)
    ErrThrow("The matrix of one subproblem is not available ");
  
  // S->remove_zeros();
  // D->remove_zeros();
  
  Output::print<1>(
           "coupled System matrix: ( S CT )\n                       ( C  D )");
  Output::print<1>(" Matrix S : Size: ", S->GetN_Rows(), " x ",
                   S->GetN_Columns(), "\t#entries: ", S->GetN_Entries(),
                   "\tnorm: ", S->GetNorm());
  Output::print<1>(" Matrix D : Size: ", D->GetN_Rows(), " x ",
                   D->GetN_Columns(), "\t#entries: ", D->GetN_Entries(),
                   "\tnorm: ", D->GetNorm());
  Output::print<1>(" Matrix CT: Size: ", CT->GetN_Rows(), " x ",
                   CT->GetN_Columns(), "\t#entries: ", CT->GetN_Entries(),
                   "\tnorm: ", CT->GetNorm());
  Output::print<1>(" Matrix C : Size: ", C->GetN_Rows(), " x ",
                   C->GetN_Columns(), "\t#entries: ", C->GetN_Entries(),
                   "\tnorm: ", C->GetNorm());
  
  std::vector<std::shared_ptr<TMatrix>> blocks = {S, CT, C, D};
  big_matrix.reset(new BlockMatrix(2, 2, blocks));
}

/** ************************************************************************ */
void StokesDarcy2D::solve_Stecklov_Poincare(InterfaceFunction& eta) 
{
  // call this routine before, so that the vector 'dofs_to_be_restricted' in
  // eta is created
  eta.restrict(*this);
  
  // compute solution with homogeneous input, so that we can easily compute the 
 
  // linear part of this operator (this operator is only affine linear)
  // eta_hom will be the right hand side of the equation
  {
    // this is essentially like this->apply(eta_hom, eta_hom); except that we
    // do not substract the homogenous solution (which would make no sense here)
    eta_hom = new InterfaceFunction(eta); // initialized with zeros
    
    // reuse eta_hom as homogeneous data to the Darcy and Stokes problem
    if(!StokesFirst)
    { // we have a Stokes-Dirichlet problem
      if(TDatabase::ParamDB->EXAMPLE == 0)
        eta_hom->set_integral(2 * TDatabase::ParamDB->SIGMA_PERM);
      else if(TDatabase::ParamDB->EXAMPLE == 1)
        eta_hom->set_integral(1.0 / 6.0);
      // if one of these cases is true, all additions to eta should have 
      // vanishing integral
    }
    
    darcy()->solve(*eta_hom); // solving Darcy with zero interface data
    stokes()->solve(*eta_hom); // solving Stokes with zero interface data
    
    eta_hom->reset();
    // at least for the Richardson iteration 'eta_hom' can be set to anything
    stokes()->map_solution_to_interface(*eta_hom, -1.0);
    // up to here eta_hom is -H_f,N(0)                (if StokesFirst==true)
    // up to here eta_hom is -H_f,D(0)                (if StokesFirst==false)
    darcy()->map_solution_to_interface(*eta_hom, -1.0);
    // up to here eta_hom is -H_f,N(0) - H_p,D(0)     (if StokesFirst==true)
    // up to here eta_hom is H_f,D(0) + H_p,N(0)      (if StokesFirst==false)
    
    // for the preconditioner we need a linear map, so we have to substract
    // the solution with zero data in every step (otherwise it would be 
    // affine linear only)
    InterfaceFunction * eta_hom_prec = NULL;
    if(StokesFirst)
    {
      InterfaceFunction *& eta_hom_p_prec = p_prec->get_homogeneous_solution();
      if(eta_hom_p_prec == NULL)
      {
        eta_hom_p_prec = new InterfaceFunction(eta); // initialized with zeros
        p_prec->solve(*eta_hom_p_prec); // solve with zero data
        p_prec->map_solution_to_interface(*eta_hom_p_prec, -1.0);
      }
      eta_hom_prec = eta_hom_p_prec;// copy pointer
    }
    else
    {
      InterfaceFunction *& eta_hom_f_prec = f_prec->get_homogeneous_solution();
      if(eta_hom_f_prec == NULL)
      {
        eta_hom_f_prec = new InterfaceFunction(eta); // initialized with zeros
        f_prec->solve(*eta_hom_f_prec); // solve with zero data
        f_prec->map_solution_to_interface(*eta_hom_f_prec, -1.0);
      }
      eta_hom_prec = eta_hom_f_prec; // copy pointer
    }
    // make sure the starting iterate is in the correct affine space
    InterfaceFunction eta_tmp(eta); // initialized with zeros
    eta_tmp = *eta_hom_prec; // copy all values from eta_hom_prec
    // remove all values _not_ corresponding to Darcy-Dirichlet dofs (if  
    // StokesFirst==true, otherwise to Stokes-Dirichlet dofs)
    eta_tmp.restrict(*this, true);
    // remove all values corresponding to Darcy-Dirichlet dofs (if 
    // StokesFirst==true, otherwise to Stokes-Dirichlet dofs)
    eta.restrict(*this, false);
    eta -= eta_tmp;
    // set the correct integral if necessary,
    if(!StokesFirst)
    { // we have a Stokes-Dirichlet problem
      if(TDatabase::ParamDB->EXAMPLE == 0)
        eta.set_integral(2 * TDatabase::ParamDB->SIGMA_PERM);
      else if(TDatabase::ParamDB->EXAMPLE == 1)
        eta.set_integral(1.0 / 6.0);
      // if one of these cases is true, all additions to eta should have 
      // vanishing integral
    }
  }
  
  // something other than 2, which would mean direct solver. We only change 
  // this parameter so that the constructor of the following 'solver' object 
  // 's' does not try to use a direct solver.
  int st = TDatabase::ParamDB->SOLVER_TYPE; // remember whatever was set
  TDatabase::ParamDB->SOLVER_TYPE = 0; 
  
  solver<StokesDarcy2D, InterfaceFunction, StokesDarcy2D> s(this, &eta, eta_hom,
                                                            this);
  switch(TDatabase::ParamDB->StoDa_updatingStrategy)
  {
    case 1: // richardson
      s.set_type(richardson);
      break;
    case 2: // cg
      s.set_type(cg);
      break;
    case 3: // cgs
      s.set_type(cgs);
      break;
    case 4: // gmres
      s.set_type(fgmres);
      break;
    case 5: // bicg
      ErrThrow("bicg solver for Stecklov-Poincare not yet implemented");
      break;
    case 6: // newton
      ErrThrow("Newton method for Stecklov-Poincare not yet implemented");
      break;
    default:
      ErrThrow("unknown updating strategy. Choose a value between 1 and 6 for ",
               "'StoDa_updatingStrategy'.");
      break;
  }
  TDatabase::ParamDB->SOLVER_TYPE = st; // reset solver type
  
  Output::print<1>("start Stecklov-Poincare iteration");
  /*
  int maxit = TDatabase::ParamDB->StoDa_nIterations;
  for(int it = 0; it < maxit; it++)
  {
    Output::print<1>("fixed point iteration ", it);
    stokes()->solve(eta);
    eta.reset();
    stokes()->map_solution_to_interface(eta, 1.0);
    //eta.PrintVals();
    p_prec->solve(eta);
    eta.reset();
    p_prec->map_solution_to_interface(eta, 1.0);
    eta.PrintVals();
  }
  eta.PrintVals("correct eta");
  darcy()->solve(eta);
  stokes()->solve(eta);
  InterfaceFunction eta2(eta);
  darcy()->map_solution_to_interface(eta2, 1.0);
  stokes()->map_solution_to_interface(eta2, 1.0);
  eta2.PrintVals("Stecklov-Poincare operator of correct eta");
  Output::print<1>("Norm of eta ", eta.norm());
  WriteVtk_and_measureErrors(*(stokes()), *(darcy()), NULL, 0);
  
  
  Output::print<1>("\nimitate preconditioned richardson iteration");
  //eta.reset();
  InterfaceFunction *& eta_hom_prec = p_prec->get_homogeneous_solution();
  if(eta_hom_prec == NULL)
    ErrThrow("homogeneous solution of preconditioner not available");
  
  
  InterfaceFunction r(eta), z(eta);
  // this->apply_scaled_add(eta, r, -1);
  darcy()->solve(eta);
  stokes()->solve(eta);
  darcy()->map_solution_to_interface(r, -1);
  stokes()->map_solution_to_interface(r, -1);
  for(int it = 0; it < 12; it++)
  {
    Output::print<1>("\n\niteration ", it);
    r.PrintVals("r");
    // this->solve(r, z);
    z.reset();
    p_prec->solve(r);
    p_prec->map_solution_to_interface(z, -1.0);
    z.PrintVals("z");
    z -= *eta_hom_prec;
    z.PrintVals("z");
    // x += z;
    eta += z;
    eta.PrintVals("eta");
    
    // r = b;
    // this->apply_scaled_add(eta, r, -1.0);
    r.reset();
    darcy()->solve(eta);
    stokes()->solve(eta);
    darcy()->map_solution_to_interface(r, -1);
    stokes()->map_solution_to_interface(r, -1);
    r.restrict(*this);
  }
  eta.PrintVals("resulting eta");
  darcy()->solve(eta);
  stokes()->solve(eta);
  WriteVtk_and_measureErrors(*(stokes()), *(darcy()), NULL, 1);
  
  exit(0);
  
  eta.restrict(*(stokes()), true);*/
  
  s.solve(TDatabase::ParamDB->StoDa_nIterations, 
          TDatabase::ParamDB->StoDa_bigResidual);
  
  Output::print<1>("finished Stecklov-Poincare iteration\n");
  
  eta.PrintGnuplotFile((char*)"hallo.gp");
  
  if(TDatabase::ParamDB->StoDa_updatingStrategy != 1)
  {
    // except for the Richardson iteration we apply the linear operator to the
    // residual, so that the eta is correct while the solution vectors for
    // velocity and pressure are not. We do another solve here, even though
    // there is probably a more elagant way to do this
    stokes()->solve(eta);
    darcy()->solve(eta);
  }
  
  if(TDatabase::ParamDB->StoDa_solutionStrategy >= 0)
  { // compute defect of big system:
    // r = ( S  CT ) ( u,p ) _ ( f_f )
    //     ( C  D  ) ( phi )   ( f_p )
    BlockVector r(*big_matrix);
    r.copy(stokes()->get_solution().get_entries(), 0);
    r.copy(darcy()->get_solution().get_entries(), 1);
    r *= *big_matrix;
    r -= *big_rhs;
    bigResidual = norm(r);
    Output::print<1>(" Norm of big system defect ", bigResidual);
  }
}

/** ************************************************************************ */
void StokesDarcy2D::apply(const InterfaceFunction & x,
                          InterfaceFunction & y) const
{
  if(y.length() == 0)
    y = *eta_hom;
  y.reset();
  
  y.add(eta_hom, 1.0);
  
  darcy()->solve(x); // Interface integrals & solving for the Darcy part
  stokes()->solve(x);// Interface integrals & solving for the Stokes part
  
  darcy()->map_solution_to_interface(y, 1.0); // restricting the solution to y
  stokes()->map_solution_to_interface(y, 1.0);// restricting the solution to y
  
  y.restrict(*this); // values at Dirichlet dofs can be set to zero
  
  //this->apply_scaled_add(x, y, 1.0);
}

/** ************************************************************************ */
void StokesDarcy2D::apply_scaled_add(const InterfaceFunction & x,
                                     InterfaceFunction & y, double a) const
{
  // y = y + a * Ax
  Output::print<1>("StokesDarcy2D::apply_scaled_add");
  // substract homogenous solution such that this is a linear operator
  y.add(eta_hom, a);
  
  darcy()->solve(x); // Interface integrals & solving for the Darcy part
  stokes()->solve(x);// Interface integrals & solving for the Stokes part
  
  darcy()->map_solution_to_interface(y, a); // restricting the solution to y
  stokes()->map_solution_to_interface(y, a);// restricting the solution to y

  y.restrict(*this); // values at Dirichlet dofs can be set to zero
}

/** ************************************************************************ */
void StokesDarcy2D::check_equalities() const
{
  double time = GetTime();
  int problem_type = TDatabase::ParamDB->StoDa_problemType;
  int solution_strategy = TDatabase::ParamDB->StoDa_solutionStrategy;
  int updating_strategy = TDatabase::ParamDB->StoDa_updatingStrategy;
  const bool C_RR = problem_type == 0 && solution_strategy == 1
                    && updating_strategy == 3;
  
  if(solution_strategy == 1 && !C_RR)
  {
    // solve the coupled problem with blockwise Gauss-Seidel or Jacobi
    // (for C_RR we can not check this)
    
    // check if the coupling matrices C and CT which are assembled during
    // the direct solution are equal to the product of sol2eta and etaToBd
    // this is true analytically and necessary so that the Jacobi and 
    // Gauss-Seidel type methods really are just that.
    
    std::shared_ptr<TMatrix> CT = c_stokes()->getCouplingMatrix();
    std::shared_ptr<TMatrix> C  = c_darcy()->getCouplingMatrix();
    std::shared_ptr<TMatrix> s_etaToBd = stokes()->get_etaToBd();
    std::shared_ptr<TMatrix> s_sol2eta = stokes()->get_map_sol2eta();
    std::shared_ptr<TMatrix> d_etaToBd = darcy()->get_etaToBd();
    std::shared_ptr<TMatrix> d_sol2eta = darcy()->get_map_sol2eta();
    
    const double tol = -1;//1e-14; // tol < 0 --> choose an appropriate value 
    
    s_etaToBd->remove_zeros(tol);
    s_sol2eta->remove_zeros(tol);
    d_etaToBd->remove_zeros(tol);
    d_sol2eta->remove_zeros(tol);
    
    
    TMatrix * CT2 = s_etaToBd->multiply(d_sol2eta.get());
    TMatrix * C2  = d_etaToBd->multiply(s_sol2eta.get());
    
    CT2->remove_zeros(tol);
    C2->remove_zeros(tol);
    
    //CT2->Printfile("CT2.txt");
    //C2->Printfile("C2.txt");
    //CT->Printfile("CT.txt");
    //C->Printfile("C.txt");
    //s_etaToBd->Printfile("s_etaToBd.txt");
    //s_sol2eta->Printfile("s_sol2eta.txt");
    //d_etaToBd->Printfile("d_etaToBd.txt");
    //d_sol2eta->Printfile("d_sol2eta.txt");
    
    //CT2->PrintFull("CT2", 2);
    //C2->PrintFull("C2", 2);
    //CT->PrintFull("CT", 2);
    //C->PrintFull("C", 2);
    //s_etaToBd->PrintFull("s_etaToBd", 2);
    //s_sol2eta->PrintFull("s_sol2eta", 2);
    //d_etaToBd->PrintFull("d_etaToBd", 2);
    //d_sol2eta->PrintFull("d_sol2eta", 2);
    
    *C2 += *C.get();
    if(C2->GetNorm() > 1e-10)
    {
      ErrThrow("Error:  C is not d_etaToBd * s_sol2eta ", C2->GetNorm());
    }
    Output::print<1>("coupling matrix C test succesfull");
    
    *CT2 += *CT;
    if(CT2->GetNorm() > 1e-10)
    {
      ErrMsg("Error: CT is not s_etaToBd * d_sol2eta " << CT2->GetNorm() );
      CT2->remove_zeros(1e-14);
      CT2->Print("CT-diff");
      ErrThrow("check_equalities failed");
    }
    Output::print<1>("coupling matrix CT test succesfull");
    
    delete CT2;
    delete C2;
  }
  if(solution_strategy == 3) //Stecklov-Poincare
  {
    // check whether solving a Neumann problem with (Neumann) data from another
    // solution yields the same solution. Same for Dirichlet. Depending on
    // 'StokesFirst' this is done for Stokes or Darcy
    InterfaceFunction eta(darcy()->getInterface(), 2);
    Output::print<1>("\n\nTesting identities");
    if(StokesFirst == 1)
    {
      eta.reset();
      // set to something other than zero
      eta = 5 / ((eta.length() - 1.0) / 2.0); 
      DarcyPrimal* d_D = darcy(); // Dirichlet Darcy problem
      DarcyPrimal* d_N = p_prec;  // Neumann Darcy problem
      StokesProblem * ns_N = stokes(); // Neumann Stokes problem
      
      d_D->solve(eta);
      eta.reset();
      d_D->map_solution_to_interface(eta, -1.0);
      //WriteVtk_and_measureErrors(*(stokes()), *d_D, NULL, 0); // dummy Stokes
      
      d_N->solve(eta);
      eta.reset();
      d_N->map_solution_to_interface(eta);
      //WriteVtk_and_measureErrors(*(stokes()), *d_N, NULL, 1); // dummy Stokes
      {
        BlockVector *D_sol = &d_D->get_solution();
        BlockVector *N_sol = &d_N->get_solution();
        *D_sol -= *N_sol;
        if(D_sol->norm() > 1e-12 * sqrt(D_sol->length()))
          Error("Difference(1) which should be zero " << D_sol->norm() << endl);
      }
      
      d_D->solve(eta);
      //WriteVtk_and_measureErrors(*(stokes()), *d_D, NULL, 2); // dummy Stokes
      {
        BlockVector *D_sol = &d_D->get_solution();
        BlockVector *N_sol = &d_N->get_solution();
        *D_sol -= *N_sol;
        if(D_sol->norm() > 1e-12 * sqrt(D_sol->length()))
          Error("Difference(2) which should be zero " << D_sol->norm() << endl);
      }
      check_linearity(d_N, eta);
      check_linearity(d_D, eta);
      check_linearity(ns_N, eta);
    }
    else
    { // Darcy first
      eta.reset();
      // we have to call restrict once before we call 'eta.set_integral' 
      eta.restrict(*this);
      
      // set to something other than zero
      eta = 5.0 / ((eta.length() - 1.0) / 2.0); 
      StokesProblem * ns_D = stokes(); // Dirichlet Stokes problem
      StokesProblem * ns_N = f_prec;   // Neumann Stokes problem
      DarcyPrimal* d_N = darcy(); // Neumann Darcy problem
      
      // for Dirichlet boundary conditions on the entire outer Stokes boundary
      // we have to be sure the input (flow) is in the correct space
      if(TDatabase::ParamDB->EXAMPLE == 0)
        eta.set_integral(2 * TDatabase::ParamDB->SIGMA_PERM);
      if(TDatabase::ParamDB->EXAMPLE == 1)
        eta.set_integral(1.0 / 6.0);
      
      ns_D->solve(eta);
      eta.reset();
      ns_D->map_solution_to_interface(eta, -1.0);
      //WriteVtk_and_measureErrors(*ns_D, *(darcy()), NULL, 0); // dummy darcy
      
      ns_N->solve(eta);
      eta.reset();
      ns_N->map_solution_to_interface(eta, 1.0);
      //WriteVtk_and_measureErrors(*ns_N, *(darcy()), NULL, 1); // dummy darcy
      {
        BlockVector *D_sol = &ns_D->get_solution();
        BlockVector *N_sol = &ns_N->get_solution();
        *D_sol -= *N_sol;
        if(D_sol->norm() > 1e-12 * sqrt(D_sol->length()))
          Error("Difference(1) which should be zero " << D_sol->norm() << endl);
      }
      
      ns_D->solve(eta);
      //WriteVtk_and_measureErrors(*ns_D, *(darcy()), NULL, 2); // dummy darcy
      {
        BlockVector *D_sol = &ns_D->get_solution();
        BlockVector *N_sol = &ns_N->get_solution();
        *D_sol -= *N_sol;
        if(D_sol->norm() > 1e-12 * sqrt(D_sol->length()))
          Error("Difference(2) which should be zero " << D_sol->norm() << endl);
      }
      /////////////////////////////////////////////////////////////////////////
      check_linearity(ns_N, eta);
      check_linearity(ns_D, eta);
      check_linearity(d_N, eta);
    }
    Output::print<1>("Testing identities finished\n");
  }
  Output::print<1>("checking equalities done within ", GetTime() - time,
                   " seconds\n");
}

/** ************************************************************************ */
double ErrorOnInterface(StokesProblem &s, DarcyPrimal &d, int it)
{
  Output::print<1>(" Compute error on interface");
  
  // normal vector to interface (pointing out of Stokes subdomain)
  double nx, ny; 
  double x0, y0, x1, y1; // coordinates of some vertices of a cell
      
  double s_u1val[3], s_u2val[3], s_pval[3]; // velocity and pressure for Stokes
  double d_pval[3]; // pressure for Darcy
  double s_u_n, d_u_n; // normal component of velocity for Stokes and Darcy
  double nTn; // normal component of normal stress (Stokes)
  double val1; // just a number to store intermediate values
  // Stokes velocity (two components)
  const TFEFunction2D *s_u1 = ((it != -1) ? s.get_velocity() : s.get_u_Direct() 
                               ).GetComponent(0);
  const TFEFunction2D *s_u2 = ((it != -1) ? s.get_velocity() : s.get_u_Direct()
                               ).GetComponent(1);
  const TFEFunction2D & s_p = (it != -1) ? s.get_pressure() : s.get_p_direct();
  const TFEFunction2D & d_p = (it != -1) ? d.getP() : d.get_p_direct();
  
  const double nu = 1 / TDatabase::ParamDB->RE_NR; // diffustion coefficient
  const double K = TDatabase::ParamDB->SIGMA_PERM; // Permeability
  double error1 = 0; // error in normal component of velocity
  double error2 = 0; // error of n.T.n + p_D
  // degree of accuracy for quadrature formula
  const double degree = 2 * TDatabase::ParamDB->VELOCITY_SPACE;
  // use same quadrature formula on every segment of the interface
  int N_LinePoints; // number of points for 1D quadrature formula
  double *LineWeights, *zeta; // weights and points of 1D quadrature formula
  QuadFormula1D LineQF = TFEDatabase2D::GetQFLineFromDegree(degree);
  TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQF);
  qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
  
  // all cells in the Stokes subdomain have the same id
  const int S_ID = s.get_pressure_space().GetCollection()->GetCell(0)
      ->GetReference_ID(); // ID of Stokes subdomain
      
  const bool writeErrorFile = true;
  std::ofstream ErrorFile;
  if(writeErrorFile)
  {
    std::ostringstream os(std::ios::trunc);
    os << "interfaceFiles/interfaceErrors" << it << ".txt" << ends;
    ErrorFile.open(os.str().c_str());
    ErrorFile << "x\t" << "DarcyFlux\t" << "StokesFlux\t" << "DarcyStress\t"
              << "-StokesStress\t" << "FluxError\t" << "StressError\t"
              << "StokesPressure\n";
  }
  for(unsigned int j = 0; j < s.getInterface().size(); j++)
  {
    // interface joints (segments)
    const TInnerInterfaceJoint *IJoint = s.getInterface()[j];
    
    // neighboring cells in Stokes and Darcy subdomain
    TBaseCell *s_cell, *d_cell; 
    s_cell = IJoint->GetNeighbour(0);
    if(s_cell->GetReference_ID() != S_ID)
    {
      d_cell = s_cell;
      s_cell = IJoint->GetNeighbour(1);
    }
    else
      d_cell = IJoint->GetNeighbour(1);
    
    // compute length of the edge
    const double hE = IJoint->GetLength();
    // normal vectors to this boundary (normalized)
    getNormal(s_cell, IJoint, nx, ny);
    {
      // index of an edge in a cell (for triangles: 0,1, or 2)
      const int iEdge = IJoint->GetIndexInNeighbor(s_cell);
      const int N_Joints = s_cell->GetN_Edges();
      s_cell->GetVertex(iEdge)->GetCoords(x0, y0);
      s_cell->GetVertex((iEdge + 1) % N_Joints)->GetCoords(x1, y1);
    }
    
    // loop of quadrature points
    for(int k = 0; k < N_LinePoints; k++)
    {
      //zeta \in [-1,1]
      const double x = x0 + ((zeta[k] + 1.0) / 2.0) * (x1 - x0);
      const double y = y0 + ((zeta[k] + 1.0) / 2.0) * (y1 - y0);
      
      // get the function values
      s_u1->FindGradientLocal(s_cell, s_cell->GetCellIndex(), x, y, s_u1val);
      s_u2->FindGradientLocal(s_cell, s_cell->GetCellIndex(), x, y, s_u2val);
      s_p.FindGradientLocal(s_cell, s_cell->GetCellIndex(), x, y, s_pval);
      d_p.FindGradientLocal(d_cell, d_cell->GetCellIndex(), x, y, d_pval);
      
      // normal component of Stokes velocity
      s_u_n = s_u1val[0] * nx + s_u2val[0] * ny;
      // normal component of Darcy velocity
      d_u_n = -K * (d_pval[1] * nx + d_pval[2] * ny);
      // difference in normal component of velocity
      val1 = s_u_n - d_u_n;
      error1 += LineWeights[k] * (hE / 2) * val1 * val1;
      
      // normal component of normal stress (Stokes)
      nTn = nx * s_u1val[1] * nx;
      nTn += nx * (s_u1val[2] + s_u2val[1]) * ny;
      nTn += ny * s_u2val[2] * ny;
      nTn *= 2 * nu; // *nu
      nTn -= s_pval[0];
      // difference to -p_D
      val1 = nTn + d_pval[0];
      error2 += LineWeights[k] * (hE / 2) * val1 * val1;
      //Output::print<1>("x ", x, " uS.n ", s_u_n, ", uD.n ", d_u_n,
      //                 ", nTn ", nTn, ", phi ", d_pval[0])
      if(writeErrorFile)
        ErrorFile << x << "\t" << d_u_n << "\t" << s_u_n << "\t" << d_pval[0]
                  << "\t" << -nTn << "\t" << s_u_n - d_u_n << "\t"
                  << nTn + d_pval[0] << "\t" << s_pval[0] << endl;
    } //for(int k=0;k<N_LinePoints;k++)
    if(writeErrorFile)
      ErrorFile << endl;
  }
  if(writeErrorFile)
    ErrorFile.close();
  
  delete s_u1;
  delete s_u2;
  
  Output::print<1>(" Flux-Error ", sqrt(error1), ", normal-stress-error ",
                   sqrt(error2), ", quotient ", sqrt(error1)/sqrt(error2));
  return sqrt(error1 + error2);
}

/** ************************************************************************ */
void WriteVtk_and_measureErrors(StokesProblem &s, DarcyPrimal &d,
                                TDomain *Domain, int it)
{
  TFEVectFunct2D & u_NSE = it != -1 ? s.get_velocity() : s.get_u_Direct();
  const TFEFunction2D *u1_NSE, *u2_NSE;
  u1_NSE = u_NSE.GetComponent(0);
  u2_NSE = u_NSE.GetComponent(1);
  const TFEFunction2D & p_NSE = it != -1 ? s.get_pressure() : s.get_p_direct();
  const TFEFunction2D & p_Darcy = (it != -1) ? (d.getP()) : (d.get_p_direct());
  /* ========================================================================
   * print minimal and maximal values of all FE - functions                
   * does not work for Raviart-Thomas elements (yet)                       */
  if(TDatabase::ParamDB->SC_VERBOSE)
  {
    p_NSE.PrintMinMax(" Navier-Stokes pressure");
    u1_NSE->PrintMinMax(" Navier-Stokes u1");
    u2_NSE->PrintMinMax(" Navier-Stokes u2");
    p_Darcy.PrintMinMax(" Darcy");
  }
  /* ==========================================================================
   *  write vtk-files, one for each subdomain                                */
  if(TDatabase::ParamDB->WRITE_VTK || TDatabase::ParamDB->WRITE_GNU)
  {
    // create one TOutput2D object for each subdomain
    TOutput2D Output_NSE(2, 1, 1, 0, Domain);
    TOutput2D Output_Darcy(1, 1, 0, 0, Domain);
    Output_NSE.AddFEFunction(&p_NSE);
    Output_NSE.AddFEVectFunct(&u_NSE);
    Output_Darcy.AddFEFunction(&p_Darcy);
    // write one vtk - file for each subdomain
    std::ostringstream os(std::ios::trunc);
    if(it != -1)
      os << TDatabase::ParamDB->OUTPUTDIR << "/" 
         << TDatabase::ParamDB->BASENAME << "NSE." << it << ".vtk" << ends;
    else
      os << TDatabase::ParamDB->OUTPUTDIR << "/" 
         << TDatabase::ParamDB->BASENAME << "NSE.direct.vtk" << ends;
    if(TDatabase::ParamDB->WRITE_VTK)
    {
      Output_NSE.WriteVtk(os.str().c_str());
    }
    os.seekp(int(os.tellp()) - 4);
    os << "gnu" << ends;
    if(TDatabase::ParamDB->WRITE_GNU)
      Output_NSE.WriteGnuplot(os.str().c_str());
    
    os.seekp(std::ios::beg);
    if(it != -1)
      os << TDatabase::ParamDB->OUTPUTDIR << "/" 
         << TDatabase::ParamDB->BASENAME << "Darcy." << it << ".vtk" << ends;
    else
      os << TDatabase::ParamDB->OUTPUTDIR << "/" 
         << TDatabase::ParamDB->BASENAME << "Darcy.direct.vtk" << ends;
    if(TDatabase::ParamDB->WRITE_VTK)
    {
      Output_Darcy.WriteVtk(os.str().c_str());
    }
    os.seekp(int(os.tellp()) - 4);
    os << "gnu" << ends;
    if(TDatabase::ParamDB->WRITE_GNU)
      Output_Darcy.WriteGnuplot(os.str().c_str());
  }
  /* ==========================================================================
   * measure errors                                                          */
  if(TDatabase::ParamDB->MEASURE_ERRORS)
  {
    TAuxParam2D aux;
    ErrorMethod2D *ErrorMeth = L2H1Errors; // in MainUtilities.C
    double err[6];
    MultiIndex2D AllDerivatives[3] = {D00, D10, D01};
    const TFESpace2D *p_space_NSE = p_NSE.GetFESpace2D();
    const TFESpace2D *v_space_NSE = u1_NSE->GetFESpace2D();
    const TFESpace2D *p_space_Darcy = d.getP().GetFESpace2D();
    
    u1_NSE->GetErrors(s.get_example().get_exact(0), 3, AllDerivatives, 2,
                      ErrorMeth, s.get_example().get_coeffs(), &aux, 1,
                      &v_space_NSE, err);
    //Output::print<1>("  L2(u1_S):      ",  err[0]);
    //Output::print<1>("  H1-semi(u1_S): ",  err[1]);
    
    u2_NSE->GetErrors(s.get_example().get_exact(1), 3, AllDerivatives, 2,
                      ErrorMeth, s.get_example().get_coeffs(), &aux, 1,
                      &v_space_NSE, err + 3);
    //Output::print<1>("  L2(u2_S):      ", err[3]);
    //Output::print<1>("  H1-semi(u2_S): ", err[4]);
    Output::print<1>(" L2(u_S):       ", sqrt(err[0]*err[0]+err[3]*err[3]));
    Output::print<1>(" H1-semi(u_S):  ", sqrt(err[1]*err[1]+err[4]*err[4]));
    
    p_NSE.GetErrors(s.get_example().get_exact(2), 3, AllDerivatives, 2,
                     ErrorMeth, s.get_example().get_coeffs(), &aux, 1,
                     &p_space_NSE, err);
    Output::print<1>(" L2(p_S):       ", err[0]);
    Output::print<1>(" H1-semi(p_S):  ", err[1]);
    
    DoubleFunct2D * const *Exact_Darcy = d.get_exact();
    p_Darcy.GetErrors(Exact_Darcy[0], 3, AllDerivatives, 2, ErrorMeth,
                      d.get_example().get_coeffs(), &aux, 1, &p_space_Darcy,
                      err);
    Output::print<1>(" L2(p_D):       ", err[0]);
    Output::print<1>(" H1-semi(p_D):  ", err[1]);
  }
  
  /* ==========================================================================
   * measure discrete errors                                                 */
  if(TDatabase::ParamDB->StoDa_solutionStrategy >=0 && it != -1)
  { // direct solve has been done before
    double relError;
    double tot_relError = 0;
    // Stokes velocity
    relError = relDiff(2 * u_NSE.GetLength(), s.get_u_Direct().GetValues(),
                       u_NSE.GetValues());
    Output::print<1>(" discrete Stokes velocity Error ", relError);
    tot_relError += relError*relError;
    
    // Stokes Pressure
    relError = relDiff(p_NSE.GetLength(), s.get_p_direct().GetValues(),
                       p_NSE.GetValues());
    Output::print<1>(" discrete Stokes pressure Error ", relError);
    tot_relError += relError*relError;
    
    // Darcy pressure
    relError = relDiff(p_Darcy.GetLength(), d.get_p_direct().GetValues(),
                       p_Darcy.GetValues());
    Output::print<1>(" discrete Darcy pressure Error  ", relError);
    tot_relError += relError*relError;
    
    Output::print<1>(" discrete total error           ", sqrt(tot_relError));
  }
  delete u1_NSE;
  delete u2_NSE;
}

/** ************************************************************************ */
void ComputeResiduals(StokesProblem &s, DarcyPrimal &d)
{
  Output::print<1>(" Computing defects");
  
  { // Stokes defects
    BlockVector r(s.get_matrix(), true);
    r = s.get_solution();
    r *= s.get_matrix(); // TODO: is matrix different??
    r -= s.get_rhs();
    
    Output::print<1>("  Norm of Stokes defect ", r.norm());
    Output::print<1>("                 velocity first component  ",
                     Dnorm(r.length(0), r.block(0)));
    Output::print<1>("                 velocity second component ",
                     Dnorm(r.length(1), r.block(1)));
    Output::print<1>("                 pressure component        ",
                     Dnorm(r.length(2), r.block(2)));
  }
  
  { // Darcy defects
    BlockVector r(d.getMat(), true);
    r = d.get_solution();
    r *= d.getMat();
    r -= d.get_rhs();
    
    Output::print<1>("  Norm of Darcy defect ", r.norm());
  }
}

