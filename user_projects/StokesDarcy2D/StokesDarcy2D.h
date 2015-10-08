#ifndef STOKESDARCY2D_H
#define STOKESDARCY2D_H

class StokesDarcy2D;

#include <StokesProblem.h>
#include <DarcyPrimal.h>
#include <Example_StokesDarcy2D.h>

#include <solver.h>


class StokesDarcy2D
{
  private:
  
  // subproblems of the coupled problem
  StokesProblem* coupled_stokes;
  DarcyPrimal* coupled_darcy;
  
  // for iteration (if not C-RR this is equal to coupled_stokes, coupled_darcy)
  StokesProblem* iteration_stokes;
  DarcyPrimal* iteration_darcy;
  
  bool StokesFirst;
  // error on interface (from previous iteration step)
  double error, initialError;
  // residual of big system using the iterative solution
  double bigResidual, initialBigResidual;
  
  // The operator (either in the fixed point or Stecklov-Poincare formulation)
  // represented by this class is affine linear. In order to use the iterative
  // solvers, we need the linear part of this operator. To achieve this we save 
  // the homogeneous (nonzero) solution. So that the linear part of this
  // operator can be computed by substracting this homogeneous solution.
  InterfaceFunction* eta_hom;
  
  /** store a preconditioner object. 
   * The preconditioner is either a free flow solve, a porous medium solve, or
   * some linear combination of both. 
   */
  StokesProblem *f_prec; // serve as preconditioner in free flow subdomain
  DarcyPrimal *p_prec; // serve as preconditioner in porous medium subdomain
  
  /** the system which we are actually trying to solve */
  std::shared_ptr<BlockMatrix> big_matrix;
  std::shared_ptr<BlockVector> big_rhs;
  std::shared_ptr<BlockVector> big_solution;
  
  // full system matrix is 
  // big_matrix = ( S  CT )
  //              ( C  D  )
  std::shared_ptr<TMatrix> C, CT;
  
  public:
  
  /** constructor */
  StokesDarcy2D(std::map<InterfaceCondition, StokesProblem*> ns_problems,
                std::map<InterfaceCondition, DarcyPrimal*> d_problems);  
  
  ~StokesDarcy2D();
  
  // appyly operator H_f(H_p(.))
  void solve_fixed_point(InterfaceFunction& eta) const;

  // apply this class as preconditioner
  void solve(const InterfaceFunction& eta_f, InterfaceFunction& eta_p) const;
  
  // apply this class as preconditioner (within fgmres)
  void solve(int i, int j, const InterfaceFunction& eta_f, 
             InterfaceFunction& eta_p) const
  { solve(eta_f, eta_p); }
  
  // solve, using two interface functions (used to solve coupled system)
  void solve_coupled_system(InterfaceFunction& eta_f,
                            InterfaceFunction& eta_p) const;

  void solveDirect();

  bool stopIteration(int it);
  
  void solve_Stecklov_Poincare(InterfaceFunction& eta);
  
  /*
   * the following functions are implemented so that the template iterative 
   * solvers will work, which is used for the Stecklov-Poincare formulation
   */
  void apply(const InterfaceFunction & x, InterfaceFunction & y) const;
  // y = y + a * Ax
  void apply_scaled_add(const InterfaceFunction & x, InterfaceFunction & y,
                        double a = 1.0) const;
  
  // depending on which problem is solved some identities are checked.
  void check_equalities() const;
  
  
  /*
   * the following function is implemented so that the compiler will not 
   * complain. This function would be needed if this operator was a matrix and
   * we wanted to use a direct solver. This is impossible here.
   */
  std::shared_ptr<TMatrix> get_combined_matrix() const
  {
    ErrThrow("The class StokesDarcy2D is not given as a matrix, therefore "
             + "get_combined_matrix() makes no sense here");
    throw; // only to avoid a compiler warning
  }
  
  // getters and setters
  
  // get individual problems for iteration
  StokesProblem* stokes() const {return iteration_stokes;};
  StokesProblem* stokes() {return iteration_stokes;};
  DarcyPrimal* darcy() const {return iteration_darcy;};
  DarcyPrimal* darcy() {return iteration_darcy;};
  
  // get individual problems for coupled system
  StokesProblem* c_stokes() const {return coupled_stokes;};
  StokesProblem* c_stokes() {return coupled_stokes;};
  DarcyPrimal* c_darcy() const {return coupled_darcy;};
  DarcyPrimal* c_darcy() {return coupled_darcy;};
  
  bool Stokes_first() const { return StokesFirst; }
  
  private:
  /**
   * Set up the combined matrix for the Stokes-Darcy system and solve
   */
  void SolveOneSystem();

  /* Return the coupling matrices which couple the Stokes and Darcy part.
   * The Matrix C has dimensions nxm where n is the number of Darcy degrees of
   * freedom and m is the number of Stokes degrees of freedom. The matrix CT has
   * dimensions mxn.
   * A call to AssembleCouplingMatrices(...) is made.
   */
  void GetCouplingMatrix();

  /**
   * Assemble the matrices coupling the Stokes to the Darcy part.
   * Rows corresponding to Dirichlet degrees of freedom are not assembled!
   */
  void AssembleCouplingMatrices() const;
  /**
   * Write the Stokes matrix 'S', the Darcy matrix 'D' and the coupling 
   * matrices 'C' and 'CT' into one big matrix which represents the discrete
   * Stokes-Darcy coupling
   *
   * (              |    )
   * (      S       | CT )
   * ( _____________|    )
   * (      C       | D  )
   */
  void CombineBigMatrix();

  /*
   * local assembling routines used to assemble the coupling matrices C, CT 
   * the first term (velocity or pressure) is for Darcy test functions and the 
   * second term (velocity or pressure) is for Stokes ansatz functions
   */
  void localAssemble_velocity_velocity(local_edge_assembling &lea) const;
  void localAssemble_velocity_pressure(local_edge_assembling &lea) const;
  void localAssemble_pressure_velocity(local_edge_assembling &lea) const;
  void localAssemble_pressure_pressure(local_edge_assembling &lea) const;

  /*
   * reorder coupling matrices (at least in the one with Stokes test space)  
   * such that for each interface dof, one row corresponds to the normal
   * component and one corresponds to the tangetial component. Otherwise 
   * interfaces not aligned with an axis won't work. To do this we will call
   * TMatrix::changeRows() on each matrix with velocity test space
   */
  void reorder_rows(std::shared_ptr<TMatrix> m) const;
};


/** Description of the function ComputeErrorOnInterface:
 * Evaluates ||u_f.n - u_p.n||_L2(\Gamma) + ||n.T.n + phi_p||_L2(\Gamma)
 * where u_f is the velocity in the Stokes part, u_p is the velocity in the 
 * Darcy part (if RT==false, u_p = -K (grad phi_p).n), phi_p is the Darcy 
 * pressure and T = 2\nu DD(u_f) - p_f I. Here p_f is the Stokes pressure, I is 
 * the identity matrix, and DD(u) = (grad u + (grad u)^T)/2 is the deformation
 * tensor.   
 */
double ErrorOnInterface(StokesProblem &s, DarcyPrimal &d, int iteration=0);

/** Description of the function WriteVtk_and_measureErrors:
 * to keep the main program short
 * the vtk files are created and errors are measured and written to console
 *
 * The parameter 'Domain' can be set to 'NULL', at least in all cases I tested.
 * The parameter 'it' should be set to -1 if only the direct solution has been
 * computed so far.
 */
void WriteVtk_and_measureErrors(StokesProblem &s, DarcyPrimal &d,
                                TDomain *Domain, int it);

/** Description of the function ComputeResiduals:
 * to keep the main program short
 * this computes residuals just for testing 
 */
void ComputeResiduals(StokesProblem &s, DarcyPrimal &d);


#endif // STOKESDARCY2D_H