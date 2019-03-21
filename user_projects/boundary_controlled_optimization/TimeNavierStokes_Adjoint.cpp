#include "TimeNavierStokes_Adjoint.hpp"
#ifdef __2D__
#include "Assemble2D.h"
#else
#include "Assemble3D.h"
#endif
#include "Hotfixglobal_AssembleNSE.h"
#include "MainUtilities.h" // BoundaryValueHomogenous
#include "Database.h" // to check TDatabase::ParamDB->NSTYPE

namespace tnse_adjoint
{
void zero_solution(double x, double y, 
#ifdef __3D__     
                   double z,
#endif
                   double *values);
template <int d>
void adjoint_assembling(double, double*, double*, double, double**, int*,
                        double***, double**);
void params_function(const double *in, double *out);
// very ugly to put this here, but I don't know how to get this into the local
// assembling
std::vector<double> cost_functional_weights;
bool restricted_curl_functional;
}


constexpr double diagonal_scaling = 1.e30;
void tnse_adjoint::zero_solution(double, double , 
#ifdef __3D__     
                   double,
#endif
                   double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
#ifdef __3D__
  values[4] = 0;
#endif
}

template<int d>
TimeNavierStokes_Adjoint<d>::TimeNavierStokes_Adjoint(const TimeNavierStokes<d>& tnse,
                             const ParameterDatabase& param_db)
 : TimeNavierStokes<d>(tnse) // copy constructor
{
  // copy remaining parts
  this->TimeNavierStokes<d>::db.merge(param_db, false);
  bool usingMultigrid = this->TimeNavierStokes<d>::solver.is_using_multigrid();
  if(usingMultigrid)
  {
    ErrThrow("Time_NSE2D_Adjoint::assemble_additional_terms not yet implemented for "
             "multigrid");
  }
  std::vector<DoubleFunction*> adjoint_solutions(d+1, tnse_adjoint::zero_solution);
  std::vector<DoubleFunction*> initial_conditions(d+1, tnse_adjoint::zero_solution);
  std::vector<BoundaryValuesFunction*> adjoint_bd(d+1, BoundaryValueHomogenous);
  this->TimeNavierStokes<d>::example = Example_TimeNSE(adjoint_solutions,
                                       this->TimeNavierStokes<d>::example.boundary_conditions,
                                       adjoint_bd,
                                       this->TimeNavierStokes<d>::example.get_coeffs(),
                                       true, true, initial_conditions);
//                                       this->TimeNavierStokes<d>::example.get_nu());
  this->TimeNavierStokes<d>::outputWriter = DataWriter<d>(param_db);
  this->TimeNavierStokes<d>::outputWriter.add_fe_vector_function(&this->get_velocity());
  this->TimeNavierStokes<d>::outputWriter.add_fe_function(&this->get_pressure());
  this->TimeNavierStokes<d>::solver = Solver<BlockFEMatrix, BlockVector>(param_db);
}

template<int d>
void TimeNavierStokes_Adjoint<d>::assemble_initial_time(const FEVectFunct& u,
                                                        const FEFunction& p)
{    
  // Preliminary verifications
  if(this->systems.size() > 1)
  {
    ErrThrow("TimeNavierStokes_Adjoint::assemble_additional_terms not yet implemented for "
             "multigrid");
  }
  if(TDatabase::ParamDB->NSTYPE != 4)
  {
    ErrThrow("The adjoint problem requires a full matrix structure. Please set "
            "'NSTYPE' to 4");
  }

  //=======================================================================
  // Assemble linear terms from TimeNavierStokes class and mass matrix
  for(auto &s : this->systems)
  {
    assemble_nse = Hotfixglobal_AssembleNSE::WITHOUT_CONVECTION;
    this->call_assembling_routine(s, LocalAssembling_type::TNSE3D_LinGAL);
    
    // copy nonactives
    s.solution.copy_nonactive(s.rhs);

  /** After copy_nonactive, the solution vectors needs to be Comm-updated in 
   * MPI-case in order to be consistently saved. It is necessary that the vector
   * is consistently saved because it is the only way to ensure that its
   * multiplication with an inconsistently saved matrix (multiplication which
   * appears in the defect and rhs computations) give the correct results. When
   * we call copy_nonactive in MPI-case, we have to remember the following: it
   * can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active. That's why we have to update the values so that the vector
   * becomes consistent again.
   */
#ifdef _MPI
  for(int i = 0; i < d; ++i)
  {
    double *ui = s.solution.block(i);
    s.velocity_space->get_communicator().consistency_update(ui, 3);
  }
  double *p = s.solution.block(d);
  s.pressure_space->get_communicator().consistency_update(p, 3);
#endif
    
    s.solution_m1 = s.solution;
    s.solution_m2 = s.solution;
  }
  
  //=======================================================================
  // assemble additional terms, which depend on the primal solution (u,p)
  // what follows is basically a wrapper to call Assemble2D or 3D
  this->assemble_adjoint_terms(u,p);
  
  //======================================================================= 
  // MANAGE OLD_RHS AND FREE SLIP
  
  // the call to assembleslip is necessary here, in order to get
  // the correct old_rhs, i.e., zeros on the slip dofs
  //if(TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=1)
  //  this->modify_slip_bc(true, true);
   
  // copy the current right hand side vector to the old_rhs
  this->old_rhs = this->systems.front().rhs;
}


template<int d>
void TimeNavierStokes_Adjoint<d>::assemble_matrices_rhs(const FEVectFunct& u,
                                                        const FEFunction& p)
{

}

template<int d>
void TimeNavierStokes_Adjoint<d>::solve()
{
  this->TimeNavierStokes<d>::solve();
  
  //System_per_grid& s = this->TimeNavierStokes<d>.systems.front();
  std::vector<std::shared_ptr<FEMatrix>> blocks = this->TimeNavierStokes<d>::systems.front().matrix.get_blocks_uniquely(
    {{0,0},{1,1}});
  for(auto mat : blocks)
  {
    mat->scale_non_active_diagonals(1./diagonal_scaling);
  }
}

template<int d>
void TimeNavierStokes_Adjoint<d>::assemble_adjoint_terms(const FEVectFunct& u,
                                                        const FEFunction& p)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
    
  const FESpace * v_space = this->TimeNavierStokes<d>::systems.front().velocity_space.get();
  const FESpace * p_space = this->TimeNavierStokes<d>::systems.front().pressure_space.get();
#ifdef __2D__  
  if(u.GetFESpace2D() != v_space || p.GetFESpace2D() != p_space)
#else
  if(u.GetFESpace3D() != v_space || p.GetFESpace3D() != p_space)
#endif
  {
    ErrThrow("primal and adjoint solutions should be defined on the same FE "
             "Space");
  }
    
  //=======================================================================
  // assemble additional terms, which depend on the primal solution (u,p)
  // what follows is basically a wrapper to call Assemble2D or 3D
  this->systems.front().rhs.reset();
  auto n_fe_spaces = 1;
  const FESpace* fe_spaces[1]{v_space};
  
  std::vector<std::shared_ptr<FEMatrix>> blocks = this->TimeNavierStokes<d>::systems.front().matrix.get_blocks_uniquely(
#ifdef __2D__
    {{0,0},{0,1},{1,0},{1,1}})
#else //__3D__
    {{0,0},{0,1},{0,2},{1,0},{1,1},{1,2},{2,0},{2,1},{2,2}})
#endif
  ;
  auto n_sq_mat = d*d;
  SquareMatrixD* sq_mat[n_sq_mat];
  for(int i = 0, j = 0; i < d*d; ++i, ++j)
  {
    if(i%d == 0 && i > 0)
      j++;
    sq_mat[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
  }
  
  auto n_rect_mat = 0;
  MatrixD** rect_mat = nullptr;
  
  auto n_rhs = d; // the velocity components on the right-hand side
  double *rhs[d] = {this->TimeNavierStokes<d>::systems.front().rhs.block(0), 
                    this->TimeNavierStokes<d>::systems.front().rhs.block(1)
#ifdef __3D__
                   ,this->TimeNavierStokes<d>::systems.front().rhs.block(2)
#endif
  };
  const FESpace *fe_rhs[d+1] = {v_space, v_space
#ifdef __3D__
    ,v_space
#endif
  };
  
  BoundaryConditionFunction * boundary_conditions[d] = {
    v_space->get_boundary_condition(), v_space->get_boundary_condition() 
#ifdef __3D__
    ,v_space->get_boundary_condition()
#endif
  };
  BoundaryValuesFunction* non_const_bound_values[d] = 
    { BoundaryValueHomogenous, BoundaryValueHomogenous 
#ifdef __3D__
    ,BoundaryValueHomogenous
#endif
    };

  //======================================================================= 
  // set up custom LocalAssembling object and assemble
  int n_terms = d+1;
  int N_FEValues = d*d+d;
  
#ifdef __2D__
  std::vector<MultiIndex2D> derivatives{D00, D10, D01};
  std::vector<int> FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure  
  std::vector<int> row_space = {0, 0, 0, 0};
  std::vector<int> column_space = {0, 0, 0, 0};
  std::vector<int> rhs_space = {0, 0};
  FEFunction *FEFunctions[d] = { u.GetComponent(0), u.GetComponent(1)};
  std::vector<int> FEValue_FctIndex{0, 1, 0, 1, 0, 1};
  std::vector<MultiIndex2D> FEValue_MultiIndex{D00, D00, D10, D10, D01, D01};
#else
  std::vector<MultiIndex3D> derivatives{D000, D100, D010, D001};
  std::vector<int> FESpaceNumber = { 0, 0, 0, 0 }; // 0: velocity, 1: pressure  
  std::vector<int> row_space = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<int> column_space = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<int> rhs_space = {0, 0, 0};
  FEFunction *FEFunctions[d] = { u.GetComponent(0), u.GetComponent(1), u.GetComponent(2)};
  std::vector<int> FEValue_FctIndex{0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
  std::vector<MultiIndex3D> FEValue_MultiIndex{D000, D000, D000, D100, D100, D100,
                                               D010, D010, D010, D001, D001, D001};
#endif
  auto coeffs = TimeNavierStokes<d>::example.get_coeffs();
  AssembleFctParam local_assembling_routine = tnse_adjoint::adjoint_assembling<d>;
  int n_matrices = d*d;
  int N_ParamFct = 1;
  std::vector<ParamFct*> ParameterFct{tnse_adjoint::params_function};
  std::vector<int> BeginParameter{0};
  int N_Parameters = 2*d+d*d; //u0,u1,u2 + x,y,z + derivatives u0,u1,u2
  
  LocalAssembling<d> la(n_terms, derivatives, FESpaceNumber, row_space, 
                       column_space, rhs_space, coeffs, 
                       local_assembling_routine, nullptr, 
                       n_matrices, n_rhs, N_ParamFct, ParameterFct,
                       BeginParameter, N_Parameters, FEFunctions, N_FEValues, 
                       FEValue_FctIndex, FEValue_MultiIndex);
  
  //do the actual assembling
#ifdef __2D__ 
  Assemble2D(
#else //__3D__
  Assemble3D(
#endif
    n_fe_spaces, fe_spaces, n_sq_mat, sq_mat, n_rect_mat, rect_mat,
             n_rhs, rhs, fe_rhs, boundary_conditions,
             non_const_bound_values, la);
  
  for (int i=0;i<d;i++)
    delete FEFunctions[i];
  
  // now make sure the Dirichlet rows are correct. Because we still need the 
  // entries in these rows after the solving step to return the Neumann data
  // on the control boundary, we can not simple set these rows to the identity.
  // Instead we multiply the diagonal by some vary large number, effectively
  // setting these dofs to zero.
#ifdef __2D__ 
  for(auto mat : {blocks[0], blocks[3]})
#else //__3D__
  for(auto mat : {blocks[0], blocks[4], blocks[8]})
#endif
  {
    mat->scale_non_active_diagonals(diagonal_scaling);
  }     
}


template <int d>
void tnse_adjoint::adjoint_assembling(double Mult, double *coeff, double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts, 
                        double ***LocMatrices, double **LocRhs)
{
  const int N_U = N_BaseFuncts[0];
  const double *U = OrigValues[0];
  const double *U_x = OrigValues[1];
  const double *U_y = OrigValues[2];
  const double *U_z = d == 2 ? nullptr : OrigValues[3];

  // solution of primal problem
  const double u1 = param[0];
  const double u2 = param[1];
  const double u3 = d == 2 ? 0. : param[2];
  const double u1x = param[d];
  const double u2x = param[d+1];
  const double u3x = d == 2 ? 0. : param[d+2];
  const double u1y = param[2*d];
  const double u2y = param[2*d+1];
  const double u3y = d == 2 ? 0. : param[2*d+2];
  const double u1z = d == 2 ? 0. : param[3*d];
  const double u2z = d == 2 ? 0. : param[3*d+1];
  const double u3z = d == 2 ? 0. : param[3*d+2];
  const double x = param[d+d*d];
  //const double y = param[d+d*d+1];
  //const double z = d == 2 ? 0. : param[d+d*d+2];
  double curl_u_1 = (u3y - u2z);
  double curl_u_2 = (u1z - u3x);
  double curl_u_3 = (u2x - u1y);
//   double norm_curl_squared = curl_u_3*curl_u_3;
//   if(d == 3)
//     norm_curl_squared += curl_u_1*curl_u_1 + curl_u_2*curl_u_2;
  bool restricted_curl = restricted_curl_functional && (x >= 8. || x <= 4.);
  double factor = restricted_curl ? 0. : tnse_adjoint::cost_functional_weights[0] * Mult;
//  factor *= (m-norm_curl_squared);
  curl_u_1 *= factor;
  curl_u_2 *= factor;
  curl_u_3 *= factor;
  
  
  const double min_u1_0 = (u1 < 0. ? u1 : 0.) * tnse_adjoint::cost_functional_weights[1] * Mult;
  const double max_u2_0 = (u2 > 0. ? u2 : 0.) * tnse_adjoint::cost_functional_weights[1] * Mult;
  
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA12 = LocMatrices[1];
  double ** MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double ** MatrixA21 = LocMatrices[d];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA23 = d == 2 ? nullptr : LocMatrices[d+2];
  double ** MatrixA31 = d == 2 ? nullptr : LocMatrices[2*d];
  double ** MatrixA32 = d == 2 ? nullptr : LocMatrices[2*d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[2*d+2];
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs3 = d == 2 ? nullptr : LocRhs[2];
  
  for(int i = 0; i < N_U; i++)
  {
    const double test = U[i];
    const double test_x = U_x[i];
    const double test_y = U_y[i];
    const double test_z = d == 2 ? 0. : U_z[i];
    
    Rhs1[i] -= d == 2 ? curl_u_3 * (-test_y) : curl_u_2*test_z - curl_u_3*test_y;
    Rhs2[i] -= d == 2 ? curl_u_3 * (-test_x) : curl_u_3*test_x - curl_u_1*test_z;
    if(d == 3)
      Rhs3[i] -= curl_u_1*test_y - curl_u_2*test_x;
    
    
    Rhs1[i] -= min_u1_0 * test;
    Rhs2[i] -= max_u2_0 * test;
    
    for(int j = 0; j < N_U; j++)
    {
      const double ansatz = U[j];
      double val;
      val = test * u1x * ansatz;
      val += (u1 * test_x + u2 * test_y + u3 * test_z) * ansatz;
      MatrixA11[i][j] += Mult * val;

      val = test * u2x * ansatz;
      MatrixA12[i][j] += Mult * val;
      
      if(d == 3)
      {
        val = test * u3x * ansatz;
        MatrixA13[i][j] += Mult * val;
      }
      
      val = test * u1y * ansatz;
      MatrixA21[i][j] += Mult * val;
      
      val = test * u2y * ansatz;
      val += (u1 * test_x + u2 * test_y + u3 * test_z) * ansatz;
      MatrixA22[i][j] += Mult * val;
      
      if(d == 3)
      {
        val = test * u3y * ansatz;
        MatrixA23[i][j] += Mult * val;
        
        val = test * u1z * ansatz;
        MatrixA31[i][j] += Mult * val;
        
        val = test * u2z * ansatz;
        MatrixA32[i][j] += Mult * val;
        
        val = test * u3z * ansatz;
        val += (u1 * test_x + u2 * test_y + u3 * test_z) * ansatz;
        MatrixA33[i][j] += Mult * val;
      }
    }                            // endfor j
  }                              // endfor i
}


void tnse_adjoint::params_function(const double *in, double *out)
{
#ifdef __2D__
  out[0] = in[2]; // u1old
  out[1] = in[3]; // u2old
  out[2] = in[4]; // D10(u1old)
  out[3] = in[5]; // D10(u2old)
  out[4] = in[6]; // D01(u1old)
  out[5] = in[7]; // D01(u2old)
  out[6] = in[0]; // x
  out[7] = in[1]; // y
#else // __3D__
  out[0] = in[3];   // u1old
  out[1] = in[4];   // u2old
  out[2] = in[5];   // u3old
  out[3] = in[6];   // D100(u1old)
  out[4] = in[7];   // D100(u2old)
  out[5] = in[8];   // D100(u3old)
  out[6] = in[9];   // D010(u1old)
  out[7] = in[10];  // D010(u2old)
  out[8] = in[11];  // D010(u3old)
  out[9] = in[12];  // D001(u1old)
  out[10] = in[13]; // D001(u2old)
  out[11] = in[14]; // D001(u3old)
  out[12] = in[0];  // x
  out[13] = in[1];  // y
  out[14] = in[2];  // y
#endif
}


#ifdef __3D__
template class TimeNavierStokes_Adjoint<3>;
#else
template class TimeNavierStokes_Adjoint<2>;
#endif
