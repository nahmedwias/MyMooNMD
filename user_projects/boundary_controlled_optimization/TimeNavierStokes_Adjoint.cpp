#include "TimeNavierStokes_Adjoint.hpp"
#ifdef ___2D__
#include "Assemble2D.h"
#else
#include "Assemble3D.h"
#endif
#include "MainUtilities.h" // BoundaryValueHomogenous
#include "Database.h" // to check TDatabase::ParamDB->NSTYPE

namespace tnse_adjoint
{
void zero_solution(double x, double y, 
#ifdef __3D__     
                   double z,
#endif
                   double *values);
void adjoint_assembling(double, double*, double*, double, double**, int*,
                        double***, double**);
void params_function(double *in, double *out);
// very ugly to put this here, but I don't know how to get this into the local
// assembling
std::vector<double> cost_functional_weights;
bool restricted_curl_functional;
}


constexpr double diagonal_scaling = 1.e30;
void tnse_adjoint::zero_solution(double x, double y, 
#ifdef __3D__     
                   double z,
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
TimeNavierStokes_Adjoint<d>::TimeNavierStokes_Adjoint(const TimeNavierStokes<d>& tnse2d,
                             const ParameterDatabase& param_db)
 : TimeNavierStokes<d>(tnse2d) // copy constructor
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
                                       false, true,initial_conditions);
//                                       this->TimeNavierStokes<d>::example.get_nu());
  this->TimeNavierStokes<d>::outputWriter = DataWriter<d>(param_db);
  this->TimeNavierStokes<d>::outputWriter.add_fe_vector_function(&this->get_velocity());
  this->TimeNavierStokes<d>::outputWriter.add_fe_function(&this->get_pressure());
  this->TimeNavierStokes<d>::solver = Solver<BlockFEMatrix, BlockVector>(param_db);
}

template<int d>
void TimeNavierStokes_Adjoint<d>::assemble(const FEVectFunct& u, const FEFunction& p,
                             const FEVectFunct& stokes_u, 
                             std::vector<double> weights,
                             bool restricted_curl)
{
//  if(systems.size() > 1)
//  {
//    ErrThrow("Time_NSE2D_Adjoint::assemble_additional_terms not yet implemented for "
//             "multigrid");
//  }
//  if(TDatabase::ParamDB->NSTYPE != 4)
//  {
//    ErrThrow("The adjoint problem requires a full matrix structure. Please set "
//            "'NSTYPE' to 4");
//  }
//  System_per_grid& s = this->systems.front();
//  const TFESpace2D * v_space = s.velocity_space.get();
//  const TFESpace2D * p_space = s.pressure_space.get();
//  if(u.GetFESpace2D() != v_space || p.GetFESpace2D() != p_space)
//  {
//    ErrThrow("primal and adjoint solutions should be defined on the same FE "
//             "Space");
//  }
//  tnse_adjoint::cost_functional_weights = weights;
//  tnse_adjoint::restricted_curl_functional = restricted_curl;
//  // delete the right-hand side and matrix
//  s.rhs.reset();
//  s.matrix.reset();
//  // assemble all terms also on the Dirichlet boundary entries.
//  TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 1;
//  this->Time_NSE2D::assemble(); // the Stokes parts
//  s.rhs.reset();
//
//  // assemble additional terms, which depend on the primal solution (u,p)
//  // what follows is basically a wrapper to call Assemble2D
//  auto n_fe_spaces = 1;
//  const TFESpace2D* fe_spaces[1]{v_space};
//
//  std::vector<std::shared_ptr<FEMatrix>> blocks = s.matrix.get_blocks_uniquely(
//    {{0,0},{0,1},{1,0},{1,1}});
//  auto n_sq_mat = 4;
//  TSquareMatrix2D* sq_mat[n_sq_mat];
//  sq_mat[0] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(0).get());
//  sq_mat[1] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());
//  sq_mat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(2).get());
//  sq_mat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());
//
//  auto n_rect_mat = 0;
//  TMatrix2D** rect_mat = nullptr;
//
//  auto n_rhs = 2; // the two velocity components on the right-hand side
//  double *rhs[2] = {s.rhs.block(0), s.rhs.block(1)};
//  const TFESpace2D *fe_rhs[3] = {v_space, v_space};
//
//  BoundCondFunct2D * boundary_conditions[2] = {
//    v_space->get_boundary_condtion(), v_space->get_boundary_condition() };
//  BoundValueFunct2D* non_const_bound_values[2] =
//    { BoundaryValueHomogenous, BoundaryValueHomogenous };
//
//  // set up custom LocalAssembling2D object:
//  int n_terms = 3;
//  std::vector<MultiIndex2D> derivatives{D00, D10, D01};
//  std::vector<int> FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//  std::vector<int> row_space = {0, 0, 0, 0};
//  std::vector<int> column_space = {0, 0, 0, 0};
//  std::vector<int> rhs_space = {0, 0};
//  CoeffFct2D coeffs = example.get_coeffs();
//  AssembleFctParam2D* local_assembling_routine = tnse_adjoint::adjoint_assembling;
//  ManipulateFct2D* manipulate_function = nullptr;
//  int n_matrices = 4;
//  int N_ParamFct = 1;
//  std::vector<ParamFct*> ParameterFct{tnse_adjoint::params_function};
//  std::vector<int> BeginParameter{0};
//  int N_Parameters = 10;
//  TFEFunction2D *FEFunctions2D[4] = { u.GetComponent(0), u.GetComponent(1),
//                                      stokes_u.GetComponent(0),
//                                      stokes_u.GetComponent(1) };
//  int N_FEValues = 8;
//  std::vector<int> FEValue_FctIndex{0, 1, 0, 1, 0, 1, 2, 3};
//  std::vector<MultiIndex2D> FEValue_MultiIndex{D00, D00, D10, D10, D01, D01,
//                                               D00, D00};
//  LocalAssembling2D la(n_terms, derivatives, FESpaceNumber, row_space,
//                       column_space, rhs_space, coeffs,
//                       local_assembling_routine, manipulate_function,
//                       n_matrices, n_rhs, N_ParamFct, ParameterFct,
//                       BeginParameter, N_Parameters, FEFunctions2D, N_FEValues,
//                       FEValue_FctIndex, FEValue_MultiIndex);
//
//  //do the actual assembling
//  Assemble2D(n_fe_spaces, fe_spaces, n_sq_mat, sq_mat, n_rect_mat, rect_mat,
//             n_rhs, rhs, fe_rhs, boundary_conditions,
//             non_const_bound_values, la, -1, true);
//
//  delete FEFunctions2D[0];
//  delete FEFunctions2D[1];
//  delete FEFunctions2D[2];
//  delete FEFunctions2D[3];
//
//  // now make sure the Dirichlet rows are correct. Because we still need the
//  // entries in these rows after the solving step to return the Neumann data
//  // on the control boundary, we can not simple set these rows to the identity.
//  // Instead we multiply the diagonal by some vary large number, effectively
//  // setting these dofs to zero.
//  for(auto mat : {blocks[0], blocks[3]})
//  {
//    mat->scale_non_active_diagonals(diagonal_scaling);
//  }
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


void tnse_adjoint::adjoint_assembling(double Mult, double *coeff, double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts, 
                        double ***LocMatrices, double **LocRhs)
{
  const int N_U = N_BaseFuncts[0];
  const double *Orig0 = OrigValues[0];        // u_x
  const double *Orig1 = OrigValues[1];        // u_y
  const double *Orig2 = OrigValues[2];        // u
  
  // solution of primal problem
  const double u1 = param[0];
  const double u2 = param[1];
  const double u1x = param[2];
  const double u2x = param[3];
  const double u1y = param[4];
  const double u2y = param[5];
  const double u1_stokes = param[6];
  const double u2_stokes = param[7];
  const double x = param[8];
  bool restricted_curl = restricted_curl_functional && (x >= 8. || x <= 4.);
  const double curl_u = restricted_curl ? 0 : (u2x - u1y) * tnse_adjoint::cost_functional_weights[0] * Mult;
  const double min_u1_0 = (u1 < 0. ? u1 : 0.) * tnse_adjoint::cost_functional_weights[1] * Mult;
  const double max_u2_0 = (u2 > 0. ? u2 : 0.) * tnse_adjoint::cost_functional_weights[1] * Mult;
  const double u1_diff = (u1 - u1_stokes) * tnse_adjoint::cost_functional_weights[2] * Mult;
  const double u2_diff = (u2 - u2_stokes) * tnse_adjoint::cost_functional_weights[2] * Mult;
  
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA12 = LocMatrices[1];
  double ** MatrixA21 = LocMatrices[2];
  double ** MatrixA22 = LocMatrices[3];
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  
  for(int i = 0; i < N_U; i++)
  {
    const double test00 = Orig0[i];
    const double test10 = Orig1[i];
    const double test01 = Orig2[i];
    
    Rhs1[i] -= curl_u * (-test01);
    Rhs2[i] -= curl_u * test10;
    Rhs1[i] -= min_u1_0 * test00;
    Rhs2[i] -= max_u2_0 * test00;
    Rhs1[i] -= u1_diff * test00;
    Rhs2[i] -= u2_diff * test00;

    for(int j = 0; j < N_U; j++)
    {
      const double ansatz00 = Orig0[j];
      double val;
      val = test00 * u1x * ansatz00;
      val += (u1 * test10 + u2 * test01) * ansatz00;
      MatrixA11[i][j] += Mult * val;

      val = test00 * u2x * ansatz00;
      MatrixA12[i][j] += Mult * val;
      
      val = test00 * u1y * ansatz00;
      MatrixA21[i][j] += Mult * val;
      
      val = test00 * u2y * ansatz00;
      val += (u1 * test10 + u2 * test01) * ansatz00;
      MatrixA22[i][j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


void tnse_adjoint::params_function(double *in, double *out)
{
  out[0] = in[2]; // u1old
  out[1] = in[3]; // u2old
  out[2] = in[4]; // D10(u1old)
  out[3] = in[5]; // D10(u2old)
  out[4] = in[6]; // D01(u1old)
  out[5] = in[7]; // D01(u2old)
  out[6] = in[8]; // u1_stokes
  out[7] = in[9]; // u2_stokes
  out[8] = in[0]; // x
  out[9] = in[1]; // y
}


#ifdef __3D__
template class TimeNavierStokes_Adjoint<3>;
#else
template class TimeNavierStokes_Adjoint<2>;
#endif
