#include "PresRobustNavierStokes.h"
#include "../include/General/templateNames.h"
#include "../include/AssembleRoutines/Assemble2D.h"
#include "../include/Matrix/SquareMatrix2D.h"
#include "vector"




template <int d>
PresRobustNavierStokes<d>::PresRobustNavierStokes(const TDomain& domain,
    const ParameterDatabase& param_db, const Example_NSPR_NSE2D& ex, 
    std::tuple<int,int,int> order) 
: TimeNavierStokes<d>(domain, param_db),
  projecti_space(this->get_velocity_space().GetCollection(), "q", "Navier--Stokes pressure", ex.get_bc(d),std::get<2>(order)),
  pr_mat({&this->get_velocity_space()},{&projecti_space, &projecti_space, &this->get_pressure_space()}),  
  rhs_xh(projecti_space.GetN_DegreesOfFreedom()),
  db(param_db), example(ex)
{
  db.merge(param_db);
  
  auto s = this->TimeNavierStokes<d>::systems.front();
  const FESpace *vsp = s.velocity_space.get();
  const FESpace *psp = s.pressure_space.get();
  
  pr_mat = BlockFEMatrixPr::Projection_NSE2D(*vsp, projecti_space,*psp);
  pr_mat.print_coloring_pattern("p");
  
  mass_reconst_=BlockFEMatrix::Mass_NSE2D_Type4(*vsp, *psp);
  
  //
  mass_reconst_.print_coloring_pattern("mass");
  Output::dash("dof of RT/BDM : ",  setw(10), this->projecti_space.GetN_DegreesOfFreedom());
  Output::dash("active dof    : ",  setw(10), this->projecti_space.GetActiveBound());
  
}
/* ************************************************************************* */
template <int d>
void PresRobustNavierStokes<d>::assemble_matrices()
{
  // assemble matrices from Navier-Stokes part
  this->assemble_initial_time();
  
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
    
  int n_fe_spaces = 2; // spaces used for assembling matrices
  const FESpace *velo_space = &this->get_velocity_space();
  const FESpace *proj_space = &this->projecti_space;
  const FESpace *pres_space = &this->get_pressure_space();
  
  const FESpace *sp_array[2] = {velo_space, proj_space};
  
  // square matrices
  int n_sq_mat_ass = 3;
  int n_re_mat_ass = 4;
  int n_rhs_ass = 1;
  
  // for the assembling of matrices or right hand side, 
  // spaces are passed as array: this corresponds to the array "sp_array"
  std::vector<int> rspace={0, 0, 1, 1, 1, 0, 0};
  std::vector<int> cspace={0, 0, 1, 0, 0, 1, 1};
  std::vector<int> rspace_rhs={1};
  
  auto& s = this->TimeNavierStokes<d>::systems.front();
  //prepare everything for storing the matrices and rhs
  //this will be used latter inside the class
  const int n_sq_mat_stored=8;
  std::vector<SquareMatrixD*> sq_matrices_stored(n_sq_mat_stored, nullptr);
  auto mblocks = mass_reconst_.get_blocks_uniquely();
  sq_matrices_stored[0] = reinterpret_cast<SquareMatrixD*>(mblocks.at(0).get());
  sq_matrices_stored[1] = reinterpret_cast<SquareMatrixD*>(mblocks.at(1).get());
  sq_matrices_stored[2] = reinterpret_cast<SquareMatrixD*>(mblocks.at(2).get());
  sq_matrices_stored[3] = reinterpret_cast<SquareMatrixD*>(mblocks.at(3).get());
  
  auto blocks = s.matrix.get_blocks_uniquely();
  sq_matrices_stored[4] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
  sq_matrices_stored[5] = reinterpret_cast<SquareMatrixD*>(blocks.at(1).get());
  sq_matrices_stored[6] = reinterpret_cast<SquareMatrixD*>(blocks.at(3).get());
  sq_matrices_stored[7] = reinterpret_cast<SquareMatrixD*>(blocks.at(4).get());
  
  const int n_re_mat_stored=0;
  int n_rhs_stored=2;
  double *rhs_stored[3] = {s.rhs.block(0), s.rhs.block(1), s.rhs.block(2)};
  const FESpace *sp_array_rhs_stored[3] = {velo_space, velo_space, pres_space};
  
  std::array<BoundaryConditionFunction*, d+1> bdcond;
  std::array<BoundaryValuesFunction*, d+1> bdval;
  for(int i = 0; i < d; ++i)
    bdcond[i] = s.velocity_space->get_boundary_condition();
  bdcond[d] = s.pressure_space->get_boundary_condition();
  for(int i = 0; i < d+1; ++i)
    bdval[i] = example.get_bd(i);
  /*
  Assemble2D_MixedFEM(n_fe_spaces, sp_array.data(), n_sq_mat_ass, 
        n_re_mat_ass, rspace, cspace, n_rhs_ass, rspace_rhs, 
        n_sq_mat_stored, nullptr, n_re_mat_stored, nullptr, 
        n_rhs_stored, rhs_stored, sp_array_rhs_stored);*/
  
//   
//   int n_rectangular_matrices = 2*d; // maximum no of rectangular matrices
//   std::vector<MatrixD*> re_matrices(n_rectangular_matrices, nullptr);
//   constexpr int n_rhs = d+1; // maximum number of right hand sides
//   std::array<double*, n_rhs> rhs_array; // right hand side 
//   // finite element function used for nonlinear term
//   std::array<FEFunction*, d+1> feFunction;
//   // boundary conditions and boundary values
//   std::array<BoundaryConditionFunction*, d+1> boundCondition;
//   std::array<BoundaryValuesFunction*, d+1> boundValues;
//   for(int i = 0; i < d+1; ++i)
//     boundValues[i] = example.get_bd()[i];
//   
//   std::array<const FESpace*, d+1> rhs_spaces;
}
/* ************************************************************************* */
#ifdef __3D__
template class PresRobustNavierStokes<3>;
#else
template class PresRobustNavierStokes<2>;
#endif

