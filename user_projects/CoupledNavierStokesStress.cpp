#include "CoupledNavierStokesStress.h"
#include "Database.h"
#include "LocalAssembleCoupledStressNS.h"
#ifdef __2D__
#include "Upwind.h"
#include "Matrix2D.h"
#include "SquareMatrix2D.h"
#include "Assemble2D.h"
#include "AuxParam2D.h"
#else
#include "Upwind3D.h"
#include "Matrix3D.h"
#include "SquareMatrix3D.h"
#include "Assemble3D.h"
#include "AuxParam3D.h"
#endif
#include "LocalProjection.h"
#include "Hotfixglobal_AssembleNSE.h"
#include "GridTransfer.h"
#include "Multigrid.h"
#ifdef _MPI
#include "ParFECommunicator3D.h"
#endif

template <int d>
ParameterDatabase CoupledNavierStokesStress<d>::default_coupled_database()
{
  Output::print<5>("creating a default TNSE parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default tnse database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TimeNavierStokes parameter database");
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(ParameterDatabase::default_solution_in_out_database());
  db.merge(LocalAssembling<d>::default_local_assembling_database());
  return db;
}

template <int d>
CoupledNavierStokesStress<d>::System_per_grid::System_per_grid(
  const Example_CoupledNS_Stress& example, TCollection& coll, std::tuple<int, int, int> order)
 : stress_space(new FESpace(&coll, "s", "stress space", example.get_bc(0), std::get<0>(order))),
   velocity_space(new FESpace(&coll, "u", "velocity space", example.get_bc(3),
                              std::get<1>(order))),
   pressure_space(new FESpace(&coll, "p", "pressure space", example.get_bc(5),
                              std::get<2>(order)))
{
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
    case 2:
    case 3:
      Output::print("Only implemented for NSTYPE 4 Laplace Type 1 ");
      exit(1);
      break;
    case 4:
      matrix = BlockFEMatrix::Stress_NSE2D_Type4(*velocity_space, *pressure_space, *stress_space);
      break;
#ifdef __2D__
#else
    Output::print("Only implemented for NSTYPE 4 Laplace Type 1 ");
    exit(1);
#endif
  }
  rhs = BlockVector(matrix, true);
  solution = BlockVector(matrix, false);
  stress1 = FEFunction(stress_space.get(), "s1", "s1", this->solution.block(0),
                 solution.length(0));
  stress2 = FEFunction(stress_space.get(), "s2", "s2", this->solution.block(1),
                 solution.length(1));
  stress3 = FEFunction(stress_space.get(), "s3", "s3", this->solution.block(2),
                 solution.length(2));
#ifdef __2D__
  u = FEVectFunct(velocity_space.get(), "u", "u", solution.block(3),
                  solution.length(3), d);
  p = FEFunction(pressure_space.get(), "p", "p", this->solution.block(5),
                 solution.length(5));
#else
  stress4 = FEFunction(stress_space.get(), "s1", "s1", this->solution.block(3),
                 solution.length(3));
  stress5 = FEFunction(stress_space.get(), "s2", "s2", this->solution.block(4),
                 solution.length(4));
  stress6 = FEFunction(stress_space.get(), "s3", "s3", this->solution.block(5),
                 solution.length(5));
  u = FEVectFunct(velocity_space.get(), "u", "u", solution.block(6),
                  solution.length(6), d);
  p = FEFunction(pressure_space.get(), "p", "p", this->solution.block(9),
                 solution.length(9));
#endif
//   for(int i=0; i<36; i++)
//   {
//     matrix.get_blocks().at(i)->info(2);
//     matrix.get_blocks()[i]->GetStructure().draw(std::string("M") + std::to_string(i));
//   }
}

template <int d>
CoupledNavierStokesStress<d>::CoupledNavierStokesStress(const TDomain& domain, 
    const ParameterDatabase& param_db, const Example_CoupledNS_Stress& ex)
: db(default_coupled_database()), systems(), example(ex), solver(param_db),
  outputWriter(param_db)
{
  db.merge(param_db);
  int s = db["stress_order"];
  int v = db["velocity_order"];
  int p = db["pressure_order"];
  auto _spaces_orders_ = std::make_tuple(s, v, p);
  this->get_velocity_pressure_orders(_spaces_orders_);
  
  bool usingMultigrid = this->solver.is_using_multigrid();
  auto collections = domain.get_grid_collections();
  TCollection *coll = collections.front(); // finest grid collection
  // create finite element space, functions, matrices, rhs and solution
  // at the finest grid
  this->systems.emplace_back(example, *coll, _spaces_orders_);
  if(usingMultigrid)
  {
    Output::print(" multigrid is not supported yet!! CoupledNavierStokesStress.cpp ");
    exit(0);
  }
  
  outputWriter.add_fe_function(&this->get_stress_xx());
  outputWriter.add_fe_function(&this->get_stress_xy());
  outputWriter.add_fe_function(&this->get_stress_yy());
  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  
  this->output_problem_size_info();
}

template<int d>
void CoupledNavierStokesStress<d>::check_and_set_parameters()
{
  if(!db["problem_type"].is(6))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 6;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to "
                      "TimeNavierStokes. It is now reset to the correct value "
                      "for TimeNavierStokes (=6).");
      db["problem_type"] = 6;
    }
  }
  
  // set the discretization parameters
  // standard Galerkin
  if(db["space_discretization_type"].is("galerkin"))
    space_disc_global = 1;
}

template<int d>
void CoupledNavierStokesStress<d>::output_problem_size_info() const
{
  int my_rank = 0;
#ifndef _MPI
  auto & stress_space = *this->systems.front().stress_space;
  auto & velocity_space = *this->systems.front().velocity_space;
  auto & pressure_space = *this->systems.front().pressure_space;

  size_t nDofs  = stress_space.GetN_DegreesOfFreedom();
  size_t nDofu  = velocity_space.GetN_DegreesOfFreedom();
  size_t nDofp  = pressure_space.GetN_DegreesOfFreedom();
  size_t nTotal = nDofs + d*nDofu + nDofp;
  size_t nActive= d*velocity_space.GetActiveBound();

  TCollection* coll = velocity_space.GetCollection();

  double hmin, hmax;
  coll->GetHminHmax(&hmin, &hmax);
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  auto velocity_comm = systems.front().velocity_space->get_communicator();
  auto pressure_comm = systems.front().pressure_space->get_communicator();
  int nDofu  = velocity_comm.get_n_global_dof();
  int nDofp  = pressure_comm.get_n_global_dof();
  int nTotal = d*nDofu + nDofp;
#endif
  if(my_rank ==0)
  {
    Output::stat("NavierStokes", "Mesh data and problem size on finest grid");
#ifndef _MPI
    Output::dash("N_Cells            :  ", setw(10), coll->GetN_Cells());
    Output::dash("h(min, max)        :  ", setw(10), hmin, setw(10), " ", hmax);
#endif
    Output::dash("dof stress         :  ", setw(10), nDofs);
    Output::dash("dof velocity       :  ", setw(10), d*nDofu);
#ifndef _MPI
    Output::dash("dof velocity active:  ", setw(10), nActive);
#endif
    Output::dash("dof pressure       :  ", setw(10), nDofp);
    Output::dash("dof total          :  ", setw(10), nTotal);
  }
}

template<int d>
void CoupledNavierStokesStress<d>::get_velocity_pressure_orders(
  std::tuple<int,int,int> &velocity_pressure_orders)
{
  int stress_order = std::get<0>(velocity_pressure_orders);
  int velocity_order = std::get<1>(velocity_pressure_orders);
  int pressure_order = std::get<2>(velocity_pressure_orders);
  int order;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5:
    case 12: case 13: case 14: case 15:
      if(velocity_order > 10)
        order = velocity_order-10;
      else
        order = velocity_order;
      break;
    case -1: case -2: case -3: case -4: case -5: case -101:
      order = velocity_order;
      break;
    case 100: case 201: case 302: case 403: case 504:
      if(d == 3)
        ErrThrow("velocity_order ", velocity_order, " not supported in 3D");
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
  
  switch(pressure_order)
  {
    case -4711:
    {
      switch(velocity_order)
      {
        case -1: case -2: case -3: case -4:
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
        case -11: case -12: case -13: case -14:
          pressure_order = -(velocity_order-1)*10;
          break;
        case 22: case 23: case 24:
          pressure_order = -(velocity_order-11)*10;
          break;
        case 100: case 201: case 302: case 403: case 504:
          pressure_order = -(velocity_order%100 + 10)*10;
          break; 
      }
      break;
    }
    case 1: case 2: case 3: case 4: case 5:
      // pressure order is chosen correctly
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
    case 100: case 201: case 302: case 403: case 504:
      if(d == 3)
        ErrThrow("pressure_order ", pressure_order, " not supported in 3D");
      // pressure order is chosen correctly
      break;
    default:
      ErrThrow("pressure space is not chosen properly ", pressure_order);
  }
  TDatabase::ParamDB->VELOCITY_SPACE = order;
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  std::get<0>(velocity_pressure_orders) = stress_order;
  std::get<1>(velocity_pressure_orders) = velocity_order;
  std::get<2>(velocity_pressure_orders) = pressure_order;

  Output::print("stress   space ", setw(6), std::get<0>(velocity_pressure_orders));
  Output::print("velocity space ", setw(6), std::get<1>(velocity_pressure_orders));
  Output::print("pressure space ", setw(6), std::get<2>(velocity_pressure_orders));
}

template<int d>
void CoupledNavierStokesStress<d>::assemble_linear_terms()
{
  for(System_per_grid & s : systems)
  {
    // assemble the stress matrices
    this->call_assembling_routine_stress(s);
    // assembling all matrices from NavierStokes model
    this->call_assembling_routine_ns(s, LocalAssembling_type::NSE3D_Linear);
  }
  
  
  // s.matrix.get_combined_matrix()->GetStructure().draw("A1");
}


template<int d>
void CoupledNavierStokesStress<d>::call_assembling_routine_stress(CoupledNavierStokesStress<d>::System_per_grid& s)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  // set arrays of spaces for matrices and rhs
  std::vector<const FESpace*> spaces_mat;
  std::vector<const FESpace*> spaces_rhs;
  std::vector<FEFunction*> fefunctions;
  // call to routine to set arrays
  set_arrays_stress(s, spaces_mat, spaces_rhs, fefunctions);
  
  // prepare matrices and rhs for assembling
  std::vector<SquareMatrixD*> sqMatrices;
  std::vector<MatrixD*> rectMatrices;
  std::vector<double*> rhs_array;
  // call the routine to prepare the matrices
  set_matrices_rhs_stress(s, LocalAssembling_type::CoupledStressNS, 
                   sqMatrices, rectMatrices, rhs_array);
  std::array<BoundaryConditionFunction*, 3> boundCondition;
  std::array<BoundaryValuesFunction*, 3> boundValues;
  for(int i = 0; i < 3; ++i)
    boundCondition[i] = s.velocity_space->get_boundary_condition();
  
  for(int i = 0; i < 3; ++i)
    boundValues[i] = example.get_bd(i);
  
  // local assembling settings
  LocalAssembleCoupledStressNS la(db, example.get_coeffs());
  
  Assemble2D(spaces_mat.size(), spaces_mat.data(), sqMatrices.size(),
               sqMatrices.data(), rectMatrices.size(), rectMatrices.data(),
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
               boundCondition.data(), boundValues.data(), la);
}

template<int d>
void CoupledNavierStokesStress<d>::call_assembling_routine_ns(
  CoupledNavierStokesStress<d>::System_per_grid& s, LocalAssembling_type type)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  // set arrays of spaces for matrices and rhs
  std::vector<const FESpace*> spaces_mat;
  std::vector<const FESpace*> spaces_rhs;
  std::vector<FEFunction*> fefunctions;
  
  // call routine to set arratys
  this->set_arrays_ns(s, spaces_mat, spaces_rhs, fefunctions);
  
  
  // prepare matrices and rhs for assembling
  std::vector<SquareMatrixD*> sqMatrices;
  std::vector<MatrixD*> rectMatrices;
  std::vector<double*> rhs_array;
  // call the routine to prepare the matrices
  this->set_matrices_rhs_ns(s, type, sqMatrices, rectMatrices, rhs_array);
  
  std::array<BoundaryConditionFunction*, d+1> boundCondition;
  std::array<BoundaryValuesFunction*, d+1> boundValues;
  for(int i = 0; i < d; ++i)
    boundCondition[i] = s.velocity_space->get_boundary_condition();
  boundCondition[d] = s.pressure_space->get_boundary_condition();
  for(int i = 0; i < d+1; ++i)
    boundValues[i] = example.get_bd(i);
  
  // local assembling settings
  LocalAssembling<d> la(this->db, type, fefunctions.data(),
                        this->example.get_coeffs(), space_disc_global);

  Assemble2D(spaces_mat.size(), spaces_mat.data(), sqMatrices.size(),
               sqMatrices.data(), rectMatrices.size(), rectMatrices.data(),
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
               boundCondition.data(), boundValues.data(), la);
}

template<int d>
void CoupledNavierStokesStress<d>::set_arrays_stress(CoupledNavierStokesStress<d>::System_per_grid& s, 
 std::vector< const FESpace* >& spaces, std::vector< const FESpace* >& spaces_rhs, 
 std::vector< FEFunction* >& functions)
{
  spaces.resize(3);
  spaces[0] = s.stress_space.get();
  spaces[1] = s.velocity_space.get();
  spaces[2] = s.pressure_space.get();
#ifdef __2D__
  spaces_rhs.resize(3); 
  spaces_rhs[0] = s.stress_space.get();
  spaces_rhs[1] = s.stress_space.get();
  spaces_rhs[2] = s.stress_space.get();  
#else
  Output::print("  ");
#endif
  // simplest model where the nonlinearity is in 
  // the convective part of the velocity
  functions.resize(d+1);
  for(int i=0; i<d; ++i)
    functions[i] = s.u.GetComponent(i);
  functions[d] = &s.p;
}

template<int d>
void CoupledNavierStokesStress<d>::set_arrays_ns(CoupledNavierStokesStress<d>::System_per_grid& s, 
 std::vector< const FESpace* >& spaces, std::vector< const FESpace* >& spaces_rhs, 
 std::vector< FEFunction* >& functions)
{
  spaces.resize(2);  
  spaces[0] = s.velocity_space.get();
  spaces[1] = s.pressure_space.get();
  
  spaces_rhs.resize(d+1);
  for(int i = 0; i < d; ++i)
    spaces_rhs[i] = s.velocity_space.get();
  spaces_rhs[d] = s.pressure_space.get();
  
  // simplest model where the nonlinearity is in 
  // the convective part of the velocity
  functions.resize(d+1);
  for(int i=0; i<d; ++i)
    functions[i] = s.u.GetComponent(i);
  functions[d] = &s.p;
}
template <int d>
void CoupledNavierStokesStress<d>::set_matrices_rhs_stress(CoupledNavierStokesStress<d>::System_per_grid& s, 
  LocalAssembling_type type, std::vector< SquareMatrixD* >& sqMat, std::vector< MatrixD* >& reMat, 
  std::vector< double* >& rhs_array)
{
  sqMat.resize(0);
  reMat.resize(0);
  // right hand side: for NSTYPE: 1,2 and 3, size is 2
  rhs_array.resize(3, 0);
  auto blocks = s.matrix.get_blocks_uniquely();

  sqMat.resize(3);
  sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
  sqMat[1] = reinterpret_cast<SquareMatrixD*>(blocks.at(3).get());
  sqMat[2] = reinterpret_cast<SquareMatrixD*>(blocks.at(6).get());
  
  reMat.resize(12);
  reMat[0] = reinterpret_cast<MatrixD*>(blocks.at(1).get());
  reMat[1] = reinterpret_cast<MatrixD*>(blocks.at(2).get());
  reMat[2] = reinterpret_cast<MatrixD*>(blocks.at(4).get());
  reMat[3] = reinterpret_cast<MatrixD*>(blocks.at(5).get());
  reMat[4] = reinterpret_cast<MatrixD*>(blocks.at(7).get());
  reMat[5] = reinterpret_cast<MatrixD*>(blocks.at(8).get());
  
  reMat[6] = reinterpret_cast<MatrixD*>(blocks.at(9).get());
  reMat[7] = reinterpret_cast<MatrixD*>(blocks.at(10).get());
  reMat[8] = reinterpret_cast<MatrixD*>(blocks.at(11).get());
  reMat[9] = reinterpret_cast<MatrixD*>(blocks.at(15).get());
  reMat[10] = reinterpret_cast<MatrixD*>(blocks.at(16).get());
  reMat[11] = reinterpret_cast<MatrixD*>(blocks.at(17).get());
  
  for(auto sm : sqMat)
    sm->reset();
  
  for(auto rm : reMat)
    rm->reset();
  
  rhs_array[0] = s.rhs.block(0);
  rhs_array[1] = s.rhs.block(1);
  rhs_array[2] = s.rhs.block(2);
  
  s.rhs.reset();
}

template<int d>
void CoupledNavierStokesStress<d>::set_matrices_rhs_ns(CoupledNavierStokesStress<d>::System_per_grid& s,
 LocalAssembling_type type, std::vector< SquareMatrixD* >& sqMat, std::vector< MatrixD* >& reMat, 
 std::vector< double* >& rhs_array)
{
  int nstype = TDatabase::ParamDB->NSTYPE;
  if(nstype != 4 )
    ErrThrow("Nstype: ", nstype, " is not supported ");
  
  sqMat.resize(0);
  reMat.resize(0);
  // right hand side: for NSTYPE: 1,2 and 3, size is 2
  rhs_array.resize(3, 0);
  rhs_array[0] = s.rhs.block(3);
  rhs_array[1] = s.rhs.block(4);
  rhs_array[2] = s.rhs.block(5);
  
  auto blocks = s.matrix.get_blocks_uniquely();

  sqMat.resize(5);
  sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(12).get());
  sqMat[1] = reinterpret_cast<SquareMatrixD*>(blocks.at(13).get());
  sqMat[2] = reinterpret_cast<SquareMatrixD*>(blocks.at(18).get());
  sqMat[2] = reinterpret_cast<SquareMatrixD*>(blocks.at(19).get());
  
  reMat.resize(4);
  reMat[0] = reinterpret_cast<MatrixD*>(blocks.at(21).get());
  reMat[1] = reinterpret_cast<MatrixD*>(blocks.at(22).get());
  reMat[2] = reinterpret_cast<MatrixD*>(blocks.at(14).get());
  reMat[3] = reinterpret_cast<MatrixD*>(blocks.at(20).get());
  
  for(auto* mat : sqMat)
  {
    if(mat != nullptr)
      mat->reset();    
  }
  for(auto* mat : reMat)
  {
    if(mat != nullptr)
      mat->reset();    
  }
}

#ifdef __3D__
template class CoupledNavierStokesStress<3>;
#else
template class CoupledNavierStokesStress<2>;
#endif