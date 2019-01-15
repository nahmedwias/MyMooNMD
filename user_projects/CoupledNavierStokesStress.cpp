#include "CoupledNavierStokesStress.h"
#include "Database.h"
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
  const Example_TimeNSE& example, TCollection& coll, std::tuple<int, int, int> order)
 : stress_space(new FESpace(&coll, "s", "stress space", example.get_bc(d), std::get<0>(order))),
   velocity_space(new FESpace(&coll, "u", "velocity space", example.get_bc(0),
                              std::get<1>(order))),
   pressure_space(new FESpace(&coll, "p", "pressure space", example.get_bc(d),
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
  p = FEFunction(pressure_space.get(), "p", "p", this->solution.block(4),
                 solution.length(4));
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
}

template <int d>
CoupledNavierStokesStress<d>::CoupledNavierStokesStress(const TDomain& domain, 
    const ParameterDatabase& param_db)
: CoupledNavierStokesStress<d>(domain, param_db, Example_TimeNSE(param_db))
{
}

template <int d>
CoupledNavierStokesStress<d>::CoupledNavierStokesStress(const TDomain& domain, 
    const ParameterDatabase& param_db, const Example_TimeNSE& ex)
: db(default_coupled_database()), systems()
{
  db.merge(param_db);
  auto t = std::make_tuple(-11, 2, 1);
}

template<int d>
void CoupledNavierStokesStress<d>::check_and_set_parameters()
{

}

template<int d>
void CoupledNavierStokesStress<d>::get_velocity_pressure_orders(
  std::tuple<int,int,int> &velocity_pressure_orders)
{
  int velocity_order = std::get<0>(velocity_pressure_orders);
  int pressure_order = std::get<1>(velocity_pressure_orders);
  int stress_order = std::get<2>(velocity_pressure_orders);
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
  std::get<0>(velocity_pressure_orders) = velocity_order;
  std::get<1>(velocity_pressure_orders) = pressure_order;

  Output::print("velocity space ", setw(6), std::get<0>(velocity_pressure_orders));
  Output::print("pressure space ", setw(6), std::get<1>(velocity_pressure_orders));
  Output::print("stress   space ", setw(6), std::get<2>(velocity_pressure_orders));
}

template<int d>
void CoupledNavierStokesStress<d>::set_arrays(CoupledNavierStokesStress< d >::System_per_grid& s, 
 std::vector< const FESpace* >& spaces, std::vector< const FESpace* >& spaces_rhs, 
 std::vector< FEFunction* >& functions)
{
  spaces.resize(3);
  spaces[0] = s.stress_space.get();
  spaces[1] = s.velocity_space.get();
  spaces[2] = s.pressure_space.get();
#ifdef __2D__
  spaces_rhs.resize(6); // 3 stress, 2 velocity, 1 pressure
  for(int i=0; i<3; i++)
    spaces_rhs[0] = s.stress_space.get();
  
  spaces_rhs[3] = s.velocity_space.get();
  spaces_rhs[4] = s.velocity_space.get();
  
  spaces_rhs[5] = s.pressure_space.get();
#else
#endif
  // simplest model where the nonlinearity is in 
  // the convective part of the velocity
  functions.resize(d);
  for(int i=0; i<d; ++i)
    functions[i] = s.u.GetComponent(i);
  functions[d] = &s.p;
}

#ifdef __3D__
template class CoupledNavierStokesStress<3>;
#else
template class CoupledNavierStokesStress<2>;
#endif