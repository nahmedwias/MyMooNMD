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
  u = FEVectFunct(velocity_space.get(), "u", "u", solution.block(0),
                  solution.length(0), d);
  p = FEFunction(pressure_space.get(), "p", "p", this->solution.block(d),
                 solution.length(d));
  
}

template <int d>
CoupledNavierStokesStress<d>::CoupledNavierStokesStress(const TDomain& domain, 
    const ParameterDatabase& param_db)
: CoupledNavierStokesStress<d>(domain, param_db, Example_TimeNSE(param_db))
{
  auto t = std::make_tuple(-11, 2, 1);
}

template <int d>
CoupledNavierStokesStress<d>::CoupledNavierStokesStress(const TDomain& domain, 
    const ParameterDatabase& param_db, const Example_TimeNSE& ex)
: db(default_coupled_database()), systems()
{
}

#ifdef __3D__
template class CoupledNavierStokesStress<3>;
#else
template class CoupledNavierStokesStress<2>;
#endif