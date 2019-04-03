#include "VolumeOfFluid.h"
#include "Database.h"


template<int d> 
VolumeOfFluid<d>::System_per_grid::System_per_grid(const VolumeOfFluid::Example_TimeNSE& example, TCollection& coll, std::tuple<int, int, int> order) 
 : velocity_space(new FESpace(&coll, "u", "u", example.get_bc(0), std::get<0>(order))),
   pressure_space(new FESpace(&coll, "p", "p", example.get_bc(d), std::get<1>(order))),
   temp_space(new FESpace(&coll, "c", "c", example.get_bc(0), std::get<2>(order)))
{
  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 4:
      // matrices for navier stoeks
      break;
    case 14:
      break;
    default:
      ErrThrow("only supported for NSTYPE 4 and 14");
  }
  
  /// convection diffusion matrix 
}


template<int d> 
VolumeOfFluid<d>::VolumeOfFluid(const TDomain& domain, const ParameterDatabase& param_db)
{
}

template<int d> 
VolumeOfFluid<d>::VolumeOfFluid(const TDomain& domain, const ParameterDatabase& param_db, const VolumeOfFluid::Example_TimeNSE& ex)
{
}
