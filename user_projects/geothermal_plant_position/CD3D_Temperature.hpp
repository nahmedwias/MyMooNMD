#ifndef CD3D_TEMPERATURE_H
#define CD3D_TEMPERATURE_H
#include "Time_CD3D.h"
#include "FEVectFunct3D.h"
#include "Domain.h"

class CD3D_Temperature : public Time_CD3D
{
  public: 
  
#ifdef _MPI
CD3D_Temperature( TDomain& domain,
                 const ParameterDatabase& param_db,
                 const Example_TimeCD3D& example, 
                 int maxSubDomainPerDof);
#else
CD3D_Temperature( TDomain& domain,
                 const ParameterDatabase& param_db,
                 const Example_TimeCD3D& example);
#endif


  
#ifdef _MPI
  CD3D_Temperature(TDomain &domain, const ParameterDatabase& param_db, int maxSubDomainPerDof);
#else
  CD3D_Temperature(TDomain &domain, const ParameterDatabase& param_db);
#endif
  
  
  void assemble(const TFEVectFunct3D& convection, const double * x, double nu);
  
  void reset_for_output();
};


#endif // CD3D_TEMPERATURE_H
