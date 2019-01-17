#include "LocalAssembleCoupledStressNS.h"
#include "Coupled_NS_Stress_Local_routines.h"
#include "Database.h"
#include "NSE_local_assembling_routines.h"
#include <numeric> // std::iota

LocalAssembleCoupledStressNS::LocalAssembleCoupledStressNS(ParameterDatabase db, CoeffFct coeffs)
: LocalAssembling<2>(db)
{
  this->N_Parameters = 0;
  this->N_ParamFct = 0;
  this->ParameterFct = {};
  this->N_FEValues = 0;
  this->FEValue_FctIndex = {};
  this->FEValue_MultiIndex = {};
  this->BeginParameter = {};
  
  N_Rhs = 3;
  RhsSpace = {0, 0, 0};
  N_Terms = 9;
  //                     s,   u,   p,   s_x, s_y, u_x, u_y, p_x, p_y
  Derivatives   = {D00, D00, D00, D10, D01, D10, D01, D10, D01}; //, D20, D02};
  FESpaceNumber = {0,   1,   2,   0,   0,   1,   1,   2,   2};
  this->Needs2ndDerivatives = new bool[3];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->Needs2ndDerivatives[2] = false;
  this->Manipulate = nullptr;
  
  N_Parameters = 0;
  N_ParamFct = 0;
  ParameterFct = {};
  N_FEValues = 0;
  FEValue_FctIndex = {};
  FEValue_MultiIndex = {};
  BeginParameter = {};
  N_Matrices = 15;
  
  RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
  ColumnSpace = { 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0 };
  using namespace std::placeholders;
  double eta = db["eta_non_newt"];
  this->local_assemblings_routines.push_back(
        std::bind(Stress_Stress<2>, _1, _2, _3, _4, _5, _6, _7, _8, eta));
  this->local_assemblings_routines.push_back(Stress_velocity<2>);
  this->local_assemblings_routines.push_back(Velocity_stress<2>);
  this->local_assemblings_routines.push_back(Stress_rhs<2>);
  
  AllOrigValues = new double** [N_Terms];
  OrigValues = new double* [N_Terms];
}
