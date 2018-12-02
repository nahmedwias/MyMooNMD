#include <LocalAssembling.h>

#include <Database.h>
#include <MainUtilities.h> 
#include <FEDatabase3D.h>
#include <FEFunction3D.h>
#include <FEFunction2D.h>
#include <MooNMD_Io.h>
#include <string.h>

#include <ConvDiff.h>
#ifdef __3D__
#include <ConvDiff3D.h>
#endif
#ifdef __2D__
#include <ConvDiff2D.h>
#include <FEDatabase2D.h>
#endif
#include "NSE_local_assembling_routines.h"
#include "CD_local_assembling_routines.h"
#include <DarcyMixed.h>

#include <Brinkman3D_Mixed.h>
#include <TNSE3D_FixPo.h>
#include <TNSE3D_ParamRout.h>
#include <TNSE3DSmagorinsky.h>
#include <TNSE3D_ParamRout.h>
#include <TNSE2DGalerkin.h>
#include <Brinkman2D_Mixed.h>

#include <numeric> // std::iota
#include <algorithm>

std::ostream& operator<<(std::ostream& out, const LocalAssembling_type value)
{
  const char* s = 0;
#define PROCESS_VAL(p) case(LocalAssembling_type::p): s = #p; break;
  switch(value)
  {
    PROCESS_VAL(Brinkman3D_Galerkin);
    PROCESS_VAL(ResidualStabPkPk_for_Brinkman3D_Galerkin1);
    PROCESS_VAL(GradDivStab_for_Brinkman3D_Galerkin1);
    PROCESS_VAL(Brinkman2D_Galerkin1);
    PROCESS_VAL(Brinkman2D_Galerkin2);
    PROCESS_VAL(Brinkman2D_Galerkin1ResidualStabP1);
    PROCESS_VAL(Brinkman2D_Galerkin1ResidualStabP2);
    PROCESS_VAL(Brinkman2D_GradDivStabilization);
    PROCESS_VAL(ConvDiff);
    PROCESS_VAL( Darcy );
    PROCESS_VAL(TCDStiffMassRhs);
    PROCESS_VAL(NSE3D_Linear);
    PROCESS_VAL(NSE3D_NonLinear);
    PROCESS_VAL(TNSE3D_LinGAL);
    PROCESS_VAL(TNSE3D_NLGAL);
    PROCESS_VAL(TNSE3D_Rhs);
    PROCESS_VAL(Custom);
    default: s = "unknown LocalAssembling_type type"; break;
  }
#undef PROCESS_VAL
  return out << s;
}

template<int d>
typename Template_names<d>::MultiIndex_vector indices_up_to_order(int order);
template<>
Template_names<2>::MultiIndex_vector indices_up_to_order<2>(int order)
{
  switch(order)
  {
    case 0: return { D00 }; break;
    case 1: return { D00, D10, D01 }; break;
    case 2: return { D00, D10, D01, D20, D11, D02 }; break;
    default:
      ErrThrow("Multi-indices only available for orders 0,1, and 2");
      break;
  }
}
template<>
Template_names<3>::MultiIndex_vector indices_up_to_order<3>(int order)
{ 
  switch(order)
  {
    case 0: return { D000 }; break;
    case 1: return { D000, D100, D010, D001 };
    case 2:
      return { D000, D100, D010, D001, D200, D110, D101, D020, D011, D002 };
      break;
    default:
      ErrThrow("Multi-indices only available for orders 0,1, and 2");
      break;
  }
}

template<int d>
LocalAssembling<d>::LocalAssembling(ParameterDatabase param_db,
                                    LocalAssembling_type t,
                                    FEFunction **fefunctions3d,
                                    CoeffFct coeffs, int disctype)
 : db(default_local_assembling_database()), type(t),
   discretization_type(disctype), Coeffs(coeffs), FEFunctions3D(fefunctions3d)
{
  Output::print<5>("Constructor of LocalAssembling3D: using type ", type);
  db.merge(param_db, false);
  Parameter disc_type{this->db["space_discretization_type"]};

  // the values below only matter if you need an existing finite element
  // function during your assembly. Change them in such a case
  this->N_Parameters = 0;
  this->N_ParamFct = 0;
  this->ParameterFct = {};
  this->N_FEValues = 0;
  this->FEValue_FctIndex = {};
  this->FEValue_MultiIndex = {};
  this->BeginParameter = {};

  this->N_Spaces=0; // is unused for built-in discretization types anyway

  // set all member variables according to type
  switch(this->type)
  {
    ///////////////////////////////////////////////////////////////////////////
    // Brinkman3D: problems and Brinkman problem
    case LocalAssembling_type::Brinkman3D_Galerkin:
    {
      if(d == 2) ErrThrow("Brinkman needs refactoring (unify 2D and 3D)");
      this->N_Terms = 5;
      //this->Derivatives = {D100, D010, D001, D000, D000};
      this->Derivatives = indices_up_to_order<d>(1);
      this->Derivatives.erase(this->Derivatives.begin());
      this->Derivatives.insert(this->Derivatives.end(), 2,
                               indices_up_to_order<d>(0)[0]);
      this->Needs2ndDerivatives = new bool[2];
      this->Needs2ndDerivatives[0] = false;
      this->Needs2ndDerivatives[1] = false;
      this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
      this->N_Rhs = 4;
      this->RhsSpace = { 0, 0, 0, 1 };
      this->Manipulate = nullptr;
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 4:
          this->N_Matrices = 15;
          this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
          this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
          this->local_assemblings_routines.push_back(
            Brinkman3DType2Galerkin);
          break;
        case 14:
          //Matrix Type 14
          this->N_Matrices = 16;
          this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0}; //u: A11,A12,A13,A21,A22,A23,A31,A32,A33,C,B1T,B2T,B3T,B1,B2,B3 (here the lying B-Blocks come first)
          this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1}; //p; (here the standing B-Blocks come first)
          this->local_assemblings_routines.push_back(
            Brinkman3DType1Galerkin);
          break;
        default:
          ErrThrow("Unknown parameter TDatabase::ParamDB->NSTYPE in LocalAssembling3D_type::Brinkman3D_Galerkin case.");
      }
      break;
    }
    case LocalAssembling_type::ResidualStabPkPk_for_Brinkman3D_Galerkin1:
    {
      if(d == 2) ErrThrow("Brinkman needs refactoring (unify 2D and 3D)");
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 14:
        {
          //Matrix Type 14
          this->N_Terms = 11;                                                                             // = #(Derivatives)
          //this->Derivatives = { D100, D010, D001, D000, D000, D100, D010, D001, D200, D020, D002};        // u_x, u_y, u_z, u, p, p_x, p_y, p_z, u_xx, u_yy, u_zz
          auto sot = indices_up_to_order<d>(2);
          this->Derivatives = sot;
          this->Derivatives.erase(this->Derivatives.begin());
          this->Derivatives.insert(this->Derivatives.end(), 2,
                                   indices_up_to_order<d>(0)[0]);
          this->Derivatives.insert(this->Derivatives.end(), 
                                   {sot[1], sot[2], sot[3], sot[4], sot[7],
                                    sot[9]});
          this->Needs2ndDerivatives = new bool[2];                                                        // usually 2nd derivatives are not needed
          this->Needs2ndDerivatives[0] = true;
          this->Needs2ndDerivatives[1] = true;
          this->FESpaceNumber = { 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 };                                      // 0: velocity space, 1: pressure space
          this->N_Matrices = 16;                                                                          // here some stabilization is allowed in the matrix C
          // in the lower right corner
          this->RowSpace =    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};                          //u: A11,A12,A13,A21,A22,A23,A31,A32,A33,C,B1T,B2T,B3T,B1,B2,B3 (here the lying B-Blocks come first)
          this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1};                          //p; (here the standing B-Blocks come first)
          this->N_Rhs = 4;                                                                                // f1, f2, g
          this->RhsSpace = { 0, 0, 0, 1 };                                                                // corresp. to velocity testspace = 0 / pressure = 1
          this->local_assemblings_routines.push_back(
            ResidualStabPkPk_for_Brinkman3DType1Galerkin);
          this->Manipulate = nullptr;
          break;
        }
        default:
          ErrThrow("Unknown parameter TDatabase::ParamDB->NSTYPE in LocalAssembling3D_type::ResidualStabPkPk_for_Brinkman3D_Galerkin1 case.");
      }
      break;

    }
    case LocalAssembling_type::GradDivStab_for_Brinkman3D_Galerkin1:
    {
      if(d == 2) ErrThrow("Brinkman needs refactoring (unify 2D and 3D)");
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 14:
        {
          //Matrix Type 14
          this->N_Terms = 11;                                                                             // = #(Derivatives)
          //this->Derivatives = { D100, D010, D001, D000, D000, D100, D010, D001, D200, D020, D002};        // u_x, u_y, u_z, u, p, p_x, p_y, p_z, u_xx, u_yy, u_zz
          auto sot = indices_up_to_order<d>(2);
          this->Derivatives = sot;
          this->Derivatives.erase(this->Derivatives.begin());
          this->Derivatives.insert(this->Derivatives.end(), 2,
                                   indices_up_to_order<d>(0)[0]);
          this->Derivatives.insert(this->Derivatives.end(), 
                                   {sot[1], sot[2], sot[3], sot[4], sot[7],
                                    sot[9]});
          this->Needs2ndDerivatives = new bool[2];                                                        // usually 2nd derivatives are not needed
          this->Needs2ndDerivatives[0] = true;
          this->Needs2ndDerivatives[1] = true;
          this->FESpaceNumber = { 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 };                                      // 0: velocity space, 1: pressure space
          this->N_Matrices = 16;                                                                          // here some stabilization is allowed in the matrix C
          // in the lower right corner
          this->RowSpace =    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};                          //u: A11,A12,A13,A21,A22,A23,A31,A32,A33,C,B1T,B2T,B3T,B1,B2,B3 (here the lying B-Blocks come first)
          this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1};                          //p; (here the standing B-Blocks come first)
          this->N_Rhs = 4;                                                                                // f1, f2, g
          this->RhsSpace = { 0, 0, 0, 1 };                                                                // corresp. to velocity testspace = 0 / pressure = 1
          this->local_assemblings_routines.push_back(
            GradDivStab_for_Brinkman3DType1Galerkin);
          this->Manipulate = nullptr;
          break;
        }
        default:
          ErrThrow("Unknown parameter TDatabase::ParamDB->NSTYPE in LocalAssembling3D_type::GradDivStab_for_Brinkman3D_Galerkin1 case.");
      }
      break;
    }
    case LocalAssembling_type::Brinkman2D_Galerkin1:
    {
      if(d == 3) ErrThrow("Brinkman needs refactoring (unify 2D and 3D)");
      //Matrix Type 14
      this->N_Terms = 4;
      //this->Derivatives = { D10, D01, D00, D00 };
      this->Derivatives = indices_up_to_order<d>(1);
      this->Derivatives.erase(this->Derivatives.begin());
      this->Derivatives.insert(this->Derivatives.end(), 2,
                              indices_up_to_order<d>(0)[0]);
      this->Needs2ndDerivatives = new bool[2];
      this->Needs2ndDerivatives[1] = false;
      this->Needs2ndDerivatives[0] = false;
      this->FESpaceNumber = { 0, 0, 0, 1 };                               // 0: velocity, 1: pressure
      this->N_Matrices = 9;
      this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0}; // A11, A12, A21, A22, C, B1, B2, B1T, B2T ; corresp. to u
      this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1}; // A11, A12, A21, A22, C, B1, B2, B1T, B2T ; corresp. to p
      this->N_Rhs = 3;
      this->RhsSpace = { 0, 0, 1 };
      this->local_assemblings_routines.push_back(BrinkmanType1Galerkin);
      this->Manipulate = nullptr;
      break;
    }
    case LocalAssembling_type::Brinkman2D_Galerkin2:
    {
      if(d == 3) ErrThrow("Brinkman needs refactoring (unify 2D and 3D)");
      //Matrix Type 14
      this->N_Terms = 6;                                                  // = #(Derivatives)
      //this->Derivatives = { D10, D01, D00, D00, D10, D01};                // u_x, u_y, u, p, p_x, p_y
      auto first_order_indices = indices_up_to_order<d>(1);
      this->Derivatives = {first_order_indices[1], first_order_indices[2]};
      this->Derivatives.insert(this->Derivatives.end(), 2,
                              indices_up_to_order<d>(0)[0]);
      this->Derivatives.insert(this->Derivatives.end(),
                               first_order_indices.begin(),
                               first_order_indices.end());
      this->Needs2ndDerivatives = new bool[2];                            // usually 2nd derivatives are not needed
      this->Needs2ndDerivatives[0] = false;
      this->Needs2ndDerivatives[1] = false;
      this->FESpaceNumber = { 0, 0, 0, 1, 1, 1 };                         // 0: velocity space, 1: pressure space
      this->N_Matrices = 9;                                               // here some stabilization is allowed in the matrix C
      // in the lower right corner
      this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0};
      this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1};
      this->N_Rhs = 3;                                                    // f1, f2, g
      this->RhsSpace = { 0, 0, 1 };                                       // corresp. to velocity testspace = 0 / pressure = 1
      this->local_assemblings_routines.push_back(BrinkmanType2Galerkin);
      this->Manipulate = nullptr;
      break;
    }
    case LocalAssembling_type::Brinkman2D_Galerkin1ResidualStabP1:
    {
      if(d == 3) ErrThrow("Brinkman needs refactoring (unify 2D and 3D)");
      //Matrix Type 14
      this->N_Terms = 6;                                                  // = #(Derivatives)
      //this->Derivatives = { D10, D01, D00, D00, D10, D01};                // u_x, u_y, u, p, p_x, p_y
      auto first_order_indices = indices_up_to_order<d>(1);
      this->Derivatives = {first_order_indices[1], first_order_indices[2]};
      this->Derivatives.insert(this->Derivatives.end(), 1,
                              indices_up_to_order<d>(0)[0]);
      this->Derivatives.insert(this->Derivatives.end(),
                               first_order_indices.begin(),
                               first_order_indices.end());
      this->Needs2ndDerivatives = new bool[2];                            // usually 2nd derivatives are not needed
      this->Needs2ndDerivatives[0] = false;
      this->Needs2ndDerivatives[1] = false;
      this->FESpaceNumber = { 0, 0, 0, 1, 1, 1 };                         // 0: velocity space, 1: pressure space
      this->N_Matrices = 9;                                               // here some stabilization is allowed in the matrix C
      // in the lower right corner
      this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0};
      this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1};
      this->N_Rhs = 3;                                                    // f1, f2, g
      this->RhsSpace = { 0, 0, 1 };                                       // corresp. to velocity testspace = 0 / pressure = 1
      this->local_assemblings_routines.push_back(
        BrinkmanType1GalerkinResidualStabP1);
      this->Manipulate = nullptr;
      break;
    }
    case LocalAssembling_type::Brinkman2D_Galerkin1ResidualStabP2:
    {
      if(d == 3) ErrThrow("Brinkman needs refactoring (unify 2D and 3D)");
      //Matrix Type 14
      this->N_Terms = 8;                                                  // = #(Derivatives)
      //this->Derivatives = { D10, D01, D00, D00, D10, D01, D20, D02};      // u_x, u_y, u, p, p_x, p_y, u_xx, u_yy
      auto first_order_indices = indices_up_to_order<d>(1);
      auto second_order_indices = indices_up_to_order<d>(2);
      this->Derivatives = {first_order_indices[1], first_order_indices[2]};
      this->Derivatives.insert(this->Derivatives.end(), 1,
                              indices_up_to_order<d>(0)[0]);
      this->Derivatives.insert(this->Derivatives.end(),
                               first_order_indices.begin(),
                               first_order_indices.end());
      this->Derivatives.insert(this->Derivatives.end(),
                               {second_order_indices[3],
                                second_order_indices[5]});
      this->Needs2ndDerivatives = new bool[2];                            // usually 2nd derivatives are not needed
      this->Needs2ndDerivatives[0] = true;
      this->Needs2ndDerivatives[1] = true;
      this->FESpaceNumber = { 0, 0, 0, 1, 1, 1, 0, 0 };                   // 0: velocity space, 1: pressure space
      this->N_Matrices = 9;                                               // here some stabilization is allowed in the matrix C
      // in the lower right corner
      this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0};
      this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1};
      this->N_Rhs = 3;                                                    // f1, f2, g
      this->RhsSpace = { 0, 0, 1 };                                       // corresp. to velocity testspace = 0 / pressure = 1
      this->local_assemblings_routines.push_back(
        BrinkmanType1GalerkinResidualStabP2);
      this->Manipulate = nullptr;
      break;
    }
    case LocalAssembling_type::Brinkman2D_GradDivStabilization:
    {
      if(d == 3) ErrThrow("Brinkman needs refactoring (unify 2D and 3D)");
      this->N_Terms = 2;                                                  // = #(Derivatives)
      //this->Derivatives = { D10, D01 };                                   // u_x, u_y, u, p, p_x, p_y, u_xx, u_yy
      this->Derivatives = indices_up_to_order<d>(1);
      this->Derivatives.erase(this->Derivatives.begin());
      this->Needs2ndDerivatives = new bool[2];                            // usually 2nd derivatives are not needed
      this->Needs2ndDerivatives[0] = false;
      this->Needs2ndDerivatives[1] = false;
      this->FESpaceNumber = { 0, 0 };                                     // 0: velocity space, 1: pressure space
      this->N_Matrices = 9;                                               // here some stabilization is allowed in the matrix C
      // in the lower right corner
      this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0};
      this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1};
      this->N_Rhs = 3;                                                    // f1, f2, g
      this->RhsSpace = { 0, 0, 1 };                                       // corresp. to velocity testspace = 0 / pressure = 1
      this->local_assemblings_routines.push_back(BrinkmanGradDivStab);
      this->Manipulate = nullptr;
      break;
    }
    ///////////////////////////////////////////////////////////////////////////
    // stationary convection diffusion problems
    case LocalAssembling_type::ConvDiff:
    {
      this->N_Matrices = 1;
      this->RowSpace = { 0 };
      this->ColumnSpace = { 0 };
      this->N_Rhs = 1;
      this->RhsSpace = { 0 };
      this->N_Terms = d+1;
      //this->Derivatives = { D000, D100, D010, D001 }; // or {D00, D10, D01}
      this->Derivatives = indices_up_to_order<d>(1);
      this->Needs2ndDerivatives = new bool[1];
      this->Needs2ndDerivatives[0] = false;
      this->FESpaceNumber = std::vector<int>(d+1, 0);
      this->Manipulate = nullptr;
      this->local_assemblings_routines.push_back(BilinearAssembleGalerkin<d>);
      if((disc_type.is("supg") || disc_type.is("gls")))
      {
        this->Derivatives = indices_up_to_order<d>(2);
        this->N_Terms = this->Derivatives.size();
        this->Needs2ndDerivatives[0] = true;
        this->FESpaceNumber = std::vector<int>(d+d+1, 0);
        if(disc_type.is("supg"))
          this->local_assemblings_routines.push_back(
                    BilinearAssemble_SD<d>);
        else
          this->local_assemblings_routines.push_back(
                    BilinearAssemble_GLS<d>);
      }
      else if(!disc_type.is("galerkin"))
      {
        ErrThrow("currently the discretization type ", disc_type,
                 " is not supported by the class CD3D");
      }
      break; // break for the type LocalAssembling3D_type::CD3D 
    }
    case LocalAssembling_type::TCDStiffMassRhs:
    case LocalAssembling_type::TCDStiffRhs:
      this->set_parameters_for_tcd(type);
     break;    
    ///////////////////////////////////////////////////////////////////////////
    case LocalAssembling_type::Darcy:
    {
      if(!disc_type.is("galerkin"))
      {
        ErrThrow("Darcy only supports a galerkin discretization currently, ",
                 disc_type);
      }
      // ( A B' )   ( 0 2 )
      // ( B C  )   ( 3 1 )
      this->N_Terms = 2*(d+1);
      //this->Derivatives = {D000, D000, D100, D010, D001, D100, D010, D001};
      auto fot = indices_up_to_order<d>(1);
      this->Derivatives = {fot[0]};
      this->Derivatives.insert(this->Derivatives.end(), fot.begin(), fot.end());
      this->Derivatives.insert(this->Derivatives.end(), fot.begin()+1,
                               fot.end());
      this->Needs2ndDerivatives = new bool[2];
      this->Needs2ndDerivatives[0] = false;
      this->Needs2ndDerivatives[1] = false;
      this->N_Matrices = 4;
      this->RowSpace = { 0, 1, 0, 1 };
      this->ColumnSpace = { 0, 1, 1, 0 };
      this->N_Rhs = 2;
      this->RhsSpace = { 0, 1 };
      this->local_assemblings_routines.push_back(
              BilinearAssembleDarcyGalerkin<d>);
      if(d == 3)
        this->FESpaceNumber = { 0, 1, 0, 0, 0, 1, 1, 1 }; // 0: velocity, 1: pressure
      else
        this->FESpaceNumber = { 0, 1, 0, 0, 1, 1};
      this->Manipulate = nullptr;
      break;
    }
    ///////////////////////////////////////////////////////////////////////////
    // NSE3D: stationary Navier-Stokes problems
    case LocalAssembling_type :: NSE3D_Linear:
    case LocalAssembling_type :: NSE3D_NonLinear:
      this->set_parameters_for_nse(type);
      break;
    ////////////////////////////////////////////////////////////////////////////
    // TNSE3D: nonstationary Navier-Stokes problems
    case LocalAssembling_type::TNSE3D_LinGAL:
    case LocalAssembling_type::TNSE3D_NLGAL:
    case LocalAssembling_type::TNSE3D_Rhs:
      switch(this->discretization_type)
      {
        case GALERKIN:
          this->set_parameters_for_tnse(type);
          break;
        case SMAGORINSKY:
          this->set_parameters_for_tnse_smagorinsky(type);
          break;
        default:
          ErrThrow("DISCTYPE", this->discretization_type , "is not supported yet!!");
      }
      break;
    default:
      ErrThrow("Unknown or unhandled LocalAssembling_type case. ", type);
  }
  
  AllOrigValues = new double** [N_Terms];
  OrigValues = new double* [N_Terms];

  // some consistency checks
  if(Coeffs == nullptr)
  {
    ErrThrow("You need to specify a valid function for the coefficients");
  }
  if(local_assemblings_routines.empty() 
     || local_assemblings_routines[0] == nullptr)
  {
    ErrThrow("A local assembling routine was not set! ", this->type);
  }
}

//========================================================================
template<int d>
LocalAssembling<d>::LocalAssembling(
  int myN_Terms, MultiIndex_vector myDerivatives,
  std::vector<int> myFESpaceNumber, std::vector<int> myRowSpace,
  std::vector<int> myColumnSpace, std::vector<int> myRhsSpace,
  CoeffFct myCoeffs, AssembleFctParam myAssembleParam,
  ManipulateFct* myManipulate, int myN_Matrices, int myN_Rhs,
  int myN_ParamFct, std::vector<ParamFct*> myParameterFct,
  std::vector<int> myBeginParameter, int myN_Parameters,
  FEFunction** myFEFunctions3D, int myN_FEValues,
  std::vector<int> myFEValue_FctIndex,
  MultiIndex_vector myFEValue_MultiIndex,
  int discretization_type_in)
 : db(default_local_assembling_database()), type{LocalAssembling_type::Custom},
   discretization_type{discretization_type_in}, N_Terms(myN_Terms),
   Derivatives(myDerivatives), FESpaceNumber(myFESpaceNumber),
   RowSpace(myRowSpace), ColumnSpace(myColumnSpace), RhsSpace(myRhsSpace),
   Coeffs(myCoeffs), local_assemblings_routines({myAssembleParam}),
   Manipulate(myManipulate),
   N_Matrices(myN_Matrices), N_Rhs(myN_Rhs), N_ParamFct(myN_ParamFct),
   ParameterFct(myParameterFct), BeginParameter(myBeginParameter),
   N_Parameters(myN_Parameters), N_FEValues(myN_FEValues),
   FEFunctions3D(myFEFunctions3D), FEValue_FctIndex(myFEValue_FctIndex),
   FEValue_MultiIndex(myFEValue_MultiIndex)
{
  // Some data members get an extra treatment - 
  // The auxiliary arrays (All)OrigValues are dynamically allocated with size
  // "N_Terms". "N_Spaces" is determined by finding the max in "FESpaceNumber"
  // (+1). "Needs2ndDerivative" is dynamically allocated to the size "N_Spaces"
  // and then filled according to the appearance of "D200", "D020", "D002", 
  // "D110", "D101" or "D011" in "Derivatives".
  
  //Catch some things which might cause trouble.
  if((int)myDerivatives.size() != N_Terms)
  {
    Output::print("Error: myDerivatives.size() != N_Terms.");
  }
  if((int)myFESpaceNumber.size() != N_Terms)
  {
    Output::print("Error: myFESpaceNumber.size() != N_Terms.");
  }
  if((int)myParameterFct.size() != N_ParamFct)
  {
    Output::print("Error: myParameterFct.size() != myN_ParamFct.");
  }
  if((int)myBeginParameter.size() != N_ParamFct)
  {
    Output::print("Error: myBeginParameter.size() != myN_ParamFct.");
  }

  //Inform the world of what's going on.
  Output::print<5>("Constructor of LocalAssembling3D: using type ", type);

  //Dynamically allocate space for auxiliary arrays
  AllOrigValues = new double** [N_Terms];
  OrigValues = new double* [N_Terms];

  //CODE taken from TDiscretForm2D::TDiscreteForm2D(...)
  // find number of spaces
  int max = -1;
  for(int i=0;i<N_Terms;i++)
  {
    int j = FESpaceNumber[i];
    if(j > max) max = j;
  }
  N_Spaces = max+1;

  //Fill the array Needs2ndDerivatives from the vector myNeeds2ndDerivatives
  Needs2ndDerivatives = new bool[N_Spaces];
  for(int i=0;i<N_Spaces;i++){
    Needs2ndDerivatives[i] = false;
  }
  for(int i=0;i<N_Terms;i++)
  {
    auto alpha = Derivatives[i];
    auto sot = indices_up_to_order<d>(2);
    int j = FESpaceNumber[i];
    if(std::find_if(sot.begin()+d+1, sot.end(),
                    [alpha](decltype(alpha) mi){ return alpha == mi;})
      != sot.end())
      Needs2ndDerivatives[j] = true;
  }
}


//========================================================================
template<int d>
LocalAssembling<d>::~LocalAssembling()
{
  delete [] AllOrigValues;
  delete [] OrigValues;
  delete [] Needs2ndDerivatives;
}


//========================================================================
template<int d>
ParameterDatabase LocalAssembling<d>::default_local_assembling_database()
{
  ParameterDatabase db("default local assembling database");
  
  db.add("with_coriolis_force", false,
         "include the coriolis force for (Navier--)Stokes in 3D. This requires "
         "special local assemblings and the (pde-) coefficients of the example "
         "must include the coriolis force vector Omega.");
  
  db.add("laplace_type_deformation", false, 
         "determine the way the laplacian is discretized.");
  
  db.add("nse_nonlinear_form", "convective",
         "Determine how the nonlinear term for Navier--Stokes is assembled. "
          "convective means ( (u.nabla) u, v), "
          "skew_symmetric means (1/2) [((u.nabla) u, v) - ((u.nabla) v, u)], "
          "rotational means ((nabla x u) x u, v), "
          "emac means (D(u)u + div(u)u, v)",
         {"convective", "skew_symmetric", "rotational", "divergence", "emac"});
  
  db.add("space_discretization_type", "galerkin",
         "The type of discretization. Note that not all types are possible for "
         "all problem classes.",
         {"galerkin", "supg", "upwind", "smagorinsky", "cip", "dg", "symm_gls",
          "nonsymm_gls", "pspg", "brezzi_pitkaeranta", "vms_projection",
          "vms_projection_expl", "local_projection", "local_projection_2_level",
          "residual_based_vms"}); 
  
  db.add("pspg_delta0", 0.1, 
         "the stabilization parameter for pspg (Pressure Stabilization Petrov "
         "Galerkin) is delta0 * h^2 /nu, where h is a cell measure (e.g. "
         "diameter), nu is the inverse of the reynolds number, and delta0 is "
         " this parameter. This parameter is also used for Brezzi-Pitkaeranta",
         0., 10.);

  ///@todo add a parameter for the characteristic length L_0 (Brinkman case)  
  db.add("graddiv_stab", 0., 
         "the stabilization parameter for Grad-Div is delta0 (nu + sigma L_0^2), "
	 " where L is a characteristic length (taken equal to 1), nu is the "
	 "inverse of the reynolds number,"
	 " sigma is the inverse of permeability, and delta0 is this parameter.", 0., 10000.);

  db.add("gls_stab", 0.,
         "the stabilization parameter for GLS stabilization is "
	 " delta0 h^2/(nu + sigma L_0^2), "
	 " where L_0 is a characteristic length (taken equal to 1), nu is the "
	 "inverse of the reynolds number,"
	 " sigma is the inverse of permeability, and delta0 is this parameter.", 0., 10000.);

  db.add("L_0", 1.,
         "the parameter for relating Stokes with Darcy terms in Brinkman problem "
   " which is a characteristic length (by default equal to 1).", 0., 10000.);

  db.add("corner_stab", 0.,
         "the stabilization parameter needed for the Darcy limit of Brinkman in order"
         "to restrict normal jumps of the velocity across corners of the domain"
         "(by default equal to 0)", 0., 10000.);
  // the following 'local projection stabilization' (lps) parameters are not 
  // used in this class, but are somewhat similar to other parameters here.
  db.add("lps_delta0", 0.1,
         "The stabilization parameter for local projection stabilization (lps) "
         "is delta0 * h^2 / nu, where h is a cell measure (e.g. diameter), nu "
         "is the inverse of the reynolds number, and delta0 is this parameter. "
         "Sometimes (in time-dependent problems) it is set to be "
         "delta0 * hK / nu.");
  db.add("lps_delta1", 0.1,
         "The stabilization parameter for local projection stabilization (lps) "
         "is delta0 * h^2 / nu, where h is a cell measure (e.g. diameter), nu "
         "is the inverse of the reynolds number, and delta0 is 'lps_delta0'. "
         "Sometimes (in time-dependent problems) it is set to be "
         "delta1 * hK / nu, where delta1 is this parameter.");
  db.add("lps_coeff_type", 0u, "Determine the way the local projection "
         "stabilization (lps) parameter is computed.", 0u, 5u);
  
  db.merge(ParameterDatabase::parmoon_default_database());
  
  return db;
}


//========================================================================
template<int d>
void LocalAssembling<d>::GetLocalForms(int N_Points, const double *weights,
                                       double *AbsDetjk,
                                       std::array<double*, d> coordinates,
                                       int *N_BaseFuncts, BaseFunct *BaseFuncts,
                                       TBaseCell *Cell, int cell_num,
                                       int N_Matrices, int N_Rhs,
                                       double ***LocMatrix, double **LocRhs,
                                       double factor)
{
  this->GetParameters(N_Points, Cell, cell_num, coordinates);
  this->local_coefficients.resize(N_Points);
  int i,j, N_Rows, N_Columns;
  double **CurrentMatrix, *MatrixRow;
  double Mult, *Coeff, *Param;
  const double hK = Cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);
  // the following is only used in the call of Coeffs and Manipulate
  double *parameters[N_Points];
  double *coefficients[N_Points];
  for(int i=0; i<N_Points; i++)
  {
    if(N_ParamFct == 0)
      parameters[i] = nullptr;
    else
      parameters[i] = &this->parameter_functions_values[i*N_Parameters];
    coefficients[i] = local_coefficients[i].data();
  }

  for(i=0; i<N_Matrices; ++i)
  {
    CurrentMatrix = LocMatrix[i];
    N_Rows = N_BaseFuncts[RowSpace[i]];
    N_Columns = N_BaseFuncts[ColumnSpace[i]];
    for(j=0;j<N_Rows;j++)
    {
      MatrixRow = CurrentMatrix[j];
      memset(MatrixRow, 0, SizeOfDouble*N_Columns);
    } // endfor j
  } // endfor i

  for(i=0; i<N_Rhs; ++i)
  {
    N_Rows = N_BaseFuncts[RhsSpace[i]];
    memset(LocRhs[i], 0, SizeOfDouble*N_Rows);
  }

  // *****************************************************
  // for 2Phase flow problems (Sashikumaar Ganesan)
  coefficients[0][0] = Cell->GetPhase_ID();
  coefficients[0][1] = Cell->GetRegionID();
  coefficients[0][2] = hK;
  // *****************************************************

  if(Coeffs)
  {
    double * X = coordinates[0];
    double * Y = coordinates[1];
#ifdef __3D__
    double * Z = coordinates[2];
    Coeffs(N_Points, X, Y, Z, parameters, coefficients);
#else
    Coeffs(N_Points, X, Y, parameters, coefficients);
#endif
  }

  if(Manipulate)
    Manipulate(N_Points, coefficients, parameters, Cell);

  for(i=0; i<N_Terms; ++i)
  {
#ifdef __3D__
    AllOrigValues[i] = 
      TFEDatabase3D::GetOrigElementValues(BaseFuncts[FESpaceNumber[i]], 
                                          Derivatives[i]);
#else
    AllOrigValues[i] =
      TFEDatabase2D::GetOrigElementValues(BaseFuncts[FESpaceNumber[i]],
                                          Derivatives[i]);
#endif
      
  }

  for(i=0; i<N_Points; ++i)
  {
    Mult = weights[i] * AbsDetjk[i] * factor;
    Coeff = local_coefficients[i].data();
    Coeff[19] = AbsDetjk[i];
   
    Param = parameters[i];

    for(j=0; j<N_Terms; j++)
      OrigValues[j] = AllOrigValues[j][i];

    for(auto& lar : local_assemblings_routines)
      lar(Mult, Coeff, Param, hK, OrigValues, N_BaseFuncts, LocMatrix, LocRhs);
  } // end loop over quadrature points 
}
//========================================================================
template<int d>
void LocalAssembling<d>::GetParameters(int n_points,
                                      TBaseCell *cell, int cellnum,
                                      std::array<double*, d> coordinates)
{
  if(N_ParamFct == 0)
    return;
  this->parameter_functions_values.resize(n_points * this->N_Parameters, 0.0);
  
  double *CurrValues, *CurrOrigValues;
  int *CurrIndex;
  int N_BaseFunct[N_FEValues];
  double *Values[N_FEValues];
  double **orig_values[N_FEValues];
  int *Index[N_FEValues];
  double Temp[d + N_FEValues];

   // collect information
  for(int j=0;j<N_FEValues;j++)
  {
    auto fefunction = FEFunctions3D[FEValue_FctIndex[j]];
    Values[j] = fefunction->GetValues();
#ifdef __2D__
    auto fespace = fefunction->GetFESpace2D();
#else
    auto fespace = fefunction->GetFESpace3D();
#endif
    auto& fe = fespace->get_fe(cellnum);
    auto BaseFunct_Id = fe.get_id();
    N_BaseFunct[j] = fe.GetN_DOF();
#ifdef __3D__
    orig_values[j] = TFEDatabase3D::GetOrigElementValues(BaseFunct_Id,
                                                         FEValue_MultiIndex[j]);
#else
    orig_values[j] = TFEDatabase2D::GetOrigElementValues(BaseFunct_Id, 
                                                         FEValue_MultiIndex[j]);
#endif
 
    Index[j] = fespace->GetGlobalDOF(cellnum);
  } // endfor j


  // loop over all quadrature points
  for(int i=0;i<n_points;i++)
  {
    // first three parameters are the coordinates
    Temp[0] = coordinates[0][i]; // x
    Temp[1] = coordinates[1][i]; // y
    if(d == 3)
      Temp[2] = coordinates[2][i]; // z
      
    // loop to calculate all FE values
    for(int k=d,j=0;j<N_FEValues;j++,k++)
    {
      double s = 0;
      int n = N_BaseFunct[j];
      CurrValues = Values[j];
      CurrOrigValues = orig_values[j][i];
      CurrIndex = Index[j];
      for(int l=0;l<n;l++)
        s += CurrValues[CurrIndex[l]]*CurrOrigValues[l];
      Temp[k] = s;
    }  // endfor j
      
    // loop to calculate all parameters
    for(int j=0;j<N_ParamFct;j++)
    {
      //currparam = param + BeginParameter[j];
      double * currparam = &this->parameter_functions_values[
        i*N_Parameters + this->BeginParameter[j]];
      ParameterFct[j](Temp, currparam);
     
    } // endfor j
  } // endfor i
}
//========================================================================
template<int d>
void LocalAssembling<d>::set_parameters_for_tcd(LocalAssembling_type type)
{
  this->N_Matrices = 2;
  this->RowSpace = { 0, 0 };
  this->ColumnSpace = { 0, 0 };
  this->N_Rhs = 1;
  this->RhsSpace = { 0 };
  this->N_Terms = d+1;
  //this->Derivatives = { D000, D100, D010, D001 }; // or {D00, D10, D01}
  this->Derivatives = indices_up_to_order<d>(1);
  this->Needs2ndDerivatives = new bool[1];
  this->Needs2ndDerivatives[0] = false;
  this->FESpaceNumber = std::vector<int>(d+1, 0);
  this->Manipulate = nullptr;
  
  Parameter disc_type{this->db["space_discretization_type"]};
  switch(type)
  {
    case LocalAssembling_type::TCDStiffMassRhs:
      // stiffness matrix and rhs 
      this->local_assemblings_routines.push_back(TCDStiff<d>);
      // mass matrix 
      this->local_assemblings_routines.push_back(TCDMass<d>);
      // rhs 
      this->local_assemblings_routines.push_back(TCDRhs<d>);
      if(disc_type.is("supg"))
      {
        this->Derivatives = indices_up_to_order<d>(2);
        this->N_Terms = this->Derivatives.size();
        this->Needs2ndDerivatives[0] = true;
        this->Needs2ndDerivatives[1] = true;
        this->FESpaceNumber = std::vector<int>(N_Terms, 0);
        // stiffness matrix and rhs 
        this->local_assemblings_routines.push_back(TCDStiffSUPG<d>);
        // mass matrix 
        this->local_assemblings_routines.push_back(TCDMassSUPG<d>);
        // rhs 
        this->local_assemblings_routines.push_back(TCDRhsSUPG<d>);
      }
      else if(!disc_type.is("galerkin"))
      {
        ErrThrow("currently the discretization type ", disc_type,
                 " is not supported by the class Time_CD2D");
      }
      break;
    case LocalAssembling_type::TCDStiffRhs:
      // stiff matrix, rhs 
      // stiffness matrix and rhs 
      this->local_assemblings_routines.push_back(TCDStiff<d>);
      // rhs 
      this->local_assemblings_routines.push_back(TCDRhs<d>);
      
      if(disc_type.is("supg"))
      {
      }
      break;
      default:
      ErrThrow("unknown LocalAssembling_type ", this->type);
      break;
  }
}
//========================================================================
template<int d>
void LocalAssembling<d>::set_parameters_for_nse( LocalAssembling_type type)
{
  //bool with_coriolis = db["with_coriolis_force"];
  //bool laplace_type_deformation = (TDatabase::ParamDB->LAPLACETYPE == 1);
  bool laplace_type_deformation = this->db["laplace_type_deformation"];
  std::string disc_type = this->db["space_discretization_type"];
  Parameter nonlin_form(db["nse_nonlinear_form"]);
  bool galerkin = (disc_type == std::string("galerkin"));
  bool pspg = (disc_type == std::string("pspg"));
  bool symm_gls = (disc_type == std::string("symm_gls"));
  bool nonsymm_gls = (disc_type == std::string("nonsymm_gls"));
  bool brezzi_pitkaeranta = (disc_type == std::string("brezzi_pitkaeranta"));
  bool local_projection = (disc_type == std::string("local_projection"));
  int nstype = TDatabase::ParamDB->NSTYPE;
  int problem_type = db["problem_type"];
  
  if(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==1)
  {
    ErrThrow("Newton method is not supported yet");
  }
  if(TDatabase::ParamDB->LAPLACETYPE == 1)
  {
    if((nstype==1) || nstype==2)
    {
      ErrThrow("LAPLACETYPE is only supported for NSTYPE 3, 4, and 14");
    }
  }
  if(nonlin_form.is("divergence"))
  {
    ErrThrow("nse_nonlinear_form '", nonlin_form, "' is not yet implemented");
  }
  if(nonlin_form.is("rotational") || nonlin_form.is("emac"))
  {
    if((nstype==1) || nstype==2)
    {
      ErrThrow("nse_nonlinear_form ", nonlin_form,
               " is only supported for NSTYPE 3, 4, and 14");
    }
  }
  if(!galerkin && !pspg && !symm_gls && !nonsymm_gls && !brezzi_pitkaeranta
     && !local_projection)
  {
    ErrThrow("unsupported space_discretization_type for NSE", d, "D: ",
             disc_type);
  }
  if((pspg || symm_gls || nonsymm_gls) && nstype != 14)
  {
    ErrThrow("for PSPG, symmetric GLS and non-symmetric GLS stabilization we "
             "need separate B and BT blocks as well as a C block, i.e., "
             "nstype 14");
  }
  // common for all NSTYPE, Discrete forms, etc
  this->N_Rhs = d+1;
  this->RhsSpace = { 0, 0, 0, 1 };
  this->N_Terms = d+2;
  auto foi = indices_up_to_order<d>(1); // first_order_index
  this->Derivatives = indices_up_to_order<d>(0);
  this->Derivatives.insert(this->Derivatives.end(), foi.begin(), foi.end());
  // Derivatives = { D000, D000, D100, D010, D001 } or { D00, D00, D10, D01}
  this->FESpaceNumber = { 0, 1, 0, 0 }; // 0: velocity, 1: pressure
  if(d == 3)
    this->FESpaceNumber.push_back(0);
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->Manipulate = nullptr;
  this->N_Parameters = d;
  this->N_ParamFct = 1;
  this->ParameterFct =  { NSParamsVelocity<d> };
  this->N_FEValues = d;
  this->FEValue_FctIndex = std::vector<int>(d);
  std::iota(this->FEValue_FctIndex.begin(), this->FEValue_FctIndex.end(), 0);
  this->FEValue_MultiIndex = MultiIndex_vector(d, indices_up_to_order<d>(0)[0]);
  this->BeginParameter = { 0 };
  this->N_Matrices = (d+1)*(d+1); // 9, 16
  if(d==3)
  {
    this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
    this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1 };
  }
  else
  {
    this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0 };
    this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1 };
  }
  if(laplace_type_deformation)
  {
    this->local_assemblings_routines.push_back(NSLaplaceDeformation<d>);
  }
  else
  {
    if(nstype == 1 || nstype == 2)
    {
      this->local_assemblings_routines.push_back(NSLaplaceGradGradSingle<d>);
    }
    else
    {
      this->local_assemblings_routines.push_back(NSLaplaceGradGrad<d>);
    }
  }

  // add Darcy (resistance) term for Brinkman problem
  if (problem_type == 7)
  {
    if(nstype == 1 || nstype == 2)
    {
      this->local_assemblings_routines.push_back(NSResistanceMassMatrixSingle<d>);
    }
    else
    {
      this->local_assemblings_routines.push_back(NSResistanceMassMatrix<d>);
    }
  }

  // stabilization
  using namespace std::placeholders;
  double pspg_delta0 = db["pspg_delta0"];
  double gls_stab = db["gls_stab"];
 
  
  if(pspg || symm_gls || nonsymm_gls) // need second derivatives
  {
    this->N_Terms = 2*d+2+d*(d+1)/2;
    auto soi = indices_up_to_order<d>(2);
    this->Derivatives.insert(this->Derivatives.end(), soi.begin()+1, soi.end());
    // Derivatives = { D000, D000, D100, D010, D001, D100, D010, D001, 
    //                 D200, D110, D101, D020, D011, D002 }
    // or            { D00, D00, D10, D01, D10, D01, D20, D11, D02}
    for(int i = 0; i < d; ++i)
      this->FESpaceNumber.push_back(1);
    for(int i = 0; i < d*(d+1)/2; ++i)
      this->FESpaceNumber.push_back(0);
    // FESpaceNumber = {0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0}
    // or              {0, 1, 0, 0, 1, 1, 0, 0, 0}
    this->Needs2ndDerivatives[0] = true;
    if(pspg)
      this->local_assemblings_routines.push_back(
        std::bind(NSPSPG<d>, _1, _2, _3, _4, _5, _6, _7, _8, pspg_delta0));
    else if(symm_gls)
      this->local_assemblings_routines.push_back(
        std::bind(NSsymmGLS<d>, _1, _2, _3, _4, _5, _6, _7, _8, gls_stab));
    else if(nonsymm_gls)
      this->local_assemblings_routines.push_back(
        std::bind(NSnonsymmGLS<d>, _1, _2, _3, _4, _5, _6, _7, _8, gls_stab));
  }
  else if(brezzi_pitkaeranta)
  {
    this->N_Terms += d; // the pressure derivatives
    this->Derivatives.insert(this->Derivatives.end(), foi.begin()+1, foi.end());
    for(int i = 0; i < d; ++i)
      this->FESpaceNumber.push_back(1);
    this->local_assemblings_routines.push_back(
      std::bind(NS_BrezziPitkaeranta<d>, _1, _2, _3, _4, _5, _6, _7, _8,
                pspg_delta0));
  }

  // grad-div stabilization
  double graddiv_stab = db["graddiv_stab"];
  if(std::abs(graddiv_stab) > 1e-10) 
  {
    this->local_assemblings_routines.push_back(
      std::bind(NSGradDiv<d>, _1, _2, _3, _4, _5, _6, _7, _8, graddiv_stab));

    this->local_assemblings_routines.push_back(
        std::bind(NSGradDiv_RightHandSide<d>,
            _1, _2, _3, _4, _5, _6, _7, _8,
            graddiv_stab));
  }

  switch(type)
  {
  case LocalAssembling_type::NSE3D_Linear:
    //this->local_assemblings_routines.push_back(NSDivergenceBlocks<d>);
    if(nonsymm_gls)
    {
      this->local_assemblings_routines.push_back(
          std::bind(NSDivergenceBlocks<d>,
              _1, _2, _3, _4, _5, _6,
              _7, _8, -1));
      this->local_assemblings_routines.push_back(std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6,
							   _7, _8, -1));
    }
    else
    {
      this->local_assemblings_routines.push_back(
          std::bind(NSDivergenceBlocks<d>, _1, _2, _3, _4, _5, _6, _7, _8, 1));
      this->local_assemblings_routines.push_back(std::bind(NSRightHandSide<d>, _1, _2, _3, _4, _5, _6,
							   _7, _8, 1));
    }
    
    if(nstype == 2 || nstype == 4 || nstype == 14)
    {
      this->local_assemblings_routines.push_back(NSGradientBlocks<d>);
    }
    if(pspg)
      this->local_assemblings_routines.push_back(
          std::bind(NSPSPG_RightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8,
              pspg_delta0));
    else if(symm_gls)
      this->local_assemblings_routines.push_back(
          std::bind(NSsymmGLS_RightHandSide<d>, _1, _2, _3, _4, _5, _6, _7, _8,
              gls_stab));
    else if(nonsymm_gls)
      this->local_assemblings_routines.push_back(
          std::bind(NSnonsymmGLS_RightHandSide<d>, _1, _2, _3, _4, _5, _6, _7,
              _8, gls_stab));
    break;
  case LocalAssembling_type::NSE3D_NonLinear:
    {
      if(nonlin_form.is("convective"))
      {
        if(nstype == 1 || nstype == 2)
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_convective_Single<d>);
        }
        else
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_convective<d>);
        }
      }
      else if(nonlin_form.is("skew_symmetric"))
      {
        if(nstype == 1 || nstype == 2)
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_skew_symmetric_Single<d>);
        }
        else
        {
          this->local_assemblings_routines.push_back(
            NSNonlinearTerm_skew_symmetric<d>);
        }
      }
      else if(nonlin_form.is("rotational"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_rotational<d>);
      }
      else if(nonlin_form.is("emac"))
      {
        this->local_assemblings_routines.push_back(
          NSNonlinearTerm_emac<d>);
      }
      else
      {
        ErrThrow("unknown type for nse_nonlinear_form ", nonlin_form);
      }
      break;
    }
    default:
      ErrThrow("unknown LocalAssembling3D_type ", this->type);
      break;
  }
}
//========================================================================
template<int d>
void LocalAssembling<d>::set_parameters_for_tnse( LocalAssembling_type la_type)
{
  unsigned int nstype = TDatabase::ParamDB->NSTYPE;
  unsigned int laplace_type = TDatabase::ParamDB->LAPLACETYPE;
  if(laplace_type == 1 && (nstype==1 || nstype==2))
  {
    ErrThrow("LAPLACETYPE ", laplace_type, " is supported only for "
               "NSTYPE 3 and 4, not for NSTYPE 1 and 2.");
  }
  // same for all nstypes; 
  //NOTE: change according to the discretization schemes used
  // changing needed for turbulent models and for the newton method
  this->N_Parameters = d;
  this->N_ParamFct = 1;
  this->ParameterFct =  { NSParamsVelocity<d> };
  this->N_FEValues = d;
  this->FEValue_FctIndex = std::vector<int>(d);
  std::iota(this->FEValue_FctIndex.begin(), this->FEValue_FctIndex.end(), 0);
  this->FEValue_MultiIndex = MultiIndex_vector(d, indices_up_to_order<d>(0)[0]);
  this->BeginParameter = { 0 };
  
  this->N_Terms = d+2;
  auto foi = indices_up_to_order<d>(1); // first_order_index
  this->Derivatives = indices_up_to_order<d>(0);
  this->Derivatives.insert(this->Derivatives.end(), foi.begin(), foi.end());
  if(d == 3) // this will be removed soon ...
  {
    //this->Derivatives = {D100, D010, D001, D000, D000};
    this->Derivatives.erase(this->Derivatives.begin());
    this->Derivatives.erase(this->Derivatives.begin());
    this->Derivatives.insert(this->Derivatives.end(), 2,
                             indices_up_to_order<d>(0)[0]);
  }
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  if(d == 3)
  {
    this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
  }
  else
  {
    this->FESpaceNumber = { 0, 1, 0, 0 }; // 0: velocity, 1: pressure
  }
  this->Manipulate = nullptr;
  this->N_Rhs = d+1;
  this->RhsSpace = { 0, 0, 0, 1};
  switch(la_type)
  {
    case LocalAssembling_type::TNSE3D_LinGAL:
      // case 0: fixed point, case 1: newton iteration
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) 
      {
        case 0: // fixed point iteration
          switch(nstype)
          {
            case 1:
            {
              this->N_Matrices = 2 + d;
              if(d == 3)
              {
                this->RowSpace    = { 0, 0, 1, 1, 1 };
                this->ColumnSpace = { 0, 0, 0, 0, 0 };
                this->local_assemblings_routines.push_back(TimeNSType1Galerkin3D);
              }
              else
              {
                this->RowSpace      = { 0, 0, 1, 1 };
                this->ColumnSpace   = { 0, 0, 0, 0 };
                this->local_assemblings_routines.push_back(TimeNSType1Galerkin);
              }
            }
              break;
            case 2:
            {
              this->N_Matrices = 2 + 2*d;
              if(d == 3)
              {
                this->RowSpace    = { 0, 0, 1, 1, 1, 0, 0, 0 };
                this->ColumnSpace = { 0, 0, 0, 0, 0, 1, 1, 1 };
                this->local_assemblings_routines.push_back(TimeNSType2Galerkin3D);
              }
              else
              {
                this->RowSpace      = { 0, 0, 1, 1, 0, 0 };
                this->ColumnSpace   = { 0, 0, 0, 0, 1, 1 };
                this->local_assemblings_routines.push_back(TimeNSType2Galerkin);
              }
            }
              break;
            case 3:
              this->N_Matrices = d + d*d + d; // mass matrix, A, B
              if(d == 3)
              {
                this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
                this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                if(laplace_type==0)
                {
                  this->local_assemblings_routines.push_back(
                    TimeNSType3Galerkin3D);
                }
                else
                {
                  this->local_assemblings_routines.push_back(
                    TimeNSType3GalerkinDD3D);
                }
              }
              else
              {
                this->RowSpace      = { 0, 0, 0, 0, 0, 0, 1, 1 };
                this->ColumnSpace   = { 0, 0, 0, 0, 0, 0, 0, 0 };
                if(laplace_type == 0)
                  this->local_assemblings_routines.push_back(TimeNSType3Galerkin);
                else
                  this->local_assemblings_routines.push_back(TimeNSType3GalerkinDD);
              }
              break;
            case 4:
              this->N_Matrices = d + d*d + 2*d; // mass matrix, A, B, BT
              if(d == 3)
              {
                this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
                this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
                if(laplace_type==0)
                {
                  this->local_assemblings_routines.push_back(
                    TimeNSType4Galerkin3D);
                }
                else
                {
                  this->local_assemblings_routines.push_back(
                    TimeNSType4GalerkinDD3D);
                }
              }
              else
              {
                this->RowSpace      = { 0, 0, 0, 0, 0, 0, 1, 1, 0, 0 };
                this->ColumnSpace   = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 };
                if(laplace_type == 0)
                  this->local_assemblings_routines.push_back(TimeNSType4Galerkin);
                else
                  this->local_assemblings_routines.push_back(TimeNSType4GalerkinDD);
              }
              break;
            case 14:
              // I have to do that 
              break;
          }
          break;
        case 1: // newton iteration
          ErrThrow("Newton method is not supported yet");
          break;
        default:
          ErrThrow("SC_NONLIN_ITE_TYPE_SADDLE ",
                   TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE,
                   "not supported");
      }
      break;
    // local assembling of nonlinear term
    case LocalAssembling_type::TNSE3D_NLGAL:
      if(d == 3) // to be removed
      {
        this->N_Terms = 4;
        this->N_Rhs = 0;
        this->RhsSpace = { };
      }
      // case 0: fixed point, case 1: newton iteration
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) 
      {
        case 0: // fixed point iteration
          switch(nstype)
          {
            case 1:
            case 2:
              this->N_Matrices = 1;
              this->RowSpace    = { 0};
              this->ColumnSpace = { 0};              
              if(d == 3)
                this->local_assemblings_routines.push_back(
                  TimeNSType1_2NLGalerkin3D);
              else
                this->local_assemblings_routines.push_back(
                  TimeNSType1_2NLGalerkin);
              break;
            case 3:
            case 4:
              this->N_Matrices = d;
              this->RowSpace    = std::vector<int>(d, 0);
              this->ColumnSpace = std::vector<int>(d, 0);
              if(d == 3)
              {
                if(laplace_type==0)
                {
                  this->local_assemblings_routines.push_back(
                    TimeNSType3_4NLGalerkin3D);
                }
                else
                {
                  this->local_assemblings_routines.push_back(
                    TimeNSType3_4NLGalerkinDD3D);
                }
              }
              else
              {
                if(laplace_type == 0)
                  this->local_assemblings_routines.push_back(TimeNSType3_4NLGalerkin);
                else
                  this->local_assemblings_routines.push_back(TimeNSType3_4NLGalerkinDD);
              }
              break;
            case 14:
              ErrThrow("NSTYPE 14 is not supported yet");
              break;
          }
          break;
        case 1: // newton iteration
          ErrThrow("Newton method is not supported yet");
          break;
        default:
          ErrThrow("SC_NONLIN_ITE_TYPE_SADDLE ",
                   TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE,
                   "not supported");
      }
      break;
    // local assembling of right hand side
    case LocalAssembling_type::TNSE3D_Rhs:
      // case 0: fixed point iteration, case 1: Newton iteration
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE)
      {
        case 0:  // fixed point iteration
          if(d == 3) // to be removed
          {
            this->N_Terms = 1;
            this->Derivatives = indices_up_to_order<d>(0);
            this->FESpaceNumber = { 0 }; // 0: velocity, 1: pressure
          }
          this->N_Matrices = 0;
          this->RowSpace = { };
          this->ColumnSpace = { };
          if(d == 3)
            this->local_assemblings_routines.push_back(TimeNSRHS3D);
          else
            this->local_assemblings_routines.push_back(TimeNSRHS);
          break;
        case 1: // Newton iteration
          ErrThrow("Newton iteration is not supported yet.");
          break;
        default:
          ErrThrow("SC_NONLIN_ITE_TYPE_SADDLE ", TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE,
                   " not supported.");
      }
      break;
    default:
      ErrThrow("Wrong LocalAssembling3D_type for set_parameters_for_tnse.");
  }
}
//========================================================================
template<>
void LocalAssembling<2>::set_parameters_for_tnse_smagorinsky( LocalAssembling_type type)
{
  ErrThrow("smagorinsky is not supported in 2D yet");
}
template<>
void LocalAssembling<3>::set_parameters_for_tnse_smagorinsky( LocalAssembling_type type)
{
  //NOTE: change according to the discretization schemes used
  // changing needed for turbulent models and for the newton method
  this->N_Parameters = 15;
  this->N_ParamFct = 1;
  this->ParameterFct =  { TimeNSParamsVelo_GradVelo3D };
  this->N_FEValues = 12;
  this->FEValue_FctIndex = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2 };
  // u1old, u2old, u3old, all derivatives of u1, u2, u3
  this->FEValue_MultiIndex = { D000, D000, D000, 
                               D100, D100, D100,
                               D010, D010, D010,
                               D001, D001, D001 };
  this->BeginParameter = { 0 };
  

  switch(type)
  {
    case LocalAssembling_type::TNSE3D_LinGAL:
      // case 0: fixed point, case 1: newton iteration
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) 
      {
        // fixed point 
        case 0:
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
              this->N_Terms = 5;
              this->Derivatives = {D100, D010, D001, D000, D000};
              this->Needs2ndDerivatives = new bool[2];
              this->Needs2ndDerivatives[0] = false;
              this->Needs2ndDerivatives[1] = false;
              this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
              this->N_Matrices = 5;
              this->RowSpace    = { 0, 0, 1, 1, 1 };
              this->ColumnSpace = { 0, 0, 0, 0, 0 };
              this->N_Rhs = 4;
              this->RhsSpace = { 0, 0, 0, 0 };
              this->local_assemblings_routines.push_back(
                TimeNSType1Smagorinsky3D);
              this->Manipulate = nullptr;              
              break;
            case 2:
              this->N_Terms = 5;
              this->Derivatives = {D100, D010, D001, D000, D000};
              this->Needs2ndDerivatives = new bool[2];
              this->Needs2ndDerivatives[0] = false;
              this->Needs2ndDerivatives[1] = false;
              this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
              this->N_Matrices = 8;
              this->RowSpace    = { 0, 0, 1, 1, 1, 0, 0, 0 };
              this->ColumnSpace = { 0, 0, 0, 0, 0, 1, 1, 1 };
              this->N_Rhs = 4;
              this->RhsSpace = { 0, 0, 0, 0 };
              this->local_assemblings_routines.push_back(
                TimeNSType2Smagorinsky3D);
              this->Manipulate = nullptr;
              break;
            case 3:
              this->N_Terms = 5;
              this->Derivatives = {D100, D010, D001, D000, D000};
              this->Needs2ndDerivatives = new bool[2];
              this->Needs2ndDerivatives[0] = false;
              this->Needs2ndDerivatives[1] = false;
              this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
              this->N_Matrices = 15;
              this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
              this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
              this->N_Rhs = 4;
              this->RhsSpace = { 0, 0, 0, 0 };

              this->Manipulate = nullptr;
              if(TDatabase::ParamDB->LAPLACETYPE==0)
              {
                this->local_assemblings_routines.push_back(
                  TimeNSType3Smagorinsky3D);
              }
              else
              {
                this->local_assemblings_routines.push_back(
                  TimeNSType3SmagorinskyDD3D);
              }
              ErrThrow("not tested and adjusted yet: ");
              break;
            case 4:
              this->N_Terms = 5;
              this->Derivatives = {D100, D010, D001, D000, D000};
              this->Needs2ndDerivatives = new bool[2];
              this->Needs2ndDerivatives[0] = false;
              this->Needs2ndDerivatives[1] = false;
              this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
              this->N_Matrices = 18;
              this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
              this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
              this->N_Rhs = 4;
              this->RhsSpace = { 0, 0, 0, 1};

              this->Manipulate = nullptr;
              if(TDatabase::ParamDB->LAPLACETYPE==0)
              {
                this->local_assemblings_routines.push_back(
                  TimeNSType4Smagorinsky3D);
              }
              else
              {
                this->local_assemblings_routines.push_back(
                  TimeNSType4SmagorinskyDD3D);
              }
              break;
          }
          break;
          // newton iteration 
        case 1:
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              ErrThrow("Newton iteration is only supported for NSTYE 3, and 4");
              break;
            case 3:
              break;
            case 4:
              break;
          }
          break;
      }// endswitch SC_NONLIN_ITE_TYPE_SADDLE
       break;
    case LocalAssembling_type::TNSE3D_NLGAL:
    {
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) 
      {
        case 0: // fixed point iteration
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              this->N_Terms = 4;
              this->Derivatives = {D100, D010, D001, D000};
              this->Needs2ndDerivatives = new bool[2];
              this->Needs2ndDerivatives[0] = false;
              this->Needs2ndDerivatives[1] = false;
              this->FESpaceNumber = { 0, 0, 0, 0 }; // 0: velocity, 1: pressure
              this->N_Matrices = 1;
              this->RowSpace    = { 0};
              this->ColumnSpace = { 0};
              this->N_Rhs = 0;
              this->RhsSpace = { };
              this->local_assemblings_routines.push_back(
                TimeNSType1_2NLSmagorinsky3D);
              this->Manipulate = nullptr;    
              break;
            case 3:
            case 4:
              this->N_Terms = 4;
              this->Derivatives = {D100, D010, D001, D000};
              this->Needs2ndDerivatives = new bool[2];
              this->Needs2ndDerivatives[0] = false;
              this->Needs2ndDerivatives[1] = false;
              this->FESpaceNumber = { 0, 0, 0, 0 }; // 0: velocity, 1: pressure
              this->N_Matrices = 3;
              this->RowSpace    = { 0, 0, 0};
              this->ColumnSpace = { 0, 0, 0};
              this->N_Rhs = 0;
              this->RhsSpace = { };
              if(TDatabase::ParamDB->LAPLACETYPE==0)
              {
                this->local_assemblings_routines.push_back(
                  TimeNSType3_4NLSmagorinsky3D);
              }
              else
              {
                this->N_Matrices = 9;
                this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
                this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
                this->local_assemblings_routines.push_back(
                  TimeNSType3_4NLSmagorinskyDD3D);
              }
              
              this->Manipulate = nullptr;    
              break;
          }
        break;
        case 1:// Newton iteration
          ErrThrow("Newton method is not yet supported");
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              break;
            case 3:
            case 4:
              break;
          }
        break;
        default:
          ErrThrow("SC_NONLIN_ITE_TYPE_SADDLE: ", TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE,
                   " is not implemented")
      }
    }
    break;
   case LocalAssembling_type::TNSE3D_Rhs:
     // case 0: fixed point iteration, case 1: Newton iteration
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE)
      {
        case 0:  // fixed point iteration
          this->N_Terms = 1;
          this->Derivatives = { D000 };
          this->Needs2ndDerivatives = new bool[2];
          this->Needs2ndDerivatives[0] = false;
          this->Needs2ndDerivatives[1] = false;
          this->FESpaceNumber = { 0 }; // 0: velocity, 1: pressure
          this->N_Matrices = 0;
          this->RowSpace = { };
          this->ColumnSpace = { };
          this->N_Rhs = 4 ; // TODO The case NSTYPE4 has to be implemented
          this->RhsSpace = {0, 0, 0, 0};
          this->local_assemblings_routines.push_back(TimeNSRHS3D);
          this->Manipulate = nullptr;
          break;
        case 1: // Newton iteration
          ErrThrow("Newton iteration is not supported yet.");
          break;
        default:
          ErrThrow("SC_NONLIN_ITE_TYPE_SADDLE ", TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE,
                   " not supported.");
      }
     break;
        default:
          ErrThrow("Unknown LocalAssembling3D_type");
  }
}

#ifdef __3D__
template class LocalAssembling<3>;
#endif
#ifdef __2D__
template class LocalAssembling<2>;
#endif
