#include <Database.h>
#include <MainUtilities.h> 
#include <FEDatabase3D.h>
#include <FEFunction3D.h>
#include <LocalAssembling3D.h>
#include <ConvDiff.h>
#include <ConvDiff3D.h>
#include <NSE3D_FixPo.h>
#include <NSE3D_FixPoSkew.h>
#include <NSE3D_Param.h>
#include <NSE3D_ParamRout.h>


#include <MooNMD_Io.h>
#include <string.h>

#include <DiscreteForm3D.h>

/** @brief a helper function returning a string with for the name of the 
 *         LocalAssembling3D_type. This returns an empty string in case the type
 *         is not known. */

std::string LocalAssembling3D_type_to_string(LocalAssembling3D_type type)
{
  switch(type)
  {
    case LocalAssembling3D_type :: CD3D:
      switch(TDatabase::ParamDB->DISCTYPE)
      {
        case GALERKIN:
          return std::string("CD3D_Galerkin");
          break;
      }
      break;
    default:
      return std::string();
  }
}

LocalAssembling3D::LocalAssembling3D(LocalAssembling3D_type type, 
                                     TFEFunction3D **fefunctions3d,
                                     CoeffFct3D *coeffs)
 : type(type), name(LocalAssembling3D_type_to_string(type)), Coeffs(coeffs),
   FEFunctions3D(fefunctions3d)
{

  Output::print<2>("Constructor of LocalAssembling3D: using type ", name);
  
  this->N_Parameters = 0;
  this->N_ParamFct = 0;
  this->ParameterFct = {};
  this->N_FEValues = 0;
  this->FEValue_FctIndex = {};
  this->BeginParameter = {};
  
  // set all member variables according to type
  switch(this->type)
  {
    case LocalAssembling3D_type::CD3D:
      switch(TDatabase::ParamDB->DISCTYPE)
      {
        case GALERKIN:
          this->N_Terms = 4;
          this->Derivatives = { D100, D010, D001, D000 };
          this->Needs2ndDerivatives = new bool[1];
          this->Needs2ndDerivatives[0] = false;
          this->FESpaceNumber = { 0, 0, 0, 0};
          this->N_Matrices = 1;
          this->RowSpace = { 0 };
          this->ColumnSpace = { 0 };
          this->N_Rhs = 1;
          this->RhsSpace = { 0 };
          this->AssembleParam = BilinearAssembleGalerkin;
          this->Manipulate = NULL;
          break;
        case SUPG:
        	ErrThrow("currently DISCTYPE ", TDatabase::ParamDB->DISCTYPE,
        	               " (SUPG) is not supported by the class CD3D");
          break;
        default:
          ErrThrow("currently DISCTYPE ", TDatabase::ParamDB->DISCTYPE,
               " is not supported by the class CD3D");
      }// endswitch TDatabase::ParamDB->DISCTYPE
      break; // break for the type LocalAssembling3D_type::CD3D 
    case LocalAssembling3D_type :: NSE3D_Linear:
    case LocalAssembling3D_type :: NSE3D_NonLinear:
      this->set_parameters_for_nse(type);
      break;
    default:
      ErrThrow("Unknown or unhandled LocalAssembling3D_type case.");
  }
  
  AllOrigValues = new double** [N_Terms];
  OrigValues = new double* [N_Terms];

  // some consistency checks
  if(Coeffs == NULL)
  {
    ErrThrow("You need to specify a valid function for the coefficients");
  }
  if(AssembleParam == NULL)
  {
    ErrThrow("S local assmebling routine was not set");
  }
}
//========================================================================
LocalAssembling3D::LocalAssembling3D(LocalAssembling3D_type type, TAuxParam3D& aux,
                                     TDiscreteForm3D& df)
  :type(type),
   name(df.getName()), N_Terms(df.getNTerms()), N_Spaces(df.getNSpaces()),
   Needs2ndDerivatives(nullptr), Derivatives(this->N_Terms, D000), 
   FESpaceNumber(this->N_Terms, 0), RowSpace(df.getNMatrices(), 0),
   ColumnSpace(df.getNMatrices(), 0), RhsSpace(df.getNRhs(), 0),
   Coeffs(df.getCoeffs()), AssembleParam(df.getAssembleParam()),
   Manipulate(df.getManipulate()), AllOrigValues(new double** [N_Terms]),
   OrigValues(new double* [N_Terms]), N_Matrices(df.getNMatrices()),
   N_Rhs(df.getNRhs()), N_ParamFct(aux.getNParamFct()),
   ParameterFct(this->N_ParamFct, nullptr), BeginParameter(this->N_ParamFct, 0),
   N_Parameters(aux.getNParameters()), N_FEValues(aux.getNFeValues()),
   FEFunctions3D(aux.getFeFunctions3D()), FEValue_FctIndex(this->N_FEValues,0),
   FEValue_MultiIndex(this->N_FEValues, D000)
{
  this->Needs2ndDerivatives = new bool[this->N_Spaces];
  for(int i = 0; i < this->N_Spaces; ++i)
    this->Needs2ndDerivatives[i] = df.GetNeeds2ndDerivatives()[i];
  
  for(int i = 0; i < this->N_Terms; ++i)
  {
    this->Derivatives.at(i) = df.getDerivative(i);
    this->FESpaceNumber.at(i) = df.getFeSpaceNumber(i);
  }
  
  for(int i = 0; i < this->N_Matrices; ++i)
  {
    this->RowSpace.at(i) = df.rowSpaceOfMat(i);
    this->ColumnSpace.at(i) = df.colSpaceOfMat(i);
  }
  
  for(int i = 0; i < this->N_Rhs; ++i)
    this->RhsSpace.at(i) = df.getRhsSpace(i);
  
  for(int i = 0; i < this->N_ParamFct; ++i)
  {
    this->ParameterFct.at(i) = aux.getParameterFct(i);
    this->BeginParameter.at(i) = aux.getBeginParameter(i);
  }
  
  for(int i = 0; i < this->N_FEValues; ++i)
  {
    this->FEValue_FctIndex.at(i) = aux.getFeValueFctIndex(i);
    this->FEValue_MultiIndex.at(i) = aux.getFeValueMultiIndex(i);
  }
  
  // some consistency checks
  if(Coeffs == NULL)
  {
    ErrThrow("You need to specify a valid function for the coefficients");
  }
  if(this->AssembleParam == NULL)
  {
    // this means in the discrete form there was only a pointer to a
    // AssembleFct3D rather than a AssembleFctParam3D.
    ErrThrow("can't create LocalAssembling3D object, missing AssembleFctParam3D");
  }
}
//========================================================================
LocalAssembling3D::~LocalAssembling3D()
{
  delete [] AllOrigValues;
  delete [] OrigValues;
  delete [] Needs2ndDerivatives;
}
//========================================================================
void LocalAssembling3D::GetLocalForms(int N_Points, double *weights,  double *AbsDetjk,
                       double *X, double *Y, double *Z,
                       int *N_BaseFuncts,
                       BaseFunct3D *BaseFuncts, 
                       double **Parameters, double **AuxArray,
                       TBaseCell *Cell, int N_Matrices,
                       int N_Rhs,
                       double ***LocMatrix, double **LocRhs,
                       double factor) const
{
  int i,j, N_Rows, N_Columns;
  double **CurrentMatrix, *MatrixRow;
  double Mult, *Coeff, *Param;
  const double hK = Cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);


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
  AuxArray[0][0] = Cell->GetPhase_ID();
  AuxArray[0][1] = Cell->GetRegionID();
  AuxArray[0][2] = hK;
  // *****************************************************

  if(Coeffs)
    Coeffs(N_Points, X, Y, Z, Parameters, AuxArray);

  if(Manipulate)
    Manipulate(N_Points, AuxArray, Parameters, Cell);

  for(i=0; i<N_Terms; ++i)
  {
    AllOrigValues[i] = 
      TFEDatabase3D::GetOrigElementValues(BaseFuncts[FESpaceNumber[i]], 
                                          Derivatives[i]);
  }

  for(i=0; i<N_Points; ++i)
  {
    Mult = weights[i] * AbsDetjk[i] * factor;
    Coeff = AuxArray[i];
    Coeff[19] = AbsDetjk[i];
   
    Param = Parameters[i];

    for(j=0; j<N_Terms; j++)
      OrigValues[j] = AllOrigValues[j][i];

    AssembleParam(Mult, Coeff, Param, hK, OrigValues, N_BaseFuncts, LocMatrix,
                  LocRhs);
  } // end loop over quadrature points 
}
//========================================================================
void LocalAssembling3D::GetParameters(int n_points, TCollection *Coll,
                                      TBaseCell *cell, int cellnum,
                                      double *x, double *y, double *z, double **Parameters) const
{
  double *param, *currparam, s;
  const TFESpace3D *fespace;
  TFEFunction3D *fefunction;
  FE3D FE_Id;
  BaseFunct3D BaseFunct_Id;
  int *GlobalNumbers, *BeginIndex;
  
  double *CurrValues, *CurrOrigValues;
  int *CurrIndex;
  int *N_BaseFunct = new int[N_FEValues];
  double **Values = new double* [N_FEValues];
  double ***orig_values = new double** [N_FEValues];
  int **Index = new int* [N_FEValues];
  double Temp[2 + N_FEValues];
  int n;

   // collect information
  for(int j=0;j<N_FEValues;j++)
  {
    fefunction = FEFunctions3D[FEValue_FctIndex[j]];
    Values[j] = fefunction->GetValues();

    fespace = fefunction->GetFESpace3D();
    FE_Id = fespace->GetFE3D(cellnum, cell);
    BaseFunct_Id = TFEDatabase3D::GetFE3D(FE_Id)->GetBaseFunct3D_ID();

    N_BaseFunct[j]=TFEDatabase3D::GetBaseFunct3D(BaseFunct_Id)->GetDimension();
    orig_values[j] = TFEDatabase3D::GetOrigElementValues
                      (BaseFunct_Id, FEValue_MultiIndex[j]);

    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    Index[j] = GlobalNumbers + BeginIndex[cellnum];
  } // endfor j


  // loop over all quadrature points
  for(int i=0;i<n_points;i++)
  {
    // parameters in this quadrature point  
    param = Parameters[i];
    // first three parameters are the coordinates
    Temp[0] = x[i];
    Temp[1] = y[i];
    Temp[2] = z[i];

    // loop to calculate all FE values
    for(int k=3,j=0;j<N_FEValues;j++,k++)
    {
      s = 0;
      n = N_BaseFunct[j];
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
      currparam = param + BeginParameter[j];
      
      /// change for twophase
      currparam[0] = cell->GetPhase_ID();
      
      ParameterFct[j](Temp, currparam);
    } // endfor j
  } // endfor i
  
  delete [] N_BaseFunct;
  delete [] Values;
  delete [] orig_values;
  delete [] Index;
}
//========================================================================
void LocalAssembling3D::set_parameters_for_nse(LocalAssembling3D_type type)
{
  switch(type)
  {
    case LocalAssembling3D_type::NSE3D_Linear:
      switch(TDatabase::ParamDB->DISCTYPE)
      {
        case GALERKIN: // GALERKIN 
          switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) // nonlinear iteration type (fixed point or Newton's type)
          {
            case 0: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0
              switch(TDatabase::ParamDB->NSTYPE)
              {
                case 1: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0, NSTYPE 1
                  if(TDatabase::ParamDB->LAPLACETYPE != 0)
                  {
                    ErrThrow("LAPLACETYPE must be set to 0 in case of NSTYPE 1");
                  }
                  this->N_Terms = 5;
                  this->Derivatives = {D100, D010, D001, D000, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 4;
                  this->RowSpace = { 0, 1, 1, 1 };
                  this->ColumnSpace = { 0, 0, 0, 0 };
                  this->N_Rhs = 3;
                  this->RhsSpace = { 0, 0, 0 };
                   if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
                     this->AssembleParam = NSType1Galerkin3D;
                   else
                     this->AssembleParam = NSType1GalerkinSkew3D;
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
                case 2: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0, NSTYPE 2
                  if(TDatabase::ParamDB->LAPLACETYPE != 0)
                  {
                    ErrThrow("LAPLACETYPE must be set to 0 in case of NSTYPE 2");
                  }
                  this->N_Terms = 5;
                  this->Derivatives = {D100, D010, D001, D000, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 7;
                  this->RowSpace    = { 0, 1, 1, 1, 0, 0, 0 };
                  this->ColumnSpace = { 0, 0, 0, 0, 1, 1, 1 };
                  this->N_Rhs = 3;
                  this->RhsSpace = { 0, 0, 0 };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
                    this->AssembleParam = NSType2Galerkin3D;
                  else
                    this->AssembleParam = NSType2GalerkinSkew3D;
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
                case 3: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0, NSTYPE 3
                  this->N_Terms = 5;
                  this->Derivatives = {D100, D010, D001, D000, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 12;
                  this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
                  this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                  this->N_Rhs = 3;
                  this->RhsSpace = { 0, 0, 0 };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM ==0)
                  {
                    if(TDatabase::ParamDB->LAPLACETYPE == 0)
                      this->AssembleParam = NSType3Galerkin3D;
                    else
                      this->AssembleParam = NSType3GalerkinDD3D;
                  }
                  else // Skew 
                  {
                    if(TDatabase::ParamDB->LAPLACETYPE == 0)
                      this->AssembleParam = NSType3GalerkinSkew3D;
                    else
                      this->AssembleParam = NSType3GalerkinDDSkew3D;
                  }
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
                case 4: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0, NSTYPE 4
                  this->N_Terms = 5;
                  this->Derivatives = {D100, D010, D001, D000, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 15;
                  this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
                  this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
                  this->N_Rhs = 3;
                  this->RhsSpace = { 0, 0, 0 };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM ==0)
                  {
                    if(TDatabase::ParamDB->LAPLACETYPE == 0)
                      this->AssembleParam = NSType4Galerkin3D;
                    else
                      this->AssembleParam = NSType4GalerkinDD3D;
                  }
                  else // Skew 
                  {
                    if(TDatabase::ParamDB->LAPLACETYPE == 0)
                      this->AssembleParam = NSType4GalerkinSkew3D;
                    else
                      this->AssembleParam = NSType4GalerkinDDSkew3D;
                  }
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
              }// endswitch for NSTYPE
              break;
            case 1: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=1
              switch(TDatabase::ParamDB->NSTYPE)
              {
                case 1:// GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=1, NSTYPE=1
                case 2:// GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=1, NSTYPE=2
                  ErrThrow("Wrong NSTYPE ", TDatabase::ParamDB->NSTYPE,
                           " for Newton's method !!! Use NSTYPE 3 or 4 !!!");
                  break;
                case 3: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=3, NSTYPE=3
                  this->N_Terms = 5;
                  this->Derivatives = {D100, D010, D001, D000, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 12;
                  this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
                  this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                  this->N_Rhs = 3;
                  this->RhsSpace = { 0, 0, 0 };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
                  {
                    // if(TDatabase::ParamDB->LAPLACETYPE==0)
                    //   this->AssembleParam = ;
                    // else
                    //  this->AssembleParam = ;
                  }
                  else // Skew 
                  {
                    // if(TDatabase::ParamDB->LAPLACETYPE==0)
                    //   this->AssembleParam = ;
                    // else
                    //   this->AssembleParam = ;
                  }
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  // this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
                case 4: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=3, NSTYPE=4
                  this->N_Terms = 5;
                  this->Derivatives = {D100, D010, D001, D000, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 15;
                  this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
                  this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
                  this->N_Rhs = 3;
                  this->RhsSpace = { 0, 0, 0 };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM ==0)
                  {
                    // if(TDatabase::ParamDB->LAPLACETYPE == 0)
                    //   this->AssembleParam = ;
                    // else
                    //  this->AssembleParam = ;
                  }
                  else // Skew 
                  {
                    // if(TDatabase::ParamDB->LAPLACETYPE == 0)
                    //   this->AssembleParam = ;
                    // else
                    //   this->AssembleParam = ;
                  }
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  // this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
              }// endswitch for NSTYPE
              break;
          }// endswitch of the iteration type (FIXED POINT OR NEWTON'S)
          break;
          default:
            ErrThrow("currently DISCTYPE ", TDatabase::ParamDB->DISCTYPE,
                     " is not supported by the class NSE3D");
      }// endswitch for the DISCTYPE 
      break; // break for the LocalAssembling3D_type NSE3D_Linear
    case LocalAssembling3D_type :: NSE3D_NonLinear:
      switch(TDatabase::ParamDB->DISCTYPE)
      {
        case GALERKIN:
          switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) // nonlinear iteration type (fixed point or Newton's type)
          {
            case 0: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0
              switch(TDatabase::ParamDB->NSTYPE)
              {
                case 1: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0, NSTYPE 1
                  if(TDatabase::ParamDB->LAPLACETYPE != 0)
                  {
                    ErrThrow("LAPLACETYPE must be set to 0 in case of NSTYPE 1");
                  }
                  this->N_Terms = 4;
                  this->Derivatives = {D100, D010, D001, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 1;
                  this->RowSpace = { 0 };
                  this->ColumnSpace = { 0 };
                  this->N_Rhs = 0;
                  this->RhsSpace = {};
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
                    this->AssembleParam = NSType1_2NLGalerkin3D;
                  else
                    this->AssembleParam = NSType1_2NLGalerkinSkew3D;
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
                case 2: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0, NSTYPE 2
                  if(TDatabase::ParamDB->LAPLACETYPE != 0)
                  {
                    ErrThrow("LAPLACETYPE must be set to 0 in case of NSTYPE 1");
                  }
                  this->N_Terms = 4;
                  this->Derivatives = {D100, D010, D001, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 1;
                  this->RowSpace = { 0 };
                  this->ColumnSpace = { 0 };
                  this->N_Rhs = 0;
                  this->RhsSpace = { };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
                    this->AssembleParam = NSType1_2NLGalerkin3D;
                  else
                    this->AssembleParam = NSType1_2NLGalerkinSkew3D;
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
                case 3: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0, NSTYPE 3
                  this->N_Terms = 4;
                  this->Derivatives = {D100, D010, D001, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0}; // 0: velocity, 1: pressure
                  this->N_Matrices = 3;
                  this->RowSpace    = { 0, 0, 0 };
                  this->ColumnSpace = { 0, 0, 0 };
                  this->N_Rhs = 0;
                  this->RhsSpace = { };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM ==0)
                  {
                    if(TDatabase::ParamDB->LAPLACETYPE == 0)
                      this->AssembleParam = NSType3_4NLGalerkin3D;
                    else
                     this->AssembleParam = NSType3_4NLGalerkinDD3D;
                  }
                  else // Skew 
                  {
                    if(TDatabase::ParamDB->LAPLACETYPE == 0)
                      this->AssembleParam = NSType3_4NLGalerkinSkew3D;
                    else
                      this->AssembleParam = NSType3_4NLGalerkinDDSkew3D;
                  }
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
                case 4: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=0, NSTYPE 4
                  this->N_Terms = 4;
                  this->Derivatives = {D100, D010, D001, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 3;
                  this->RowSpace    = { 0, 0, 0};
                  this->ColumnSpace = { 0, 0, 0};
                  this->N_Rhs = 0;
                  this->RhsSpace = { };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM ==0)
                  {
                    if(TDatabase::ParamDB->LAPLACETYPE == 0)
                      this->AssembleParam = NSType3_4NLGalerkin3D;
                    else
                     this->AssembleParam = NSType3_4NLGalerkinDD3D;
                  }
                  else // Skew 
                  {
                    if(TDatabase::ParamDB->LAPLACETYPE == 0)
                      this->AssembleParam = NSType3_4NLGalerkinSkew3D;
                    else
                      this->AssembleParam = NSType3_4NLGalerkinDDSkew3D;
                  }
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
              }// endswitch for NSTYPE
              break;
            case 1: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=1
              switch(TDatabase::ParamDB->NSTYPE)
              {
                case 1:// GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=1, NSTYPE=1
                case 2:// GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=1, NSTYPE=2
                  ErrThrow("Wrong NSTYPE ", TDatabase::ParamDB->NSTYPE, " for Newton's method !!!", "Use NSTYPE 3 or 4 !!!");
                  break;
                case 3: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=3, NSTYPE=3
                  this->N_Terms = 5;
                  this->Derivatives = {D100, D010, D001, D000, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 9;
                  this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                  this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                  this->N_Rhs = 0;
                  this->RhsSpace = { };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM==0)
                  {
                    // if(TDatabase::ParamDB->LAPLACETYPE==0)
                    //   this->AssembleParam = ;
                    // else
                    //  this->AssembleParam = ;
                  }
                  else // Skew 
                  {
                    // if(TDatabase::ParamDB->LAPLACETYPE==0)
                    //   this->AssembleParam = ;
                    // else
                    //   this->AssembleParam = ;
                  }
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  // this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
                case 4: // GALERKIN , SC_NONLIN_ITE_TYPE_SADDLE=3, NSTYPE=4
                  this->N_Terms = 5;
                  this->Derivatives = {D100, D010, D001, D000, D000};
                  this->Needs2ndDerivatives = new bool[1];
                  this->Needs2ndDerivatives[0] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 9;
                  this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
                  this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
                  this->N_Rhs = 0;
                  this->RhsSpace = { };
                  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM ==0)
                  {
                    // if(TDatabase::ParamDB->LAPLACETYPE == 0)
                    //   this->AssembleParam = ;
                    // else
                    //  this->AssembleParam = ;
                  }
                  else // Skew 
                  {
                    // if(TDatabase::ParamDB->LAPLACETYPE == 0)
                    //   this->AssembleParam = ;
                    // else
                    //   this->AssembleParam = ;
                  }
                  this->Manipulate = NULL;
                  
                  this->N_Parameters = 3;
                  this->N_ParamFct = 1;
                  // this->ParameterFct =  { NSParamsVelo3D };
                  this->N_FEValues = 3;
                  this->FEValue_FctIndex = { 0, 1, 2 };
                  this->FEValue_MultiIndex = { D000, D000, D000 };
                  this->BeginParameter = { 0 };
                  break;
              }// endswitch for NSTYPE
              break;
          }// endswitch of the iteration type (FIXED POINT OR NEWTON'S)          
          break; // break for the GALERKIN
          default:
             ErrThrow("currently DISCTYPE ", TDatabase::ParamDB->DISCTYPE,
                     " is not supported by the class NSE3D");
      } // endswitch for the DISCTYPE
      break; // endswitch for the LocalAssembling3D_type NSE3D_NonLinear
  } // endswitch (type)
}
//========================================================================
