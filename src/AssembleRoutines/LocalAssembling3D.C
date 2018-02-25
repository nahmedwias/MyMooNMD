#include <LocalAssembling3D.h>

#include <Database.h>
#include <MainUtilities.h> 
#include <FEDatabase3D.h>
#include <FEFunction3D.h>
#include <MooNMD_Io.h>
#include <string.h>

#include <DiscreteForm3D.h>
#include <ConvDiff.h>
#include <ConvDiff3D.h>
#include <NSE3D_FixPo.h>
#include <NSE3D_FixPoSkew.h>
#include <NSE3D_Param.h>
#include <NSE3D_ParamRout.h>

#include <Brinkman3D_Mixed.h>
#include <TCD3D.h> // local routines for time convection-diffusion-reaction
#include <TNSE3D_FixPo.h>
#include <TNSE3D_ParamRout.h>

#include <TNSE3DSmagorinsky.h>
#include <TNSE3DProjBasedVMS.h>
#include <TNSE3DSUPG.h>
/**
 * TODO There is still a lot of cases where the array "Needs2ndDerivatives" is
 * constructed with length 1 only although there are two spaces available
 * (Navier--Stokes, velo and pressure).
 * This will produce valgrind errors and might lead to even worse things.
 * Be prepared! And fix it eventually:
 *
 *            this->Needs2ndDerivatives = new bool[2];
 *            this->Needs2ndDerivatives[0] = false;
 *            this->Needs2ndDerivatives[1] = false;
 */

/** @brief a helper function returning a string with for the name of the 
 *         LocalAssembling3D_type. This returns an empty string in case the type
 *         is not known. */

std::string LocalAssembling3D_type_to_string(LocalAssembling3D_type type, int disctype)
{
  switch(type)
  {
    case LocalAssembling3D_type :: CD3D:
      switch(disctype)
      {
        case GALERKIN:
          return std::string("CD3D_Galerkin");
          break;
      }
      break;
    //////////////////////////////////////////////////
    case LocalAssembling3D_type::TCD3DStiffRhs:
      switch(disctype)
      {
	case GALERKIN:
	  return std::string("TCD3D_Stiff_Rhs");	  
	  break;
	case SUPG:
	  return std::string("TCD3D_Stiff_Rhs_SUPG");
	  break;
      }
      break;
    case LocalAssembling3D_type::TCD3D:
      switch(disctype)
      {
	case GALERKIN:
	  return std::string("TCD3D_AllGalerkin");
	case SUPG:
	  return std::string("TCD3D_AllSUPG");
      }
      break;
      ///////////////////////////////////////////////////////////////////////////
      // Brinkman3D: Brinkman problems
  case LocalAssembling3D_type::Brinkman3D_Galerkin:
      return std::string("Brinkman3D_Galerkin");
      break;
  case LocalAssembling3D_type::ResidualStabPkPk_for_Brinkman3D_Galerkin1:
      return std::string("ResidualStabPkPk_for_Brinkman3D_Galerkin1");
      break;
  case LocalAssembling3D_type::GradDivStab_for_Brinkman3D_Galerkin1:
      return std::string("GradDivStab_for_Brinkman3D_Galerkin1");
      break;

    default:
      return std::string();
  }
  return std::string(); //avoid compiler warning
}

LocalAssembling3D::LocalAssembling3D(LocalAssembling3D_type type, 
                                     TFEFunction3D **fefunctions3d,
                                     CoeffFct3D coeffs,
                                     int disctype)
 : type(type), discretization_type(disctype),
   name(LocalAssembling3D_type_to_string(type,disctype)), Coeffs(coeffs),
   FEFunctions3D(fefunctions3d)
{

  Output::print<5>("Constructor of LocalAssembling3D: using type ", name);
  
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
        case LocalAssembling3D_type::Brinkman3D_Galerkin:
            switch(TDatabase::ParamDB->NSTYPE)
        {
            case 4:
                this->N_Terms = 5;
                this->Derivatives = {D100, D010, D001, D000, D000};
                this->Needs2ndDerivatives = new bool[2];
                this->Needs2ndDerivatives[0] = false;
                this->Needs2ndDerivatives[1] = false;
                this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                this->N_Matrices = 15;
                this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
                this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
                //this->N_Rhs = 3;
                //CB DEBUG
                this->N_Rhs = 4;
                //END DEBUG
                this->RhsSpace = { 0, 0, 0, 1 };
                this->AssembleParam = Brinkman3DType2Galerkin;
                this->Manipulate = nullptr;
                break;
            case 14:
                //Matrix Type 14
                this->N_Terms = 5;
                this->Derivatives = { D100, D010, D001, D000, D000 };
                this->Needs2ndDerivatives = new bool[2];
                this->Needs2ndDerivatives[0] = false;
                this->Needs2ndDerivatives[1] = false;
                this->FESpaceNumber = { 0, 0, 0, 0, 1 };                               // 0: velocity, 1: pressure
                this->N_Matrices = 16;
                this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0}; //u: A11,A12,A13,A21,A22,A23,A31,A32,A33,C,B1T,B2T,B3T,B1,B2,B3 (here the lying B-Blocks come first)
                this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1}; //p; (here the standing B-Blocks come first)
                this->N_Rhs = 4;
                this->RhsSpace = { 0, 0, 0, 1 };
                this->AssembleParam = Brinkman3DType1Galerkin;
                this->Manipulate = nullptr;
                break;
        }
            break;

        case LocalAssembling3D_type::ResidualStabPkPk_for_Brinkman3D_Galerkin1:
            switch(TDatabase::ParamDB->NSTYPE)
        {
            case 14:
                //Matrix Type 14
                this->N_Terms = 11;                                                                             // = #(Derivatives)
                this->Derivatives = { D100, D010, D001, D000, D000, D100, D010, D001, D200, D020, D002};        // u_x, u_y, u_z, u, p, p_x, p_y, p_z, u_xx, u_yy, u_zz
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
                this->AssembleParam = ResidualStabPkPk_for_Brinkman3DType1Galerkin;
                this->Manipulate = nullptr;
                break;
        }
            break;

        case LocalAssembling3D_type::GradDivStab_for_Brinkman3D_Galerkin1:
            switch(TDatabase::ParamDB->NSTYPE)
        {
            case 14:
                //Matrix Type 14
                this->N_Terms = 11;                                                                             // = #(Derivatives)
                this->Derivatives = { D100, D010, D001, D000, D000, D100, D010, D001, D200, D020, D002};        // u_x, u_y, u_z, u, p, p_x, p_y, p_z, u_xx, u_yy, u_zz
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
                this->AssembleParam = GradDivStab_for_Brinkman3DType1Galerkin;
                this->Manipulate = nullptr;
                break;
        }
            break;


    ///////////////////////////////////////////////////////////////////////////
    // CD3D: stationary convection diffusion problems
    case LocalAssembling3D_type::CD3D:
      switch(this->discretization_type)
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
          this->Manipulate = nullptr;
          break;
        case SUPG:
          // second derivatives are not supported yet
          // the method is used for the convection dominant
          // => coefficient of laplace is smaller
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
          this->AssembleParam = BilinearAssemble_SD;
          this->Manipulate = nullptr;
          break;
        default:
          ErrThrow("currently DISCTYPE ", this->discretization_type,
               " is not supported by the class CD3D");
      }// endswitch this->discretization_type
      break; // break for the type LocalAssembling3D_type::CD3D 
    ///////////////////////////////////////////////////////////////////////////
    // TCD3D: nonstationary convection-diffusion-reaction problems
    case LocalAssembling3D_type::TCD3D:
      switch(this->discretization_type)
      {
	case GALERKIN:
	  this->N_Terms = 4;
          this->Derivatives = { D100, D010, D001, D000 };
          this->Needs2ndDerivatives = new bool[1];
          this->Needs2ndDerivatives[0] = false;
          this->FESpaceNumber = { 0, 0, 0, 0};
          this->N_Matrices = 2; // Mass and Stiffness Matrices
          this->RowSpace = { 0, 0 };
          this->ColumnSpace = { 0, 0 };
          this->N_Rhs = 1;
          this->RhsSpace = { 0 };
          this->AssembleParam = MatrixMARhsAssemble;
          this->Manipulate = nullptr;
	  break;
	case SUPG:
	  this->N_Terms = 4;
          this->Derivatives = { D100, D010, D001, D000 };
          this->Needs2ndDerivatives = new bool[1];
          this->Needs2ndDerivatives[0] = false;
          this->FESpaceNumber = { 0, 0, 0, 0};
          this->N_Matrices = 2; // Mass and Stiffness matrices, NOTE: M = (u, v + delta * bgradv)
          this->RowSpace = { 0, 0 };
          this->ColumnSpace = { 0, 0 };
          this->N_Rhs = 1;
          this->RhsSpace = { 0 };
          this->AssembleParam = MatricesMARhsAssemble_SUPG; 
          this->Manipulate = nullptr;
	  break;
      }
      break;
    case LocalAssembling3D_type::TCD3DStiffRhs:      
      switch(this->discretization_type)
      {
	case GALERKIN:
	  this->N_Terms = 4;
	  this->Derivatives = { D100, D010, D001, D000 };
	  this->Needs2ndDerivatives = new bool[1];
	  this->Needs2ndDerivatives[0] = false;
	  this->FESpaceNumber = { 0, 0, 0, 0 };
	  this->N_Matrices = 1;
	  this->RowSpace = { 0 };
	  this->ColumnSpace = { 0 };
	  this->N_Rhs = 1;
	  this->RhsSpace = { 0 };
	  this->Manipulate = nullptr;
	  this->Manipulate = nullptr;
	  this->AssembleParam = MatrixARhsAssemble;
	  break;
	case SUPG:
	  this->N_Terms = 4;
	  this->Derivatives = { D100, D010, D001, D000 };
	  this->Needs2ndDerivatives = new bool[1];
	  this->Needs2ndDerivatives[0] = false;
	  this->FESpaceNumber = { 0, 0, 0, 0 };
	  this->N_Matrices = 2;
	  this->RowSpace = { 0, 0 };
	  this->ColumnSpace = { 0, 0 };
	  this->N_Rhs = 1;
	  this->RhsSpace = { 0 };
	  this->Manipulate = nullptr;
	  this->Manipulate = nullptr;
	  this->AssembleParam = MatricesMARhsAssemble_SUPG;
	  break;
      }
      break;
    ///////////////////////////////////////////////////////////////////////////
    // NSE3D: stationary Navier-Stokes problems
    case LocalAssembling3D_type :: NSE3D_Linear:
    case LocalAssembling3D_type :: NSE3D_NonLinear:
      this->set_parameters_for_nse(type);
      break;
    ////////////////////////////////////////////////////////////////////////////
    // TNSE3D: nonstationary Navier-Stokes problems
    case LocalAssembling3D_type::TNSE3D_LinGAL:
    case LocalAssembling3D_type::TNSE3D_NLGAL:
    case LocalAssembling3D_type::TNSE3D_Rhs:
      switch(this->discretization_type)
      {
        case GALERKIN:
          this->set_parameters_for_tnse(type);
          break;
        case SMAGORINSKY:
	case SMAGORINSKY_COARSE:
          this->set_parameters_for_tnse_smagorinsky(type);
          break;
	case VMS_PROJECTION:
	  this->set_parameters_for_tnse_vms_projection(type);
	  break;
	case SUPG:
	  this->set_parameters_for_tnse_supg(type);
	  break;
        default:
          ErrThrow("DISCTYPE", this->discretization_type , "is not supported yet!!");
      }
      break;
    default:
      ErrThrow("Unknown or unhandled LocalAssembling3D_type case.");
  }
  
  AllOrigValues = new double** [N_Terms];
  OrigValues = new double* [N_Terms];

  // some consistency checks
  if(Coeffs == nullptr)
  {
    ErrThrow("You need to specify a valid function for the coefficients");
  }
  if(AssembleParam == nullptr)
  {
    ErrThrow("A local assembling routine was not set!");
  }
}
//========================================================================
LocalAssembling3D::LocalAssembling3D(LocalAssembling3D_type type, TAuxParam3D& aux,
                                     TDiscreteForm3D& df)
  :type(type), discretization_type(0), //default value for custom constructor
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
  if(Coeffs == nullptr)
  {
    ErrThrow("You need to specify a valid function for the coefficients");
  }
  if(this->AssembleParam == nullptr)
  {
    // this means in the discrete form there was only a pointer to a
    // AssembleFct3D rather than a AssembleFctParam3D.
    ErrThrow("can't create LocalAssembling3D object, missing AssembleFctParam3D");
  }
}
//========================================================================
LocalAssembling3D::LocalAssembling3D(
  int myN_Terms, std::vector<MultiIndex3D> myDerivatives,
  std::vector<int> myFESpaceNumber, std::vector<int> myRowSpace,
  std::vector<int> myColumnSpace, std::vector<int> myRhsSpace,
  CoeffFct3D myCoeffs, AssembleFctParam3D* myAssembleParam,
  ManipulateFct3D* myManipulate, int myN_Matrices, int myN_Rhs,
  int myN_ParamFct, std::vector<ParamFct*> myParameterFct,
  std::vector<int> myBeginParameter, int myN_Parameters,
  TFEFunction3D** myFEFunctions3D, int myN_FEValues,
  std::vector<int> myFEValue_FctIndex,
  std::vector<MultiIndex3D> myFEValue_MultiIndex,
  int discretization_type_in)
 : type{LocalAssembling3D_type::Custom}, discretization_type{discretization_type_in},
   N_Terms(myN_Terms),
   Derivatives(myDerivatives), FESpaceNumber(myFESpaceNumber),
   RowSpace(myRowSpace), ColumnSpace(myColumnSpace), RhsSpace(myRhsSpace),
   Coeffs(myCoeffs), AssembleParam(myAssembleParam), Manipulate(myManipulate),
   N_Matrices(myN_Matrices), N_Rhs(myN_Rhs), N_ParamFct(myN_ParamFct),
   ParameterFct(myParameterFct), BeginParameter(myBeginParameter),
   N_Parameters(myN_Parameters), N_FEValues(myN_FEValues),
   FEFunctions3D(myFEFunctions3D), FEValue_FctIndex(myFEValue_FctIndex),
   FEValue_MultiIndex(myFEValue_MultiIndex)
{
  // Some data members get an extra treatment - "name" is set to CUSTOMIZED,
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

  name = std::string("CUSTOMIZED");
  //Inform the world of what's going on.
  Output::print<5>("Constructor of LocalAssembling3D: using type ", name);

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
    MultiIndex3D alpha = Derivatives[i];
    int j = FESpaceNumber[i];
    if(alpha == D200 || alpha == D020 || alpha == D002 || alpha == D110 
       || alpha == D101 || alpha == D011)
      Needs2ndDerivatives[j] = true;
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
  double Temp[3 + N_FEValues];
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
    {
      switch(this->discretization_type)
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 15;
                  this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
                  this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
                  //this->N_Rhs = 3;
                  //CB DEBUG
                  this->N_Rhs = 4;
                  //END DEBUG
                  this->RhsSpace = { 0, 0, 0, 1 };
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
                  this->FESpaceNumber = { 0, 0, 0, 0, 1 }; // 0: velocity, 1: pressure
                  this->N_Matrices = 15;
                  this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 };
                  this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 };
                  //this->N_Rhs = 3;
                  //CB DEBUG
                  this->N_Rhs = 4;
                  //END DEBUG
                  this->RhsSpace = { 0, 0, 0, 1 };
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
                  this->Manipulate = nullptr;
                  
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
            ErrThrow("currently DISCTYPE ", this->discretization_type,
                     " is not supported by the class NSE3D");
      }// endswitch for the DISCTYPE 
      break; // break for the LocalAssembling3D_type NSE3D_Linear
  }
    case LocalAssembling3D_type :: NSE3D_NonLinear:
      switch(this->discretization_type)
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
                  this->Needs2ndDerivatives = new bool[2];
                  this->Needs2ndDerivatives[0] = false;
                  this->Needs2ndDerivatives[1] = false;
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
                  this->Manipulate = nullptr;
                  
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
             ErrThrow("currently DISCTYPE ", this->discretization_type,
                     " is not supported by the class NSE3D");
      } // endswitch for the DISCTYPE
      break; // endswitch for the LocalAssembling3D_type NSE3D_NonLinear
      default:
        ErrThrow("Wrong LocalAssembling3D_type for set_parameters_for_nse.");
  } // endswitch (type)
}
//========================================================================
void LocalAssembling3D::set_parameters_for_tnse(LocalAssembling3D_type la_type)
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
  this->N_Parameters = 3;
  this->N_ParamFct = 1;
  this->ParameterFct =  { TimeNSParamsVelo3D };
  this->N_FEValues = 3;
  this->FEValue_FctIndex = { 0, 1, 2 };
  this->FEValue_MultiIndex = { D000, D000, D000 };
  this->BeginParameter = { 0 };
  
  switch(la_type)
  {
    case LocalAssembling3D_type::TNSE3D_LinGAL:
      // case 0: fixed point, case 1: newton iteration
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) 
      {
        case 0: // fixed point iteration
          switch(nstype)
          {
            case 1:
            {
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
              this->AssembleParam = TimeNSType1Galerkin3D;
              this->Manipulate = nullptr;              
            }
              break;
            case 2:
            {
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
              this->AssembleParam = TimeNSType2Galerkin3D;
              this->Manipulate = nullptr;
            }
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
              
              if(laplace_type==0)
                this->AssembleParam = TimeNSType3Galerkin3D;
              else 
                this->AssembleParam = TimeNSType3GalerkinDD3D;
              
              this->Manipulate = nullptr;
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
              
              if(laplace_type==0)
                this->AssembleParam = TimeNSType4Galerkin3D;
              else 
                this->AssembleParam = TimeNSType4GalerkinDD3D;
              
              this->Manipulate = nullptr;
              break;
            case 14:
              // I have to do that 
              break;
          }
          break;
        case 1: // newton iteration
          switch(nstype)
          {
            case 1:
            case 2:
              ErrThrow("NEWTON method is only supported for NSTYPE 3 and 4");
              break;
            case 3:
              if(laplace_type==1)
                ErrThrow("Newton method only for LAPLACETYPE 0 is supported");
              break;
            case 4:
              ErrThrow("Newton method is not supported yet");
              break;
          }
          break;
        default:
          ErrThrow("SC_NONLIN_ITE_TYPE_SADDLE ",  TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE, 
                   "not supported");
      }
      break;
    // local assembling of nonlinear term
    case LocalAssembling3D_type::TNSE3D_NLGAL:
      // case 0: fixed point, case 1: newton iteration
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) 
      {
        case 0: // fixed point iteration
          switch(nstype)
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
              this->AssembleParam = TimeNSType1_2NLGalerkin3D;
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
              if(laplace_type==0)
                this->AssembleParam = TimeNSType3_4NLGalerkin3D;
              else
                this->AssembleParam = TimeNSType3_4NLGalerkinDD3D;
              this->Manipulate = nullptr;    
              break;
            case 14:
              ErrThrow("NSTYPE 14 is not supported yet");
              break;
          }
          break;
        case 1: // newton iteration
          switch(nstype)
          {
            case 1:
            case 2:
              ErrThrow("NEWTON iteration is only supported for NSTYPE 3 and 4");
              break;
            case 3:
              ErrThrow("Newton method is not supported yet");
              break;
            case 4:
              ErrThrow("Newton method is not supported yet");
              break;
          }
          break;
        default:
          ErrThrow("SC_NONLIN_ITE_TYPE_SADDLE ",  TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE, 
                   "not supported");
      }
      break;
    // local assembling of right hand side
    case LocalAssembling3D_type::TNSE3D_Rhs:
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
          this->AssembleParam =TimeNSRHS3D;
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
      ErrThrow("Wrong LocalAssembling3D_type for set_parameters_for_tnse.");
  }
}
//========================================================================
void LocalAssembling3D::set_parameters_for_tnse_smagorinsky(LocalAssembling3D_type type)
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
    case LocalAssembling3D_type::TNSE3D_LinGAL:
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
              this->AssembleParam=TimeNSType1Smagorinsky3D;
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

              this->AssembleParam=TimeNSType2Smagorinsky3D;
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
                this->AssembleParam=TimeNSType3Smagorinsky3D;
              else 
                this->AssembleParam=TimeNSType3SmagorinskyDD3D;
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
              this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
		                    0, 0, 0, 
				    1, 1, 1, 0, 0, 0 };
              this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		                    0, 0, 0, 
				    0, 0, 0, 1, 1, 1 };
              this->N_Rhs = 4;
              this->RhsSpace = { 0, 0, 0, 1};

              this->Manipulate = nullptr;
              if(TDatabase::ParamDB->LAPLACETYPE==0)
                this->AssembleParam=TimeNSType4Smagorinsky3D;
              else 
                this->AssembleParam=TimeNSType4SmagorinskyDD3D;
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
    case LocalAssembling3D_type::TNSE3D_NLGAL:
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
              this->AssembleParam = TimeNSType1_2NLSmagorinsky3D;
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
                this->AssembleParam = TimeNSType3_4NLSmagorinsky3D;
              else
              {
                this->N_Matrices = 9;
                this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
                this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
                this->AssembleParam = TimeNSType3_4NLSmagorinskyDD3D;
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
   case LocalAssembling3D_type::TNSE3D_Rhs:
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
          this->AssembleParam = TimeNSRHS3D;
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

void LocalAssembling3D::set_parameters_for_tnse_vms_projection(LocalAssembling3D_type type)
{
  int nstype = TDatabase::ParamDB->NSTYPE;
  if(nstype !=3 && nstype !=4)
  {
    ErrThrow(" VMS is only supported for NSTYE 3 and 4");
  }
  if(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE)
  {
    ErrThrow(" VMS in only supported for SC_NONLIN_ITE_TYPE_SADDLE ",
    "fixed point iteration, i.e., SC_NONLIN_ITE_TYPE_SADDLE = 0");
  }
  
  this->N_Parameters = 16;
  this->N_ParamFct = 1;
  this->ParameterFct =  { TimeNSParamsVelo_GradVelo_LargeScale3D };
  this->N_FEValues = 13;
  this->FEValue_FctIndex = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 3};

  this->FEValue_MultiIndex = { D000, D000, D000, // u1old, u2old, u3old
                               D100, D100, D100, // u1old_x, u2old_x, u3old_x
                               D010, D010, D010, // u1old_y, u2old_y, u3old_y
                               D001, D001, D001, // u1old_z, u2old_z, u3old_z
                               D000 // pold
  };
  this->BeginParameter = { 0 };
 // switch over all (linear + nonlinear ) and nonlinear (only ) types
  switch(type)
  {
    case LocalAssembling3D_type::TNSE3D_LinGAL:
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) 
      {
        case 0: // fixed point iteration
          switch(nstype)
          {
            case 3:
	      ErrThrow("NSTYPE 3 is not yet supported for VMS_PROJECTION");
              break;
            case 4:
              this->N_Terms = 6;
              this->Derivatives = {D100, D010, D001, D000, D000, D000};
              this->Needs2ndDerivatives = new bool[3];
              this->Needs2ndDerivatives[0] = false;
              this->Needs2ndDerivatives[1] = false;
              this->Needs2ndDerivatives[2] = false;
              this->FESpaceNumber = { 0, 0, 0, 0, 1, 2 }; // 0: velocity, 1: pressure, 2: projection
              this->N_Matrices = 25;
              this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, // A-block
                                    0, 0, 0,// Mass 
                                    2, // L-matrix (vms)
                                    1, 1, 1, 0, 0, 0, // B-block
                                    0, 0, 0, 2, 2, 2 }; // vms matrices
              this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                    0, 0, 0, // Mass 
                                    2, // L-matrix (vms)
                                    0, 0, 0, 1, 1, 1, // B-block
                                    2, 2, 2, 0, 0, 0}; // vms matrices
              this->N_Rhs = 3;
              this->RhsSpace = { 0, 0, 0};
              this->Manipulate = NULL;
              this->AssembleParam=TimeNSType4VMS_ProjectionDD3D;
              break;
          }
          break; // fixed point iteration
        case 1: // newton iteration
          ErrThrow("Set parameters for newton iteration");
          break; // newton iteration
      }
      break; // LocalAssembling3D_type::TNSE3D_LinGAL:
    case LocalAssembling3D_type::TNSE3D_NLGAL:
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) 
      {
        case 0:// fixed point iteration
          switch(nstype)
          {
            case 3: case 4:
              this->N_Terms = 5;
              this->Derivatives = {D100, D010, D001, D000, D000};
              this->Needs2ndDerivatives = new bool[2];
              this->Needs2ndDerivatives[0] = false;
              this->Needs2ndDerivatives[1] = false;
              this->FESpaceNumber = { 0, 0, 0, 0, 2 }; // 0: velocity, 2: projection
              this->N_Matrices = 12;
              this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, // A - block 
                                    0, 0, 0 }; // vms matrices
              this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                    2, 2, 2}; // vms matrices
              this->N_Rhs = 0;
              this->RhsSpace = { };
          
              this->Manipulate = NULL;
              this->AssembleParam=TimeNSType3_4NLVMS_ProjectionDD3D;
              break;
          }
          break;// fixed point iteration 
        case 1: // newton iteration 
          ErrThrow("Set parameters for newton iteration");
          break; // newton iteration
      }
      break; // LocalAssembling3D_type::TNSE3D_NLGAL:
    case LocalAssembling3D_type::TNSE3D_Rhs:
      switch(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE) 
      {
        case 0: // fixed point iteration 
          this->N_Terms = 1;
          this->Derivatives = { D000 };
          this->Needs2ndDerivatives = new bool[2];
          this->Needs2ndDerivatives[0] = false;
          this->Needs2ndDerivatives[1] = false;
          this->FESpaceNumber = { 0 }; // 0: velocity, 1: pressure
          this->N_Matrices = 0;
          this->RowSpace = { };
          this->ColumnSpace = { };
          this->N_Rhs = 3 ; 
          this->RhsSpace = {0, 0, 0};
          this->AssembleParam = TimeNSRHS3D;
          this->Manipulate = NULL;
         break;
        case 1: // newton iteration 
          ErrThrow("Set parameters for newton iteration");
          break; // newton iteration
      }
      break;
    default:
      ErrThrow("Unknown LocalAssembling3D_type");
  }
}

void LocalAssembling3D::set_parameters_for_tnse_supg(LocalAssembling3D_type type)
{
  
  if(TDatabase::ParamDB->NSTYPE < 4 )
  { 
    ErrThrow("SUPG method is only supported for NSTYPE 4 and 14 ", 
             TDatabase::ParamDB->NSTYPE);
  }
  if(TDatabase::ParamDB->NSTYPE==4)
  {
    this->N_Parameters = 5;
    this->N_ParamFct = 1;
    this->ParameterFct =  { TimeNSType4Params_SUPG };
    this->N_FEValues = 3;
    this->FEValue_MultiIndex = { D000, D000, D000};
    this->FEValue_FctIndex = { 0, 1, 2};
    this->BeginParameter = { 0 };
  }
  if(TDatabase::ParamDB->NSTYPE==14)
  {
    this->N_Parameters = 9;
    this->N_ParamFct = 1;
    this->ParameterFct =  { TimeNSType14Params_SUPG };
    this->N_FEValues = 6;
    this->FEValue_MultiIndex = { D000, D000, D000, 
                                 D000, D000, D000 };
    this->FEValue_FctIndex = { 0, 1, 2, 
                               3, 4, 5 };
    this->BeginParameter = { 0 };
  }
  
  switch(type)
  {
    case LocalAssembling3D_type::TNSE3D_LinGAL:
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 4:
          this->N_Terms = 11;
          this->Derivatives = {D100, D010, D001, D000, // u_x, u_y, u_z, u
                               D000, D100, D010, D001, // p, p_x, p_y, p_z
                               D200, D020, D002};      // u_xx, u_yy, u_zz
          this->Needs2ndDerivatives = new bool[3];
          this->Needs2ndDerivatives[0] = true;
          this->Needs2ndDerivatives[1] = true;
          this->Needs2ndDerivatives[2] = true;
          this->FESpaceNumber = { 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 }; // 0: velocity, 1: pressure
          this->N_Matrices = 18;
          this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, // A-Block
	                        0, 0, 0, // Mass matrices
		                1, 1, 1, // BT Blocks
				0, 0, 0 };
          this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, // A blocks
		                0, 0, 0, // Mass matrices
				0, 0, 0, // BT Blocks
				1, 1, 1};
          this->N_Rhs = 3;
          this->RhsSpace = { 0, 0, 0 };
          this->AssembleParam = TimeNSType4SUPGDD3D;
          this->Manipulate = nullptr;
          break; // NSTYPE4
        case 14:
          // TODO 
          ErrThrow("NSTYPE 14 for SUPG is not supported yet");
          break; // NSTYPE14 
      }
      break; // LocalAssembling3D_type::TNSE3D_LinGAL:
    case LocalAssembling3D_type::TNSE3D_NLGAL:
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 4:
          this->N_Terms = 11;
          this->Derivatives = {D100, D010, D001, D000, // u_x, u_y, u_z, u
                               D000, D100, D010, D001, // p, p_x, p_y, p_z
	                      D200, D020, D002};      // u_xx, u_yy, u_zz
          this->Needs2ndDerivatives = new bool[3];
          this->Needs2ndDerivatives[0] = true;
          this->Needs2ndDerivatives[1] = true;
          this->Needs2ndDerivatives[2] = true;
          this->FESpaceNumber = { 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 }; // 0: velocity, 1: pressure
          this->N_Matrices = 15;
          this->RowSpace    = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
          this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1};
          this->N_Rhs = 3;
          this->RhsSpace = {0, 0, 0 };
          this->AssembleParam = TimeNSType4NLSUPGDD3D;
          this->Manipulate = nullptr;
          break; // NSTYPE4
        case 14:
          // TODO 
          ErrThrow("NSTYPE 14 for SUPG is not supported yet");
          break; // NSTYPE14 
      }
      break; // LocalAssembling3D_type::TNSE3D_NLGAL:
    case LocalAssembling3D_type::TNSE3D_Rhs:
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 4:
          this->N_Terms = 4;
          this->Derivatives = {D100, D010, D001, D000}; // u_x, u_y, u_z, u
          this->Needs2ndDerivatives = new bool[3];
          this->Needs2ndDerivatives[0] = false;
          this->Needs2ndDerivatives[1] = false;
          this->Needs2ndDerivatives[2] = false;
          this->FESpaceNumber = {0, 0, 0, 0 }; // 0: velocity, 1: pressure
          this->N_Matrices = 0; // 3BT matrices and 1 mass matrix
          this->RowSpace    = { };
          this->ColumnSpace = { };
          this->N_Rhs = 3;
          this->RhsSpace = { 0, 0, 0 };
          this->AssembleParam = TimeNSType4RHSSUPG3D;
          this->Manipulate = nullptr;
          break; // NSTYPE4
        case 14:
          // TODO 
          ErrThrow("NSTYPE 14 for SUPG is not supported yet");
          break; // NSTYPE14 
      }
      break; // LocalAssembling3D_type::TNSE3D_Rhs:
      default:
      ErrThrow("Wrong LocalAssembling3D_type for set_parameters_for_tnse.");
  }
}
