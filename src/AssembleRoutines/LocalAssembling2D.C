#include <LocalAssembling2D.h>
#include <Database.h>
#include <MainUtilities.h> // linfb, ave_l2b_quad_points
#include <FEDatabase2D.h>
#include <FEFunction2D.h>
#include <MooNMD_Io.h>
#include <string.h>
#include <DiscreteForm2D.h> // to be removed
#include <Brinkman2D_Mixed.h>// local assembling routines for 2D Navier-Stokes
#include <ConvDiff.h>
#include <ConvDiff2D.h> // local assembling routines for 2D convection-diffusion
#include <DarcyMixed.h> // local assembling routines for 2D Darcy problems
#include <NSE2D_EquOrd_FixPo.h> // local assembling routines for equal order elements
#include <NSE2DGalerkin.h>// local assembling routines for 2D Navier-Stokes
#include <NSE2DSUPG.h>// local assembling routines for 2D Navier-Stokes
#include <NSE2D_Newton.h>

#include "TNSE2DGalerkin.h"
#include "NSE2DGalerkin.h"
#include "TNSE2DSUPG.h"
#include "TNSE2DResBasedVMS.h"

#include "CommonRoutineTNSE2D.h"

//==============================================================================
/** @brief a helper function returning a string with the name of the
 *         LocalAssembling2D_type. This returns an empty string in case the type
 *         is not known. */

std::string LocalAssembling2D_type_to_string(LocalAssembling2D_type type, int disctype)
{
    switch(type)
    {
            ///////////////////////////////////////////////////////////////////////////
            // CD2D: stationary convection diffusion problems
        case LocalAssembling2D_type::ConvDiff:
            switch(disctype)
        {
            case GALERKIN:
                if(TDatabase::ParamDB->Axial3D)
                    return std::string("CD2D_Axiax3D_Galerkin");
                else
                    return std::string("CD2D_Galerkin");
            case SUPG:
                return std::string("CD2D_SUPG");
            case GLS:
                return std::string("CD2D_GLS");
        }
            break;
            ///////////////////////////////////////////////////////////////////////////
            // TCD2D: time dependent convection diffusion problems
        case LocalAssembling2D_type::TCD2D:
            switch(disctype)
        {
            case GALERKIN:
                return std::string("TCD2D_Stiff_Rhs");
            case SUPG:
                return std::string("TCD2D_Stiff_Rhs_SUPG");
        }
        case LocalAssembling2D_type::TCD2D_Mass:
            return std::string("TCD2D_Mass");
            ///////////////////////////////////////////////////////////////////////////
            // NSE2D: stationary Navier-Stokes problems
        case NSE2D_All:
            return std::string("NSE2D_All");
        case NSE2D_NL:
            return std::string("NSE2D_NL");
            ///////////////////////////////////////////////////////////////////////////
        case NSE2D_SUPG:
            return std::string("NSE2D_SUPG");
        case NSE2D_SUPG_NL:
            return std::string("NSE2D_SUPG_NL");
            ///////////////////////////////////////////////////////////////////////////
            // Darcy2D: stationary Darcy problems
        case Darcy2D_Galerkin:
            return std::string("Darcy2D_Galerkin");
            ///////////////////////////////////////////////////////////////////////////
            // Brinkman2D: Brinkman problems
        case Brinkman2D_Galerkin1:
            return std::string("Brinkman2D_Galerkin1");
            
        case Brinkman2D_Galerkin2:
            return std::string("Brinkman2D_Galerkin2");
            
        case Brinkman2D_Galerkin1ResidualStabP1:
            return std::string("Brinkman2D_Galerkin1ResidualStabP1");
            
        case Brinkman2D_Galerkin1ResidualStabP2:
            return std::string("Brinkman2D_Galerkin1ResidualStabP2");
        
        case Brinkman2D_GradDivStabilization:
            return std::string("Brinkman2D_GradDivStabilization");

            ///////////////////////////////////////////////////////////////////////////
            // TNSE2D: nonstationary Navier-Stokes
        case LocalAssembling2D_type::TNSE2D:
            switch(disctype)
        {
            case GALERKIN:
                return std::string("TNSE2D");
            case SUPG:
                return std::string("TNSE2D_SUPG");
        }
            break;
        case LocalAssembling2D_type::TNSE2D_NL:
            switch(disctype)
        {
            case GALERKIN:
                return std::string("TNSE2D_NLGalerkin");
            case SUPG:
                return std::string("TNSE2D_NLSUPG");
        }
            break;
        case LocalAssembling2D_type::TNSE2D_Rhs:
            switch(disctype)
        {
            case GALERKIN:
                return std::string("TNSE2D_Rhs");
            case SUPG:
                return std::string("TNSE2D_RhsSUPG");
        }
            break;
        case LocalAssembling2D_type::Custom:
            return std::string("customized");
    }
    return std::string();
    
    
}

//==============================================================================
LocalAssembling2D::LocalAssembling2D(LocalAssembling2D_type type, 
                                     TFEFunction2D **fefunctions2d,
                                     CoeffFct2D coeffs,
                                     int disctype)
 : db(default_local_assembling_database()),type(type),
   discretization_type(disctype),
   name(LocalAssembling2D_type_to_string(type, disctype)), Coeffs(coeffs),
   FEFunctions2D(fefunctions2d)
{
    Output::print<5>("Constructor of LocalAssembling2D: using type ", name);
    
    // the values below only matter if you need an existing finite element
    // function during your assembly. Change them in such a case
    this->N_Parameters = 0;
    this->N_ParamFct = 0;
    this->ParameterFct = {};
    this->N_FEValues = 0;
    this->FEValue_FctIndex = {};
    this->FEValue_MultiIndex = {};
    this->BeginParameter = {};

    // set all member variables according to the LocalAssembling2D_type
    switch(type)
    {
            ///////////////////////////////////////////////////////////////////////////
            // CD2D: stationary convection diffusion problems
        case LocalAssembling2D_type::ConvDiff:
            this->N_Matrices = 1;
            this->RowSpace = { 0 };
            this->ColumnSpace = { 0 };
            this->N_Rhs = 1;
            this->RhsSpace = { 0 };
            switch(this->discretization_type)
        {
            case GALERKIN:
                this->N_Terms = 3;
                this->Derivatives = { D10, D01, D00 };
                this->Needs2ndDerivatives = new bool[1];
                this->Needs2ndDerivatives[0] = false;
                this->FESpaceNumber = { 0, 0, 0 };
                
                if(TDatabase::ParamDB->Axial3D)
                {
                  this->local_assemblings_routines.push_back(
                      BilinearAssemble_Axial3D);
                }
                else
                {
                  this->local_assemblings_routines.push_back(
                    BilinearAssembleGalerkin);
                }
                this->Manipulate = nullptr;
                break;
            case SUPG:
            case GLS:
                this->N_Terms = 5;
                this->Derivatives = { D10, D01, D00, D20, D02 };
                this->Needs2ndDerivatives = new bool[1];
                this->Needs2ndDerivatives[0] = true;
                this->FESpaceNumber = { 0, 0, 0, 0, 0 };
                this->local_assemblings_routines.push_back(
                    BilinearAssembleGalerkin);
                if(this->discretization_type==SUPG)
                {
                  this->local_assemblings_routines.push_back(
                    BilinearAssemble_SD);
                }
                else
                {
                  this->local_assemblings_routines.push_back(
                    BilinearAssemble_GLS);
                }
                
                if(TDatabase::ParamDB->SDFEM_NORM_B==0)
                    this->Manipulate = linfb;
                else
                    this->Manipulate = ave_l2b_quad_points;
                
                break;
            default:
                ErrMsg("currently DISCTYPE " << this->discretization_type <<
                       " is not supported by the class CD2D");
                throw("unsupported DISCTYPE");
        }
            break;
            ///////////////////////////////////////////////////////////////////////////
            // TCD2D: time dependent convection diffusion problems
        case LocalAssembling2D_type::TCD2D:
            this->N_Matrices = 1;
            this->RowSpace = { 0 };
            this->ColumnSpace = { 0 };
            this->N_Rhs = 1;
            this->RhsSpace = { 0 };
            this->Manipulate = nullptr;
            
            switch(this->discretization_type)
        {
            case GALERKIN:
                this->N_Terms = 3;
                this->Derivatives = { D10, D01, D00 };
                this->Needs2ndDerivatives = new bool[1];
                this->Needs2ndDerivatives[0] = false;
                this->FESpaceNumber = { 0, 0, 0 };
                this->local_assemblings_routines.push_back(LocalMatrixARhs);
                break;
            case SUPG:
            case GLS:
                this->N_Terms = 5;
                this->Derivatives = { D10, D01, D00, D20, D02 };
                this->Needs2ndDerivatives = new bool[1];
                this->Needs2ndDerivatives[0] = true;
                this->FESpaceNumber = { 0, 0, 0, 0, 0 }; // number of terms = 5
                
                if(this->discretization_type==SUPG)
                {
                  this->local_assemblings_routines.push_back(
                    LocalMatrixARhs_SUPG);
                }
                else
                {
                    ErrMsg("currently DISCTYPE " << this->discretization_type <<
                           " is not supported by the class CD2D");
                    throw("unsupported DISCTYPE");
                }
                break;
        }
            break;// case LocalAssembling2D_type::TCD2D:
        case LocalAssembling2D_type::TCD2D_Mass:
            this->N_Matrices = 1;
            this->RowSpace = { 0 };
            this->ColumnSpace = { 0 };
            this->N_Rhs = 0;
            this->RhsSpace = { 0 };
            this->Manipulate = nullptr;
            switch(this->discretization_type)
        {
            case GALERKIN:
                this->N_Terms = 1;
                this->Derivatives = { D00 };
                this->Needs2ndDerivatives = new bool[1];
                this->Needs2ndDerivatives[0] = false;
                this->FESpaceNumber = { 0 };
                this->local_assemblings_routines.push_back(LocalMatrixM);
                break;
            case SUPG:
                this->N_Terms = 3;
                this->Derivatives = { D10, D01, D00 };
                this->Needs2ndDerivatives = new bool[1];
                this->Needs2ndDerivatives[0] = false;
                this->FESpaceNumber = { 0, 0, 0 };
                this->local_assemblings_routines.push_back(LocalMatrixM_SUPG);
                break;
        }
            break;  //LocalAssembling2D_type::TCD2D_Mass
            ///////////////////////////////////////////////////////////////////////////
            // Brinkman2D: problems and Brinkman problem
            
        case LocalAssembling2D_type::Brinkman2D_Galerkin1:
            //Matrix Type 14
           
            this->N_Terms = 4;
            this->Derivatives = { D10, D01, D00, D00 };
            this->Needs2ndDerivatives = new bool[2];
            this->Needs2ndDerivatives[0] = false;
            this->Needs2ndDerivatives[1] = false;
            this->FESpaceNumber = { 0, 0, 0, 1 };                               // 0: velocity, 1: pressure
            this->N_Matrices = 9;
            this->RowSpace =    { 0, 0, 0, 0, 1, 1, 1, 0, 0}; // A11, A12, A21, A22, C, B1, B2, B1T, B2T ; corresp. to u
            this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1}; // A11, A12, A21, A22, C, B1, B2, B1T, B2T ; corresp. to p
            this->N_Rhs = 3;
            this->RhsSpace = { 0, 0, 1 };
            this->local_assemblings_routines.push_back(BrinkmanType1Galerkin);
            this->Manipulate = nullptr;
            break;
            
        case LocalAssembling2D_type::Brinkman2D_Galerkin2:
            //Matrix Type 14
            this->N_Terms = 6;                                                  // = #(Derivatives)
            this->Derivatives = { D10, D01, D00, D00, D10, D01};                // u_x, u_y, u, p, p_x, p_y
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
            
        case LocalAssembling2D_type::Brinkman2D_Galerkin1ResidualStabP1:
            //Matrix Type 14
            this->N_Terms = 6;                                                  // = #(Derivatives)
            this->Derivatives = { D10, D01, D00, D00, D10, D01};                // u_x, u_y, u, p, p_x, p_y
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
            
        case LocalAssembling2D_type::Brinkman2D_Galerkin1ResidualStabP2:
            //Matrix Type 14
            this->N_Terms = 8;                                                  // = #(Derivatives)
            this->Derivatives = { D10, D01, D00, D00, D10, D01, D20, D02};      // u_x, u_y, u, p, p_x, p_y, u_xx, u_yy
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

        case LocalAssembling2D_type::Brinkman2D_GradDivStabilization:
            this->N_Terms = 2;                                                  // = #(Derivatives)
            this->Derivatives = { D10, D01 };                                   // u_x, u_y, u, p, p_x, p_y, u_xx, u_yy
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
        


            ///////////////////////////////////////////////////////////////////////////
            // NSE2D: stationary Navier-Stokes problems problem
        case NSE2D_All:
        case NSE2D_NL:
            this->set_parameters_for_nseGalerkin(type);
            break;
            ///////////////////////////////////////////////////////////////////////////
        case Darcy2D_Galerkin:
            this->N_Terms = 6;
            this->Derivatives = { D00, D00, D10, D01, D10, D01 };
            this->Needs2ndDerivatives = new bool[2];
            this->Needs2ndDerivatives[0] = false;
            this->Needs2ndDerivatives[1] = false;
            this->FESpaceNumber = { 0, 1, 0, 0, 1, 1};
            this->N_Matrices = 4;
            this->RowSpace = {0, 1, 0, 1};
            this->ColumnSpace = { 0, 1, 1, 0};
            this->N_Rhs = 2;
            this->RhsSpace = { 0, 1 };
            this->local_assemblings_routines.push_back(
              BilinearAssembleDarcyGalerkin<2>);
            this->Manipulate = nullptr;
            break;
            ////////////////////////////////////////////////////////////////////////////
            // TNSE2D: nonstationary Navier-Stokes problems
        case LocalAssembling2D_type::TNSE2D:
        case LocalAssembling2D_type::TNSE2D_NL:
        case LocalAssembling2D_type::TNSE2D_Rhs:
	  switch(discretization_type){
	    case GALERKIN:
	    case LOCAL_PROJECTION:
	      this->set_parameters_for_tnseGalerkin(type);
	      break;
	    case SUPG:
	      this->set_parameters_for_tnseSUPG(type);
	      break;
	    case SMAGORINSKY:
	      this->set_parameters_for_tnseSmagorinsky(type);
	      break;
	    case RBVMS:
	      this->set_parameters_for_tnseRBVMS(type);
	      break;
	    default:
	      ErrThrow("space discretization type: ", discretization_type, " is not supported");
	  }
	  // this->set_parameters_for_tnse(type);
	  break;
            
        default:
            ErrMsg("unknown LocalAssembling2D_type " << type << " " << this->name);
            throw("unknown LocalAssembling2D_type");
    }
    
    AllOrigValues = new double** [N_Terms];
    OrigValues = new double* [N_Terms];
    
    // some consistency checks
    if(Coeffs == nullptr)
    {
        ErrThrow("You need to specify a valid function for the coefficients");
    }
    if(this->local_assemblings_routines.empty() 
     || this->local_assemblings_routines[0] == nullptr)
    {
      ErrThrow("A local assembling routine was not set!");
    }
    parameter_functions_values.resize(0);
}

//==============================================================================
/*! @brief Customized constructor. */
LocalAssembling2D::LocalAssembling2D(int myN_Terms,
		std::vector<MultiIndex2D> myDerivatives, std::vector<int> myFESpaceNumber,
		std::vector<int> myRowSpace, std::vector<int> myColumnSpace, std::vector<int> myRhsSpace,
		CoeffFct2D myCoeffs, AssembleFctParam2D* myAssembleParam, ManipulateFct2D* myManipulate,
		int myN_Matrices, int myN_Rhs,
		int myN_ParamFct, std::vector<ParamFct*> myParameterFct, std::vector<int> myBeginParameter, int myN_Parameters,
		TFEFunction2D **myFEFunctions2D,  int myN_FEValues,
		std::vector<int> myFEValue_FctIndex, std::vector<MultiIndex2D> myFEValue_MultiIndex)

: db(default_local_assembling_database()), type{LocalAssembling2D_type::Custom},
  N_Terms(myN_Terms), Derivatives(myDerivatives), FESpaceNumber(myFESpaceNumber),
  RowSpace(myRowSpace), ColumnSpace(myColumnSpace), RhsSpace(myRhsSpace),
  Coeffs(myCoeffs), local_assemblings_routines({myAssembleParam}),
  Manipulate(myManipulate),
  N_Matrices(myN_Matrices), N_Rhs(myN_Rhs),
  N_ParamFct(myN_ParamFct), ParameterFct(myParameterFct), BeginParameter(myBeginParameter), N_Parameters(myN_Parameters),
  N_FEValues(myN_FEValues), FEFunctions2D(myFEFunctions2D),
  FEValue_FctIndex(myFEValue_FctIndex), FEValue_MultiIndex(myFEValue_MultiIndex)

{
  // Some data members get an extra treatment - "name" is set to CUSTOMIZED,
  // The auxiliary arrays (All)OrigValues are dynamically allocated with size "N_Terms".
  // "N_Spaces" is determined by finding the max in "FESpaceNumber" (+1).
  // "Needs2ndDerivative" is dynamically allocated to the size "N_Spaces" and then filled
  // according to the appearance of "D20", "D11" or "D02" in "Derivatives".
  
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
  Output::print<5>("Constructor of LocalAssembling2D: using type ", name);

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
		MultiIndex2D alpha = Derivatives[i];
		int j = FESpaceNumber[i];
		if(alpha == D20 || alpha == D11 || alpha == D02)
			Needs2ndDerivatives[j] = true;
	}
	//END code taken from TDiscretForm2D::TDiscreteForm2D(...)

	parameter_functions_values.resize(0);

}

//==============================================================================
LocalAssembling2D::~LocalAssembling2D()
{
  delete [] AllOrigValues;
  delete [] OrigValues;
  delete [] Needs2ndDerivatives;
}

//==============================================================================
ParameterDatabase LocalAssembling2D::default_local_assembling_database()
{
  ParameterDatabase db("default local assembling database");
  
  db.add("with_coriolis_force", false,
         "include the coriolis force for (Navier--)Stokes in 3D. This requires "
         "special local assemblings and the (pde-) coefficients of the example "
         "must include the coriolis force vector Omega.");
  
  db.add("laplace_type_deformation", false, 
         "determine the way the laplacian is discretized.");
  
  db.add("nse_nonlinear_form", "convective",
         "Determine how the nonlinear term for Navier--Stokes is assembled.",
         {"convective", "skew-symmetric", "rotational"});
  
  return db;
}


//==============================================================================
void LocalAssembling2D::GetLocalForms(int N_Points,
				      double *weights, 
                                      double *AbsDetjk,
				      double *X, double *Y,
                                      int *N_BaseFuncts,
                                      BaseFunct2D *BaseFuncts,
                                      double **Parameters,
				      double **AuxArray,
                                      TBaseCell *Cell, int N_Matrices,
                                      int N_Rhs,
                                      double ***LocMatrix, double **LocRhs,
                                      double factor)
{
    const double hK = Cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);
    
    //this->GetParameters(N_Points, nullptr, Cell, ??Cell->GetCellIndex(), X, Y,
    //                    Parameters);
    
    for(int i=0; i<N_Matrices; ++i)
    {
        double **CurrentMatrix = LocMatrix[i];
        int N_Rows = N_BaseFuncts[RowSpace[i]];
        int N_Columns = N_BaseFuncts[ColumnSpace[i]];
        for(int j=0;j<N_Rows;j++)
        {
            double *MatrixRow = CurrentMatrix[j];
            memset(MatrixRow, 0, SizeOfDouble*N_Columns);
        } // endfor j
    } // endfor i
    
    for(int i=0; i<N_Rhs; ++i)
    {
        int N_Rows = N_BaseFuncts[RhsSpace[i]];
        memset(LocRhs[i], 0, SizeOfDouble*N_Rows);
    }
    
    // *****************************************************
    // for 2Phase flow problems (Sashikumaar Ganesan)
    AuxArray[0][0] = Cell->GetPhase_ID();
    // *****************************************************
    
    if(Coeffs)
        Coeffs(N_Points, X, Y, Parameters, AuxArray);
    
    if(Manipulate)
        Manipulate(N_Points, AuxArray, Parameters, Cell);
    
    for(int i=0; i<N_Terms; ++i)
    {
        AllOrigValues[i] =
        TFEDatabase2D::GetOrigElementValues(BaseFuncts[FESpaceNumber[i]],
                                            Derivatives[i]);
    }
    
    for(int i=0; i<N_Points; ++i)
    {
        double Mult = weights[i] * AbsDetjk[i] * factor;
        double *Coeff = AuxArray[i];
        Coeff[19] = AbsDetjk[i];
        
        if(TDatabase::ParamDB->Axial3DAxis == 1)
        {
            // r in axial3D (X: symmetric) problems (Sashikumaar Ganesan)
            Coeff[20] = Y[i];
        }
        else
        {
            // r in axial3D (Y: symmetric) problems (Sashikumaar Ganesan)
            Coeff[20] = X[i];
        }
        
        double *Param = Parameters[i];
        
        for(int j=0; j<N_Terms; j++)
            OrigValues[j] = AllOrigValues[j][i];
        
        for(auto lar : local_assemblings_routines)
        {
          lar(Mult, Coeff, Param, hK, OrigValues, N_BaseFuncts, LocMatrix,
              LocRhs);
        }
    } // end loop over quadrature points 
}



//HIER/////////////////////////////////////////////////////////////////////////
//==============================================================================
//==============================================================================
void LocalAssembling2D::get_local_forms(int N_Points,
                                        double *weights,
                                        double *AbsDetjk,
                                        double *X,
                                        double *Y,
                                        int *N_BaseFuncts,
                                        BaseFunct2D *BaseFuncts,
                                        //double **AuxArray,
                                        TBaseCell *Cell,
                                        int N_Matrices,
                                        int N_Rhs,
                                        double ***LocMatrix,
                                        double **LocRhs,
                                        double factor)
{
  // get cell measure (depending on the input parameter CELL_MEASURE)
  const double hK = Cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);
    
  ///@todo this is only a temporary implementation  (double** Parameters) to improve the old code
  double **Parameters;
  Parameters = new double*[this->parameter_functions_values.size()];
  for(unsigned int i=0; i<this->parameter_functions_values.size(); i++)
  {
    Parameters[i] = &this->parameter_functions_values[i][0];
  }
  

    // allocate the memory for matrices
  for(int i = 0; i < N_Matrices; i++)
  {
   // CurrenMatrix: a pointer to i-th local matrix
    double **CurrentMatrix = LocMatrix[i];
     
      int N_Rows = N_BaseFuncts[RowSpace[i]];
   int N_Columns = N_BaseFuncts[ColumnSpace[i]];
   for(int j=0;j<N_Rows;j++)
    {
      // allocate the memory of the Row pointer
      double *MatrixRow = CurrentMatrix[j];
      memset(MatrixRow, 0, SizeOfDouble*N_Columns);
    } // endfor j
  } // endfor i

  // allocate mmemory for rhs
  for(int i=0; i<N_Rhs; i++)
  {
    int N_Rows = N_BaseFuncts[RhsSpace[i]];
    memset(LocRhs[i], 0, SizeOfDouble*N_Rows);
  }

  // get the values (and derivatives) of basis functions
  // the type of values needed are specified in Derivatives[i], e.g. D00,D01,...
  ///@todo can this be done somewhere else? (it does not depend on function input)
  for(int i=0; i<this->N_Terms; ++i)
  {
    this->AllOrigValues[i] =
      TFEDatabase2D::GetOrigElementValues(BaseFuncts[this->FESpaceNumber[i]],
					  this->Derivatives[i]);
  }

  // *****************************************************
  // for 2Phase flow problems (Sashikumaar Ganesan)
  ///@attention Alfonso commented out this part.
  //AuxArray[0][0] = Cell->GetPhase_ID();
  // *****************************************************

  // the double** coefficient_values replaces the double** AuxArray, with
  // AuxArray[i] = values of coefficients (or parameters) at Gauss point i
  // note: the (used) length of AuxArray[i] changes depending on the
  // example. In the original implementation AuxArray[i] was allocated as 40*double
  // the function Coeffs (depending on Example) fills AuxArray
  double **coefficient_values;
  coefficient_values = new double*[N_Points];
  for(int i=0; i<N_Points; ++i)
  {
    // do not allocate less than 20 (see below)
    coefficient_values[i] = new double[30];
  }
  
  if(Coeffs)
    Coeffs(N_Points, X, Y, Parameters, coefficient_values);

  ///@attention the Manipulate function appears not to be used in the current implementation
  // Alfonso commented out this part.
  //if(Manipulate)
  //  Manipulate(N_Points, AuxArray, Parameters, Cell);
  
  // loop over Gauss points
  for(int i=0; i<N_Points; ++i)
  {
    double Mult = weights[i] * AbsDetjk[i] * factor;
    // coefficient_values[i][19] is set by hand
    coefficient_values[i][19] = AbsDetjk[i];
    
    if(TDatabase::ParamDB->Axial3DAxis == 1)
      coefficient_values[i][20] = Y[i]; // r in axial3D (X: symmetric) problems (Sashikumaar Ganesan)
    else
      coefficient_values[i][20] = X[i]; // r in axial3D (Y: symmetric) problems (Sashikumaar Ganesan)

    // OrigValues: values of basis functions (e.g., D00, D01, D11, ...) at Gauss(i)
    for(int j=0; j<N_Terms; j++)
      OrigValues[j] = AllOrigValues[j][i];
    
    for(auto lar : local_assemblings_routines)
    {
      lar(Mult, coefficient_values[i], Parameters[i], hK, OrigValues,
          N_BaseFuncts, LocMatrix, LocRhs);
    }
  } // end loop over quadrature points
 
  delete [] Parameters;
}


//==============================================================================
void LocalAssembling2D::GetLocalForms(int N_Points,
                                      double *weights,
                                      double *AbsDetjk,
                                      double *X,
                                      double *Y,
                                      int *N_BaseFuncts,
                                      BaseFunct2D *BaseFuncts,
                                      TBaseCell *Cell,
                                      double ***LocMatrix,
                                      double **LocRhs,
                                      double factor)
{
    const double hK = Cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);
    double *Coefficients[N_Points];
    double *aux = new double [N_Points*20]; // do not change below 20
    
    for(int j=0;j<N_Points;j++)
        Coefficients[j] = aux + j*20;
    
    if(Coeffs)
        Coeffs(N_Points, X, Y, nullptr, Coefficients);
    
    if(Manipulate)
        Manipulate(N_Points, Coefficients, nullptr, Cell);
    for(int j=0;j<N_Terms;j++)
    {
        AllOrigValues[j] =
        TFEDatabase2D::GetOrigElementValues(BaseFuncts[FESpaceNumber[j]],
                                            Derivatives[j]);
    }
    
    for(int i=0;i<N_Points;i++)
    {
        double Mult = weights[i]*AbsDetjk[i];
        Coefficients[i][19] = AbsDetjk[i];
        
        for(int j=0;j<N_Terms;j++) {
            OrigValues[j] = AllOrigValues[j][i];
        }
        
        for(auto lar : local_assemblings_routines)
        {
          lar(Mult, Coefficients[i], nullptr, hK, OrigValues, N_BaseFuncts,
              LocMatrix, LocRhs);
        }
    } // endfor i
    delete [] aux;
}

//==============================================================================
void LocalAssembling2D::GetParameters(int n_points,
                                      TCollection *Coll,
                                      TBaseCell *cell,
                                      int cellnum,
                                      double *x,
                                      double *y,
                                      double **Parameters)
{
  int *N_BaseFunct = new int[N_FEValues];
  double **Values = new double* [N_FEValues];
  double ***orig_values = new double** [N_FEValues];
  int **Index = new int* [N_FEValues];
  double Temp[2 + N_FEValues];
  // collect information
  for(int j=0; j<this->N_FEValues; j++)
  {
    TFEFunction2D *fefunction = this->FEFunctions2D[this->FEValue_FctIndex[j]];
    
    Values[j] = fefunction->GetValues();
    
    const TFESpace2D *fespace = fefunction->GetFESpace2D();
    FE2D FE_Id = fespace->GetFE2D(cellnum, cell);
    BaseFunct2D BaseFunct_Id = TFEDatabase2D::GetFE2D(FE_Id)->GetBaseFunct2D_ID();

    N_BaseFunct[j]=TFEDatabase2D::GetBaseFunct2D(BaseFunct_Id)->GetDimension();
    
    orig_values[j] = TFEDatabase2D::GetOrigElementValues(BaseFunct_Id, 
                                                         FEValue_MultiIndex[j]);
    Index[j] = fespace->GetGlobalDOF(cellnum);
  } // endfor j

  // loop over all quadrature points
  if(N_ParamFct != 0)
  {
    for(int i=0; i<n_points; ++i)
    {
      double *param = Parameters[i];

      Temp[0] = x[i];
      Temp[1] = y[i];

      // loop to calculate all FE values
      for(int k=2,j=0; j<N_FEValues; j++,k++)
      {
        double s = 0;
        int n = N_BaseFunct[j];
        double  *CurrValues = Values[j];
        double  *CurrOrigValues = orig_values[j][i];
        int *CurrIndex = Index[j];
        for(int l=0;l<n;l++)
          s += CurrValues[CurrIndex[l]]*CurrOrigValues[l];
        Temp[k] = s;
      }  // endfor j

      // loop to calculate all parameters
      for(int j=0; j<N_ParamFct; j++)
      {
        double *currparam = param + this->BeginParameter[j];
        this->ParameterFct[j](Temp, currparam);
      } // endfor j
    } // endfor i
  }
  
  delete [] N_BaseFunct;
  delete [] Values;
  delete [] orig_values;
  delete [] Index;
}

//-------------------------------------------------------------------------------
void LocalAssembling2D::compute_parameters(int n_points,
    TCollection *Coll,
    TBaseCell *cell,
    int cellnum,
    double *x,
    double *y)
{
  //std::vector<std::vector<double>> parameter_functions_values;
  this->parameter_functions_values.resize(n_points);

  for (unsigned int i = 0; i < parameter_functions_values.size(); i++)
  {
    parameter_functions_values[i].resize(this->N_FEValues);
  }

  ///@todo check the case N_ParamFct > 1
  if (N_ParamFct > 1)
  {
    Output::print(" ** WARNING: the function should not work for N_ParamFct > 1");
  }

  int *N_BaseFunct = new int[this->N_FEValues];
  double **Values = new double* [this->N_FEValues];
  double ***orig_values = new double** [this->N_FEValues]; // create an array of doublepointers (pointers to pointers) to doubles

  int **Index = new int* [this->N_FEValues];
  double Temp[this->N_FEValues + 2];

  // collect information
  for(int j = 0; j < this->N_FEValues; j++)
  {
   TFEFunction2D *fefunction = this->FEFunctions2D[this->FEValue_FctIndex[j]];

    // get all values (and derivatives) of basis function at given Gauss point
    Values[j] = fefunction->GetValues();
    const TFESpace2D *fespace = fefunction->GetFESpace2D();
    FE2D FE_Id = fespace->GetFE2D(cellnum, cell);
    BaseFunct2D BaseFunct_Id = TFEDatabase2D::GetFE2D(FE_Id)->GetBaseFunct2D_ID();
    
    N_BaseFunct[j]=TFEDatabase2D::GetBaseFunct2D(BaseFunct_Id)->GetDimension();

    orig_values[j] = TFEDatabase2D::GetOrigElementValues(BaseFunct_Id, FEValue_MultiIndex[j]);

    Index[j] = fespace->GetGlobalDOF(cellnum);
  } // endfor j

  // loop over all quadrature points
  if(N_ParamFct != 0)
  {
    for(int i = 0; i < n_points; i++)
    {
      // all values at quadrature point i
      //double *param = Parameters[i];
      Temp[0] = x[i];
      Temp[1] = y[i];
      // loop to calculate all FE values
      for(int j = 0; j < N_FEValues; j++)
      {
        double s = 0;
        int n = N_BaseFunct[j];
        double  *CurrValues = Values[j];
        double  *CurrOrigValues = orig_values[j][i];
        int *CurrIndex = Index[j];
        for(int l = 0; l < n; l++)
          s += CurrValues[CurrIndex[l]] * CurrOrigValues[l];
        Temp[2+j] = s;
      }  // endfor j

      // loop to calculate all parameters
      for(int j = 0; j < N_ParamFct; j++)
      {
        double *currparam = new double[N_FEValues];// = param + this->BeginParameter[j];
        this->ParameterFct[j](Temp, currparam);

        for (unsigned int l = 0; l < parameter_functions_values[i].size(); l++)
        {
          parameter_functions_values[i][l] = currparam[l];
        } 
      }// endfor j
    } // endfor i
  }
  delete [] N_BaseFunct;
  delete [] Values;
  delete [] orig_values;
  delete [] Index;
}

void LocalAssembling2D::set_parameters_for_nseGalerkin(LocalAssembling2D_type type)
{
  int nstype = TDatabase::ParamDB->NSTYPE;
  bool laplace_type_deformation = (TDatabase::ParamDB->LAPLACETYPE == 1);
  //bool laplace_type_deformation = this->db["laplace_type_deformation"];
  if(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==1)
  {
    ErrThrow("Newton method is not supported yet");
  }
  if(TDatabase::ParamDB->LAPLACETYPE == 1)
  {
    if((nstype==1) || nstype==2)
    {
      ErrThrow("LAPLACETYPE is only supported for NSTYPE 3, and 4");
    }
  }
  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM>0)
  {
    ErrThrow("Skew symmetric case is not implemented for all NSTYPE");
  }
  // common for all NSTYPE, Discrete forms, etc
  this->N_Rhs = 3;
  this->RhsSpace = { 0, 0, 1 };
  
  this->N_Terms = 4;
  this->Derivatives = { D00, D00, D10, D01};
  this->FESpaceNumber = { 0, 1, 0, 0 }; // 0: velocity, 1: pressure
  
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->Manipulate = nullptr;
  
  this->N_Parameters = 2;
  this->N_ParamFct = 1;
  this->ParameterFct =  { NSParamsVelo };
  this->N_FEValues = 2;
  this->FEValue_FctIndex = { 0, 1 };
  this->FEValue_MultiIndex = { D00, D00 };
  this->BeginParameter = { 0 };
  
  this->N_Matrices = 9;
  this->RowSpace = { 0, 0, 0, 0, 1, 1, 1, 0, 0 };
  this->ColumnSpace = { 0, 0, 0, 0, 1, 0, 0, 1, 1 };
  
  if(laplace_type_deformation)
  {
    this->local_assemblings_routines.push_back(NSLaplaceDeformation);
  }
  else
  {
    if(nstype == 1 || nstype == 2)
    {
      this->local_assemblings_routines.push_back(NSLaplaceGradGradSingle);
    }
    else
    {
      this->local_assemblings_routines.push_back(NSLaplaceGradGrad);
    }
  }
  switch(type)
  {
    case NSE2D_All:
      this->local_assemblings_routines.push_back(NSDivergenceBlocks);
      this->local_assemblings_routines.push_back(NSRightHandSide);
      if(nstype == 2 || nstype == 4 || nstype == 14)
      {
        this->local_assemblings_routines.push_back(NSGradientBlocks);
      }
      break; // NSE2D_ALL
    case NSE2D_NL:
      if(nstype == 1 || nstype == 2)
      {
        this->local_assemblings_routines.push_back(NSNonlinearTermSingle);
      }
      else
      {
        this->local_assemblings_routines.push_back(NSNonlinearTerm);
      }
      break; // NSE2D_NL
    default:
      ErrThrow("unknown LocalAssembling2D_type ", type, "  ", this->name);
      break;
  }
}

//=========================================================================
void LocalAssembling2D::set_parameters_for_tnseGalerkin(LocalAssembling2D_type type)
{
  int nstype = TDatabase::ParamDB->NSTYPE;
  // few checks
  if(TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==1)
  {
    ErrMsg("Newton method is not supported yet");
    exit(1);
  }
  if(TDatabase::ParamDB->LAPLACETYPE == 1)
  {
    if((nstype==1) || nstype==2)
    {
      ErrMsg("LAPLACETYPE is only supported for NSTYPE 3, and 4");
      exit(1);
    }
  }
  
  if(TDatabase::ParamDB->NSE_NONLINEAR_FORM>0)
  {
    ErrMsg("Skew symmetric case is not implemented for all NSTYPE");
    exit(1);
  }
  
  // common for all NSTYPE, Discrete forms, etc
  this->N_Rhs = 3;
  this->RhsSpace = { 0, 0, 1 };
  this->N_Terms = 4;
  this->Derivatives = { D00, D00, D10, D01};
  this->FESpaceNumber = { 0, 1, 0, 0 }; // 0: velocity, 1: pressure
  
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->Manipulate = nullptr;
  
  this->N_Parameters = 2;
  this->N_ParamFct = 1;
  this->ParameterFct = { NSParamsVelo };
  this->N_FEValues = 2;
  this->BeginParameter = { 0 };
  this->FEValue_MultiIndex = { D00, D00 };
  this->FEValue_FctIndex = { 0, 1 };    
  
  switch(type)
  {
    case TNSE2D:
      switch(nstype)
      {
        case 1:
          this->N_Matrices    = 4;
          this->RowSpace      = { 0, 0, 1, 1 };
          this->ColumnSpace   = { 0, 0, 0, 0 };
          this->local_assemblings_routines.push_back(TimeNSType1Galerkin);
          break; // nstype 1
        case 2:
          this->N_Matrices    = 6;
          this->RowSpace      = { 0, 0, 1, 1, 0, 0 };
          this->ColumnSpace   = { 0, 0, 0, 0, 1, 1 };
          this->local_assemblings_routines.push_back(TimeNSType2Galerkin);
          break; // nstype 2
        case 3:
          this->N_Matrices    = 8;
          this->RowSpace      = { 0, 0, 0, 0, 0, 0, 1, 1 };
          this->ColumnSpace   = { 0, 0, 0, 0, 0, 0, 0, 0 };
          if(TDatabase::ParamDB->LAPLACETYPE == 0)
            this->local_assemblings_routines.push_back(TimeNSType3Galerkin);
          else
            this->local_assemblings_routines.push_back(TimeNSType3GalerkinDD);
          break; // nstype 3
        case 4:
          this->N_Matrices    = 10;
          this->RowSpace      = { 0, 0, 0, 0, 0, 0, 1, 1, 0, 0 };
          this->ColumnSpace   = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 };
          if(TDatabase::ParamDB->LAPLACETYPE == 0)
            this->local_assemblings_routines.push_back(TimeNSType4Galerkin);
          else
            this->local_assemblings_routines.push_back(TimeNSType4GalerkinDD);
	  break;
	case 14:
	  this->N_Matrices    = 11;
          this->RowSpace      = { 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0 };
          this->ColumnSpace   = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1 };
          if(TDatabase::ParamDB->LAPLACETYPE == 0)
            this->local_assemblings_routines.push_back(TimeNSType14LPS);
          else
            this->local_assemblings_routines.push_back(TimeNSType14LPSDD);
          break; // nstype 4
      }// switch nstype
      break;// case TNSE2D 
    //==============================
    case TNSE2D_NL:
      switch(nstype)
      {
        case 1:
        case 2:
          this->N_Matrices    = 1;
          this->RowSpace      = { 0 };
          this->ColumnSpace   = { 0 };
          this->local_assemblings_routines.push_back(TimeNSType1_2NLGalerkin);
          break; // nstype 1, 2
        case 3:
        case 4:
        case 14:
          this->N_Matrices    = 2;
          this->RowSpace      = { 0, 0 };
          this->ColumnSpace   = { 0, 0 };
          if(TDatabase::ParamDB->LAPLACETYPE==0)
            this->local_assemblings_routines.push_back(TimeNSType3_4NLGalerkin);
          else
            this->local_assemblings_routines.push_back(TimeNSType3_4NLGalerkinDD);
          break; // nstype 3, 4
      }// switch nstype
      break; // case TNSE2D_NL
    //==============================
    case TNSE2D_Rhs:
      this->N_Matrices = 0;
      this->RowSpace = {};
      this->ColumnSpace = { };
      this->local_assemblings_routines.push_back(TimeNSRHS);
      this->Manipulate = nullptr;
      break; // case TNSE2D_Rhs
   //==============================
    default:
      ErrThrow("That's the wrong LocalAssembling2D_type ", type, " to come here.");
  }
}


void LocalAssembling2D::set_parameters_for_tnseSUPG(LocalAssembling2D_type type)
{
  if(TDatabase::ParamDB->NSTYPE < 4 )
  { 
    ErrThrow("SUPG method is only supported for NSTYPE 4 and 14 ", TDatabase::ParamDB->NSTYPE);
  }
  if(TDatabase::ParamDB->NSTYPE == 4)
  {
    this->N_Parameters = 2;
    this->N_ParamFct = 1;
    this->ParameterFct = {TimeNSParams2};
    this->N_FEValues = 2;
    this->BeginParameter = { 0 };
    this->FEValue_MultiIndex = { D00, D00 };
    this->FEValue_FctIndex = { 0, 1 };    
  }
  if(TDatabase::ParamDB->NSTYPE == 14)
  {
    this->N_Parameters = 4;
    this->N_ParamFct = 1;
    this->ParameterFct = {TimeNSType14ParamsSUPG};
    this->N_FEValues = 4;
    this->BeginParameter = { 0 };
    this->FEValue_MultiIndex = { D00, D00, D00, D00 };
    this->FEValue_FctIndex = { 0, 1, 2, 3 };
  }
  
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = true;
  this->Needs2ndDerivatives[1] = true;
  
  switch(type)
  {
    case TNSE2D:
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 4:
          this->N_Terms = 8; // Derivatives
                          // u_x, u_y, u, p, p_x, p_y, u_xx, u_yy
          this->Derivatives = { D10, D01, D00, D00, D10, D01, D20, D02};
          // 0: velocity space, 1: pressure space
          this->FESpaceNumber = { 0, 0, 0, 1, 1, 1, 0, 0 };
          // total number of matrices 
          this->N_Matrices = 10;
          // in the lower right corner
          this->RowSpace =    { 0, 0, 0, 0, 0, 0, 1, 1, 0, 0};
          this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1};
          this->N_Rhs = 2; // f1, f2, 
          this->RhsSpace = { 0, 0 };
          this->AssembleParam = TimeNSType4SUPG;
          this->Manipulate = NULL;
          break;
        case 14:
          this->N_Terms = 8; // Derivatives
                          // u_x, u_y, u, p, p_x, p_y, u_xx, u_yy
          this->Derivatives = { D10, D01, D00, D00, D10, D01, D20, D02};
          // 0: velocity space, 1: pressure space
          this->FESpaceNumber = { 0, 0, 0, 1, 1, 1, 0, 0 };
          // total number of matrices 
          this->N_Matrices = 11;
          // in the lower right corner
          this->RowSpace =    { 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0};
          this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1};
          this->N_Rhs = 3; // f1, f2, f3 (pressure part)
          this->RhsSpace = { 0, 0, 1 };
          this->AssembleParam = TimeNSType14SUPG;
          this->Manipulate = NULL;
         break;
      }
      break;
   //==============================
    case TNSE2D_NL:
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 4:
          this->N_Terms = 8; // Derivatives
                          // u_x, u_y, u, p, p_x, p_y, u_xx, u_yy
          this->Derivatives = { D10, D01, D00, D00, D10, D01, D20, D02};
          // 0: velocity space, 1: pressure space
          this->FESpaceNumber = { 0, 0, 0, 1, 1, 1, 0, 0 };
          // total number of matrices 
          this->N_Matrices = 8;
          // in the lower right corner
          this->RowSpace =    { 0, 0, 0, 0, 0, 0, 0, 0};
          this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 1, 1};
          this->N_Rhs = 2; // only stabilization terms 
          this->RhsSpace = {0, 0 };
          this->AssembleParam = TimeNSType4NLSUPG;
          this->Manipulate = NULL;
         break;
        case 14:
          this->N_Terms = 8; // Derivatives
                          // u_x, u_y, u, p, p_x, p_y, u_xx, u_yy
          this->Derivatives = { D10, D01, D00, D00, D10, D01, D20, D02};
          // 0: velocity space, 1: pressure space
          this->FESpaceNumber = { 0, 0, 0, 1, 1, 1, 0, 0 };
          // total number of matrices 
          this->N_Matrices = 11;
          // in the lower right corner
          this->RowSpace =    { 0, 0, 0, 0, // A-Blocks
	                        0, 0, // M-Blocks
				1, // C-Block
				1, 1, // B-Blocks
				0, 0 // BT-Blocks
				};
          this->ColumnSpace = { 0, 0, 0, 0, // A-Blocks
	                        0, 0, // M-Blocks
				1, // C Block
				0, 0, // B-Blocks
				1, 1  // BT-Blocks
				};
          this->N_Rhs = 3; // only stabilization terms 
          this->RhsSpace = {0, 0, 1 };
          this->AssembleParam = TimeNSType14NLSUPG;
          this->Manipulate = NULL;
          break;
      }
      break;
     //==============================
    case TNSE2D_Rhs:
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 4:
          this->N_Terms = 3;
          this->Derivatives = { D10, D01, D00};
          this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
          this->N_Matrices = 0;
          this->RowSpace = { };
          this->ColumnSpace = { };
          this->N_Rhs = 2 ;
          this->RhsSpace = {0, 0 };
          this->AssembleParam =TimeNSType4RHSSUPG;
          this->Manipulate = NULL;
          break;
        case 14:
          this->N_Terms = 5;
          this->Derivatives = { D10, D01, D00, D10, D01};
          this->FESpaceNumber = { 0, 0, 0, 1, 1 }; // 0: velocity, 1: pressure
          this->N_Matrices = 0;
          this->RowSpace = {};
          this->ColumnSpace = { };
          this->N_Rhs = 3 ;
          this->RhsSpace = {0, 0, 1};
          this->AssembleParam =TimeNSType14RHSSUPG; 
          this->Manipulate = NULL;
        break;
      }
      break;
     //==============================
    default:
      ErrThrow("That's the wrong LocalAssembling2D_type ", type, " to come here.");
  }
}

void LocalAssembling2D::set_parameters_for_tnseSmagorinsky(LocalAssembling2D_type type)
{
  this->N_Parameters = 6;
  this->N_ParamFct = 1;
  this->ParameterFct =  { TimeNSParamsVelo_GradVelo };
  this->N_FEValues = 6;
  this->FEValue_FctIndex = { 0, 1, 0, 1, 0, 1 };
  // u1old, u2old, u3old, all derivatives of u1, u2, u3
  this->FEValue_MultiIndex = { D00, D00, 
                               D10, D10,
                               D01, D01};     
  this->BeginParameter = { 0 };
  // common for all NSTYPE, Discrete forms, etc
  if(type==TNSE2D)
  {
    this->N_Terms = 4;
    this->Derivatives = { D10, D01, D00, D00 };  
    this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure    
    this->N_Rhs = 3; // NOTE: check why is this always three??
    this->RhsSpace = { 0, 0, 0 };
  }
  else if(type==TNSE2D_NL)
  {
    this->N_Terms = 3;
    this->Derivatives = { D10, D01, D00 };    
    this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure        
    this->N_Rhs = 0;
    this->RhsSpace = {};
  }
  
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = false;
  this->Needs2ndDerivatives[1] = false;
  this->Manipulate = NULL;
  int nstype = TDatabase::ParamDB->NSTYPE;
  switch(type)
  {
    case TNSE2D:
      switch(nstype)
      {
        case 1:
          this->N_Matrices    = 4;
          this->RowSpace      = { 0, 0, 1, 1 };
          this->ColumnSpace   = { 0, 0, 0, 0 };
          this->AssembleParam = TimeNSType1Smagorinsky;
          break; // nstype 1
        case 2:
          this->N_Matrices    = 6;
          this->RowSpace      = { 0, 0, 1, 1, 0, 0 };
          this->ColumnSpace   = { 0, 0, 0, 0, 1, 1 };
          this->AssembleParam = TimeNSType2Smagorinsky;
          break; // nstype 2
        case 3:
          this->N_Matrices    = 8;
          this->RowSpace      = { 0, 0, 0, 0, 0, 0, 1, 1 };
          this->ColumnSpace   = { 0, 0, 0, 0, 0, 0, 0, 0 };
          if(TDatabase::ParamDB->LAPLACETYPE == 0)
            this->AssembleParam = TimeNSType3Smagorinsky;
          else
            this->AssembleParam = TimeNSType3SmagorinskyDD;
          break; // nstype 3
        case 4:
          this->N_Matrices    = 10;
          this->RowSpace      = { 0, 0, 0, 0, 0, 0, 1, 1, 0, 0 };
          this->ColumnSpace   = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 };
          if(TDatabase::ParamDB->LAPLACETYPE == 0)
            this->AssembleParam = TimeNSType4Smagorinsky;
          else
            this->AssembleParam = TimeNSType4SmagorinskyDD;
          break; // nstype 4
      }// switch nstype
      break;// case TNSE2D 
  //==============================
    case TNSE2D_NL:
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          this->N_Matrices    = 1;
          this->RowSpace      = { 0 };
          this->ColumnSpace   = { 0 };
          this->AssembleParam = TimeNSType1_2NLSmagorinsky;
          break; // nstype 1, 2
        case 3:
        case 4:
          if(TDatabase::ParamDB->LAPLACETYPE==0)
	  {
	    this->N_Matrices    = 4;
	    this->RowSpace      = { 0, 0, 0, 0 };
	    this->ColumnSpace   = { 0, 0, 0, 0 };
	    this->AssembleParam = TimeNSType3_4NLSmagorinsky;
	  }
          else
	  {
	    this->N_Matrices    = 4;
	    this->RowSpace      = { 0, 0, 0, 0 };
	    this->ColumnSpace   = { 0, 0, 0, 0 };
	    this->AssembleParam = TimeNSType3_4NLSmagorinskyDD;
	  }
          break; // nstype 3, 4
      }
      break;
  //==============================
    case TNSE2D_Rhs:
      this->N_Terms = 1;
      this->Derivatives = { D00 };      
      this->FESpaceNumber = { 0 }; // 0: velocity, 1: pressure
      this->N_Matrices = 0;
      this->RowSpace = {};
      this->ColumnSpace = { };
      this->N_Rhs = 3;
      this->RhsSpace = {0, 0, 0};
      this->AssembleParam =TimeNSRHS;
      this->Manipulate = NULL;
     break;
  //==============================
    default:
      ErrThrow("That's the wrong LocalAssembling2D_type ", type, " to come here.");
  }
}

void LocalAssembling2D::set_parameters_for_tnseRBVMS(LocalAssembling2D_type type)
{
  
  if(TDatabase::ParamDB->NSTYPE < 4 )
  { 
    ErrThrow("Residual Based VMS method is only supported for NSTYPE 4 and 14 ", 
             TDatabase::ParamDB->NSTYPE);
  }
  this->N_Parameters = 17;
  this->N_ParamFct = 1;
  this->ParameterFct = {TimeNSParams_Residual_VMS};
  this->N_FEValues = 17;
  this->BeginParameter = { 0 };
  this->FEValue_MultiIndex = { D00, D00, // u1sigma, u2sigma
                               D00, D00, // combination of previous 2 solutions
                               D10, D10, D01, D01, // u1sigma_x, u2sigma_x, u1sigma_y, u2sigma_y
                               D20, D20, D02, D02,// u1sigma_xx, u2sigma_xx, u1sigma_yy, u2sigma_yy
                               D10, D01, D00, // pold_x, pold_y, pold
                               D00, D00 // time derivative
  };
  this->FEValue_FctIndex = { 0, 1,        //  
                             2, 3,        // 
                             0, 1, 0, 1,  //
                             0, 1, 0, 1,  //
                             4, 4, 4,     //
                             5, 6         //
  };

  
  this->N_Terms = 8;
  this->Derivatives = { D10, D01, D00, D00, D10, D01, D20, D02};
  this->FESpaceNumber = { 0, 0, 0, 1, 1, 1, 0, 0 };  
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives = new bool[2];
  this->Needs2ndDerivatives[0] = true;
  this->Needs2ndDerivatives[1] = true;
  
  switch(type){
    case TNSE2D:
      switch(TDatabase::ParamDB->NSTYPE){
        case 4:
          this->N_Matrices = 12;
          this->RowSpace =    { 0, 0, 0, 0, // A-Blocks
	                        0, 0, 0, 0, //M-Blocks
				1, 1, //B-Blocks
				0, 0 // BT-Blocks
				};
          this->ColumnSpace = { 0, 0, 0, 0, //A-Blocks
                                0, 0, 0, 0, //M-Blocks
				0, 0, //B-Blocks
				1, 1 //BT-Blocks
				};
          this->N_Rhs = 2; 
          this->RhsSpace = { 0, 0 };
          this->AssembleParam = TimeNSType4Residual_VMS;
          this->Manipulate = NULL;
          break;
        case 14:
          this->N_Matrices = 13;
          this->RowSpace =    { 0, 0, 0, 0, // A-Blocks
	                        0, 0, 0, 0, //M-Blocks
				1, // C-Block
				1, 1, // B-Blocks
				0, 0 //BT-Blocks
				};
          this->ColumnSpace = { 0, 0, 0, 0, // A-Blocks
	                        0, 0, 0, 0, //M-Blocks
				1, // C-Block
				0, 0, // B-Blocks
				1, 1 // BT-Blocks
				};
          this->N_Rhs = 3; 
          this->RhsSpace = { 0, 0, 1 };
          this->AssembleParam = TimeNSType14Residual_VMS;
          this->Manipulate = NULL;
          break;
      }// switch nstype     
      break; // case TNSE2D
//--------------------------
    case TNSE2D_NL:
      switch(TDatabase::ParamDB->NSTYPE){
        case 4:
          this->N_Matrices = 10;
          this->RowSpace =    { 0, 0, 0, 0, // A-Blocks
	                        0, 0, 0, 0, //M-Blocks
				0, 0
				};
          this->ColumnSpace = { 0, 0, 0, 0, // A-Blocks
	                        0, 0, 0, 0, //M-Blocks
				1, 1
				};
          this->N_Rhs = 2; 
          this->RhsSpace = { 0, 0 };
          this->AssembleParam = TimeNSType4NLResidual_VMS;
          this->Manipulate = NULL;
          break;
        case 14:
          this->N_Matrices = 13;
          this->RowSpace =    { 0, 0, 0, 0, // A-Blocks
	                        0, 0, 0, 0, //M-Blocks
				1, // C-Block
				1, 1, // B-Blocks
				0, 0 // BT-Blocks
				};
          this->ColumnSpace = { 0, 0, 0, 0, // A-Blocks
	                        0, 0, 0, 0, //M-Blocks
				1, // C-Block
				0, 0, // B-Blocks
				1, 1 // BT-Blocks
				};
          this->N_Rhs = 3; 
          this->RhsSpace = { 0, 0, 1 };
          this->AssembleParam = TimeNSType14Residual_VMS;
          this->Manipulate = NULL;
          break;
      }// switch nstype   
      break;// case TNSE2D_NL
//---------------------------      
    case TNSE2D_Rhs:
      switch(TDatabase::ParamDB->NSTYPE){
        case 4:
          
          this->N_Terms = 3;
          this->Derivatives = { D10, D01, D00};
          this->Needs2ndDerivatives = new bool[1];
          this->Needs2ndDerivatives[0] = false;
          this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
          this->N_Matrices = 0;
          this->RowSpace = { };
          this->ColumnSpace = { };
          this->N_Rhs = 2 ;
          this->RhsSpace = {0, 0 };
          this->AssembleParam = TimeNSType4RHS_Residual_VMS;
          this->Manipulate = NULL;
          break;
        case 14:
          this->N_Terms = 5;
          this->Derivatives = { D10, D01, D00, D10, D01};
          this->Needs2ndDerivatives = new bool[1];
          this->Needs2ndDerivatives[0] = false;
          this->FESpaceNumber = { 0, 0, 0, 1, 1 }; // 0: velocity, 1: pressure
          this->N_Matrices = 0;
          this->RowSpace = {};
          this->ColumnSpace = { };
          this->N_Rhs = 3 ;
          this->RhsSpace = {0, 0, 1};
          this->AssembleParam = TimeNSType14RHS_Residual_VMS;
          this->Manipulate = NULL;
          break;
      }// switch nstype   
      break; // case TNSE2D_Rhs
  }
}
