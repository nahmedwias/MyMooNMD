// =======================================================================
// @(#)AuxParam2D.C        1.2 09/17/99
// 
// Class:       TAuxParam2D
// Purpose:     store parameter functions and FE functions
//
// Author:      Gunar Matthies (06.08.98)
//
// History:     start of implementation 06.08.98 (Gunar Matthies)
//
// =======================================================================

#include <AuxParam2D.h>
#include <FEDatabase2D.h>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>

TAuxParam2D::TAuxParam2D(int n_paramfct, int n_fevalues,
                         TFEFunction2D **fefunctions2d, ParamFct **parameterfct,
                         int *fevalue_fctindex,
                         MultiIndex2D *fevalue_multiindex, int n_parameters,
                         int *beginparameter)
{
  N_ParamFct = n_paramfct;
  N_FEValues = n_fevalues;

  FEFunctions2D = fefunctions2d;
  ParameterFct = parameterfct;

  FEValue_FctIndex = fevalue_fctindex;
  FEValue_MultiIndex = fevalue_multiindex;

  N_Parameters = n_parameters;
  BeginParameter = beginparameter;

  Temp = new double[2 + N_FEValues];

  Values = new const double* [N_FEValues];
  OrigValues = new double** [N_FEValues];
  Index = new int* [N_FEValues];
  N_BaseFunct = new int[N_FEValues];
}

TAuxParam2D::TAuxParam2D() 
 : TAuxParam2D(0, 0, nullptr, nullptr, nullptr, nullptr, 0, nullptr)
{
}



/** return all parameters at all quadrature points */
void TAuxParam2D::GetParameters(int N_Points, TCollection *, TBaseCell *cell,
                                int cellnum, const double *, const double *,
                                double *X, double *Y, double **Parameters)
{
  // collect information
  for(int j=0;j<N_FEValues;j++)
  {
    const TFEFunction2D *fefunction = FEFunctions2D[FEValue_FctIndex[j]];
    Values[j] = fefunction->GetValues();
    //  if (N_FEValues==8)
    //  OutPut("aac " << (int) fefunction << " " <<  Values[j][0]<< endl);

    auto fespace = fefunction->GetFESpace2D();
    FE2D FE_Id = fespace->GetFE2D(cellnum, cell);
    BaseFunct2D BaseFunct_Id = TFEDatabase2D::GetFE2D(FE_Id)->GetBaseFunct2D_ID();

    N_BaseFunct[j]=TFEDatabase2D::GetBaseFunct2D(BaseFunct_Id)->GetDimension();
    
    OrigValues[j] = TFEDatabase2D::GetOrigElementValues(BaseFunct_Id, FEValue_MultiIndex[j]);
    Index[j] = fespace->GetGlobalDOF(cellnum);
  } // endfor j

  // loop over all quadrature points
  for(int i=0;i<N_Points;i++)
  {
    double * param = Parameters[i];

    Temp[0] = X[i];
    Temp[1] = Y[i];

    // loop to calculate all FE values
    for(int k=2,j=0;j<N_FEValues;j++,k++)
    {
      double s = 0;
      int n = N_BaseFunct[j];
      const double *CurrValues = Values[j];
      const double *CurrOrigValues = OrigValues[j][i];
      const int *CurrIndex = Index[j];
      for(int l=0;l<n;l++)
        s += CurrValues[CurrIndex[l]]*CurrOrigValues[l];
      Temp[k] = s;
    }  // endfor j

    // loop to calculate all parameters
    for(int j=0;j<N_ParamFct;j++)
    {
      double * currparam = param + BeginParameter[j];
      ParameterFct[j](Temp, currparam);
    } // endfor j
  } // endfor i
}

/** return all parameters at all quadrature points on the boundary*/
void TAuxParam2D::GetParameters(int N_Points, TCollection *Coll,
                                TBaseCell *cell, int cellnum, const double *t,
                                int joint, double **Parameters)
{
  int i,j,k,l,n; //N_Cells;
  double xi, eta; // eps = 1e-20;
  double s;
  double *param, *currparam, *CurrOrigValues;
  int *CurrIndex;
  TFEFunction2D *fefunction;
  BaseFunct2D BaseFunct_Id;
  FE2D FE_ID = C_P00_2D_T_A; //avoid uninit warning
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf, **AllBaseFuncts;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D], uref[MaxN_BaseFunctions2D];
  double uxiref[MaxN_BaseFunctions2D], uetaref[MaxN_BaseFunctions2D];
  
//  int *Numbers;
//  double u, ux, uy;
//  double val;
  int *GlobalNumbers, *BeginIndex;

  double X, Y, absdetjk;

  AllBaseFuncts = new TBaseFunct2D*[N_FEValues];

  for(j=0;j<N_FEValues;j++)
  {
    fefunction = FEFunctions2D[FEValue_FctIndex[j]];
    
    Values[j] = fefunction->GetValues();

    auto fespace = fefunction->GetFESpace2D();
    FE_ID = fespace->GetFE2D(cellnum, cell);
    BaseFunct_Id = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D_ID();

    AllBaseFuncts[j] = TFEDatabase2D::GetBaseFunct2D(BaseFunct_Id);
    N_BaseFunct[j] = AllBaseFuncts[j]->GetDimension();
    
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    Index[j] = GlobalNumbers + BeginIndex[cellnum];
  } // endfor j

  FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
  RefTrans = FE_Obj->GetRefTransID();
  // set cell for reference transformation
  TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
  
  for(i=0;i<N_Points;i++)
  {
    switch(joint)
    {
    case 0: 
      xi=t[i];
      eta=-1;
      break; 
    case 1: 
      xi=1;
      eta=t[i];
      break;
    case 2: 
      xi=-t[i];
      eta=1;
      break;
    case 3: 
      xi=-1;
      eta=-t[i];
      break;
    }//switch
//     cout << "xi eta " << xi << " " << eta << endl;

    param = Parameters[i];

    TFEDatabase2D::GetOrigFromRef(RefTrans, 1, &xi, &eta, &X, &Y, &absdetjk);
    // Temp[0] = X[i];
    // Temp[1] = Y[i];

    // loop to calculate all FE values
    for(k=2,j=0;j<N_FEValues;j++,k++)
    {
      s = 0;
      n = N_BaseFunct[j];
      const double * CurrValues = Values[j];

      // get values and derivatives of basis functions on the
      // reference mesh cell
      bf = AllBaseFuncts[j];
      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);
  
      // compute values on the original mesh cell 
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
				   uref, uxiref, uetaref, 
				   uorig, uxorig, uyorig);
      switch(FEValue_MultiIndex[j])
      {
        case D00:
          CurrOrigValues = uorig;
	break;
        case D10:
          CurrOrigValues = uxorig;
	break;
        case D01:
          CurrOrigValues = uyorig;
	break;
        default:
         cerr << "Second derivatives not added, see AuxParam2D "  << endl;
          exit (-1);
         break;
      } // endswitch
      // cout << "CurrOrigValues" << *uorig << endl;
      // cout << "CurrOrigValuesx" << *uxorig << endl;
      CurrIndex = Index[j];
      for(l=0;l<n;l++)
      { 
	s += CurrValues[CurrIndex[l]]*CurrOrigValues[l];
	// cout << "Parameter   " << CurrOrigValues[l] << endl;
      }
      Temp[k] = s;
    }  // endfor j

    // loop to calculate all parameters
    for(j=0;j<N_ParamFct;j++)
    {
      currparam = param + BeginParameter[j];
      ParameterFct[j](Temp, currparam);
    } // endfor j
  } // endfor i 
  delete [] AllBaseFuncts;
}

/** destructor */
TAuxParam2D::~TAuxParam2D()
{
  delete [] Temp;
  delete [] Values;
  delete [] OrigValues;
  delete [] Index;
  delete [] N_BaseFunct;
}
