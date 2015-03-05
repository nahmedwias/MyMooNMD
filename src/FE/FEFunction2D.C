// =======================================================================
// @(#)FEFunction2D.C        1.6 04/13/00
//
// Class:       TFEFunction2D
// Purpose:     a function from a finite element space in 2D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
// =======================================================================

#include <Database.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <FEDatabase2D.h>
#include <FEFunction2D.h>
#include <string.h>
#include <AllRefTrans.h>
#include <NodalFunctional2D.h>
#include <MainUtilities.h>
#include <MooNMD_Io.h>

#include <InterfaceJoint.h>
#include <BdLine.h>
#include <LinAlg.h>

#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <MooNMD_Io.h>
// #include <malloc.h>
#include <dirent.h> 
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

void OnlyDirichlet(int i, double t, BoundCond &cond)
{
	cond = DIRICHLET;
}

/** constructor with vector initialization */
TFEFunction2D::TFEFunction2D(TFESpace2D *fespace2D, char *name,
char *description, double *values, int length)
{
  FESpace2D=fespace2D;

  Name=strdup(name);

  Description=strdup(description);

  Values=values;

  Length=length;
}


TFEFunction2D::~TFEFunction2D()
{
  delete Name;
  delete Description;
}


/** calculate errors to given function */
/** parallel with MPI, Sashi : 10.10.09 */
void TFEFunction2D::GetErrors(DoubleFunct2D *Exact, int N_Derivatives,
                              MultiIndex2D *NeededDerivatives,
                              int N_Errors, ErrorMethod2D *ErrorMeth, 
                              CoeffFct2D *Coeff, 
                              TAuxParam2D *Aux,
                              int n_fespaces, TFESpace2D **fespaces,
                              double *errors)
{
  int i,j,k,l,n,m, ij, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_, N_Edges;
  int Used[N_FEs2D], *N_BaseFunct;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  TFE2D *ele;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double *Derivatives[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK,xy_error,xy_error_max, x_error_max, y_error_max;
  bool *SecondDer;
  int out_xy_error;
  out_xy_error = 0;
  xy_error_max = 0;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  N_Parameters = Aux->GetN_Parameters();
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux + j*N_Parameters;

  aux = new double [MaxN_QuadPoints_2D*N_Derivatives];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux + j*N_Derivatives;

  aux = new double [MaxN_QuadPoints_2D * 4];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    ExactVal[j] = aux + j*4;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*20;

  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  BeginIndex = FESpace2D->GetBeginIndex();

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;
  errors[N_Errors] = 0.0;  // for L_infty-error

  // ########################################################################
  // loop over all cells
  // ########################################################################
  Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

#ifdef _MPI
      int ID, rank;
      MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
      ID  = cell->GetSubDomainNo();

      if(rank!=ID) // halo cells errors will not be calculated
       {
        continue; 
       }
#endif

    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
	case 0: // diameter
	    hK = cell->GetDiameter();
	    break;
	    //case 1: // with reference map
	    //OutPut("cell measure " << endl);
	    //hK = cell->GetLengthWithReferenceMap();
	    //break;
	case 2: // shortest edge
	    hK = cell->GetShortestEdge();
	    break;
	case 1: // with reference map
	case 3: // measure
 	    hK = cell->GetMeasure();
	    hK = sqrt(hK);
	    break;
	case 4: // mesh size in convection direction, this is just a dummy
	    hK = cell->GetDiameter();
	    break;
	case 5: // take value from an array
	    // this is in general not the diameter but a pw constant value
	    // which is needed for some reasons
	    hK = cell->GetDiameter();
	    break;
	default: // diameter
	    hK = cell->GetDiameter();
	    break;
    }

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    memset(Used, 0, N_FEs2D*SizeOfInt);
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i, cell);
      Used[CurrentElement] = 1;
    }

    N_LocalUsedElements = 0;
    memset(LocalUsedElements, 0, SizeOfInt*N_FEs2D);
    j = 0;
    for(k=0;k<N_FEs2D;k++)
      if(Used[k])
    {
      LocalUsedElements[j] = (FE2D)k;
      j++;
    }
    N_LocalUsedElements = j;
    

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);

    if(N_Parameters>0)
      Aux->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);

    if ((TDatabase::ParamDB->DISCTYPE == SDFEM)||
	(TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM))
    {
	TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
	N_Edges = cell->GetN_Edges();
	for (ij=0;ij<N_Edges;ij++)
	{
	    TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
	    TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
	}
	if (N_Edges==3)
	    TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
	TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
    }

    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace2D->GetFE2D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = Values[DOF[l]];

    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
        NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
        }                                         // endfor l
        Derivatives[j][k] = value;
      }                                           // endfor j
    }                                             // endfor k

    //errors[N_Errors] = 0.0;  //added D.Sirch
    for(j=0;j<N_Points;j++)
    {
      Exact(X[j], Y[j], ExactVal[j]);
       // D.Sirch: computation of L^\inf-error
      if(fabs(*ExactVal[j] - Derivatives[j][0]) > errors[N_Errors])
      {
	errors[N_Errors] = fabs(*ExactVal[j] - Derivatives[j][0]);
      }
    }

    /* Compute errors of the function in all quadrature points*/
    if (out_xy_error == 1)
    {
      for(j=0;j<N_Points;j++)
      {
        xy_error = fabs(Derivatives[j][0]-ExactVal[j][0]);
        //OutPut("xy_Error in: " << X[j] << " " << Y[j] << " is " << xy_error  << endl);
        OutPut(X[j] << " " <<  Y[j] << " " << xy_error << endl);
        if(xy_error > xy_error_max)
        {
          xy_error_max = xy_error;
          x_error_max = X[j];
          y_error_max = Y[j];
        }
      }
    }
    if(Coeff)
      Coeff(N_Points, X, Y, Param, AuxArray);
    // special rule for L1 error and P1 FE
#ifdef __2D__
    if ((ErrorMeth== &L1Error) && ( LocalUsedElements[0]==C_P1_2D_T_A))
    {
      // compute coordinates of vertices
      for (j=0;j<3; j++)
      {
        X[j] = cell->GetVertex(j)->GetX();
        Y[j] = cell->GetVertex(j)->GetY();
      }
      // set flag
      hK = -4711;
    }
#endif

    ErrorMeth(N_Points, X, Y, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

#ifdef __2D__
    if (!(ErrorMeth== &SDFEMErrors))
      {
       for(j=0;j<N_Errors;j++)
	  errors[j] += LocError[j];
      }
    else
      {
       for(j=0;j<N_Errors-1;j++)
        errors[j] += LocError[j];
        // L_infty error
        if (errors[N_Errors-1] <  LocError[N_Errors-1])
         errors[N_Errors-1] = LocError[N_Errors-1];
      }
#endif
  }                                               // endfor i


#ifndef _MPI // sqrt(errors[j]) in the main programm after collecting error from all subdomain
#ifdef __2D__
  if (!(ErrorMeth== &L1Error))
   {
    if (!(ErrorMeth== &SDFEMErrors))
      {
	for(j=0;j<N_Errors;j++)
	  errors[j] = sqrt(errors[j]);
      }
    else
     {
	for(j=0;j<N_Errors-1;j++)
	  errors[j] = sqrt(errors[j]);
     }
   }
#endif
#ifdef __3D__
 	for(j=0;j<N_Errors;j++)
	  errors[j] = sqrt(errors[j]);
#endif
#endif

  if (out_xy_error == 1)
    OutPut("Maximal xy_Error in: " <<  x_error_max << " " <<  y_error_max << " value: " << xy_error_max  << endl);

  errors[N_Errors] =  xy_error_max;

  delete AuxArray[0];
  delete [] SecondDer;
  delete ExactVal[0];
  delete Derivatives[0];
}                                                 // TFEFunction2D::GetErrors


#ifdef __2D__
/** calculate errors to given function */
// G. Matthies: 18.6.09
void TFEFunction2D::GetErrorsAdapt(DoubleFunct2D *Exact, int N_Derivatives,
			      MultiIndex2D *NeededDerivatives,
			      int N_Errors, ErrorMethod2D *ErrorMeth, 
			      CoeffFct2D *Coeff, 
			      TAuxParam2D *Aux,
			      int n_fespaces, TFESpace2D **fespaces,
			      double *errors)
{
	int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
	int N_Cells, N_Points, N_Parameters, N_;
	int Used[N_FEs2D], *N_BaseFunct;
	TFESpace2D *fespace;
	FE2D LocalUsedElements[N_FEs2D], CurrentElement;
	BaseFunct2D BaseFunct, *BaseFuncts;
	TCollection *Coll;
	TBaseCell *cell;
	TFE2D *ele;
	double *weights, *xi, *eta;
	double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
	double AbsDetjk[MaxN_QuadPoints_2D];
	RefTrans2D RefTrans;
	double *Param[MaxN_QuadPoints_2D], *aux;
	double *Derivatives[MaxN_QuadPoints_2D];
	double *ExactVal[MaxN_QuadPoints_2D];
	double *AuxArray[MaxN_QuadPoints_2D];
	int *DOF, ActiveBound, DirichletBound, end, last;
	double **OrigFEValues, *Orig, value;
	double FEFunctValues[MaxN_BaseFunctions2D];
	int *GlobalNumbers, *BeginIndex;
	double LocError[4];
	double hK;
	double *deriv, *exactval, w, t;
	bool *SecondDer;
	
	TGridCell *RootCell;
	TBaseCell **LocalCell;
	TCollection *LocalColl;
	FE2D *LocalFEs;
	TFESpace2D *RootSpace;
	TFEFunction2D *RootFct;
	double *RootVal;
	int RootLen;
	
	const int MaxLev = 5;
	int Lev;
	TCollection *FineColls[MaxLev];
	TFEFunction2D *FineFcts[MaxLev];
	TFESpace2D *FineSpaces[MaxLev];
	TCollection *CurrColl;
	TGridCell *CurrCell;
	int N_FineCells, N_CoarseCells, CellId;
	TBaseCell **FineCells;
	FE2D *FineFEs[MaxLev];
	TBdLine *Edges[4];
	double xv[4], yv[4];
	double *FineVals[MaxLev];
	double *AuxVals[MaxLev];
	int FineDOF;
	TAuxParam2D *FineAux;
	int *NewDOF, *NewGlobalNumbers, *NewBeginIndex;
	
	char Name[] = "name";
	char Description[] = "description";
	
	Edges[0] = new TBdLine(12345);
	Edges[1] = new TBdLine(12345);
	Edges[2] = new TBdLine(12345);
	Edges[3] = new TBdLine(12345);
  
	BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
	N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

	SecondDer = new bool[n_fespaces];
	for(i=0;i<n_fespaces;i++)
		SecondDer[i] = FALSE;

	N_Parameters = Aux->GetN_Parameters();
	aux = new double [MaxN_QuadPoints_2D*N_Parameters];
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		Param[j] = aux + j*N_Parameters;

	aux = new double [MaxN_QuadPoints_2D*N_Derivatives];
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		Derivatives[j] = aux + j*N_Derivatives;
  
	aux = new double [MaxN_QuadPoints_2D * 4];
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		ExactVal[j] = aux + j*4;

  // 20 <= number of term
	aux = new double [MaxN_QuadPoints_2D*20]; 
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		AuxArray[j] = aux + j*20;

	GlobalNumbers = FESpace2D->GetGlobalNumbers();
	BeginIndex = FESpace2D->GetBeginIndex();

	for(i=0;i<N_Errors;i++)
		errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
	Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
	N_Cells = Coll->GetN_Cells();
 
	for(i=0;i<N_Cells;i++)
	{
		cell = Coll->GetCell(i);
		
		
		hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
		N_LocalUsedElements = 1;
		LocalUsedElements[0] = fespaces[0]->GetFE2D(i, cell);
		
		CurrentElement = FESpace2D->GetFE2D(i, cell);
		BaseFunct = BaseFuncts[CurrentElement];
		N_ = N_BaseFunct[CurrentElement];
		
		DOF = GlobalNumbers + BeginIndex[i];
		
		RootCell = new TGridCell(cell->GetRefDesc(),0);
		j = cell->GetN_Vertices();
		for(k=0;k<j;k++)
			cell->GetVertex(k)->GetCoords(xv[k], yv[k]);
		
		for(k=0;k<j;k++)
		{
			RootCell->SetVertex(k, cell->GetVertex(k));
			RootCell->SetJoint(k, new TBoundEdge(Edges[k], 0, 1));
			Edges[k]->SetParams(xv[k], yv[k], xv[(k+1)%j]-xv[k],
					    yv[(k+1)%j]-yv[k]);
		}
		LocalCell = new TBaseCell*[1];
		LocalCell[0] = (TBaseCell *)RootCell;
		LocalColl = new TCollection(1, LocalCell);
		LocalFEs = new FE2D[1];
		LocalFEs[0] = LocalUsedElements[0];
		
		RootSpace = new TFESpace2D(LocalColl, Name, Description,
			               OnlyDirichlet,
			               LocalFEs, NULL);
		RootLen = N_;
		RootVal = new double[RootLen];
		NewDOF = RootSpace->GetGlobalNumbers() + RootSpace->GetBeginIndex()[0];
		for(l=0;l<N_;l++)
			RootVal[NewDOF[l]] = Values[DOF[l]];
		RootFct = new TFEFunction2D(RootSpace, Name, Description, 
					    RootVal, RootLen);
		FineAux = new TAuxParam2D(1, 0, 0, 0, &RootSpace, NULL, NULL, 
					  NULL, NULL, 0, NULL);
		RootFct->GetErrors(Exact, N_Derivatives, NeededDerivatives,
				   N_Errors, ErrorMeth, 
				   Coeff, FineAux, 1, &RootSpace, LocError);
		delete FineAux;

		FineColls[0] = LocalColl;
		FineFcts[0] = RootFct;
		FineSpaces[0] = RootSpace;
		FineFEs[0] = LocalFEs;
		FineVals[0] = RootVal;
   		
		for(Lev=1;Lev<MaxLev;Lev++)
		{
			CurrColl = FineColls[Lev-1];
			N_CoarseCells = CurrColl->GetN_Cells();
			N_FineCells = 4*N_CoarseCells;
			FineCells = new TBaseCell*[N_FineCells];
			CellId = 0;
			for(j=0;j<N_CoarseCells;j++)
			{
				CurrCell = (TGridCell*)(CurrColl->GetCell(j));
				CurrCell->SetRegRefine();
				CurrCell->Refine(Lev);
				for(k=0;k<4;k++)
				{
					FineCells[CellId] =
							(TBaseCell*)(CurrCell->GetChild(k));
					CellId++;
				}
			}
			FineColls[Lev] = new TCollection(N_FineCells, FineCells);
			FineFEs[Lev] = new FE2D[N_FineCells];
			for(j=0;j<N_FineCells;j++)
				FineFEs[Lev][j] = LocalUsedElements[0];
			FineSpaces[Lev] = new TFESpace2D(FineColls[Lev], Name, Description,
					OnlyDirichlet,
					FineFEs[Lev], NULL);
			FineDOF = FineSpaces[Lev]->GetN_DegreesOfFreedom();
			FineVals[Lev] = new double[FineDOF];
			AuxVals[Lev] = new double[FineDOF];
			FineFcts[Lev] = new TFEFunction2D(FineSpaces[Lev], 
					Name, Description, 
					FineVals[Lev], FineDOF);
			Prolongate(FineSpaces[Lev-1], FineSpaces[Lev],
				   FineVals[Lev-1], FineVals[Lev], AuxVals[Lev]);
			FineAux = new TAuxParam2D(1, 0, 0, 0, &(FineSpaces[Lev]), 
					NULL, NULL, NULL, NULL, 0, NULL);
			FineFcts[Lev]->GetErrors(Exact, N_Derivatives, NeededDerivatives,
					   N_Errors, ErrorMeth, 
					   Coeff, FineAux, 1, &(FineSpaces[Lev]), LocError);
			delete FineAux;
		} // end for Lev < MaxLev
		
		// delete all allocated memory
		for(Lev=MaxLev-1;Lev>0;Lev--)
		{
			N_CoarseCells = FineColls[Lev-1]->GetN_Cells();
			for(j=0;j<N_CoarseCells;j++)
			{
				FineColls[Lev-1]->GetCell(j)->Derefine();
			}
			delete AuxVals[Lev];
			delete FineVals[Lev];
			delete FineColls[Lev];
			delete FineFcts[Lev];
			delete FineSpaces[Lev];
		}
		delete RootFct;
		delete RootSpace;
		delete RootVal;
		delete LocalColl;
		delete RootCell;
		
		for(j=0;j<N_Errors;j++)
			errors[j] += LocError[j]*LocError[j];

	} // endfor i

	for(j=0;j<N_Errors;j++)
		errors[j] = sqrt(errors[j]);

	delete Param[0];
	delete AuxArray[0];
	delete SecondDer;
	delete ExactVal[0];
	delete Derivatives[0];
} // TFEFunction2D::GetErrorsAdapt
#endif
/** calculate errors to given function taylored for use in OPTPDE (is called from
    GetErrorsAdaptOPTPDE) */
void TFEFunction2D::GetErrorsOPTPDE(DoubleFunct2D *Exact, int N_Derivatives,
			      MultiIndex2D *NeededDerivatives,
			      int N_Errors, ErrorMethod2D *ErrorMeth, 
			      CoeffFct2D *Coeff, 
			      TAuxParam2D *Aux,
			      int n_fespaces, TFESpace2D **fespaces, int& kink,
			      double upper, double lower, double *errors)
{
	int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
	int N_Cells, N_Points, N_Parameters, N_;
	int Used[N_FEs2D], *N_BaseFunct;
	TFESpace2D *fespace;
	FE2D LocalUsedElements[N_FEs2D], CurrentElement;
	BaseFunct2D BaseFunct, *BaseFuncts;
	TCollection *Coll;
	TBaseCell *cell;
	TFE2D *ele;
	double *weights, *xi, *eta;
	double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
	double AbsDetjk[MaxN_QuadPoints_2D];
	RefTrans2D RefTrans;
	double *Param[MaxN_QuadPoints_2D], *aux;
	double *Derivatives[MaxN_QuadPoints_2D];
	double *ExactVal[MaxN_QuadPoints_2D];
	double *AuxArray[MaxN_QuadPoints_2D];
	int *DOF, ActiveBound, DirichletBound, end, last;
	double **OrigFEValues, *Orig, value;
	double FEFunctValues[MaxN_BaseFunctions2D];
	int *GlobalNumbers, *BeginIndex;
	double LocError[4];
	double hK;
	bool *SecondDer;
	TVertex *vertex;
	double loc_x, loc_y, loc_r;
  
	BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
	N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

	SecondDer = new bool[n_fespaces];
	for(i=0;i<n_fespaces;i++)
		SecondDer[i] = FALSE;

	N_Parameters = Aux->GetN_Parameters();
	aux = new double [MaxN_QuadPoints_2D*N_Parameters];
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		Param[j] = aux + j*N_Parameters;

	aux = new double [MaxN_QuadPoints_2D*N_Derivatives];
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		Derivatives[j] = aux + j*N_Derivatives;
  
	aux = new double [MaxN_QuadPoints_2D * 4];
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		ExactVal[j] = aux + j*4;

  // 20 <= number of term
	aux = new double [MaxN_QuadPoints_2D*20]; 
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		AuxArray[j] = aux + j*20;

	GlobalNumbers = FESpace2D->GetGlobalNumbers();
	BeginIndex = FESpace2D->GetBeginIndex();

	for(i=0;i<N_Errors;i++)
		errors[i] = 0.0;
	errors[N_Errors] = 0.0;  // for L_infty-error

// ########################################################################
// loop over all cells
// ########################################################################
	Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
	N_Cells = Coll->GetN_Cells();
 
	for(i=0;i<N_Cells;i++)
	{
		cell = Coll->GetCell(i);

		hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
		memset(Used, 0, N_FEs2D*SizeOfInt);
		for(j=0;j<n_fespaces;j++)
		{
			CurrentElement = fespaces[j]->GetFE2D(i, cell);
			Used[CurrentElement] = 1;
		}

		N_LocalUsedElements = 0;
		memset(LocalUsedElements, 0, SizeOfInt*N_FEs2D);
		j = 0;
		for(k=0;k<N_FEs2D;k++)
			if(Used[k])
		{
			LocalUsedElements[j] = (FE2D)k;
			j++;
		}
		N_LocalUsedElements = j;
    

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
		TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
				       Coll, cell, SecondDer,
				       N_Points, xi, eta, weights, X, Y, AbsDetjk);

		if(N_Parameters>0)
			Aux->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param); 
    
    // calculate all needed derivatives of this FE function
		CurrentElement = FESpace2D->GetFE2D(i, cell);
		BaseFunct = BaseFuncts[CurrentElement];
		N_ = N_BaseFunct[CurrentElement];

		DOF = GlobalNumbers + BeginIndex[i];
		for(l=0;l<N_;l++)
			FEFunctValues[l] = Values[DOF[l]];

		for(k=0;k<N_Derivatives;k++)
		{
			OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
					NeededDerivatives[k]);
			for(j=0;j<N_Points;j++)
			{
				Orig = OrigFEValues[j];
				value = 0;
				for(l=0;l<N_;l++)
				{
					value += FEFunctValues[l] * Orig[l];
				} // endfor l
				// here we have to perform the projection for control
				if(k==0)
				{
					if (value > upper) 
					{
						value = upper;
						kink = 1;
					}
					else if (value < lower)
					{
						value = lower;
						kink = 1;
					}
				}
				Derivatives[j][k] = value;
			} // endfor j
		} // endfor k

		//errors[N_Errors] = 0.0;  //added D.Sirch
		for(j=0;j<N_Points;j++)
		{
			Exact(X[j], Y[j], ExactVal[j]);
       // D.Sirch: computation of L^\inf-error
			if(fabs(*ExactVal[j] - Derivatives[j][0]) > errors[N_Errors])
			{
				errors[N_Errors] = fabs(*ExactVal[j] - Derivatives[j][0]);
			}
		}

		if(Coeff)
			Coeff(N_Points, X, Y, Param, AuxArray);

		ErrorMeth(N_Points, X, Y, AbsDetjk, weights, hK, Derivatives, 
			  ExactVal, AuxArray, LocError);

		for(j=0;j<N_Errors;j++)
			errors[j] += LocError[j];

	} // endfor i

	for(j=0;j<N_Errors;j++)
		errors[j] = sqrt(errors[j]);

	delete Param[0];
	delete AuxArray[0];
	delete SecondDer;
	delete ExactVal[0];
	delete Derivatives[0];
} // TFEFunction2D::GetErrorsOPTPDE

#ifdef __2D__
/** calculate errors to given function taylored for use in OPTPDE */
void TFEFunction2D::GetErrorsAdaptOPTPDE(DoubleFunct2D *Exact, int N_Derivatives,
				   MultiIndex2D *NeededDerivatives,
				   int N_Errors, ErrorMethod2D *ErrorMeth, 
				   CoeffFct2D *Coeff, 
				   TAuxParam2D *Aux,
				   int n_fespaces, TFESpace2D **fespaces,
				   double radius, double upper, double lower,double *errors)
{
	int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
	int N_Cells, N_Points, N_Parameters, N_;
	int Used[N_FEs2D], *N_BaseFunct;
	TFESpace2D *fespace;
	FE2D LocalUsedElements[N_FEs2D], CurrentElement;
	BaseFunct2D BaseFunct, *BaseFuncts;
	TCollection *Coll;
	TBaseCell *cell;
	TFE2D *ele;
	double *weights, *xi, *eta;
	double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
	double AbsDetjk[MaxN_QuadPoints_2D];
	RefTrans2D RefTrans;
	double *Param[MaxN_QuadPoints_2D], *aux;
	double *Derivatives[MaxN_QuadPoints_2D];
	double *ExactVal[MaxN_QuadPoints_2D];
	double *AuxArray[MaxN_QuadPoints_2D];
	int *DOF, ActiveBound, DirichletBound, end, last;
	double **OrigFEValues, *Orig, value;
	double FEFunctValues[MaxN_BaseFunctions2D];
	int *GlobalNumbers, *BeginIndex;
	double LocError[4];
	double hK;
	double *deriv, *exactval, w, t;
	bool *SecondDer;
	
	TGridCell *RootCell;
	TBaseCell **LocalCell;
	TCollection *LocalColl;
	FE2D *LocalFEs;
	TFESpace2D *RootSpace;
	TFEFunction2D *RootFct;
	double *RootVal, RootError;
	int RootLen;
	
	const int MaxLev = 5;
	int Lev, LastLev;
	TCollection *FineColls[MaxLev];
	TFEFunction2D *FineFcts[MaxLev];
	TFESpace2D *FineSpaces[MaxLev];
	TCollection *CurrColl;
	TGridCell *CurrCell;
	int N_FineCells, N_CoarseCells, CellId;
	TBaseCell **FineCells;
	FE2D *FineFEs[MaxLev];
	TBdLine *Edges[4];
	double xv[4], yv[4];
	double *FineVals[MaxLev];
	double *AuxVals[MaxLev];
	int FineDOF;
	TAuxParam2D *FineAux;
	int *NewDOF, *NewGlobalNumbers, *NewBeginIndex;
	int kink;
	
	char Name[] = "name";
	char Description[] = "description";
	
	Edges[0] = new TBdLine(12345);
	Edges[1] = new TBdLine(12345);
	Edges[2] = new TBdLine(12345);
	Edges[3] = new TBdLine(12345);
  
	BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
	N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

	SecondDer = new bool[n_fespaces];
	for(i=0;i<n_fespaces;i++)
		SecondDer[i] = FALSE;

	N_Parameters = Aux->GetN_Parameters();
	aux = new double [MaxN_QuadPoints_2D*N_Parameters];
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		Param[j] = aux + j*N_Parameters;

	aux = new double [MaxN_QuadPoints_2D*N_Derivatives];
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		Derivatives[j] = aux + j*N_Derivatives;
  
	aux = new double [MaxN_QuadPoints_2D * 4];
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		ExactVal[j] = aux + j*4;

  // 20 <= number of term
	aux = new double [MaxN_QuadPoints_2D*20]; 
	for(j=0;j<MaxN_QuadPoints_2D;j++)
		AuxArray[j] = aux + j*20;

	GlobalNumbers = FESpace2D->GetGlobalNumbers();
	BeginIndex = FESpace2D->GetBeginIndex();

	for(i=0;i<N_Errors;i++)
		errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
	Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
	N_Cells = Coll->GetN_Cells();
 
	for(i=0;i<N_Cells;i++)
	{
		cell = Coll->GetCell(i);
		
		
		hK = cell->GetDiameter();
		kink=0;

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
		N_LocalUsedElements = 1;
		LocalUsedElements[0] = fespaces[0]->GetFE2D(i, cell);
		
		CurrentElement = FESpace2D->GetFE2D(i, cell);
		BaseFunct = BaseFuncts[CurrentElement];
		N_ = N_BaseFunct[CurrentElement];
		
		DOF = GlobalNumbers + BeginIndex[i];
		
		RootCell = new TGridCell(cell->GetRefDesc(),0);
		j = cell->GetN_Vertices();
		for(k=0;k<j;k++)
			cell->GetVertex(k)->GetCoords(xv[k], yv[k]);
		
		for(k=0;k<j;k++)
		{
			RootCell->SetVertex(k, cell->GetVertex(k));
			RootCell->SetJoint(k, new TBoundEdge(Edges[k], 0, 1));
			Edges[k]->SetParams(xv[k], yv[k], xv[(k+1)%j]-xv[k],
					    yv[(k+1)%j]-yv[k]);
		}
		LocalCell = new TBaseCell*[1];
		LocalCell[0] = (TBaseCell *)RootCell;
		LocalColl = new TCollection(1, LocalCell);
		LocalFEs = new FE2D[1];
		LocalFEs[0] = LocalUsedElements[0];
		
		RootSpace = new TFESpace2D(LocalColl, Name, Description,
					   OnlyDirichlet,
					   LocalFEs, NULL);
		RootLen = N_;
		RootVal = new double[RootLen];
		NewDOF = RootSpace->GetGlobalNumbers() + RootSpace->GetBeginIndex()[0];
		for(l=0;l<N_;l++)
			RootVal[NewDOF[l]] = Values[DOF[l]];
		RootFct = new TFEFunction2D(RootSpace, Name, Description, 
					    RootVal, RootLen);
		FineAux = new TAuxParam2D(1, 0, 0, 0, &RootSpace, NULL, NULL, 
					  NULL, NULL, 0, NULL);
		RootFct->GetErrorsOPTPDE(Exact, N_Derivatives, NeededDerivatives,
				   N_Errors, ErrorMeth, 
				   Coeff, FineAux, 1, &RootSpace, 
				   kink, upper, lower, LocError);
		delete FineAux;
		
		  // if the element has nonempty intersection with active set, do some adaptive integration
		if(kink==1)
//		if(1)
		{	
			//cout<<"kink="<<kink<<endl;
			RootError = LocError[0];
	
			FineColls[0] = LocalColl;
			FineFcts[0] = RootFct;
			FineSpaces[0] = RootSpace;
			FineFEs[0] = LocalFEs;
			FineVals[0] = RootVal;
			
			for(Lev=1;Lev<MaxLev;Lev++)
			{
				CurrColl = FineColls[Lev-1];
				N_CoarseCells = CurrColl->GetN_Cells();
				N_FineCells = 4*N_CoarseCells;
				FineCells = new TBaseCell*[N_FineCells];
				CellId = 0;
				for(j=0;j<N_CoarseCells;j++)
				{
					CurrCell = (TGridCell*)(CurrColl->GetCell(j));
					CurrCell->SetRegRefine();
					CurrCell->Refine(Lev);
					for(k=0;k<4;k++)
					{
						FineCells[CellId] =
							(TBaseCell*)(CurrCell->GetChild(k));
						CellId++;
					}
				}
				FineColls[Lev] = new TCollection(N_FineCells, FineCells);
				FineFEs[Lev] = new FE2D[N_FineCells];
				for(j=0;j<N_FineCells;j++)
					FineFEs[Lev][j] = LocalUsedElements[0];
				FineSpaces[Lev] = new TFESpace2D(FineColls[Lev], Name,
						Description, OnlyDirichlet,
						FineFEs[Lev], NULL);
				FineDOF = FineSpaces[Lev]->GetN_DegreesOfFreedom();
				FineVals[Lev] = new double[FineDOF];
				AuxVals[Lev] = new double[FineDOF];
				FineFcts[Lev] = new TFEFunction2D(FineSpaces[Lev], 
						Name, Description, 
						FineVals[Lev], FineDOF);
				Prolongate(FineSpaces[Lev-1], FineSpaces[Lev],
					FineVals[Lev-1], FineVals[Lev], AuxVals[Lev]);
				FineAux = new TAuxParam2D(1, 0, 0, 0, &(FineSpaces[Lev]), 
						NULL, NULL, NULL, NULL, 0, NULL);
				FineFcts[Lev]->GetErrorsOPTPDE(Exact, N_Derivatives,
						NeededDerivatives, N_Errors, ErrorMeth, 
						Coeff, FineAux, 1, &(FineSpaces[Lev]), 
						kink, upper, lower, LocError);
				delete FineAux;
 				if(fabs(LocError[0] - RootError) < 1e-3*RootError) break;
 				else RootError = LocError[0];
			} // end for Lev < MaxLev
			
			// delete all allocated memory
			 if(Lev==MaxLev) LastLev = MaxLev-1;
			 else LastLev = Lev;
			//for(Lev=MaxLev-1;Lev>0;Lev--)
			for(Lev=LastLev;Lev>0;Lev--)
			{
				N_CoarseCells = FineColls[Lev-1]->GetN_Cells();
				for(j=0;j<N_CoarseCells;j++)
				{
					FineColls[Lev-1]->GetCell(j)->Derefine();
				}
				delete AuxVals[Lev];
				delete FineVals[Lev];
				delete FineColls[Lev];
				delete FineFcts[Lev];
				delete FineSpaces[Lev];
			}
		}
		
		delete RootFct;
		delete RootSpace;
		delete RootVal;
		delete LocalColl;
		delete RootCell;
		
		for(j=0;j<N_Errors;j++)
			errors[j] += LocError[j]*LocError[j];

	} // endfor i

	for(j=0;j<N_Errors;j++)
		errors[j] = sqrt(errors[j]);
	errors[N_Errors] = LocError[N_Errors];

	delete Param[0];
	delete AuxArray[0];
	delete SecondDer;
	delete ExactVal[0];
	delete Derivatives[0];
} // TFEFunction2D::GetErrorsAdaptOPTPDE
#endif

/** determine the value of function and its first derivatives at
    the given point */
void TFEFunction2D::FindGradient(double x, double y, double *values)
{
  int i,j,k, N_Cells;
  double xv, yv, xi, eta;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  int N_BaseFunct;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref;

  int *Numbers, N_Found;
  double u, ux, uy;
  double val;
  int *GlobalNumbers, *BeginIndex;

  N_Found = 0;
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;

  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->PointInCell(x,y))
    {
      N_Found++;
      // cout << "point found" << endl;
      FE_ID = FESpace2D->GetFE2D(i, cell);
      FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
      RefTrans = FE_Obj->GetRefTransID();

      // set cell for reference transformation
      TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);

      // find local coordinates of the given point
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      // cout << " xi: " << xi << endl;
      // cout << "eta: " << eta << endl;

      // get base function object
      bf = FE_Obj->GetBaseFunct2D();
      N_BaseFunct = bf->GetDimension();

      uorig = new double[N_BaseFunct];
      uxorig = new double[N_BaseFunct];
      uyorig = new double[N_BaseFunct];

      uref = new double[N_BaseFunct];
      uxiref = new double[N_BaseFunct];
      uetaref = new double[N_BaseFunct];

      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);

      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
        uref, uxiref, uetaref, uorig, uxorig, uyorig);

      u = 0;
      ux = 0;
      uy = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(j=0;j<N_BaseFunct;j++)
      {
        val = Values[Numbers[j]];      
        u += uorig[j]*val;
        ux += uxorig[j]*val;
        uy += uyorig[j]*val;
      }

      values[0] += u;
      values[1] += ux;
      values[2] += uy;

      delete uorig;
      delete uxorig;
      delete uyorig;
      delete uref;
      delete uxiref;
      delete uetaref;

    }                                             // endif
  }                                               // endfor

  if(N_Found>0)
  {
    values[0] /= N_Found;
    values[1] /= N_Found;
    values[2] /= N_Found;
  }
  else
  {
   OutPut(" Point not found " <<   "x " <<x <<   " y " << y <<  endl);
   exit(0);
  }
  
}


/** determine the value of function and its first derivatives at
    the given point which lies within !! the cell *cell (not on the boundary
    of that cell !!) */
void TFEFunction2D::FindGradientLocal(TBaseCell *cell, int cell_no,
double x, double y, double *values)
{
  int j,k, N_Cells;
  double xv, yv, xi, eta, eps = 1e-20;
  TCollection *Coll;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  int N_BaseFunct;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref;

  int *Numbers;
  double u, ux, uy;
  double val;
  int *GlobalNumbers, *BeginIndex;

  Coll = FESpace2D->GetCollection();

  // get properties of the fe space
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();

  FE_ID = FESpace2D->GetFE2D(cell_no, cell);
  FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
  RefTrans = FE_Obj->GetRefTransID();

  // set cell for reference transformation
  TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);

  // find local coordinates of the given point
  TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
  // cout << " xi: " << xi << endl;
  // cout << "eta: " << eta << endl;

  // get base function object
  bf = FE_Obj->GetBaseFunct2D();
  N_BaseFunct = bf->GetDimension();

  // allocate arrays
  uorig = new double[N_BaseFunct];
  uxorig = new double[N_BaseFunct];
  uyorig = new double[N_BaseFunct];

  uref = new double[N_BaseFunct];
  uxiref = new double[N_BaseFunct];
  uetaref = new double[N_BaseFunct];

  // get values and derivatives of basis functions on the
  // reference mesh cell
  bf->GetDerivatives(D00, xi, eta, uref);
  bf->GetDerivatives(D10, xi, eta, uxiref);
  bf->GetDerivatives(D01, xi, eta, uetaref);

  // compute values on the original mesh cell
  TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
    uref, uxiref, uetaref,
    uorig, uxorig, uyorig);

  u = 0;
  ux = 0;
  uy = 0;
  Numbers = GlobalNumbers + BeginIndex[cell_no];
  for(j=0;j<N_BaseFunct;j++)
  {
    k = Numbers[j];
    val = Values[k];
    u += uorig[j]*val;
    ux += uxorig[j]*val;
    uy += uyorig[j]*val;
  }

  if (fabs(ux)<eps)
    ux=0.0;
  if (fabs(uy)<eps)
    uy=0.0;

  values[0] = u;
  values[1] = ux;
  values[2] = uy;

  // set memory free
  delete uorig;
  delete uxorig;
  delete uyorig;
  delete uref;
  delete uxiref;
  delete uetaref;
}


/** determine the value of function at
    the given point which lies within !! the cell *cell (not on the boundary
    of that cell !!) */
void TFEFunction2D::FindValueLocal(TBaseCell *cell, int cell_no,
double x, double y, double *values)
{
  int j,k, N_Cells;
  double xv, yv, xi, eta, eps = 1e-20;
  TCollection *Coll;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  int N_BaseFunct;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref;

  int *Numbers;
  double u;
  double val;
  int *GlobalNumbers, *BeginIndex;

  // get properties of the fe space
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();

  Coll = FESpace2D->GetCollection();

  FE_ID = FESpace2D->GetFE2D(cell_no, cell);
  FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
  RefTrans = FE_Obj->GetRefTransID();

  // set cell for reference transformation
  TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);

  // find local coordinates of the given point
  TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
  // cout << " xi: " << xi << endl;
  // cout << "eta: " << eta << endl;

  // get base function object
  bf = FE_Obj->GetBaseFunct2D();
  N_BaseFunct = bf->GetDimension();

  // allocate arrays
  uorig = new double[N_BaseFunct];
  uxorig = new double[N_BaseFunct];
  uyorig = new double[N_BaseFunct];

  uref = new double[N_BaseFunct];
  uxiref = new double[N_BaseFunct];
  uetaref = new double[N_BaseFunct];

  // get values and derivatives of basis functions on the
  // reference mesh cell
  bf->GetDerivatives(D00, xi, eta, uref);
  bf->GetDerivatives(D10, xi, eta, uxiref);
  bf->GetDerivatives(D01, xi, eta, uetaref);

  // compute values on the original mesh cell
  TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
    uref, uxiref, uetaref,
    uorig, uxorig, uyorig);

  u = 0;
  Numbers = GlobalNumbers + BeginIndex[cell_no];
  for(j=0;j<N_BaseFunct;j++)
  {
    k = Numbers[j];
    val = Values[k];
    u += uorig[j]*val;
  }

  values[0] = u;

  // set memory free
  delete uorig;
  delete uxorig;
  delete uyorig;
  delete uref;
  delete uxiref;
  delete uetaref;
}


/** calculate the interpolation of an exact function */
void TFEFunction2D::Interpolate(DoubleFunct2D *Exact)
{
  int i,j,k,l;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_, N_Points;
  double s, *xi, *eta;
  double Val[MaxN_BaseFunctions2D];
  double OutVal[MaxN_BaseFunctions2D];
  int *DOF, Index;
  RefTrans2D F_K;
  TRefTrans2D *rt;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D];
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double FctVal[4];
  int PolynomialDegree, ApproxOrder;
  QuadFormula2D QuadFormula;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Edges;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans, *RefTransArray;

  // begin code

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_DOFs = FESpace2D->GetN_DegreesOfFreedom();

  memset(Values, 0, SizeOfDouble*N_DOFs);
  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ApproxOrder = TFEDatabase2D::GetAccuracyFromFE2D(FEId);

    RefElement = Element->GetBaseFunct2D()->GetRefElement();

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree
          (3*PolynomialDegree);
        N_Edges = 4;
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree
          (3*PolynomialDegree-1);
        N_Edges = 3;
        break;
    }

    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint ||
          jointtype == IsoBoundEdge)
          IsIsoparametric = TRUE;
      }
    }                                             // endif

    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          break;
      }
    }                                             // endif IsIsoparametric

    switch(RefTrans)
    {
      case QuadAffin:
        rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
        ((TQuadAffin *)rt)->SetCell(cell);
        F_K = QuadAffin;
        break;
      case QuadBilinear:
        rt = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)rt)->SetCell(cell);
        F_K = QuadBilinear;
        break;
      case QuadIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(QuadIsoparametric);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        F_K = QuadIsoparametric;
        break;
      case TriaAffin:
        rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)rt)->SetCell(cell);
        F_K = TriaAffin;
        break;
      case TriaIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(TriaIsoparametric);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        F_K = TriaIsoparametric;
        break;
      default:
        cout << "unknown reftrans id: " << RefTrans << endl;
    }
    TFEDatabase2D::GetOrigFromRef(F_K, N_Points, xi, eta,
      X, Y, AbsDetjk);

    for(j=0;j<N_Points;j++)
    {
      Exact(X[j], Y[j], FctVal);
      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues,
                          FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];
  }
}

/** calculate parameters which are connected to a mesh cell */
void TFEFunction2D::GetMeshCellParams(DoubleFunct2D *Exact, int N_Derivatives,
                                      MultiIndex2D *NeededDerivatives,
                                      int N_Errors, ErrorMethod2D *ErrorMeth,
                                      CoeffFct2D *Coeff,
                                      TAuxParam2D *Aux,
                                      int n_fespaces, TFESpace2D **fespaces,
                                      double *errors, double *cell_parameters)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs2D], *N_BaseFunct;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  TFE2D *ele;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double *Derivatives[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool *SecondDer;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  N_Parameters = Aux->GetN_Parameters();
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux + j*N_Parameters;

  aux = new double [MaxN_QuadPoints_2D*N_Derivatives];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux + j*N_Derivatives;

  aux = new double [MaxN_QuadPoints_2D * 4];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    ExactVal[j] = aux + j*4;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*20;

  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  BeginIndex = FESpace2D->GetBeginIndex();

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;

  // ########################################################################
  // loop over all cells
  // ########################################################################
  Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    memset(Used, 0, N_FEs2D*SizeOfInt);
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i, cell);
      Used[CurrentElement] = 1;
    }

    N_LocalUsedElements = 0;
    memset(LocalUsedElements, 0, SizeOfInt*N_FEs2D);
    j = 0;
    for(k=0;k<N_FEs2D;k++)
      if(Used[k])
    {
      LocalUsedElements[j] = (FE2D)k;
      j++;
    }
    N_LocalUsedElements = j;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
      Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);

    if(N_Parameters>0)
      Aux->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);

    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace2D->GetFE2D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = Values[DOF[l]];

    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
        NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
        }                                         // endfor l
        Derivatives[j][k] = value;
      }                                           // endfor j
    }                                             // endfor k

    for(j=0;j<N_Points;j++)
      Exact(X[j], Y[j], ExactVal[j]);

    if(Coeff)
      Coeff(N_Points, X, Y, Param, AuxArray);

    ErrorMeth(N_Points, X, Y, AbsDetjk, weights, hK, Derivatives,
      ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
    {
      errors[j] += LocError[j];
      cell_parameters[i + j *N_Cells] = LocError[j];
    }

  }                                               // endfor i

  for(j=0;j<N_Errors;j++)
    errors[j] = sqrt(errors[j]);

  delete Param[0];
  delete AuxArray[0];
  delete SecondDer;
  delete ExactVal[0];
  delete Derivatives[0];
}                                                 // TFEFunction2D::GetErrors


/** calculate the super-convergence interpolation of an exact function */
void TFEFunction2D::InterpolateSuper(DoubleFunct2D *Exact)
{
  int i,j,k,l;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  BaseFunct2D BF;
  TNodalFunctional2D *nf;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_, N_Points;
  double s, *xi, *eta;
  double Val[MaxN_BaseFunctions2D];
  double OutVal[MaxN_BaseFunctions2D];
  int *DOF, Index;
  RefTrans2D F_K;
  TRefTrans2D *rt;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D];
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double FctVal[4];
  int PolynomialDegree, ApproxOrder;
  QuadFormula2D QuadFormula;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Edges;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans, *RefTransArray;

  // begin code

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_DOFs = FESpace2D->GetN_DegreesOfFreedom();

  memset(Values, 0, SizeOfDouble*N_DOFs);
  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    switch(nf->GetID())
    {
      case NF_C_Q_Q2_2D:
        nf = TFEDatabase2D::GetNodalFunctional2D(NF_S_Q_Q2_2D);
        break;
      default:
        cout << "unknown reftrans id: " << RefTrans << endl;
        exit(0);
      break;
    }
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ApproxOrder = TFEDatabase2D::GetAccuracyFromFE2D(FEId);

    RefElement = Element->GetBaseFunct2D()->GetRefElement();

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree
          (3*PolynomialDegree);
        N_Edges = 4;
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree
          (3*PolynomialDegree);
        N_Edges = 3;
        break;
    }

    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
      }
    }                                             // endif

    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          break;
      }
    }                                             // endif IsIsoparametric

    switch(RefTrans)
    {
      case QuadAffin:
        rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
        ((TQuadAffin *)rt)->SetCell(cell);
        F_K = QuadAffin;
        break;
      case QuadBilinear:
        rt = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)rt)->SetCell(cell);
        F_K = QuadBilinear;
        break;
      case QuadIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(QuadIsoparametric);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        F_K = QuadIsoparametric;
        break;
      case TriaAffin:
        rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)rt)->SetCell(cell);
        F_K = TriaAffin;
        break;
      case TriaIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(TriaIsoparametric);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        F_K = TriaIsoparametric;
        break;
      default:
        cout << "unknown reftrans id: " << RefTrans << endl;
    }
    TFEDatabase2D::GetOrigFromRef(F_K, N_Points, xi, eta,
      X, Y, AbsDetjk);

    for(j=0;j<N_Points;j++)
    {
      Exact(X[j], Y[j], FctVal);
      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues, FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];
  }
}


/** write the solution into a data file - written by Sashi **/
void TFEFunction2D::WriteSol()
{
  int i, N_Joints, N_Cells;
  static int img=0;
  char *BaseName, Dquot;

  TCollection *Coll;
  TBaseCell *cell;

  Dquot = 34; //  see ASCII Chart
  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  i=0;
  cell =  Coll->GetCell(i);
  N_Joints = cell->GetN_Joints();
  BaseName = TDatabase::ParamDB->VTKBASENAME;

  std::ostringstream os;
  os << " ";

  OutPut("Writing solution into "<< BaseName << img<< ".Sol MooNMD file"<< endl);

  os.seekp(std::ios::beg);
  os << BaseName<<img<<".Sol" << ends;
  std::ofstream dat(os.str().c_str());
  img++;
  if (!dat)
   {
    cerr << "cannot open file for output" << endl;
    exit(0);
   }

    dat << "# Solution of the scalar "<<Dquot<<Name<<Dquot<<", written by MooNMD"<< endl;
    dat << "# Author: Sashikumaar Ganesan" <<  endl;
    dat << "#  N_Cells, Cell_Type, N_Dim, N_Dof" <<  endl;
    dat <<N_Cells << " " << N_Joints << " " <<  Length<< endl;
    dat <<  endl;

    dat << "# Dof "<< " Nodal Values"<< endl;

    for(i=0;i<Length;i++)
     dat << i << " " << Values[i] << endl;
}


/** Read the solution from a given data file - written by Sashi **/
void TFEFunction2D::ReadSol(char *BaseName)
{
 int i, j, N_Joints, N_Cells, N_cells, N_joints, length;
 char line[100];

 TCollection *Coll;
 TBaseCell *cell;

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  i=0;
  cell =  Coll->GetCell(i);
  N_Joints = cell->GetN_Joints();

  std::ifstream dat(BaseName);
  if (!dat)
   {
    cerr << "cannot open '" <<  BaseName << "' for input" << endl;
    exit(-1);
   }

  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat >> N_cells >> N_joints >> length;
  dat.getline (line, 99);

  if(N_cells!=N_Cells || N_joints!=N_Joints || length!=Length )
   {
    OutPut("Given data file does not match with this FE scalar function !"<<endl);
    exit(0);
   }

  dat.getline (line, 99);
  OutPut("Reading nodal values of the FE saclar function !"<<endl);

  for(i=0;i<Length;i++)
   {
    dat.getline (line, 99);
    dat >> j >> Values[i];
   }

  dat.close();
}


/** interpolate the old mesh fe function values to the new fe function */
void TFEFunction2D::Interpolate(TFEFunction2D *OldFeFunction)
{
  int i,j,k,l, N_Cells, N_Edges;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_, N_Points, *DOF, Index;
  int PolynomialDegree, ApproxOrder;
  int *IntIndex;
  
  double s, *xi, *eta;
  double Val[MaxN_BaseFunctions2D];
  double OutVal[MaxN_BaseFunctions2D];
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D];
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double FctVal[4], values[4];
  
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  RefTrans2D F_K;
  TRefTrans2D *rt;
  QuadFormula2D QuadFormula;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  BF2DRefElements RefElement;
  RefTrans2D RefTrans, *RefTransArray;

  // begin code

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_DOFs = FESpace2D->GetN_DegreesOfFreedom();

  
  IntIndex = new int[N_DOFs];
  memset(IntIndex, 0, SizeOfInt*N_DOFs);
   
  memset(Values, 0, SizeOfDouble*N_DOFs);
  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ApproxOrder = TFEDatabase2D::GetAccuracyFromFE2D(FEId);

    RefElement = Element->GetBaseFunct2D()->GetRefElement();

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree
          (3*PolynomialDegree);
        N_Edges = 4;
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree
          (3*PolynomialDegree-1);
        N_Edges = 3;
        break;
    }

    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint ||
          jointtype == IsoBoundEdge)
          IsIsoparametric = TRUE;
      }
    }                                             // endif

    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          break;
      }
    }                                             // endif IsIsoparametric

    switch(RefTrans)
    {
      case QuadAffin:
        rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
        ((TQuadAffin *)rt)->SetCell(cell);
        F_K = QuadAffin;
        break;
      case QuadBilinear:
        rt = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)rt)->SetCell(cell);
        F_K = QuadBilinear;
        break;
      case QuadIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(QuadIsoparametric);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        F_K = QuadIsoparametric;
        break;
      case TriaAffin:
        rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)rt)->SetCell(cell);
        F_K = TriaAffin;
        break;
      case TriaIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(TriaIsoparametric);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        F_K = TriaIsoparametric;
        break;
      default:
        cout << "unknown reftrans id: " << RefTrans << endl;
    }
    TFEDatabase2D::GetOrigFromRef(F_K, N_Points, xi, eta, X, Y, AbsDetjk);

    for(j=0;j<N_Points;j++)
    {
      OldFeFunction->FindGradient(X[j], Y[j], values);
      PointValues[j] = values[0]; 
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues,
                          FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
     {
      Values[DOF[j]] += FunctionalValues[j];
      IntIndex[DOF[j]] ++;
     }
  } //for(i=0;i<N_
  
  
  for(i=0;i<N_DOFs;i++)
   {
    if(IntIndex[i])
     { Values[i] /=(double)IntIndex[i];}
    else
    {     
      OutPut(" Error in interpolateion value not interpolated : " <<    i  << endl); 
      exit(1);
    }
   }
   
   
   delete [] IntIndex;
   
}

/** correct the sol to the Old_Mass (remessing, temp, surfact, psd, etc) - added by sashi */
void  TFEFunction2D::CorrectMass(double OldMass)
{
  int i;
  
  double sum=0., MeanDiff, params[3], *w, *OrigSol, WeightedVol;
  
//   N_DOFs = FESpace2D->GetN_DegreesOfFreedom();
     
  for(i=0;i<Length;i++)
   { sum +=Values[i]; }
 
  // if the mean value is zero, do nothing
  if(fabs(sum)>1e-12 )
   { 
    // add the MeanMassDiff to the sol with weights to avoid negative values
    w = new double[Length];
    OrigSol = new double[Length];
    
    // store the orignal values of the fefunction
    memcpy(OrigSol, Values, Length*SizeOfDouble);
 
    //weights, sol with 0 values will not be altered, otherwise sol vale become negative due to correction 
    for(i=0;i<Length;i++)
     { 
       w[i] = Values[i]/sum; 
//         w[i] = 1.;
     }
     
    // calculate the weitgted volume, if w=1, this step is not needed, but sol may take negative
    memcpy(Values, w, Length*SizeOfDouble);   
    this->GetMassAndMean(params); 
    WeightedVol = params[0];
    
    // copy back the orig values to fefunction
    memcpy(Values, OrigSol, Length*SizeOfDouble);   
    this->GetMassAndMean(params);
    MeanDiff = (OldMass - params[0])/WeightedVol;
      

//    cout << "MeanDiff " << MeanDiff << endl;
 
    for(i=0;i<Length;i++)
     {Values[i] += (MeanDiff*w[i]); }
     
     delete [] w;
     delete [] OrigSol;
    } // if(fabs(sum)>1e-12 )
   
} // CorrectMass


/** Retun the mass, domain volume and mean values of the function - added by Sashi*/
void  TFEFunction2D::GetMassAndMean(double *OutVal)
{
  int i,j,k,l, polydegree, N_QFPoints, ORDER;
  int N_Cells, N_Joints, N_Vertices;
  int *BeginIndex, *GlobalNumbers, *DOF, N_BF;

  double Mult, r_axial, val, mass, volume;
  double *weights, *xi, *eta;
  double values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double AbsDetjk[MaxN_QuadPoints_2D], X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];

  TJoint *joint;
  TBaseCell *cell;
  TCollection *coll;
  JointType jointtype;
  BoundTypes bdtype;
  RefTrans2D RefTrans;
  boolean IsIsoparametric;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  FE2D FEid;
  TBaseFunct2D *bf;
  TRefTrans2D *F_K;

//   FeSpace = fefunction->GetFESpace2D();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
//   U = fefunction->GetValues();
 
  coll = FESpace2D->GetCollection();
  N_Cells = coll->GetN_Cells();

  mass = 0.;
  volume = 0.;
//   Concentration = 0.;
  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    FEid = FESpace2D->GetFE2D(i, cell);
    
     RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEid);
     N_Joints = cell->GetN_Joints();
     IsIsoparametric = FALSE;
     if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
      {
       for(j=0;j<N_Joints;j++)
        {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
         {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)  IsIsoparametric = TRUE;
         }
        if(jointtype == InterfaceJoint)
         {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
         }
        if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
         IsIsoparametric = TRUE;

        } // for(j=0;j<  
       } // if(TDatabase::ParamDB->USE_ISOPARAMETRIC)   
     
   if(IsIsoparametric)
    {
      switch(N_Joints)
      {
        case 4:
          RefTrans = QuadIsoparametric;
        break;

        case 3:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEid);
    switch(RefTrans)
    {
      case TriaAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaAffin *)F_K)->SetCell(cell);
//         locvol = ((TTriaAffin *)rt)->GetVolume();
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case TriaIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)F_K)->SetCell(cell);
//         locvol = ((TTriaIsoparametric *)F_K)->GetVolume();
        ((TTriaIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadAffin *)F_K)->SetCell(cell);
//         locvol = ((TQuadAffin *)rt)->GetVolume();
        ((TQuadAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadBilinear:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)F_K)->SetCell(cell);
//         locvol = ((TQuadIsoparametric *)rt)->GetVolume();
        ((TQuadIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;
    }
       
     // find basis functions on cell i
    bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
    N_BF = bf->GetDimension();
    DOF = GlobalNumbers + BeginIndex[i];

    for(k=0;k<N_QFPoints;k++)
     {
      bf->GetDerivatives(D00, xi[k], eta[k], values[k]);
      
      if(TDatabase::ParamDB->Axial3D)
       {
        if(TDatabase::ParamDB->Axial3DAxis==1)
         {
          r_axial = Y[k];  //   (X: symmetric  problems   
         }
        else
         {
          r_axial = X[k];  // Y: symmetric problems     
         }   

      if(r_axial<=0)
       {
        cout <<"X[k] negative in GetMassAndMean, change Quad rule " <<  r_axial <<endl;
//         exit(0);
       }
       
       r_axial = fabs(r_axial);      
       }
      else
      {
       r_axial = 1.;
      }

      Mult = r_axial*weights[k]*AbsDetjk[k];
      val = 0.;
      for(l=0;l<N_BF;l++)
       {
        j = DOF[l];
        val += Values[j]*values[k][l];
       }

     mass += val*Mult;
     volume += Mult;
    } //  for(k=0;k<N_QFPoints;      
   } // for(i=0;i<N_Cells
  
   if(TDatabase::ParamDB->Axial3D)  
    {
     OutVal[0] = 2.*Pi*mass;
     OutVal[1] = 2.*Pi*volume;
    }
   else
    {
     OutVal[0] = mass;
     OutVal[1] = volume;     
    }

   OutVal[2] = mass/volume;
} // GetMassAndMean







