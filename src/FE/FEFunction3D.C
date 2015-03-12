// =======================================================================
// %W% %G%
// 
// Class:       TFEFunction3D
// Purpose:     a function from a finite element space in 3D
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
#include <BoundFace.h>
#include <FEDatabase3D.h>
#include <FEFunction3D.h>
#include <string.h>
#include <AllRefTrans3D.h>
// #include <NodalFunctional3D.h>

#include <stdlib.h>
#include <InterfaceJoint.h>
#include <BdPlane.h>
#include <LinAlg.h>

void OnlyDirichlet(double x, double y, double z, BoundCond &cond)
{
	cond = DIRICHLET;
}

/** constructor with vector initialization */
TFEFunction3D::TFEFunction3D(TFESpace3D *fespace3D, char *name, 
                             char *description, double *values, int length)
{
  FESpace3D=fespace3D;

  Name=strdup(name);

  Description=strdup(description);

  Values=values;

  Length=length;
}

TFEFunction3D::~TFEFunction3D()
{
  free(Name);
  free(Description);
}


/** calculate errors to given function */
void TFEFunction3D::GetErrors(DoubleFunct3D *Exact, int N_Derivatives,
                              MultiIndex3D *NeededDerivatives,
                              int N_Errors, ErrorMethod3D *ErrorMeth, 
                              CoeffFct3D *Coeff, 
                              TAuxParam3D *Aux,
                              int n_fespaces, TFESpace3D **fespaces,
                              double *errors)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs3D], *N_BaseFunct;
  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  RefTrans3D RefTrans;
  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool *SecondDer;
  double loc_x, loc_y, loc_z, loc_r;

  
#ifdef _MPI
   int ID, rank;
   
   MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);    
#endif
      
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  N_Parameters = Aux->GetN_Parameters();
  
  if(N_Parameters==0)
   aux = NULL;
  else
   aux = new double [MaxN_QuadPoints_3D*N_Parameters];
  
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Param[j] = aux + j*N_Parameters;

  aux = new double [MaxN_QuadPoints_3D*N_Derivatives];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux + j*N_Derivatives;
  
  aux = new double [MaxN_QuadPoints_3D * 5];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    ExactVal[j] = aux + j*5;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_3D*20]; 
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;

  GlobalNumbers = FESpace3D->GetGlobalNumbers();
  BeginIndex = FESpace3D->GetBeginIndex();

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  
  TVertex *vertex;

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

#ifdef _MPI
    ID = cell->GetSubDomainNo();
    if(rank!=ID) continue;
#endif
    
    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    memset(Used, 0, N_FEs3D*SizeOfInt);
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE3D(i, cell);
      Used[CurrentElement] = 1;
    }

    N_LocalUsedElements = 0;
    memset(LocalUsedElements, 0, SizeOfInt*N_FEs3D);
    j = 0;
    for(k=0;k<N_FEs3D;k++)
      if(Used[k])
      {
        LocalUsedElements[j] = (FE3D)k;
        j++;
      }
    N_LocalUsedElements = j;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta, weights,
                           X, Y, Z, AbsDetjk);

    if(N_Parameters>0)
      Aux->GetParameters(N_Points, cell, i, xi, eta, zeta,
                         X, Y, Z, Param); 

    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace3D->GetFE3D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = Values[DOF[l]];

    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct,
                                      NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
        } // endfor l
        Derivatives[j][k] = value;
      } // endfor j
    } // endfor k
 
    errors[N_Errors] = 0.0;  //added D.Sirch
    for(j=0;j<N_Points;j++)
    {
      Exact(X[j], Y[j], Z[j], ExactVal[j]);
      
       // D.Sirch: computation of L^\inf-error
      if(fabs(*ExactVal[j] - Derivatives[j][0]) > errors[N_Errors])
        errors[N_Errors] = fabs(*ExactVal[j] - Derivatives[j][0]);
    }

    if(Coeff)
      Coeff(N_Points, X, Y, Z, Param, AuxArray);

    ErrorMeth(N_Points, X, Y, Z, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
      errors[j] += LocError[j];

  } // endfor i

#ifndef _MPI // sqrt(errors[j]) in the main programm after collecting error from all subdomains
  for(j=0;j<N_Errors;j++)
    errors[j] = sqrt(errors[j]);
#endif
  
  delete [] AuxArray[0];
  delete [] SecondDer;
  delete [] ExactVal[0];
  delete [] Derivatives[0];
  
  if(Param[0])
   delete [] Param[0];
  
} // TFEFunction3D::GetErrors

void TFEFunction3D::GetErrorsAdapt(DoubleFunct3D *Exact, int N_Derivatives,
				   MultiIndex3D *NeededDerivatives,
				   int N_Errors, ErrorMethod3D *ErrorMeth, 
				   CoeffFct3D *Coeff, 
				   TAuxParam3D *Aux,
				   int n_fespaces, TFESpace3D **fespaces,
				   double *errors)
{
	int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
	int N_Cells, N_Points, N_Parameters, N_;
	int Used[N_FEs3D], *N_BaseFunct;
	TFESpace3D *fespace;
	FE3D LocalUsedElements[N_FEs3D], CurrentElement;
	BaseFunct3D BaseFunct, *BaseFuncts;
	TCollection *Coll;
	TBaseCell *cell;
	TFE3D *ele;
	double *weights, *xi, *eta;
	double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
	double AbsDetjk[MaxN_QuadPoints_3D];
	RefTrans3D RefTrans;
	double *Param[MaxN_QuadPoints_3D], *aux;
	double *Derivatives[MaxN_QuadPoints_3D];
	double *ExactVal[MaxN_QuadPoints_3D];
	double *AuxArray[MaxN_QuadPoints_3D];
	int *DOF, ActiveBound, DirichletBound, end, last;
	double **OrigFEValues, *Orig, value;
	double FEFunctValues[MaxN_BaseFunctions3D];
	int *GlobalNumbers, *BeginIndex;
	double LocError[4];
	double hK;
	double *deriv, *exactval, w, t;
	bool *SecondDer;
	
	TGridCell *RootCell;
	TBaseCell **LocalCell;
	TCollection *LocalColl;
	FE3D *LocalFEs;
	TFESpace3D *RootSpace;
	TFEFunction3D *RootFct;
	double *RootVal;
	int RootLen;
	
	const int MaxLev = 3;
	int Lev;
	TCollection *FineColls[MaxLev];
	TFEFunction3D *FineFcts[MaxLev];
	TFESpace3D *FineSpaces[MaxLev];
	TCollection *CurrColl;
	TGridCell *CurrCell;
	int N_FineCells, N_CoarseCells, CellId;
	TBaseCell **FineCells;
	FE3D *FineFEs[MaxLev];
	TBdPlane* Faces[6];
	double xv[8], yv[8], zv[8];
	double param1[4], param2[4];
	double a_x, a_y, a_z, b_x, b_y, b_z, norm;
	double c_x, c_y, c_z, d_x, d_y, d_z, n_x, n_y, n_z;
	double *FineVals[MaxLev];
	double *AuxVals[MaxLev];
	int FineDOF;
	TAuxParam3D *FineAux;
	int *NewDOF, *NewGlobalNumbers, *NewBeginIndex;
        int N_Faces, N_Vertices;
        const int *TmpFV, *TmpLen;
        int MaxLen, first, second, third;
	
	char Name[] = "name";
	char Description[] = "description";
	
	Faces[0] = new TBdPlane(12345);
	Faces[1] = new TBdPlane(12345);
	Faces[2] = new TBdPlane(12345);
	Faces[3] = new TBdPlane(12345);
	Faces[4] = new TBdPlane(12345);
	Faces[5] = new TBdPlane(12345);
  
	BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
	N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

	SecondDer = new bool[n_fespaces];
	for(i=0;i<n_fespaces;i++)
		SecondDer[i] = FALSE;

	N_Parameters = Aux->GetN_Parameters();
	aux = new double [MaxN_QuadPoints_3D*N_Parameters];
	for(j=0;j<MaxN_QuadPoints_3D;j++)
		Param[j] = aux + j*N_Parameters;

	aux = new double [MaxN_QuadPoints_3D*N_Derivatives];
	for(j=0;j<MaxN_QuadPoints_3D;j++)
		Derivatives[j] = aux + j*N_Derivatives;
  
	aux = new double [MaxN_QuadPoints_3D * 5];
	for(j=0;j<MaxN_QuadPoints_3D;j++)
		ExactVal[j] = aux + j*5;

  // 20 <= number of term
	aux = new double [MaxN_QuadPoints_3D*20]; 
	for(j=0;j<MaxN_QuadPoints_3D;j++)
		AuxArray[j] = aux + j*20;

	GlobalNumbers = FESpace3D->GetGlobalNumbers();
	BeginIndex = FESpace3D->GetBeginIndex();

	for(i=0;i<N_Errors;i++)
		errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
	Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
	N_Cells = Coll->GetN_Cells();
 
	for(i=0;i<N_Cells;i++)
	{
                // OutPut("memory at start of cell " << i << " " << GetMemory() << endl);
		cell = Coll->GetCell(i);
		
                if(cell->GetType() != Tetrahedron)
                {
                  Error("Only tetrahedra are allowed." << endl);
                  exit(-1);
                }
		
		hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
		N_LocalUsedElements = 1;
		LocalUsedElements[0] = fespaces[0]->GetFE3D(i, cell);
		
		CurrentElement = FESpace3D->GetFE3D(i, cell);
		BaseFunct = BaseFuncts[CurrentElement];
		N_ = N_BaseFunct[CurrentElement];
		
		DOF = GlobalNumbers + BeginIndex[i];
		
		RootCell = new TGridCell(cell->GetRefDesc(),0);
		N_Vertices = cell->GetN_Vertices();
		for(k=0;k<N_Vertices;k++)
                {
			cell->GetVertex(k)->GetCoords(xv[k], yv[k], zv[k]);
 			RootCell->SetVertex(k, cell->GetVertex(k));
		}

                RootCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

		N_Faces = cell->GetN_Faces();

		for(k=0;k<N_Faces;k++)
		{	
                  // cout << endl << "joint: " << k << endl;
                  first  = TmpFV[k*MaxLen + 0];
                  second = TmpFV[k*MaxLen + 1];
                  third  = TmpFV[k*MaxLen + 2];

                  // vector a = second vertex - first vertex
                  a_x = xv[second] - xv[first];
                  a_y = yv[second] - yv[first];
                  a_z = zv[second] - zv[first];
                  // cout << "a: " << a_x << " " << a_y << " " << a_z << endl;

                  // vector c = third vertex - first vertex
                  c_x = xv[third] - xv[first];
                  c_y = yv[third] - yv[first];
                  c_z = zv[third] - zv[first];
                  // cout << "c: " << c_x << " " << c_y << " " << c_z << endl;
			
                  // outer normal = a \times c
		  n_x = -(a_y*c_z - a_z*c_y);
		  n_y = -(a_z*c_x - a_x*c_z);
		  n_z = -(a_x*c_y - a_y*c_x); 
                  // cout << "n: " << n_x << " " << n_y << " " << n_z << endl;
                  norm = sqrt(n_x*n_x + n_y*n_y + n_z*n_z);
                  n_x /= norm;
                  n_y /= norm;
                  n_z /= norm;
                  // cout << "n: " << n_x << " " << n_y << " " << n_z << endl;
			
		  Faces[k]->SetParams(xv[first],yv[first],zv[first],
                                      a_x, a_y, a_z, n_x, n_y, n_z);
			
                  // vector b = a \times n
		  b_x = a_y*n_z - a_z*n_y;
		  b_y = a_z*n_x - a_x*n_z;
		  b_z = a_x*n_y - a_y*n_x; 
                  // cout << "b: " << b_x << " " << b_y << " " << b_z << endl;

		  param1[0] = 0.0;
		  param2[0] = 0.0;

		  param1[1] = 1.0;
		  param2[1] = 0.0;

                  norm = a_x * b_y - a_y * b_x;
                  if(fabs(norm)<1e-8)
                  {
                    norm = a_y * b_z - a_z * b_y;
                    if(fabs(norm)<1e-8)
                    {
                      norm = a_x * b_z - a_z * b_x;
		      param1[2] = ( b_z*c_x - b_x*c_z )/norm;
		      param2[2] = (-a_z*c_x + a_x*c_z )/norm; 
                    }
                    else
                    {
		      param1[2] = ( b_z*c_y - b_y*c_z )/norm;
		      param2[2] = (-a_z*c_y + a_y*c_z )/norm; 
                    }
                  }
                  else
                  {
		    param1[2] = ( b_y*c_x - b_x*c_y )/norm;
		    param2[2] = (-a_y*c_x + a_x*c_y )/norm; 
                  } 
                  // cout << "Norm: " << norm << endl;

                  // cout << "Param: " << param1[2] << "  " << param2[2] << endl;
			
		  RootCell->SetJoint(k, new TBoundFace(Faces[k], param1, param2));
 		}

		LocalCell = new TBaseCell*[1];
		LocalCell[0] = (TBaseCell *)RootCell;
		LocalColl = new TCollection(1, LocalCell);
		LocalFEs = new FE3D[1];
		LocalFEs[0] = LocalUsedElements[0];
		
		RootSpace = new TFESpace3D(LocalColl, Name, Description,
					   OnlyDirichlet,LocalFEs);
		RootLen = N_;
		RootVal = new double[RootLen];
		NewDOF = RootSpace->GetGlobalNumbers() + RootSpace->GetBeginIndex()[0];
		for(l=0;l<N_;l++)
			RootVal[NewDOF[l]] = Values[DOF[l]];
		RootFct = new TFEFunction3D(RootSpace, Name, Description, 
					    RootVal, RootLen);
		FineAux = new TAuxParam3D(1, 0, 0, 0, &RootSpace, NULL, NULL, 
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
			N_FineCells = 8*N_CoarseCells;
			FineCells = new TBaseCell*[N_FineCells];
			CellId = 0;
			for(j=0;j<N_CoarseCells;j++)
			{
				CurrCell = (TGridCell*)(CurrColl->GetCell(j));
				CurrCell->SetRegRefine();
				CurrCell->Refine(Lev);
				for(k=0;k<8;k++)
				{
					FineCells[CellId] = (TBaseCell*)(CurrCell->GetChild(k));
					CellId++;
				}
			}
			FineColls[Lev] = new TCollection(N_FineCells, FineCells);
			FineFEs[Lev] = new FE3D[N_FineCells];
			for(j=0;j<N_FineCells;j++)
				FineFEs[Lev][j] = LocalUsedElements[0];
			FineSpaces[Lev] = new TFESpace3D(FineColls[Lev], Name, Description,
					OnlyDirichlet,FineFEs[Lev]);
			FineDOF = FineSpaces[Lev]->GetN_DegreesOfFreedom();
			FineVals[Lev] = new double[FineDOF];
			AuxVals[Lev] = new double[FineDOF];
			FineFcts[Lev] = new TFEFunction3D(FineSpaces[Lev], 
					Name, Description, 
					FineVals[Lev], FineDOF);
			Prolongate(FineSpaces[Lev-1], FineSpaces[Lev],
				   FineVals[Lev-1], FineVals[Lev], AuxVals[Lev]);
			FineAux = new TAuxParam3D(1, 0, 0, 0, &(FineSpaces[Lev]), 
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

                // OutPut("memory at end of cell " << i << " " << GetMemory() << endl);
	} // endfor i

	for(j=0;j<N_Errors;j++)
		errors[j] = sqrt(errors[j]);

	delete Param[0];
	delete AuxArray[0];
	delete SecondDer;
	delete ExactVal[0];
	delete Derivatives[0];
} // TFEFunction3D::GetErrorsAdapt

/** calculate errors to given function taylored for OPTPDE*/
void TFEFunction3D::GetErrorsOPTPDE(DoubleFunct3D *Exact, int N_Derivatives,
			      MultiIndex3D *NeededDerivatives,
			      int N_Errors, ErrorMethod3D *ErrorMeth, 
			      CoeffFct3D *Coeff, 
			      TAuxParam3D *Aux,
			      int n_fespaces, TFESpace3D **fespaces,
			      double radius, double upper, double lower, double *errors)
{
	int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
	int N_Cells, N_Points, N_Parameters, N_;
	int Used[N_FEs3D], *N_BaseFunct;
	TFESpace3D *fespace;
	FE3D LocalUsedElements[N_FEs3D], CurrentElement;
	BaseFunct3D BaseFunct, *BaseFuncts;
	TCollection *Coll;
	TBaseCell *cell;
	TFE3D *ele;
	double *weights, *xi, *eta, *zeta;
	double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
	double Z[MaxN_QuadPoints_3D];
	double AbsDetjk[MaxN_QuadPoints_3D];
	RefTrans3D RefTrans;
	double *Param[MaxN_QuadPoints_3D], *aux;
	double *Derivatives[MaxN_QuadPoints_3D];
	double *ExactVal[MaxN_QuadPoints_3D];
	double *AuxArray[MaxN_QuadPoints_3D];
	int *DOF, ActiveBound, DirichletBound, end, last;
	double **OrigFEValues, *Orig, value;
	double FEFunctValues[MaxN_BaseFunctions3D];
	int *GlobalNumbers, *BeginIndex;
	double LocError[4];
	double hK;
	bool *SecondDer;

	BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
	N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

	SecondDer = new bool[n_fespaces];
	for(i=0;i<n_fespaces;i++)
		SecondDer[i] = FALSE;

	N_Parameters = Aux->GetN_Parameters();
	aux = new double [MaxN_QuadPoints_3D*N_Parameters];
	for(j=0;j<MaxN_QuadPoints_3D;j++)
		Param[j] = aux + j*N_Parameters;

	aux = new double [MaxN_QuadPoints_3D*N_Derivatives];
	for(j=0;j<MaxN_QuadPoints_3D;j++)
		Derivatives[j] = aux + j*N_Derivatives;
  
	aux = new double [MaxN_QuadPoints_3D * 5];
	for(j=0;j<MaxN_QuadPoints_3D;j++)
		ExactVal[j] = aux + j*5;

  // 20 <= number of term
	aux = new double [MaxN_QuadPoints_3D*20]; 
	for(j=0;j<MaxN_QuadPoints_3D;j++)
		AuxArray[j] = aux + j*20;

	GlobalNumbers = FESpace3D->GetGlobalNumbers();
	BeginIndex = FESpace3D->GetBeginIndex();

	for(i=0;i<N_Errors;i++)
		errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
	Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
	N_Cells = Coll->GetN_Cells();
  
	TVertex *vertex;
	double loc_x, loc_y, loc_z, loc_r;
	for(i=0;i<N_Cells;i++)
	{
		cell = Coll->GetCell(i);
    
    // D.Sirch: Compute error only in area around corner/edge
		/*----------------------*/
		loc_r=2.0;
// 		loc_x=0.0;
// 		loc_y=0.0;
// 		loc_z=0.0;
		  //compute distance of the element to the edge
		for(k=0;k<4;k++)
		{
			vertex = cell->GetVertex(k);
			vertex->GetCoords(loc_x,loc_y, loc_z);
// 			loc_x+= 0.25*loc_x;
// 			loc_y+= 0.25*loc_y;
// 			loc_z+= 0.25*loc_z;
			
			if(sqrt(loc_x*loc_x + loc_y*loc_y) < loc_r)
			{
				loc_r = sqrt(loc_x*loc_x + loc_y*loc_y);
			}
		}
// 		loc_r = sqrt(loc_x*loc_x + loc_y*loc_y);
		if(loc_r > radius) continue;
    
		/*-----------------------*/

		hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
		memset(Used, 0, N_FEs3D*SizeOfInt);
		for(j=0;j<n_fespaces;j++)
		{
			CurrentElement = fespaces[j]->GetFE3D(i, cell);
			Used[CurrentElement] = 1;
		}

		N_LocalUsedElements = 0;
		memset(LocalUsedElements, 0, SizeOfInt*N_FEs3D);
		j = 0;
		for(k=0;k<N_FEs3D;k++)
			if(Used[k])
		{
			LocalUsedElements[j] = (FE3D)k;
			j++;
		}
		N_LocalUsedElements = j;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
		TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
				       Coll, cell, SecondDer,
				       N_Points, xi, eta, zeta, weights,
				       X, Y, Z, AbsDetjk);

		if(N_Parameters>0)
			Aux->GetParameters(N_Points, cell, i, xi, eta, zeta,
					   X, Y, Z, Param); 

    // calculate all needed derivatives of this FE function
		CurrentElement = FESpace3D->GetFE3D(i, cell);
		BaseFunct = BaseFuncts[CurrentElement];
		N_ = N_BaseFunct[CurrentElement];

		DOF = GlobalNumbers + BeginIndex[i];
		for(l=0;l<N_;l++)
			FEFunctValues[l] = Values[DOF[l]];

		for(k=0;k<N_Derivatives;k++)
		{
			OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct,
					NeededDerivatives[k]);
			for(j=0;j<N_Points;j++)
			{
				Orig = OrigFEValues[j];
				value = 0;
				for(l=0;l<N_;l++)
				{
					value += FEFunctValues[l] * Orig[l];
					//cout<<Orig[l]<<endl;
				} // endfor l
				
				// here we have to perform the projection for control
				if (value > upper) value = upper;
				else if (value < lower) value = lower;
				
				Derivatives[j][k] = value;
			} // endfor j
		} // endfor k
    
		errors[N_Errors] = 0.0;  //added D.Sirch
		for(j=0;j<N_Points;j++)
		{
			Exact(X[j], Y[j], Z[j], ExactVal[j]);
       // D.Sirch: computation of L^\inf-error
			if(fabs(*ExactVal[j] - Derivatives[j][0]) > errors[N_Errors])
			{
				errors[N_Errors] = fabs(*ExactVal[j] - Derivatives[j][0]);
			}
		}

		if(Coeff)
			Coeff(N_Points, X, Y, Z, Param, AuxArray);

		ErrorMeth(N_Points, X, Y, Z, AbsDetjk, weights, hK, Derivatives, 
			  ExactVal, AuxArray, LocError);

		for(j=0;j<N_Errors;j++)
			errors[j] += LocError[j];

	} // endfor i

	for(j=0;j<N_Errors;j++)
		errors[j] = sqrt(errors[j]);

	delete AuxArray[0];
	delete SecondDer;
	delete ExactVal[0];
	delete Derivatives[0];
	delete Param[0];
  
} // TFEFunction3D::GetErrorsOPTPDE

void TFEFunction3D::GetMeshCellParams(DoubleFunct3D *Exact, int N_Derivatives,
                              MultiIndex3D *NeededDerivatives,
                              int N_Errors, ErrorMethod3D *ErrorMeth, 
                              CoeffFct3D *Coeff, 
                              TAuxParam3D *Aux,
                              int n_fespaces, TFESpace3D **fespaces,
                              double *errors, double *cell_parameters)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs3D], *N_BaseFunct;
  TFESpace3D *fespace;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  TFE3D *ele;
  double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  RefTrans3D RefTrans;
  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool *SecondDer;

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  N_Parameters = Aux->GetN_Parameters();
  aux = new double [MaxN_QuadPoints_3D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Param[j] = aux + j*N_Parameters;

  aux = new double [MaxN_QuadPoints_3D*N_Derivatives];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux + j*N_Derivatives;
  
  aux = new double [MaxN_QuadPoints_3D * 5];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    ExactVal[j] = aux + j*5;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_3D*20]; 
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    AuxArray[j] = aux + j*20;

  GlobalNumbers = FESpace3D->GetGlobalNumbers();
  BeginIndex = FESpace3D->GetBeginIndex();

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
    memset(Used, 0, N_FEs3D*SizeOfInt);
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE3D(i, cell);
      Used[CurrentElement] = 1;
    }

    N_LocalUsedElements = 0;
    memset(LocalUsedElements, 0, SizeOfInt*N_FEs3D);
    j = 0;
    for(k=0;k<N_FEs3D;k++)
      if(Used[k])
      {
        LocalUsedElements[j] = (FE3D)k;
        j++;
      }
    N_LocalUsedElements = j;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, zeta, weights,
                           X, Y, Z, AbsDetjk);

    if(N_Parameters>0)
      Aux->GetParameters(N_Points, cell, i, xi, eta, zeta,
                         X, Y, Z, Param); 

    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace3D->GetFE3D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = Values[DOF[l]];

    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = TFEDatabase3D::GetOrigElementValues(BaseFunct,
                                      NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
        } // endfor l
        Derivatives[j][k] = value;
      } // endfor j
    } // endfor k

    for(j=0;j<N_Points;j++)
      Exact(X[j], Y[j], Z[j], ExactVal[j]);

    if(Coeff)
      Coeff(N_Points, X, Y, Z, Param, AuxArray);

    ErrorMeth(N_Points, X, Y, Z, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
    {
      errors[j] += LocError[j];
      cell_parameters[i + j *N_Cells] = LocError[j];
    }
  } // endfor i

  for(j=0;j<N_Errors;j++)
    errors[j] = sqrt(errors[j]);

  delete AuxArray[0];
  delete SecondDer;
  delete ExactVal[0];
  delete Derivatives[0];
  delete Param[0];
  
} // TFEFunction3D::GetErrors

/** determine the value of function and its first derivatives at
    the given point */

void TFEFunction3D::FindGradient(double x, double y, double z, double *values)
{
  int i,j,k, N_Cells;
  double xv, yv, zv, xi, eta, zeta;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FE_ID;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  int N_BaseFunct;
  double *uorig, *uxorig, *uyorig, *uzorig, *uref, *uxiref, *uetaref, *uzetaref;
  
  int *Numbers, N_Found;
  double u, ux, uy, uz;
  double val;
  int *GlobalNumbers, *BeginIndex;
  
  N_Found = 0;
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();

  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->PointInCell(x,y,z))
    {
      N_Found++;      
      // cout << "point found in cell: " << i << endl;
      FE_ID = FESpace3D->GetFE3D(i, cell);
      FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
      RefTrans = FE_Obj->GetRefTransID();

      // set cell for reference transformation
      TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

      // find local coordinates of the given point
      TFEDatabase3D::GetRefFromOrig(RefTrans, x, y, z, xi, eta, zeta);
      // cout << " xi: " << xi << endl;
      // cout << "eta: " << eta << endl;
      // cout << "zeta: " << zeta << endl;

      // get base function object
      bf = FE_Obj->GetBaseFunct3D();
      N_BaseFunct = bf->GetDimension();

      uorig = new double[N_BaseFunct];
      uxorig = new double[N_BaseFunct];
      uyorig = new double[N_BaseFunct];
      uzorig = new double[N_BaseFunct];

      uref = new double[N_BaseFunct];
      uxiref = new double[N_BaseFunct];
      uetaref = new double[N_BaseFunct];
      uzetaref = new double[N_BaseFunct];

      bf->GetDerivatives(D000, xi, eta, zeta, uref);
      bf->GetDerivatives(D100, xi, eta, zeta, uxiref);
      bf->GetDerivatives(D010, xi, eta, zeta, uetaref);
      bf->GetDerivatives(D001, xi, eta, zeta, uzetaref);

      TFEDatabase3D::GetOrigValues(RefTrans, xi, eta, zeta,
                     bf, Coll, (TGridCell *)cell,
                     uref, uxiref, uetaref, uzetaref,
                     uorig, uxorig, uyorig, uzorig);

      u = 0;
      ux = 0;
      uy = 0;
      uz = 0;
      Numbers = GlobalNumbers + BeginIndex[i];
      for(j=0;j<N_BaseFunct;j++)
      {
        val = Values[Numbers[j]];
        // cout << j << " " << val << endl;
        u  +=  uorig[j]*val;
        ux += uxorig[j]*val;
        uy += uyorig[j]*val;
        uz += uzorig[j]*val;
        // cout << " uorig[j]: " << uorig[j] << endl;
        // cout << " uxorig[j]: " << uxorig[j]  << endl;
        // cout << " uyorig[j]: " << uyorig[j] << endl;
        // cout << " uzorig[j]: " << uzorig[j] << endl;
      } 

      values[0] += u;
      values[1] += ux;
      values[2] += uy;
      values[3] += uz;

      delete uorig;
      delete uxorig;
      delete uyorig;
      delete uzorig;
      delete uref;
      delete uxiref;
      delete uetaref;
      delete uzetaref;

    } // endif
  } // endfor

  if(N_Found>0)
  {
    values[0] /= N_Found;
    values[1] /= N_Found;
    values[2] /= N_Found;
    values[3] /= N_Found;
    // cout << " values[0]: " << values[0] << endl;
    // cout << " values[1]: " << values[1] << endl;
    // cout << " values[2]: " << values[2] << endl;
    // cout << " values[3]: " << values[3] << endl;
  } 
} 

void TFEFunction3D::FindGradientLocal(TBaseCell *cell, int cell_no, 
                                      double x, double y, double z, 
                                      double *values)
{
  int j,k, N_Cells;
  double xv, yv, zv, xi, eta, zeta;
  FE3D FE_ID;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  int N_BaseFunct;
  double *uorig, *uxorig, *uyorig, *uzorig, *uref, *uxiref, *uetaref, *uzetaref;
  
  int *Numbers, N_Found;
  double u, ux, uy, uz;
  double val;
  int *GlobalNumbers, *BeginIndex;
  TCollection *Coll;
  
  Coll = FESpace3D->GetCollection();
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();

  FE_ID = FESpace3D->GetFE3D(cell_no, cell);
  FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
  RefTrans = FE_Obj->GetRefTransID();

  // set cell for reference transformation
  TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);
  
  // find local coordinates of the given point
  TFEDatabase3D::GetRefFromOrig(RefTrans, x, y, z, xi, eta, zeta);
  // cout << " xi: " << xi << endl;
  // cout << "eta: " << eta << endl;
  // cout << "zeta: " << zeta << endl;
  
  // get base function object
  bf = FE_Obj->GetBaseFunct3D();
  N_BaseFunct = bf->GetDimension();
  
  uorig = new double[N_BaseFunct];
  uxorig = new double[N_BaseFunct];
  uyorig = new double[N_BaseFunct];
  uzorig = new double[N_BaseFunct];
  
  uref = new double[N_BaseFunct];
  uxiref = new double[N_BaseFunct];
  uetaref = new double[N_BaseFunct];
  uzetaref = new double[N_BaseFunct];
  
  bf->GetDerivatives(D000, xi, eta, zeta, uref);
  bf->GetDerivatives(D100, xi, eta, zeta, uxiref);
  bf->GetDerivatives(D010, xi, eta, zeta, uetaref);
  bf->GetDerivatives(D001, xi, eta, zeta, uzetaref);
  
  TFEDatabase3D::GetOrigValues(RefTrans, xi, eta, zeta,
                 bf, Coll, (TGridCell *)cell,
                 uref, uxiref, uetaref, uzetaref,
                 uorig, uxorig, uyorig, uzorig);
  u = 0;
  ux = 0;
  uy = 0;
  uz = 0;
  Numbers = GlobalNumbers + BeginIndex[cell_no];
  for(j=0;j<N_BaseFunct;j++)
  {
    val = Values[Numbers[j]];
    // cout << j << " " << val << endl;
    u  +=  uorig[j]*val;
    ux += uxorig[j]*val;
    uy += uyorig[j]*val;
    uz += uzorig[j]*val;
    // cout << " uorig[j]: " << uorig[j] << endl;
    // cout << " uxorig[j]: " << uxorig[j]  << endl;
    // cout << " uyorig[j]: " << uyorig[j] << endl;
    // cout << " uzorig[j]: " << uzorig[j] << endl;
  } 
  values[0] = u;
  values[1] = ux;
  values[2] = uy;
  values[3] = uz;
  
  delete uorig;
  delete uxorig;
  delete uyorig;
  delete uzorig;
  delete uref;
  delete uxiref;
  delete uetaref;
  delete uzetaref;
  
}

void TFEFunction3D::FindValueLocal(TBaseCell *cell, int cell_no, 
                                      double x, double y, double z, 
                                      double *values)
{
  int j,k, N_Cells;
  double xv, yv, zv, xi, eta, zeta;
  FE3D FE_ID;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  int N_BaseFunct;
  double *uorig, *uxorig, *uyorig, *uzorig, *uref, *uxiref, *uetaref, *uzetaref;
  
  int *Numbers, N_Found;
  double u;
  double val;
  int *GlobalNumbers, *BeginIndex;

  TCollection *Coll;

  Coll = FESpace3D->GetCollection();
  
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();

  FE_ID = FESpace3D->GetFE3D(cell_no, cell);
  FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
  RefTrans = FE_Obj->GetRefTransID();

  // set cell for reference transformation
  TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);
  
  // find local coordinates of the given point
  TFEDatabase3D::GetRefFromOrig(RefTrans, x, y, z, xi, eta, zeta);
  // cout << " xi: " << xi << endl;
  // cout << "eta: " << eta << endl;
  // cout << "zeta: " << zeta << endl;
  
  // get base function object
  bf = FE_Obj->GetBaseFunct3D();
  N_BaseFunct = bf->GetDimension();
  
  uorig = new double[N_BaseFunct];
  uxorig = new double[N_BaseFunct];
  uyorig = new double[N_BaseFunct];
  uzorig = new double[N_BaseFunct];
  
  uref = new double[N_BaseFunct];
  uxiref = new double[N_BaseFunct];
  uetaref = new double[N_BaseFunct];
  uzetaref = new double[N_BaseFunct];
  
  bf->GetDerivatives(D000, xi, eta, zeta, uref);
  bf->GetDerivatives(D100, xi, eta, zeta, uxiref);
  bf->GetDerivatives(D010, xi, eta, zeta, uetaref);
  bf->GetDerivatives(D001, xi, eta, zeta, uzetaref);
  
  TFEDatabase3D::GetOrigValues(RefTrans, xi, eta, zeta,
                 bf, Coll, cell,
                 uref, uxiref, uetaref, uzetaref,
                 uorig, uxorig, uyorig, uzorig);
  u = 0;
  Numbers = GlobalNumbers + BeginIndex[cell_no];
  for(j=0;j<N_BaseFunct;j++)
  {
    val = Values[Numbers[j]];
    u  +=  uorig[j]*val;
  } 
  values[0] = u;
  
  delete uorig;
  delete uxorig;
  delete uyorig;
  delete uzorig;
  delete uref;
  delete uxiref;
  delete uetaref;
  delete uzetaref;
  
}

/** calculate the interpolation of an exact function */
void TFEFunction3D::Interpolate(DoubleFunct3D *Exact)
{
  int i,j,k,l;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId;
  TFE3D *Element;
  BaseFunct3D BF;
  TNodalFunctional3D *nf;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_, N_Points;
  double s, *xi, *eta, *zeta;
  double Val[MaxN_BaseFunctions3D];
  double OutVal[MaxN_BaseFunctions3D];
  int *DOF, Index;
  RefTrans3D F_K;
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AbsDetjk[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  double FctVal[5];
  int PolynomialDegree, ApproxOrder;
  QuadFormula3D QuadFormula;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Faces;
  BF3DRefElements RefElement;
  RefTrans3D RefTrans, *RefTransArray;

  // begin code
  
  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();
  N_DOFs = FESpace3D->GetN_DegreesOfFreedom();

  memset(Values, 0, SizeOfDouble*N_DOFs);
  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace3D->GetFE3D(i, cell);
    Element = TFEDatabase3D::GetFE3D(FEId);
    nf = Element->GetNodalFunctional3D();
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    N_LocalDOFs = Element->GetN_DOF();

    PolynomialDegree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);
    ApproxOrder = TFEDatabase3D::GetAccuracyFromFE3D(FEId);

    RefElement = Element->GetBaseFunct3D()->GetRefElement();
    switch(RefElement)
    {
      case BFUnitHexahedron:
        QuadFormula = TFEDatabase3D::GetQFHexaFromDegree
                         (3*PolynomialDegree);
        N_Faces = 6;
      break;

      case BFUnitTetrahedron:
        QuadFormula = TFEDatabase3D::GetQFTetraFromDegree
                         (3*PolynomialDegree);
        N_Faces = 4;
      break;
    }

    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Faces;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryFace)
        {
          bdtype = ((TBoundFace *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Plane)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoBoundFace)
          IsIsoparametric = TRUE;
      }
    } // endif 
  
    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          RefTrans = HexaIsoparametric;
        break;
  
        case BFUnitTetrahedron:
          RefTrans = TetraIsoparametric;
        break;
      }
    } // endif IsIsoparametric
    // cout << "IsIsoparametric: " << IsIsoparametric << endl;
  
    switch(RefTrans)
    {
      case HexaAffin:
        rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        F_K = HexaAffin;
        break;
      case HexaTrilinear:
        rt = TFEDatabase3D::GetRefTrans3D(HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        F_K = HexaTrilinear;
        break;
      case HexaIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((THexaIsoparametric *)rt)->SetCell(cell);
        F_K = HexaIsoparametric;
        break;
      case TetraAffin:
        rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        F_K = TetraAffin;
        break;
      case TetraIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        F_K = TetraIsoparametric;
        break;
    }
    TFEDatabase3D::GetOrigFromRef(F_K, N_Points, xi, eta, zeta,
                                X, Y, Z, AbsDetjk);

    // cout << "----------------" << endl;
    for(j=0;j<N_Points;j++)
    {
      // cout << j << " ";
      // cout << "ref: " << xi[j] << " " << eta[j] << " " << zeta[j] << endl;
      // cout << "Ori: " << X[j] << " " << Y[j] << " " << Z[j] << endl;
      Exact(X[j], Y[j], Z[j], FctVal);
      // cout << FctVal[0] << endl;
      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues,FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];
  }
}

/** calculate the super-convergence interpolation of an exact function */
void TFEFunction3D::InterpolateSuper(DoubleFunct3D *Exact)
{
  int i,j,k,l;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId;
  TFE3D *Element;
  BaseFunct3D BF;
  TNodalFunctional3D *nf;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_, N_Points;
  double s, *xi, *eta, *zeta;
  double Val[MaxN_BaseFunctions3D];
  double OutVal[MaxN_BaseFunctions3D];
  int *DOF, Index;
  RefTrans3D F_K;
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AbsDetjk[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  double FctVal[5];
  int PolynomialDegree, ApproxOrder;
  QuadFormula3D QuadFormula;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Faces;
  BF3DRefElements RefElement;
  RefTrans3D RefTrans, *RefTransArray;

  // begin code
  
  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();
  N_DOFs = FESpace3D->GetN_DegreesOfFreedom();

  memset(Values, 0, SizeOfDouble*N_DOFs);
  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace3D->GetFE3D(i, cell);
    Element = TFEDatabase3D::GetFE3D(FEId);
    nf = Element->GetNodalFunctional3D();
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    switch(nf->GetID())
    {
      case NF_C_H_Q2_3D:
	nf = TFEDatabase3D::GetNodalFunctional3D(NF_S_H_Q2_3D);
      break;
    }
    N_LocalDOFs = Element->GetN_DOF();

    PolynomialDegree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);
    ApproxOrder = TFEDatabase3D::GetAccuracyFromFE3D(FEId);

    RefElement = Element->GetBaseFunct3D()->GetRefElement();
    switch(RefElement)
    {
      case BFUnitHexahedron:
        QuadFormula = TFEDatabase3D::GetQFHexaFromDegree
                         (3*PolynomialDegree);
        N_Faces = 6;
      break;

      case BFUnitTetrahedron:
        QuadFormula = TFEDatabase3D::GetQFTetraFromDegree
                         (3*PolynomialDegree);
        N_Faces = 4;
      break;
    }

    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Faces;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryFace)
        {
          bdtype = ((TBoundFace *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Plane)
            IsIsoparametric = TRUE;
        }
      }
    } // endif 
  
    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          RefTrans = HexaIsoparametric;
        break;
  
        case BFUnitTetrahedron:
          RefTrans = TetraIsoparametric;
        break;
      }
    } // endif IsIsoparametric
    // cout << "IsIsoparametric: " << IsIsoparametric << endl;
  
    switch(RefTrans)
    {
      case HexaAffin:
        rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        F_K = HexaAffin;
        break;
      case HexaTrilinear:
        rt = TFEDatabase3D::GetRefTrans3D(HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        F_K = HexaTrilinear;
        break;
      case HexaIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((THexaIsoparametric *)rt)->SetCell(cell);
        F_K = HexaIsoparametric;
        break;
      case TetraAffin:
        rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        F_K = TetraAffin;
        break;
      case TetraIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        F_K = TetraIsoparametric;
        break;
    }
    TFEDatabase3D::GetOrigFromRef(F_K, N_Points, xi, eta, zeta,
                                X, Y, Z, AbsDetjk);

    // cout << "----------------" << endl;
    for(j=0;j<N_Points;j++)
    {
      // cout << j << " ";
      // cout << "ref: " << xi[j] << " " << eta[j] << " " << zeta[j] << endl;
      // cout << "Ori: " << X[j] << " " << Y[j] << " " << Z[j] << endl;
      Exact(X[j], Y[j], Z[j], FctVal);
      // cout << FctVal[0] << endl;
      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues,FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];
  }
}



/** calculate the interpolation of an exact function */
void TFEFunction3D::Interpolate(int N_Coord, double *Coords, int N_AuxFeFcts,  TFEFunction3D **AuxFeFcts, DoubleFunctND *Exact)
{
  int i,j, jj, k,l;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId;
  TFE3D *Element;
  BaseFunct3D BF;
  TNodalFunctional3D *nf;
  int N_Cells, N_CoordAll;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_, N_Points;
  double s, *xi, *eta, *zeta;
  double Val[MaxN_BaseFunctions3D];
  double OutVal[MaxN_BaseFunctions3D];
  int *DOF, Index;
  RefTrans3D F_K;
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AbsDetjk[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  double *FctVal, *coords, **auxFeValues;
  double values[4];
    
  int PolynomialDegree, ApproxOrder;
  QuadFormula3D QuadFormula;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Faces;
  BF3DRefElements RefElement;
  RefTrans3D RefTrans, *RefTransArray;

  // begin code
  
  coords = new double[N_Coord+3];
  auxFeValues = new double*[N_AuxFeFcts];  
  if(N_AuxFeFcts<5)
   { FctVal = new double[5];}
  else
   { FctVal = new double[N_AuxFeFcts]; }
  
  for(i=0;i<N_Coord;i++)
   coords[i] = Coords[i];
  
  for(i=0;i<N_AuxFeFcts;i++)
   auxFeValues[i] = new double[4];
  
  N_CoordAll = N_Coord+3; // this is 3D functon  
    
    
  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();
  N_DOFs = FESpace3D->GetN_DegreesOfFreedom();

  memset(Values, 0, SizeOfDouble*N_DOFs);
  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace3D->GetFE3D(i, cell);
    Element = TFEDatabase3D::GetFE3D(FEId);
    nf = Element->GetNodalFunctional3D();
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    N_LocalDOFs = Element->GetN_DOF();

    PolynomialDegree = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);
    ApproxOrder = TFEDatabase3D::GetAccuracyFromFE3D(FEId);

    RefElement = Element->GetBaseFunct3D()->GetRefElement();
    switch(RefElement)
    {
      case BFUnitHexahedron:
        QuadFormula = TFEDatabase3D::GetQFHexaFromDegree
                         (3*PolynomialDegree);
        N_Faces = 6;
      break;

      case BFUnitTetrahedron:
        QuadFormula = TFEDatabase3D::GetQFTetraFromDegree
                         (3*PolynomialDegree);
        N_Faces = 4;
      break;
    }

    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Faces;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryFace)
        {
          bdtype = ((TBoundFace *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Plane)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoBoundFace)
          IsIsoparametric = TRUE;
      }
    } // endif 
  
    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          RefTrans = HexaIsoparametric;
        break;
  
        case BFUnitTetrahedron:
          RefTrans = TetraIsoparametric;
        break;
      }
    } // endif IsIsoparametric
    // cout << "IsIsoparametric: " << IsIsoparametric << endl;
  
    switch(RefTrans)
    {
      case HexaAffin:
        rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        F_K = HexaAffin;
        break;
      case HexaTrilinear:
        rt = TFEDatabase3D::GetRefTrans3D(HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        F_K = HexaTrilinear;
        break;
      case HexaIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((THexaIsoparametric *)rt)->SetCell(cell);
        F_K = HexaIsoparametric;
        break;
      case TetraAffin:
        rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        F_K = TetraAffin;
        break;
      case TetraIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        F_K = TetraIsoparametric;
        break;
    }
    TFEDatabase3D::GetOrigFromRef(F_K, N_Points, xi, eta, zeta,
                                X, Y, Z, AbsDetjk);

    // cout << "----------------" << endl;

    for(j=0;j<N_Points;j++)
    {
      // cout << j << " ";
      // cout << "ref: " << xi[j] << " " << eta[j] << " " << zeta[j] << endl;
      // cout << "Ori: " << X[j] << " " << Y[j] << " " << Z[j] << endl;
      
      //set the coordinate
       coords[N_Coord  ] = X[j];
       coords[N_Coord+1] = Y[j];      
       coords[N_Coord+2] = Z[j];      
      
      //get auxFevalues
      for(jj=0;jj<N_AuxFeFcts;jj++)
       {
        AuxFeFcts[jj]->FindValueLocal(cell, i, X[j], Y[j], Z[j], values);
        FctVal[jj] = values[0];
       }
      Exact(N_CoordAll, Coords, FctVal);      
     // Exact(X[j], Y[j], Z[j], FctVal);
      // cout << FctVal[0] << endl;
      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues,FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];
  }
}






