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
#include <AllRefTrans3D.h>
// #include <NodalFunctional3D.h>
#include <GridCell.h>

#include <stdlib.h>
#include <InterfaceJoint.h>
#include <BdPlane.h>
#include <LinAlg.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

void OnlyDirichlet(double, double, double, BoundCond &cond)
{
	cond = DIRICHLET;
}

TFEFunction3D::TFEFunction3D() :
    Name("dummy_fe_fct_3d"), Description("dummy_fe_fct_3d")
{
  FESpace3D=nullptr;
  Values=nullptr;
  Length=0;
}

/** constructor with vector initialization */
TFEFunction3D::TFEFunction3D(std::shared_ptr<const TFESpace3D> fespace3D,
                             const std::string& name,
                             const std::string& description, double *values, int length)
: Name(name), Description(description)
{
  FESpace3D=fespace3D;

  Values=values;

  Length=length;
}

TFEFunction3D& TFEFunction3D::operator=(const TFEFunction3D& other)
{
  this->Name        = other.Name;
  this->Description = other.Description;
  this->FESpace3D   = other.FESpace3D;
  this->Values      = other.Values;
  this->Length      = other.Length;

  return *this;
}

TFEFunction3D::~TFEFunction3D()
{

}





/** calculate errors to given function */
void TFEFunction3D::GetErrors(DoubleFunct3D *Exact, int N_Derivatives,
                              MultiIndex3D *NeededDerivatives,
                              int N_Errors, ErrorMethod *ErrorMeth, 
                              CoeffFct3D Coeff, 
                              TAuxParam3D *Aux,
                              int n_fespaces, const TFESpace3D **fespaces,
                              double *errors) const
{
  int i,j,k,l, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs3D], *N_BaseFunct;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  const double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF;
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions3D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4];
  double hK;
  bool *SecondDer;

  
#ifdef _MPI
   int ID, rank;
   
   MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);    
#endif
      
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = false;

  N_Parameters = Aux->GetN_Parameters();
  
  if(N_Parameters==0)
   aux = nullptr;
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
 
    for(j=0;j<N_Points;j++)
    {
      Exact(X[j], Y[j], Z[j], ExactVal[j]);
    }

    if(Coeff)
      Coeff(N_Points, X, Y, Z, Param, AuxArray);

    ErrorMeth(N_Points, {{X, Y, Z}}, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors-1;j++)
      errors[j] += LocError[j];
    // L-inf error
    if(errors[N_Errors-1] < LocError[N_Errors-1])
      errors[N_Errors-1] = LocError[N_Errors-1];
  } // endfor i

#ifndef _MPI // sqrt(errors[j]) in the main programm after collecting error from all subdomains
  for(j=0;j<N_Errors-1;j++)
    errors[j] = sqrt(errors[j]);
#endif
  
  delete [] AuxArray[0];
  delete [] SecondDer;
  delete [] ExactVal[0];
  delete [] Derivatives[0];
  
  if(Param[0])
   delete [] Param[0];
  
} // TFEFunction3D::GetErrors

void TFEFunction3D::GetErrorsForVectorValuedFunction(
    DoubleFunct3D * const * const Exact, ErrorMethod * const ErrMeth,
  double * const errors)
{
  // write zeros into the array "errors"
  memset(errors,0,3*SizeOfDouble);
  TCollection *Coll = FESpace3D->GetCollection(); 
  int N_Cells = Coll->GetN_Cells(); // number of cells
  
  // second derivatives are not needed
  bool Needs2ndDer = false;
  
  // for loop over cells
  for(int i =0; i < N_Cells; i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    FE3D CurrentElement = FESpace3D->GetFE3D(i, cell);
    // number of basis functions
    int N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D(CurrentElement);
    BaseFunct3D BaseFunct = 
                      TFEDatabase3D::GetBaseFunct3D_IDFromFE3D(CurrentElement);
    int * DOF = FESpace3D->GetGlobalDOF(i);
    RefTrans3D RefTrans = 
                       TFEDatabase3D::GetRefTrans3D_IDFromFE3D(CurrentElement);
    TRefTrans3D *F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);
    // get quadrature formula
    QuadFormula3D qf_id = BaryCenterTetra; //to avoid uninit warning
    switch(RefTrans)
    {
      case TetraAffin:
        //((TetraAffin*)F_K)->SetCell(cell);
        qf_id = TFEDatabase3D::GetQFTetraFromDegree(
            TDatabase::ParamDB->INPUT_QUAD_RULE);
        break;
      case HexaAffin:
        //((HexaAffin*)F_K)->SetCell(cell);
        qf_id = TFEDatabase3D::GetQFHexaFromDegree(
            TDatabase::ParamDB->INPUT_QUAD_RULE);
        break;
      case HexaTrilinear:
        //((HexaTrilinear*)F_K)->SetCell(cell);
        qf_id = TFEDatabase3D::GetQFHexaFromDegree(
            TDatabase::ParamDB->INPUT_QUAD_RULE);
        break;
      default:
        ErrMsg("No such reference transformation supported!");
        break;
    }
    // quadrature formula
    TQuadFormula3D *qf = TFEDatabase3D::GetQuadFormula3D(qf_id);
    int N_Points; // number of quadrature points in one cell
    const double *xi, *eta, *zeta; // quadrature points in reference cell
    const double *weights; // quadrature weights
    qf->GetFormulaData(N_Points, weights, xi, eta, zeta);
    TFEDatabase3D::GetBaseFunct3D(BaseFunct)->MakeRefElementData(qf_id);
    
    // modulus of determinant of reference transformation for 
    // each quadrature point
    double * AbsDetjk = new double[N_Points];
    double * X = new double[N_Points];
    double * Y = new double[N_Points];
    double * Z = new double[N_Points];
    // get quadrature coordinates on original cell (AbsDetjk is filled, too)
    TFEDatabase3D::GetOrigFromRef(RefTrans, N_Points, xi, eta ,zeta, X, Y, Z, 
                                  AbsDetjk);
    switch(RefTrans)
    {
      case TetraAffin:
        //((TetraAffin*)F_K)->GetOrigFromRef(N_Points, xi, eta ,zeta, X, Y, Z, 
        //                                   AbsDetjk);
        // fill arrays, so call of GetOrigElementValues is now possible
        ((TTetraAffin*)F_K)->GetOrigValues(1, &BaseFunct, N_Points, xi, eta,
                                          zeta, qf_id, &Needs2ndDer);
        break;
      case HexaAffin:
        //((HexaAffin*)F_K)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z,
        //                                  AbsDetjk);
        ((THexaAffin*)F_K)->GetOrigValues(1, &BaseFunct, N_Points, xi, eta,
                                         zeta, qf_id, &Needs2ndDer);
        break;
      case HexaTrilinear:
        //((HexaTrilinear*)F_K)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z,
        //                                      AbsDetjk);
        ((THexaTrilinear*)F_K)->GetOrigValues(1, &BaseFunct, N_Points, xi, eta, 
                                             zeta, qf_id, &Needs2ndDer);
        break;
      default:
        ErrMsg(" No such reference transformation supported!");
        break;
    }
    // fill ExactVal
    double ** ExactVal = new double*[N_Points];
    for(int j=0;j<N_Points;j++)
    {
      // 5 values for each of the three components:
      // (value, 3 derivatives, Laplace)
      ExactVal[j] = new double[15]; 
      Exact[0](X[j], Y[j], Z[j], ExactVal[j]     ); // x-component
      Exact[1](X[j], Y[j], Z[j], ExactVal[j] +  5); // y-component
      Exact[2](X[j], Y[j], Z[j], ExactVal[j] + 10); // z-component
    }
    // will store values of this FEFunction and its first derivatives at all 
    // quadrature points
    double ** Derivatives = new double*[N_Points];
    for(int j=0;j<N_Points;j++)
    {
      // the 12 means values and 3 first derivatives for all three components
      Derivatives[j] = new double[12];
    }
    // Get the function values of this FE-function at the local dofs.
    // some local dofs get a negative sign according to global orientation of 
    // the normals
    double * FEFunctValues = new double[N_BaseFuncts];
    for(int l=0;l<N_BaseFuncts;l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      // revert sign at some inner edges
      // edge number where a local DOF is on
      int face = TFEDatabase3D::GetFE3D(CurrentElement)
             ->GetFEDesc3D()->GetJointOfThisDOF(l);
      if(face != -1) // edge==-1 means inner dof
      {
        FEFunctValues[l] *= cell->GetNormalOrientation(face);
      }
    }
    
    // these arrays were created in GetOrigValues called earlier
    double **AllOrigValues[4];
    AllOrigValues[0] = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
    AllOrigValues[1] = TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
    AllOrigValues[2] = TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
    AllOrigValues[3] = TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
    
    // loop over all needed derivatives
    for(int k = 0; k < 4; k++)
    {
      // loop over all quadrature points
      for(int j = 0; j < N_Points; j++) 
      {
        double value_x = 0;
        double value_y = 0;
        double value_z = 0;
        // loop over all basis functions
        for(int l=0;l<N_BaseFuncts;l++)
        {
          value_x += FEFunctValues[l] * AllOrigValues[k][j][l                ];
          value_y += FEFunctValues[l] * AllOrigValues[k][j][l+   N_BaseFuncts];
          value_z += FEFunctValues[l] * AllOrigValues[k][j][l+ 2*N_BaseFuncts];
        }
        Derivatives[j][k    ] = value_x;
        Derivatives[j][k + 4] = value_y;
        Derivatives[j][k + 8] = value_z;
      }
    }
    // cell diameter, we set it one here since it is not needed in 
    // ErrorMeth=L2DivH1Errors
    double hK = 1;
    double LocError[3]; // L^2 error in value, divergence and first derivative
    ErrMeth(N_Points, {{X, Y, Z}}, AbsDetjk, weights, hK, Derivatives, ExactVal,
            nullptr, LocError);
    for(int j=0;j<3;j++) 
    {
      errors[j] += LocError[j];
    }
    // delete everything which was created with "new" within this loop
    // otherwise one would get (many) memory leaks
    delete [] AbsDetjk; AbsDetjk = nullptr;
    delete [] X;        X        = nullptr;
    delete [] Y;        Y        = nullptr;
    delete [] Z;        Z        = nullptr;
    for (int j=0; j<N_Points; j++)
    {
      delete [] ExactVal[j];    ExactVal[j] = nullptr;
      delete [] Derivatives[j]; Derivatives[j] = nullptr;
    }
    delete [] ExactVal;      ExactVal = nullptr;
    delete [] Derivatives;   Derivatives = nullptr;
    delete [] FEFunctValues; FEFunctValues = nullptr;
    
  } // end loop over all cells
  
  for(int j=0;j<3;j++)  
    errors[j] = sqrt(errors[j]);
}


void TFEFunction3D::GetMeshCellParams(DoubleFunct3D *Exact, int N_Derivatives,
                              MultiIndex3D *NeededDerivatives,
                              int N_Errors, ErrorMethod *ErrorMeth, 
                              CoeffFct3D Coeff, 
                              TAuxParam3D *Aux,
                              int n_fespaces, const TFESpace3D **fespaces,
                              double *errors, double *cell_parameters)
{
  int i,j,k,l, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs3D], *N_BaseFunct;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  const double *weights, *xi, *eta, *zeta;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double *Param[MaxN_QuadPoints_3D], *aux;
  double *Derivatives[MaxN_QuadPoints_3D];
  double *ExactVal[MaxN_QuadPoints_3D];
  double *AuxArray[MaxN_QuadPoints_3D];
  int *DOF;
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
    SecondDer[i] = false;

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

    ErrorMeth(N_Points, {{X, Y, Z}}, AbsDetjk, weights, hK, Derivatives, 
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


bool TFEFunction3D::FindGradient(double x, double y, double z,
                                 std::vector<double>& values) const
{

  if(values.size() != 4)
  {
    ErrThrow("TFEFunction3D::FindGradient expects vector of size 4 !=", values.size());
  }

  int i,j, N_Cells;
  double xi, eta, zeta;
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
  values.at(0) = 0;
  values.at(1) = 0;
  values.at(2) = 0;
  values.at(3) = 0;
  
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
      FE_ID = FESpace3D->GetFE3D(i, cell);
      FE_Obj = TFEDatabase3D::GetFE3D(FE_ID);
      RefTrans = FE_Obj->GetRefTransID();

      // set cell for reference transformation
      TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

      // find local coordinates of the given point
      TFEDatabase3D::GetRefFromOrig(RefTrans, x, y, z, xi, eta, zeta);

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
      } 

      values.at(0) += u;
      values.at(1) += ux;
      values.at(2) += uy;
      values.at(3) += uz;

      delete[] uorig;
      delete[] uxorig;
      delete[] uyorig;
      delete[] uzorig;
      delete[] uref;
      delete[] uxiref;
      delete[] uetaref;
      delete[] uzetaref;

    } // endif
  } // endfor

  if(N_Found>0)
  {
    values.at(0) /= N_Found;
    values.at(1) /= N_Found;
    values.at(2) /= N_Found;
    values.at(3) /= N_Found;
    return true;
  }
  return false;
} 

void TFEFunction3D::FindGradientLocal(const TBaseCell *cell, int cell_no, 
                                      double x, double y, double z, 
                                      double *values)
{
  int j;
  double xi, eta, zeta;
  FE3D FE_ID;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  int N_BaseFunct;
  double *uorig, *uxorig, *uyorig, *uzorig, *uref, *uxiref, *uetaref, *uzetaref;
  
  int *Numbers;
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

void TFEFunction3D::FindValueLocal(const TBaseCell *cell, int cell_no, 
                                   double x, double y, double z, 
                                   double *values) const
{
  double xi, eta, zeta;
  FE3D FE_ID;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  int N_BaseFunct;
  double *uorig, *uxorig, *uyorig, *uzorig, *uref, *uxiref, *uetaref, *uzetaref;
  
  TCollection *Coll;

  Coll = FESpace3D->GetCollection();
  
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
  int BaseVectDim = bf->GetBaseVectDim(); // either 1 or 3
  
  uorig = new double[N_BaseFunct*BaseVectDim];
  uxorig = new double[N_BaseFunct*BaseVectDim];
  uyorig = new double[N_BaseFunct*BaseVectDim];
  uzorig = new double[N_BaseFunct*BaseVectDim];
  
  uref = new double[N_BaseFunct*BaseVectDim];
  uxiref = new double[N_BaseFunct*BaseVectDim];
  uetaref = new double[N_BaseFunct*BaseVectDim];
  uzetaref = new double[N_BaseFunct*BaseVectDim];
  
  bf->GetDerivatives(D000, xi, eta, zeta, uref);
  bf->GetDerivatives(D100, xi, eta, zeta, uxiref);
  bf->GetDerivatives(D010, xi, eta, zeta, uetaref);
  bf->GetDerivatives(D001, xi, eta, zeta, uzetaref);
  
  TFEDatabase3D::GetOrigValues(RefTrans, xi, eta, zeta,
                 bf, Coll, cell,
                 uref, uxiref, uetaref, uzetaref,
                 uorig, uxorig, uyorig, uzorig);
  // for vector valued basis functions (e.g. Raviart-Thomas elements), some 
  // signs must be changed
  if(BaseVectDim>1)
  {
    for (int i = 0; i < BaseVectDim; i++)
    {
      for(int j = 0 ; j < N_BaseFunct; j++)
      {
        int face = FE_Obj->GetFEDesc3D()->GetJointOfThisDOF(j);
        if(face!=-1) // face ==-1 means inner dof
          uorig[j+i*N_BaseFunct] *= cell->GetNormalOrientation(face);
      }
    }
  }
  
  int *Numbers = FESpace3D->GetGlobalDOF(cell_no);
  for (int i = 0; i < BaseVectDim; i++)
  {
    double u = 0;
    for(int j = 0; j < N_BaseFunct; j++)
    {
      double val = Values[Numbers[j]];
      u  +=  uorig[j + i*N_BaseFunct]*val;
    }
    values[i] = u;
  }

  
  delete [] uorig;
  delete [] uxorig;
  delete [] uyorig;
  delete [] uzorig;
  delete [] uref;
  delete [] uxiref;
  delete [] uetaref;
  delete [] uzetaref;
  
}

/** calculate the interpolation of an exact function */
void TFEFunction3D::Interpolate(DoubleFunct3D *Exact)
{
  int i,j;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId;
  TFE3D *Element;
  TNodalFunctional3D *nf;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  const double *xi, *eta, *zeta;
  int *DOF;
  RefTrans3D F_K = TetraAffin; //avoid uninit warning
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AbsDetjk[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  double FctVal[5];
  int PolynomialDegree, ApproxOrder;
  QuadFormula3D QuadFormula = BaryCenterTetra;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Faces =0;
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

    IsIsoparametric = false;
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
            IsIsoparametric = true;
        }
        if(jointtype == IsoBoundFace)
          IsIsoparametric = true;
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
    {
      Values[DOF[j]] = FunctionalValues[j];
    }
  }
}

void TFEFunction3D::Interpolate_vector_valued_function(
    const std::vector<DoubleFunct3D*>& Exact)
{
  if(Exact.size() != 3)
  {
    ErrMsg("TFEFunction3D::Interpolate_vector_valued_function You have to " <<
           "provide three functions describing the exact solution " <<
           Exact.size());
    exit(1);
  }

  // reset function values
  int N_DOFs = this->FESpace3D->GetN_DegreesOfFreedom();
  memset(this->Values, 0, SizeOfDouble*N_DOFs);


  // loop over all cells
  TCollection * Coll = this->FESpace3D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  for(int i = 0; i < N_Cells; i++)
  {
    TBaseCell* cell = Coll->GetCell(i);
    FE3D FEId = FESpace3D->GetFE3D(i, cell);
    TFE3D *Element = TFEDatabase3D::GetFE3D(FEId);
    TNodalFunctional3D * nf = Element->GetNodalFunctional3D();
    int n_functionals = nf->n_functionals();
    // number of points needed to evaluate nodal functionals
    int N_Points = nf->n_points();
    // coordinates of points on reference element to evaluate nodal
    // functionals
    const double *xi, *eta, *zeta;
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    int N_LocalDOFs = Element->GetN_DOF();

    RefTrans3D RefTrans = TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FEId);

    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      ErrMsg("Interpolate_vector_valued_function with isoparametric elements"
             << " not yet implemented");
      exit(1);
    }

    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

    // coordinates of points needed to evaluate the nodal functionals on
    // cell in grid.
    double X[N_Points], Y[N_Points], Z[N_Points], AbsDetjk[N_Points];
    TFEDatabase3D::GetOrigFromRef(RefTrans, N_Points, xi, eta, zeta,
                                  X, Y, Z, AbsDetjk);

    for(int markus = 0; markus < N_Points; markus++)
    {
        //OutPut("AbsDetjk(" << markus << "): " << AbsDetjk[markus] << endl);
    }

    double PointValues[Exact.size()*N_Points]; // Exact.size() == 3
    for(int j = 0; j < N_Points; j++)
    {
/*       cout << j << " ";
       cout << "ref: " << xi[j] << " " << eta[j] << " " << zeta[j] << endl;
       cout << "Ori: " << X[j] << " " << Y[j] << " " << Z[j] << endl;*/
      for(unsigned int dim = 0; dim < Exact.size(); dim++) // Exact.size() == 3
      {
        double FctVal[5];
        Exact[dim](X[j], Y[j], Z[j], FctVal);

        // cout << FctVal[0] << endl;
        PointValues[j + dim*N_Points] = FctVal[0];
      }
    }

    double FunctionalValues[n_functionals];
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues, FunctionalValues);

    int *DOF = this->GetFESpace3D()->GetGlobalDOF(i);

    for(int j = 0; j < N_LocalDOFs; j++)
    {
      this->Values[DOF[j]] = FunctionalValues[j];
      // consider sign on an inner face
      int face = TFEDatabase3D::GetFE3D(FEId)
                   ->GetFEDesc3D()->GetJointOfThisDOF(j);
      if(face != -1) // face==-1 means inner dof
      {
        this->Values[DOF[j]] *= cell->GetNormalOrientation(face);
      }
    }
  }
}


/** calculate the super-convergence interpolation of an exact function */
void TFEFunction3D::InterpolateSuper(DoubleFunct3D *Exact)
{
  int i,j;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId;
  TFE3D *Element;
  TNodalFunctional3D *nf;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  const double *xi, *eta, *zeta;
  int *DOF;
  RefTrans3D F_K = TetraAffin;
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AbsDetjk[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  double FctVal[5];
  int PolynomialDegree, ApproxOrder;
  QuadFormula3D QuadFormula = BaryCenterTetra;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Faces = 0;
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
	  default:
	    
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

    IsIsoparametric = false;
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
            IsIsoparametric = true;
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
  int i,j, jj;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId;
  TFE3D *Element;
  TNodalFunctional3D *nf;
  int N_Cells, N_CoordAll;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  const double *xi, *eta, *zeta;
  int *DOF;
  RefTrans3D F_K = TetraAffin; //to avoid uninit warning
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AbsDetjk[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];
  double *FctVal, *coords, **auxFeValues;
  double values[4];
    
  int PolynomialDegree, ApproxOrder;
  QuadFormula3D QuadFormula = BaryCenterTetra;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Faces = 0;
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

    IsIsoparametric = false;
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
            IsIsoparametric = true;
        }
        if(jointtype == IsoBoundFace)
          IsIsoparametric = true;
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

void TFEFunction3D::add(AnalyticFunction)
{
  ErrThrow("TFEFunction3D::add is not yet implemented");
}


/** compute integral and measure */
void TFEFunction3D::compute_integral_and_measure(double& integral,
                                                 double& measure) const
{

  TCollection *coll = FESpace3D->GetCollection();

  integral = 0.0; // variable to store integral value of this TFEFunction3D
  measure = 0.0; // variable to store the measure of the domain

  // loop over all cells, find out integral value of this FEFunction3D and the`
  // measure of its domain
  const int n_cells = coll->GetN_Cells();
  for(int i = 0; i < n_cells; i++)
  {
    TBaseCell *cell = coll->GetCell(i); // current cell
#ifdef _MPI // skip halo cells
    if (cell->IsHaloCell())
    {
      continue;
    }
#endif
    FE3D feID = FESpace3D->GetFE3D(i, cell); // id of finite element

    // calculate values on original element (i.e. prepare reference
    // transformation)
    bool SecondDer = false; // defined in include/General/Constants.h
    const double *weights, *xi, *eta, *zeta;//quadrature weights and points in reference cell
    // ugly, we need to change GetOrig!!
    double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D]; // quadrature points
    double AbsDetjk[MaxN_QuadPoints_3D]; // determinant of transformation
    int n_points = 0;
    TFEDatabase3D::GetOrig(1, &feID, coll, cell, &SecondDer,
                           n_points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);

    // finite element on the current cell
    TFE3D *fe = TFEDatabase3D::GetFE3D(feID);
    const int n_loc_dof = fe->GetN_DOF(); // number of local dofs
    int * DOF = FESpace3D->GetGlobalDOF(i);

    // id of the local basis functions
    BaseFunct3D base_fc_id = fe->GetBaseFunct3D()->GetID();
    // transformed values of basis functions
    double **orig_values = TFEDatabase3D::GetOrigElementValues(base_fc_id, D000);
    // local integration (loop over all quadrature points)
    for(int j = 0; j < n_points; j++)
    {
      // local transformed values on this quadrature point
      double * orig = orig_values[j];
      double value = 0; // value of this TFEFunction3D at this quadrature point
      for(int l = 0; l < n_loc_dof; l++)
      {
        // entry in the vector of this TFEFunction3D times basis function
        value += Values[DOF[l]] * orig[l];
      } // endfor l

      const double w = weights[j] * AbsDetjk[j];
      integral += w * value;
      measure += w;
    } // endfor j
  }

#ifdef _MPI //communicate the results and add up over all processes
  double sendbuf[2] = {0.0, 0.0};
  double recvbuf[2] = {0.0, 0.0};
  sendbuf[0] = integral; //partial values
  sendbuf[1] = measure;
  MPI_Allreduce(sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  integral = recvbuf[0]; //fill in the summed up values
  measure = recvbuf[1];
#endif
}


/** project function into the space L20 (having zero mean value, or in general a mean value) */
void TFEFunction3D::project_into_L20(double a)
{
  // compute current integral and measure of the domain:
  double integral, measure;
  this->compute_integral_and_measure(integral, measure);
  double new_mean = (integral - a)/measure;

  // vector of the same length as this TFEFunction3D. It represents a function
  // which has the constant value 'mean' for all nodal functionals. The last
  // step in this projection will be to substract this vector from the vector of
  // this TFEFunction3D
  // for standard P_k or Q_k finite elements this is a constant function
  double *interpol = new double[Length];

  TCollection *coll = FESpace3D->GetCollection();
  const int n_cells = coll->GetN_Cells();
  for(int i = 0; i < n_cells; i++)
  {
    TBaseCell *cell = coll->GetCell(i); // current cell
    FE3D feID = FESpace3D->GetFE3D(i, cell); // id of finite element
    // finite element on the current cell
    TFE3D *fe = TFEDatabase3D::GetFE3D(feID);
    const int n_loc_dof = fe->GetN_DOF(); // number of local dofs
    int * DOF = FESpace3D->GetGlobalDOF(i);
    TNodalFunctional3D *nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(feID);
    int n_points; // number of evaluation points to compute nodal functionals
    const double *xi, *eta, *zeta; // coordinates of evaluation points in reference cell
    nf->GetPointsForAll(n_points, xi, eta, zeta);
    double *point_values = new double[n_points];
    for(int j = 0; j < n_points; j++)
      point_values[j] = new_mean;
    // evaluate nodal functionals
    double *functional_values = new double[n_loc_dof];
    nf->GetAllFunctionals(coll, cell, point_values, functional_values);
    for(int j = 0; j < n_loc_dof; j++)
      interpol[DOF[j]] = functional_values[j];

    delete [] point_values;
    delete [] functional_values;
  }

  // the vector 'interpol' is now complete
  // substract it from the vector of this TFEFunction3D
  for(int i = 0; i < Length; i++)
    Values[i] -= interpol[i];

  delete [] interpol;
}

/**Set Dirichlet values according to boundary conditions*/
void TFEFunction3D::SetDirichletBC(BoundCondFunct3D *BoundaryCondition,
                                   BoundValueFunct3D *BoudaryValue)
{
  int i,j, m;
  TBaseCell *cell;
  TCollection *Coll;
  FE3D FEId;
  TFE3D *Element;
  TNodalFunctional3D *nf;  
  int *BeginIndex, *GlobalNumbers;
  int N_Cells, N_Points;
  const double *xi, *eta, *zeta;
  int *DOF;
  RefTrans3D F_K = TetraAffin;
  TRefTrans3D *rt;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AbsDetjk[MaxN_PointsForNodal3D];
  double PointValues[MaxN_PointsForNodal3D];
  double FunctionalValues[MaxN_PointsForNodal3D];  
  int ApproxOrder;
  QuadFormula3D QuadFormula = BaryCenterTetra; //avoid uninit warning
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;  
  BF3DRefElements RefElement;
  RefTrans3D RefTrans, *RefTransArray;
  TFEDesc3D *FEDesc_Obj;
  int N_Faces=0, N_EdgeDOF, N_Joints;
  bool InnerBoundary, OuterBoundary;
  BoundCond Cond0;
  const int *TmpFV, *TmpLen;
  int MaxLen;
  double t0, xf, yf, zf;
  const double *t, *s;
  double LinComb[4] = {0,0,0,0};
  int *EdgeDOF;
  
  
  Coll = FESpace3D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  BeginIndex = FESpace3D->GetBeginIndex();
  GlobalNumbers = FESpace3D->GetGlobalNumbers();  
  
  
  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace3D->GetFE3D(i, cell);
    Element = TFEDatabase3D::GetFE3D(FEId);
    nf = Element->GetNodalFunctional3D();
    nf->GetPointsForAll(N_Points, xi, eta, zeta);
    
    DOF = GlobalNumbers+BeginIndex[i];
    
    ApproxOrder = TFEDatabase3D::GetAccuracyFromFE3D(FEId);

    RefElement = Element->GetBaseFunct3D()->GetRefElement();
    switch(RefElement)
    {
      case BFUnitHexahedron:
        QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(TDatabase::ParamDB->INTERNAL_QUAD_HEXA);
        N_Faces = 6;
      break;

      case BFUnitTetrahedron:
        QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(TDatabase::ParamDB->INTERNAL_QUAD_TETRA);
        N_Faces = 4;
      break;
    }

    RefTrans = RefTransArray[FEId];
    
    IsIsoparametric = false;
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
            IsIsoparametric = true;
        }
        if(jointtype == IsoBoundFace)
          IsIsoparametric = true;
      }
    }
  
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
    } 
    
  
    switch(RefTrans)
    {
      case HexaAffin:
        rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
        ((THexaAffin *)rt)->SetCell(cell);
        F_K = HexaAffin;
        // cout << "HexaAffin: " << endl;
        break;
      case HexaTrilinear:
        rt = TFEDatabase3D::GetRefTrans3D(HexaTrilinear);
        ((THexaTrilinear *)rt)->SetCell(cell);
        F_K = HexaTrilinear;
        // cout << "HexaTrilinear: " << endl;
        break;
      case HexaIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(HexaIsoparametric);
        ((THexaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((THexaIsoparametric *)rt)->SetCell(cell);
        F_K = HexaIsoparametric;
        // cout << "HexaIsoparametric: " << endl;
        break;
      case TetraAffin:
        rt = TFEDatabase3D::GetRefTrans3D(TetraAffin);
        ((TTetraAffin *)rt)->SetCell(cell);
        F_K = TetraAffin;
        // cout << "TetraAffin: " << endl;
        break;
      case TetraIsoparametric:
        rt = TFEDatabase3D::GetRefTrans3D(TetraIsoparametric);
        ((TTetraIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        F_K = TetraIsoparametric;
        // cout << "TetraIsoparametric: " << endl;
        break;
    }
    
    FEDesc_Obj = Element->GetFEDesc3D();
    N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
    N_Joints = cell->GetN_Faces();
    
    for(m=0;m<N_Joints;m++)
    {
      joint = cell->GetJoint(m);
      InnerBoundary = false;
      OuterBoundary = false;
      
      if(joint->GetType() == BoundaryFace ||
           joint->GetType() == IsoBoundFace)
          OuterBoundary = true;
          
       if(InnerBoundary || OuterBoundary)
       {
         cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
         
         t0 = 1./TmpLen[m];
         
         xf =0; yf =0; zf =0;
         for(j=0;j<TmpLen[m];j++)
         {
           cell->GetVertex(TmpFV[m*MaxLen+j])->GetCoords(X[j], Y[j], Z[j]);
           
           xf += t0*X[j];
           yf += t0*Y[j];
           zf += t0*Z[j];                    
         }
         BoundaryCondition(xf, yf, zf, Cond0);
         switch(Cond0)
         {
           case DIRICHLET:
             if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
             {
               /// cout << "USE_ISOPARAMETRIC: " << endl;
               nf->GetPointsForFace(N_Points, t, s);
               for(j=0;j<N_Points;j++)
               {
                 switch(TmpLen[m])
                 {
                   case 4:
                     LinComb[0] = (1-t[j])*(1-s[j]);
                     LinComb[1] = t[j]*(1-s[j]);
                     LinComb[2] = t[j]*s[j];
                     LinComb[3] = (1-t[j])*s[j];
                      
                     xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                          +LinComb[2]*X[2] + LinComb[3]*X[3];
                     yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                          +LinComb[2]*Y[2] + LinComb[3]*Y[3];
                     zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                          +LinComb[2]*Z[2] + LinComb[3]*Z[3];
                     break;
                   case 3:
                     LinComb[0] = 1-t[j]-s[j];
                     LinComb[0] = t[j];
                     LinComb[0] = s[j];
                      
                     xf = LinComb[0]*X[0] + LinComb[1]*X[1]
                          +LinComb[2]*X[2];
                     yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]
                          +LinComb[2]*Y[2];
                     zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]
                          +LinComb[2]*Z[2];
                     break;
                 }
                 BoudaryValue(xf, yf, zf, PointValues[j]);
               }
             }
             else
             {
               /// cout << "Non_ISOPARAMETRIC: " << endl;
               nf->GetPointsForFace(m, N_Points, xi, eta, zeta);
                TFEDatabase3D::GetOrigFromRef(F_K, N_Points, xi, eta, zeta,
                                                  X, Y, Z, AbsDetjk);

              
                for(j=0;j<N_Points;j++)             
                  BoudaryValue(X[j], Y[j], Z[j], PointValues[j]);
             }
             nf->GetFaceFunctionals(Coll, cell, m, PointValues, FunctionalValues);
              
             EdgeDOF = FEDesc_Obj->GetJointDOF(m);
             N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
              
            for(j=0;j<N_EdgeDOF;j++)
            {
               Values[DOF[EdgeDOF[j]]] = FunctionalValues[j];
               // cout << i << setw(20) << Values[DOF[EdgeDOF[j]]] << endl;
             }
             break;
	  default:
	    
	  break;	     
	     
         }         
       }///endif
    }///for m<N_Joints
  }///endfor i<N_Cells
}

void TFEFunction3D::MinMax(double & min, double & max) const
{
#ifdef _MPI
  //in MPI case, compare only master values
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  const int* masters = this->FESpace3D->get_communicator().GetMaster();
#endif

  double val;
  max = -1e100, min = 1e100;
  
  for(int i = 0; i<Length; i++)
  {
#ifdef _MPI
    //skip slave values
    if(masters[i] != my_rank)
      continue;
#endif
    val = Values[i];
    if(val>max) max = val;
    if(val<min) min = val;
  }
}

void TFEFunction3D::PrintMinMax(const std::string& name) const
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif

  double min, max;
  this->MinMax(min, max);

#ifdef _MPI
  // reduce min and max in root
  double sendbuf_min = min;
  double sendbuf_max = max;
  MPI_Reduce(&sendbuf_min, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
  MPI_Reduce(&sendbuf_max, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
#endif

  if(my_rank ==0) //only root has results - only root prints
  {
    if( min <= max )
    {
      if(name.empty())
        Output::stat("MinMax", this->Name, " min ", min, ", max ", max);
      else
        Output::stat("MinMax", name, " min ", min, ", max ", max);
    }
    else
    {
        Output::warn("TFEFunction3D::MinMax was not successful!");
    }
  }
}





