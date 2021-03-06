// =======================================================================
// @(#)FEVectFunct2D.C        1.3 07/22/99
// 
// Class:       TFEVectFunct2D
// Purpose:     a function from a finite element space in 2D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
//              WriteSol/ReadSol    08.09.09 (Sashikumaar Ganesan)
// =======================================================================
#ifdef __2D__

#ifdef _MPI
# include "mpi.h"
#endif

#include <FEVectFunct2D.h>
#include <FEDatabase2D.h>
#include <NodalFunctional2D.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <Database.h>
#include <AuxParam2D.h>
#include <GridCell.h>

#include <Joint.h>
#include <BoundEdge.h>
#include <InterfaceJoint.h>

#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <MooNMD_Io.h>
// #include <malloc.h>
#include <dirent.h> 
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <BoundaryAssembling2D.h>

/// Default constructor. Constructs an empty object.
TFEVectFunct2D::TFEVectFunct2D()
{
  N_Components = 0;
}


/** constructor with vector initialization */
TFEVectFunct2D::TFEVectFunct2D(std::shared_ptr<const TFESpace2D> fespace2D,
                               const std::string& name,
                               const std::string& description, double *values,
                               int length, int n_components)
  : TFEFunction2D(fespace2D, name, description, values, length)
{
  N_Components = n_components;
}

TFEVectFunct2D& TFEVectFunct2D::operator=( const TFEVectFunct2D & other)
{
  //call base class copy assignment
  TFEFunction2D::operator=(other);

  this->N_Components  = other.N_Components;

  return *this;
}

//====================================================================
/** calculate errors to given vector function */
void TFEVectFunct2D::GetDeformationTensorErrors( 
  DoubleFunct2D *Exact, DoubleFunct2D *Exact1,
  int N_Derivatives,
  MultiIndex2D *NeededDerivatives,
  int N_Errors, TFEFunction2D::ErrorMethod *ErrorMeth, 
  const CoeffFct2D& Coeff, 
  TAuxParam2D *Aux,
  int n_fespaces, TFESpace2D **fespaces,
  double *errors)
{
  int i,j,k,l,N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;  //N_U;
  int Used[N_FEs2D], *N_BaseFunct;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  const double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Param[MaxN_QuadPoints_2D], *aux, *aux1, *aux2, *aux3;
  double *Derivatives[2*MaxN_QuadPoints_2D];
  double *ExactVal[2*MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF;
  double **OrigFEValues, *Orig, value, value1;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4], *Values0,*Values1;
  double hK;
  bool *SecondDer;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = false;

  N_Parameters = Aux->GetN_Parameters();
  aux1 = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux1 + j*N_Parameters;

  aux2 = new double [2*MaxN_QuadPoints_2D*N_Derivatives];
  for(j=0;j<2*MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux2 + j*N_Derivatives;
  
  aux3 = new double [2*MaxN_QuadPoints_2D * 4];
  for(j=0;j<2*MaxN_QuadPoints_2D;j++)
    ExactVal[j] = aux3 + j*4;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_2D*20]; 
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*20;

  fespace = fespaces[0];
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
//  N_U = Length;
  Values0 = Values;
  Values1 = Values+Length;

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  auto Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);

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
    CurrentElement = fespace->GetFE2D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values0[DOF[l]];
      FEFunctValues1[l] = Values1[DOF[l]];
    }

    // for all needed derivatives
    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                      NeededDerivatives[k]);
      // for all quadrature points
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        value1 = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
        } // endfor l
        Derivatives[j][k] = value;
        Derivatives[j+N_Points][k] = value1;
      } // endfor j
    } // endfor k

    // exact value for first component
    for(j=0;j<N_Points;j++)
      Exact(X[j], Y[j], ExactVal[j]);

    // exact value for second component
    for(j=0;j<N_Points;j++)
      Exact1(X[j], Y[j], ExactVal[j+N_Points]);

    if(Coeff)
      Coeff(N_Points, X, Y, Param, AuxArray);      

    ErrorMeth(N_Points, {{X, Y}}, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
      errors[j] += LocError[j];

  } // endfor i

  for(j=0;j<N_Errors;j++)
  {
    if (errors[j] > 0)
    errors[j] = sqrt(errors[j]);
  }

  delete aux;
  delete aux1;
  delete aux2;
  delete aux3;
  delete SecondDer;

} // TFEVectFunct2D::GetDeformationTensorErrors

//================================================================
std::pair<double, double> TFEVectFunct2D::get_L2_norm_divergence_curl() const
{
  std::vector<double> values(2);
  auto f = [](std::vector<double>& v, std::array<double, 8> e)
           {
             v[0] += e[4] + e[7]; // divergence
             v[1] += e[5] - e[6]; // curl
           };
  this->get_functional_value(values, f);
  return {values[0], values[1]};

} // TFEVectFunct2D::GetL2NormDivergence

//================================================================
void TFEVectFunct2D::get_functional_value(std::vector<double>& values,
  const std::function<void(std::vector<double>&, std::array<double, 8>)>&
    functional)
const
{
  BaseFunct2D *BaseFuncts;
  int *N_BaseFunct;
  bool SecondDer[1] = { false };
  int N_Points;
  const double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetJK[MaxN_QuadPoints_2D];
  double FEValues0[MaxN_BaseFunctions2D];
  double FEValues1[MaxN_BaseFunctions2D];
  
  auto n_values = values.size();
  std::fill(values.begin(), values.end(), 0.0);

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  double *Values0 = Values;
  double *Values1 = Values+Length;

  auto Coll = FESpace2D->GetCollection();
  auto N_Cells = Coll->GetN_Cells();

  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
    auto FEid = FESpace2D->GetFE2D(i, cell);
    TFEDatabase2D::GetOrig(1, &FEid, Coll, cell, SecondDer, N_Points, xi, eta,
                           weights, X, Y, AbsDetJK);

    // calculate all needed derivatives of this FE function
    auto BaseFunct = BaseFuncts[FEid];
    auto N_Bf = N_BaseFunct[FEid];

    auto DOF = FESpace2D->GetGlobalDOF(i);
    
    for(int j = 0; j < N_Bf; j++)
    {
      int k = DOF[j];
      FEValues0[j] = Values0[k];
      FEValues1[j] = Values1[k];
    }

    auto OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    auto OrigFEValuesX = TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
    auto OrigFEValuesY = TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);

    // for all quadrature points
    for(int j = 0; j < N_Points; j++)
    {
      double * Orig = OrigFEValues[j];
      double * OrigX = OrigFEValuesX[j];
      double * OrigY = OrigFEValuesY[j];
      std::vector<double> local_values(n_values, 0.0);
      std::array<double, 8> evaluations{{X[j], Y[j], 0., 0., 0., 0., 0., 0.}};
      for(int l = 0; l < N_Bf; l++)
      {
        evaluations[2] += Orig[l] * FEValues0[l];
        evaluations[3] += Orig[l] * FEValues1[l];
        evaluations[4] += OrigX[l] * FEValues0[l];
        evaluations[5] += OrigX[l] * FEValues1[l];
        evaluations[6] += OrigY[l] * FEValues0[l];
        evaluations[7] += OrigY[l] * FEValues1[l];
      } // endfor l
      functional(local_values, evaluations);
      double local_weight = AbsDetJK[j] * weights[j];
      for(auto k = 0ul; k < n_values; ++k)
      {
        values[k] += local_weight * (local_values[k] * local_values[k]);
      }
    } // endfor j
  } // endfor i
  for(auto k = 0ul; k < n_values; ++k)
  {
    values[k] = std::sqrt(values[k]);
  }
} // TFEVectFunct2D::get_functional_value


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** calculate L2-norm of (u-u_h).n-error - written by Laura Blank 01.03.18*/
double TFEVectFunct2D::GetL2NormNormalComponentError(BoundValueFunct2D *Exact_u1, BoundValueFunct2D *Exact_u2, bool rescale_by_h_E)
{
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  double *Values0 = Values;
  double *Values1 = Values+Length;
  int *GlobalNumbers = FESpace2D->GetGlobalNumbers();
  int *BeginIndex = FESpace2D->GetBeginIndex();

  auto Coll = FESpace2D->GetCollection();

  // Initialize Pointer (necessary: Resource allocation is initialization Paradigm)
  double *ExactVal_u1[MaxN_QuadPoints_2D], *ExactVal_u2[MaxN_QuadPoints_2D];
  double *aux1 = new double [MaxN_QuadPoints_2D * 4]{0.};
  double *aux2 = new double [MaxN_QuadPoints_2D * 4]{0.};
  for(int ii = 0; ii < MaxN_QuadPoints_2D; ii++)
  {
    ExactVal_u1[ii] = aux1 + ii*4;
    ExactVal_u2[ii] = aux2 + ii*4;
  }

  double FEFunctValues0[MaxN_BaseFunctions2D], FEFunctValues1[MaxN_BaseFunctions2D]; 
  double final_boundary_error_l2[1];
  final_boundary_error_l2[0] = 0;

  for(int k = 0; k < TDatabase::ParamDB->n_nitsche_boundary; k++)
  {
    int boundary_component_id = TDatabase::ParamDB->nitsche_boundary_id[k];
    // Create a list of those boundary edges that are on the boundary component with given ID
    std::vector<TBoundEdge*> boundaryEdgeList;
    Coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

    double boundary_error_l2_on_Component = 0;

    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
      TBoundEdge *boundedge = boundaryEdgeList[m];
      TBaseCell *cell = boundedge->GetNeighbour(0);
      FE2D CurrentElement_FEID = FESpace2D->GetFE2D(0, cell);
 
      int joint_id = boundedge->get_index_in_neighbour(cell); 
      // get a quadrature formula good enough for the argument of the integral
      int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement_FEID);
      QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree((TDatabase::ParamDB->INPUT_QUAD_RULE < 2*fe_degree)? 2*fe_degree : TDatabase::ParamDB->INPUT_QUAD_RULE);
      std::vector<double> quadWeights, quadPoints;
      BoundaryAssembling2D::get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);
      // compute values of all basis functions and their first partial derivatives at all quadrature points
      std::vector< std::vector<double> > uorig, u_dx_orig ,u_dy_orig;
      int BaseVectDim = 1;
      BoundaryAssembling2D::get_original_values(CurrentElement_FEID, joint_id, cell, quadPoints, BaseVectDim, uorig, u_dx_orig, u_dy_orig, LineQuadFormula);

      // calculate all needed derivatives of this FE function                 
      int N_Bf = N_BaseFunct[CurrentElement_FEID];

      int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()];
      for(int l = 0; l < N_Bf; l++)
      {
        FEFunctValues0[l] = Values0[DOF[l]];
        FEFunctValues1[l] = Values1[DOF[l]];
      }

      double summing_boundary_error_l2_on_edge = 0;
      double edge_length = boundedge->get_length();
      double reference_edge_length = 2; // [-1,1] is the reference edge here↲   
      // normal vector to this boundary (normalized)
      double n1, n2;
      boundedge->get_normal(n1, n2);

      for(size_t j = 0; j < quadPoints.size(); j++)
      {
        double value = 0;
        double value_u1 = 0;
        double value_u2 = 0;
        double exact_val = 0;
        for(int l = 0; l < N_Bf; l++)
        { // compute u_h|T = \sum \alpha_i \Phi_i 
          value_u1 += FEFunctValues0[l] * uorig[j][l];
          value_u2 += FEFunctValues1[l] * uorig[j][l];
        } 
        value += value_u1 * n1 + value_u2 * n2;

        auto comp = boundedge->GetBoundComp();
        int comp_ID = comp->GetID();
        double t0, t1;
        boundedge->GetParameters(t0, t1);
        double t = t0 + 0.5 * (t1-t0) * (quadPoints[j]+1); 
        Exact_u1(comp_ID, t, ExactVal_u1[j][0]);
        Exact_u2(comp_ID, t, ExactVal_u2[j][0]);
        exact_val += ExactVal_u1[j][0] * n1 + ExactVal_u2[j][0] * n2;

        if(rescale_by_h_E)
        { 
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * 1/edge_length  * (exact_val - value) * (exact_val - value);
        }
        else
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * (exact_val - value) * (exact_val - value);
        }
      }
      boundary_error_l2_on_Component += summing_boundary_error_l2_on_edge;
    } 
    final_boundary_error_l2[0] +=  boundary_error_l2_on_Component;
  }
  final_boundary_error_l2[0] = sqrt(final_boundary_error_l2[0]);
  return final_boundary_error_l2[0];
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** calculate L2-norm of (u-u_h).n-error without the global database */
double TFEVectFunct2D::GetL2NormNormalComponentError(BoundValueFunct2D *Exact_u1, BoundValueFunct2D *Exact_u2, int boundary_component_id, bool rescale_by_h_E)
{
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  double *Values0 = Values;
  double *Values1 = Values+Length;
  int *GlobalNumbers = FESpace2D->GetGlobalNumbers();
  int *BeginIndex = FESpace2D->GetBeginIndex();

  auto Coll = FESpace2D->GetCollection();

  // Initialize Pointer (necessary: Resource allocation is initialization Paradigm)
  double *ExactVal_u1[MaxN_QuadPoints_2D], *ExactVal_u2[MaxN_QuadPoints_2D];
  double *aux1 = new double [MaxN_QuadPoints_2D * 4]{0.};
  double *aux2 = new double [MaxN_QuadPoints_2D * 4]{0.};
  for(int ii = 0; ii < MaxN_QuadPoints_2D; ii++)
  {
    ExactVal_u1[ii] = aux1 + ii*4;
    ExactVal_u2[ii] = aux2 + ii*4;
  }

  double FEFunctValues0[MaxN_BaseFunctions2D], FEFunctValues1[MaxN_BaseFunctions2D];
  double final_boundary_error_l2[1];
  final_boundary_error_l2[0] = 0;

    // Create a list of those boundary edges that are on the boundary component with given ID
    std::vector<TBoundEdge*> boundaryEdgeList;
    Coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

    double boundary_error_l2_on_Component = 0;

    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
      TBoundEdge *boundedge = boundaryEdgeList[m];
      TBaseCell *cell = boundedge->GetNeighbour(0);
      FE2D CurrentElement_FEID = FESpace2D->GetFE2D(0, cell);

      int joint_id = boundedge->get_index_in_neighbour(cell);
      // get a quadrature formula good enough for the argument of the integral
      int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement_FEID);
      QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree((TDatabase::ParamDB->INPUT_QUAD_RULE < 2*fe_degree)? 2*fe_degree : TDatabase::ParamDB->INPUT_QUAD_RULE);
      std::vector<double> quadWeights, quadPoints;
      BoundaryAssembling2D::get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);
      // compute values of all basis functions and their first partial derivatives at all quadrature points
      std::vector< std::vector<double> > uorig, u_dx_orig ,u_dy_orig;
      int BaseVectDim = 1;
      BoundaryAssembling2D::get_original_values(CurrentElement_FEID, joint_id, cell, quadPoints, BaseVectDim, uorig, u_dx_orig, u_dy_orig, LineQuadFormula);

      // calculate all needed derivatives of this FE function
      int N_Bf = N_BaseFunct[CurrentElement_FEID];

      int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()];
      for(int l = 0; l < N_Bf; l++)
      {
        FEFunctValues0[l] = Values0[DOF[l]];
        FEFunctValues1[l] = Values1[DOF[l]];
      }

      double summing_boundary_error_l2_on_edge = 0;
      double edge_length = boundedge->get_length();
      double reference_edge_length = 2; // [-1,1] is the reference edge here↲
      // normal vector to this boundary (normalized)
      double n1, n2;
      boundedge->get_normal(n1, n2);

      for(size_t j = 0; j < quadPoints.size(); j++)
      {
        double value = 0;
        double value_u1 = 0;
        double value_u2 = 0;
        double exact_val = 0;
        for(int l = 0; l < N_Bf; l++)
        { // compute u_h|T = \sum \alpha_i \Phi_i
          value_u1 += FEFunctValues0[l] * uorig[j][l];
          value_u2 += FEFunctValues1[l] * uorig[j][l];
        }
        value += value_u1 * n1 + value_u2 * n2;

        auto comp = boundedge->GetBoundComp();
        int comp_ID = comp->GetID();
        double t0, t1;
        boundedge->GetParameters(t0, t1);
        double t = t0 + 0.5 * (t1-t0) * (quadPoints[j]+1);
        Exact_u1(comp_ID, t, ExactVal_u1[j][0]);
        Exact_u2(comp_ID, t, ExactVal_u2[j][0]);
        exact_val += ExactVal_u1[j][0] * n1 + ExactVal_u2[j][0] * n2;

        if(rescale_by_h_E)
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * 1/edge_length  * (exact_val - value) * (exact_val - value);
        }
        else
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * (exact_val - value) * (exact_val - value);
        }
      }
      boundary_error_l2_on_Component += summing_boundary_error_l2_on_edge;
    }
    final_boundary_error_l2[0] =  boundary_error_l2_on_Component;

  return final_boundary_error_l2[0];
}

//==========================================================================
/** calculate L2-norm of divergence error - written by Laura Blank 03.01.18*/
double TFEVectFunct2D::GetL2NormDivergenceError(DoubleFunct2D *Exact_u1,DoubleFunct2D *Exact_u2)
{
  BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  double *Values0 = Values;
  double *Values1 = Values+Length;

  int *GlobalNumbers = FESpace2D->GetGlobalNumbers();
  int *BeginIndex = FESpace2D->GetBeginIndex();

  auto Coll = FESpace2D->GetCollection();
  int N_Cells = Coll->GetN_Cells();

  // Initialize Pointer (necessary: Resource allocation is initialization Paradigm)
  double *ExactVal_u1[MaxN_QuadPoints_2D];
  double *ExactVal_u2[MaxN_QuadPoints_2D];
  double *aux1 = new double [MaxN_QuadPoints_2D * 4]{0.};
  double *aux2 = new double [MaxN_QuadPoints_2D * 4]{0.};
  for(int ii = 0; ii < MaxN_QuadPoints_2D; ii++)
  {
    ExactVal_u1[ii] = aux1 + ii*4;
    ExactVal_u2[ii] = aux2 + ii*4;
  }

  double loc_div_values_exact_solution = 0.;
  double loc_div_values_FEsolution = 0.;
  int N_UsedElements = 1;
  FE2D UsedElements[1];
  double FEFunctValues0[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double **OrigFEValuesX, *OrigX, value;
  double **OrigFEValuesY, *OrigY;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  bool SecondDer[1] = { false };
  double AbsDetJK[MaxN_QuadPoints_2D];

  double div_error = 0.;
  // loop over all cells
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
    FE2D FEid = FESpace2D->GetFE2D(i, cell);
    UsedElements[0] = FEid;

    // compute transformation to reference cell
    const double *xi, *eta, *weights;
    int N_Points;
    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements, 
        Coll, cell, SecondDer,
        N_Points, xi, eta, weights, X, Y, AbsDetJK);

    // calculate all needed derivatives of this FE function
    BaseFunct2D BaseFunct = BaseFuncts[FEid];
    int N_Bf = N_BaseFunct[FEid];

    int *DOF = GlobalNumbers + BeginIndex[i];
    for(int jj = 0; jj < N_Bf; jj++)
    {
      int k = DOF[jj];
      FEFunctValues0[jj] = Values0[k];
      FEFunctValues1[jj] = Values1[k];
    }

    OrigFEValuesX = TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
    OrigFEValuesY = TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);

    // loop over all quadrature points
    for(int j = 0; j < N_Points; j++)
    {
      OrigX = OrigFEValuesX[j];
      OrigY = OrigFEValuesY[j];
      value = 0;
      for(int l = 0; l < N_Bf; l++)
      {
        value += FEFunctValues0[l] * OrigX[l] + FEFunctValues1[l] * OrigY[l];
      }
      loc_div_values_FEsolution = value;

      Exact_u1(X[j], Y[j], ExactVal_u1[j]);
      Exact_u2(X[j], Y[j], ExactVal_u2[j]);
      loc_div_values_exact_solution = ExactVal_u1[j][1] + ExactVal_u2[j][2];

      div_error += AbsDetJK[j] * weights[j] * fabs(loc_div_values_FEsolution - loc_div_values_exact_solution) * fabs(loc_div_values_FEsolution - loc_div_values_exact_solution);
    } // endfor j
  } // endfor i
  delete [] aux1;
  delete [] aux2;
  div_error = sqrt(div_error);
  return div_error;
} // TFEVectFunct2D::GetL2NormDivergenceError


//==========================================================================
/** write the solution into a data file - written by Sashi **/
void TFEVectFunct2D::WriteSol(double t, const std::string& directory,
                              const std::string& basename)
{
  int i, N_Joints, N_Cells;

  char Dquot;

  #ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  #endif

  Dquot = 34; //  see ASCII Chart
  auto Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  i=0;
  auto cell =  Coll->GetCell(i);
  N_Joints = cell->GetN_Joints();
  const char* BaseName = basename.c_str();
  const char* output_directory = directory.c_str();

  std::ostringstream os;
  os << " ";

  #ifdef _MPI
  OutPut("Writing solution into "<< output_directory << "/" << BaseName << rank
         << ".Sol MooNMD file"<< endl);
  os.seekp(std::ios::beg);
  os << output_directory << "/" << BaseName<<rank<<".Sol" << ends;
  #else
  OutPut("Writing solution into "<< output_directory << "/" << BaseName << t
         << ".Sol MooNMD file"<< endl);
  os.seekp(std::ios::beg);
  os << output_directory << "/" << BaseName << t<<".Sol" << ends;
  #endif

  std::ofstream dat(os.str().c_str());

  if (!dat)
   {
    cerr << "cannot open file for output" << endl;
    exit(0);
   }

    dat << "# Solution of the vector "<<Dquot<<Name<<Dquot<<", written by MooNMD"<< endl;
    dat << "# N_Cells, Cell_Type, N_Dim, N_Dof" <<  endl;
    dat <<N_Cells << " " << N_Joints << " " << N_Components << " " << Length<< endl;
    dat <<  endl;

    dat << "# Dof "<< " Nodal Values"<< endl;

    for(i=0;i<Length;i++)
     dat << i << " " << Values[i] << " " << Values[Length + i] << endl;
  dat.close();
}


/** Read the solution from a given data file - written by Sashi **/
void TFEVectFunct2D::ReadSol(const std::string& BaseName)
{
 int i, j, N_Joints, N_Cells, N_cells, N_joints, N_components, length;
 char line[100];

  auto Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  i=0;
  auto cell =  Coll->GetCell(i);
  N_Joints = cell->GetN_Joints();

  std::ifstream dat(BaseName);
  if (!dat)
   {
    cerr << "cannot open '" <<  BaseName << "' for input" << endl;
    exit(-1);
   }

  dat.getline (line, 99);
  dat.getline (line, 99);
  dat >> N_cells >> N_joints >> N_components >> length;
  dat.getline (line, 99);

  if(N_cells!=N_Cells || N_joints!=N_Joints || N_components!=N_Components || length!=Length )
   {
    OutPut("Given data file does not match with this FE Vector function !"<<endl);
    OutPut(N_cells <<", "<< N_joints<< ", " << N_components<< " , "<< length<< " , " <<endl);
    OutPut(N_Cells <<", "<< N_Joints<< ", " << N_Components<< " , "<< Length<< " , " <<endl);
    exit(0);
   }

  dat.getline (line, 99);
  OutPut("Reading nodal values of the FE Vector function !"<<endl);

  for(i=0;i<Length;i++)
   {
    dat.getline (line, 99);
    dat >> j >> Values[i] >> Values[Length + i];
   }

  dat.close();
}

/** interpolate the old vect value to the new function */
void TFEVectFunct2D::Interpolate(TFEVectFunct2D *OldVectFunct)
{
 int i, j, N_Cells, N_DOFs, N_LocalDOFs, N_Points, N_Edges;
 int *BeginIndex, *GlobalNumbers, *DOF;
 int PolynomialDegree, ApproxOrder;

 const double *xi, *eta;
 double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
 double AbsDetjk[MaxN_PointsForNodal2D];
 double PointValues1[MaxN_PointsForNodal2D];
 double PointValues2[MaxN_PointsForNodal2D];
 double FunctionalValues[MaxN_PointsForNodal2D];
 double val1[3], val2[3];
  
 bool IsIsoparametric;
   
 RefTrans2D RefTrans, *RefTransArray;
 FE2D FEId;
 TFE2D *Element;
 TNodalFunctional2D *nf;
 QuadFormula2D QuadFormula; 
 BF2DRefElements RefElement; 
 TJoint *joint;
 JointType jointtype;
 BoundTypes bdtype;
 RefTrans2D F_K;
 TRefTrans2D *rt;

  auto Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_DOFs = FESpace2D->GetN_DegreesOfFreedom();
  
  memset(Values, 0, SizeOfDouble*N_Components*N_DOFs);
  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

    
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
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

    IsIsoparametric = false;
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
            IsIsoparametric = true;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == IsoInterfaceJoint ||
           jointtype == IsoBoundEdge)
          IsIsoparametric = true;
      }
    } // endif    
 
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
    } // endif IsIsoparametric
 

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
     OldVectFunct->FindVectGradient(X[j], Y[j], val1, val2);
     PointValues1[j] = val1[0];
     PointValues2[j] = val2[0];
    }// for(j=0;j<N_Points;j++)

    
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues1, FunctionalValues);    

    DOF = GlobalNumbers+BeginIndex[i];
    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];
    
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues2, FunctionalValues);     

    for(j=0;j<N_LocalDOFs;j++)
      Values[N_DOFs + DOF[j]] = FunctionalValues[j];

      
  } // for(i=0;i<N_Cells;i++)  
} // Interpolate
   
   
/** determine the value of a vect function and its first derivatives at
 the given point */
void TFEVectFunct2D::FindVectGradient(double x, double y, double *val1, double *val2)
{
  int i, j, N_Cells, N_Found, N_BaseFunct, N_DOFs;
  int *BeginIndex, *GlobalNumbers, *Numbers, N_Edges;
  
  double xi, eta, U1, U2;
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref;
  double u1, u1x, u1y, u2, u2x, u2y;
  
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  JointType jointtype;
  TJoint *joint;
  BoundTypes bdtype;
  BF2DRefElements RefElement; 
 
  bool IsIsoparametric;
  
  N_Found = 0;
  
  val1[0] = 0;
  val1[1] = 0;
  val1[2] = 0;

  val2[0] = 0;
  val2[1] = 0;
  val2[2] = 0;  

  
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_DOFs = FESpace2D->GetN_DegreesOfFreedom();
  
  auto Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();  

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    
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
      // mainly used for remeshing, so Affine trans os enough
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      
      // values are calculated in oldfunction, so isoparam        
      IsIsoparametric = false;           
      if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
      {
       RefElement = FE_Obj->GetBaseFunct2D()->GetRefElement();     

       switch(RefElement)
       {
         case BFUnitSquare:    
          N_Edges = 4;
         break;

         case BFUnitTriangle:
          N_Edges = 3;
         break;

         default:
         break;
        }      

       for(j=0;j<N_Edges;j++)
       {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
         bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
         if(bdtype != Line)
          IsIsoparametric = true;
        }
       if(jointtype == InterfaceJoint)
       {
        bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
        if(bdtype != Line)
          IsIsoparametric = true;
       }
       
       if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
        IsIsoparametric = true;
       }
      }// endif 

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
     } // endif IsIsoparametric           
           
      // set cell for reference transformation
      TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
      
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

      u1 = 0;
      u1x = 0;
      u1y = 0;
      
      u2 = 0;
      u2x = 0;
      u2y = 0;

      Numbers = GlobalNumbers + BeginIndex[i];
      for(j=0;j<N_BaseFunct;j++)
      {
        U1 = Values[Numbers[j]];
        u1 += uorig[j]*U1;
        u1x += uxorig[j]*U1;
        u1y += uyorig[j]*U1;

        U2 = Values[N_DOFs + Numbers[j]];        
        u2 += uorig[j]*U2;
        u2x += uxorig[j]*U2;
        u2y += uyorig[j]*U2;        
      }
  
      val1[0] += u1;
      val1[1] += u1x;
      val1[2] += u1y;
      
      val2[0] += u2;
      val2[1] += u2x;
      val2[2] += u2y;     

      delete [] uorig;
      delete [] uxorig;
      delete [] uyorig;
      delete [] uref;
      delete [] uxiref;
      delete [] uetaref;      
      
    } // if(cell->PointInCell(x,y))
  } // for(i=0;i<N_Cells;i+



  if(N_Found>0)
   {
    val1[0] /= (double)N_Found;
    val1[1] /= (double)N_Found;
    val1[2] /= (double)N_Found;
    
    val2[0] /= (double)N_Found;
    val2[1] /= (double)N_Found;
    val2[2] /= (double)N_Found;
   }
  else 
   {
    cout<<"("<<x<<" , " <<y<<" ) Point not found !!!!!"<<endl;
    exit(0);
   }
}

void TFEVectFunct2D::FindValueLocal(const TBaseCell* cell, int cell_no, 
				    double x, double y, double* values) const
{
 this->TFEFunction2D::FindValueLocal(cell, cell_no, x, y, values);
 auto u2 = this->GetComponent(1);
 u2->FindValueLocal(cell, cell_no, x, y, values+1);
 delete u2;
}


TFEVectFunct2D& TFEVectFunct2D::operator*=(double alpha)
{
  int N_Active = FESpace2D->GetActiveBound();
  for(int n=0; n<N_Components; n++)
  {
    for (int i=0; i<N_Active; i++)
    {
      Values[i+n*Length] *= alpha;
    }
  }
  return *this;
}

TFEVectFunct2D & TFEVectFunct2D::operator+=(const TFEVectFunct2D & rhs)
{
  if(FESpace2D != rhs.FESpace2D)
  {
    OutPut("ERROR: TFEVectFunct2D::operator+=() The two arguments "
           << "have different fe spaces. Exiting" << endl);
    exit(1);
  }
  if(Length != rhs.Length)
  {
    OutPut("ERROR: TFEVectFunct2D::operator+=() The two arguments "
           << "have different lengths. Exiting" << endl);
    exit(1);
  }
  if(Values == rhs.Values)
  {
    OutPut("ERROR: TFEVectFunct2D::operator+=() The two arguments "
           << "have the same solution vector. This operation would be "
           << "equivalent to a multiplication by 2! Exiting" << endl);
    exit(1);
  }
  if(N_Components != rhs.N_Components)
  {
    OutPut("ERROR: TFEVectFunct2D::operator+=() The two arguments "
           << "have different number of components. You are trying to add a "
           << rhs.N_Components << "-dimensional vector to a " << N_Components
           << "-dimensional vector! Exiting" << endl);
    exit(1);
  }
  int N_Active = FESpace2D->GetActiveBound();
  for(int n=0; n<N_Components; n++)
  {
    for (int i=0; i<N_Active; i++)
    {
      Values[i+n*Length] += rhs.Values[i+n*Length];
    }
  }
  return *this;
}

#endif // #ifdef __2D__
