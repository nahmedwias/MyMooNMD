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
#include <AuxParam2D.h>
#include <GridCell.h>
#include <ConvDiff.h>

#include <InterfaceJoint.h>
#include <IsoBoundEdge.h>
#include <BdLine.h>
#include <LinAlg.h>

#include <fstream>
#include <sstream>
#include <dirent.h> 
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

#include <BaseCell.h>
#include <BoundaryAssembling2D.h>

void OnlyDirichlet(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

TFEFunction2D::TFEFunction2D() :
  Name("dummy_fe_fct_2d"), Description("dummy_fe_fct_2d")
{
  FESpace2D=nullptr;
  Values=nullptr;
  Length=0;
}

/** constructor with vector initialization */
TFEFunction2D::TFEFunction2D(const TFESpace2D *fespace2D, std::string name,
    std::string description, double *values, int length)
: Name(name), Description(description)
{
  Output::print<5>("Constructor of TFEFunction2D");

  FESpace2D=fespace2D;

  Values=values;

  Length=length;
}


TFEFunction2D::~TFEFunction2D()
{
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** Calculate L2-errors to a given function at the boundary */
/** written by Laura Blank 27.02.18 */
void TFEFunction2D::GetL2BoundaryError(BoundValueFunct2D *Exact,
                                       TAuxParam2D *, int,
                                       const TFESpace2D **fespaces,
                                       double *final_boundary_error_l2,
                                       bool rescale_by_h_E)
{
  int *GlobalNumbers = FESpace2D->GetGlobalNumbers();
  int *BeginIndex = FESpace2D->GetBeginIndex();
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // ########################################################################
  // loop over all Nitsche edges
  // ########################################################################
  TCollection *Coll = fespaces[0]->GetCollection();  // all spaces use the same Coll

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
      // get basis dimension and FE space data of the current cell
      // Here it is assumed that only one local FE space is used for the hole mesh!
      FE2D CurrentElement_FEID = FESpace2D->GetFE2D(0, cell);

      BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
      // calculate all needed values and derivatives of the basis functions
      BaseFunct2D BaseFunct = BaseFuncts[CurrentElement_FEID];
      if(TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim() > 1)
      {
        ErrMsg("The method "
            << "TFEFunction2D::GetL2BoundaryErrors"
            << "is not designed for vector valued functions");
        OutPut("No errors were computed\n");
        return;
      }

      int BaseVectDim = 1; 
      int joint_id = boundedge->get_index_in_neighbour(cell);
      // get a quadrature formula good enough for the argument of the integral
      int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement_FEID);
      QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree((TDatabase::ParamDB->INPUT_QUAD_RULE < 2*fe_degree)? 2*fe_degree : TDatabase::ParamDB->INPUT_QUAD_RULE);
      std::vector<double> quadWeights, quadPoints;
      BoundaryAssembling2D::get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);
      // compute values of all basis functions and their first partial derivatives at all quadrature points
      std::vector< std::vector<double> > uorig, u_dx_orig ,u_dy_orig;
      BoundaryAssembling2D::get_original_values(CurrentElement_FEID, joint_id, cell, quadPoints, BaseVectDim, uorig, u_dx_orig, u_dy_orig, LineQuadFormula);

      int N_ = N_BaseFunct[CurrentElement_FEID];

      double FEFunctValues[MaxN_BaseFunctions2D];
      int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()];
      for(int l = 0; l < N_; l++)
      {
        FEFunctValues[l] = Values[DOF[l]];
      }

      double summing_boundary_error_l2_on_edge = 0;
      double edge_length = boundedge->get_length();
      double reference_edge_length = 2; // [-1,1] is the reference edge here

      for(size_t j = 0; j < quadPoints.size(); j++)
      {
        double value = 0;
        for(int l = 0; l < N_; l++)
        { 
          value += FEFunctValues[l] * uorig[j][l]; // compute u_h|T = \sum \alpha_i \Phi_i
        }
        auto comp = boundedge->GetBoundComp();
        int comp_ID = comp->GetID();
        double t0, t1;
        boundedge->GetParameters(t0, t1);
        double t = t0 + 0.5 * (t1-t0) * (quadPoints[j]+1);
        double ExactVal, Diff;
        Exact(comp_ID, t, ExactVal);
        Diff = ExactVal - value;

        if(rescale_by_h_E)
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * 1/edge_length * Diff * Diff; 
        }
        else
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * Diff * Diff; 
        }
      }
      boundary_error_l2_on_Component += summing_boundary_error_l2_on_edge;
    }
    final_boundary_error_l2[0] +=  boundary_error_l2_on_Component;
  }
 }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** Calculate L2-errors to a given function at the boundary (without the global database)*/
void TFEFunction2D::GetL2BoundaryError(BoundValueFunct2D *Exact,
                                       TAuxParam2D *, int,
                                       const TFESpace2D **fespaces,
                                       double *final_boundary_error_l2,
                                       int boundary_component_id,
                                       bool rescale_by_h_E)
{
  int *GlobalNumbers = FESpace2D->GetGlobalNumbers();
  int *BeginIndex = FESpace2D->GetBeginIndex();
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // ########################################################################
  // loop over all Nitsche edges
  // ########################################################################
  TCollection *Coll = fespaces[0]->GetCollection();  // all spaces use the same Coll

  final_boundary_error_l2[0] = 0;

    // Create a list of those boundary edges that are on the boundary component with given ID
    std::vector<TBoundEdge*> boundaryEdgeList;
    Coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

    double boundary_error_l2_on_Component = 0;

    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
      TBoundEdge *boundedge = boundaryEdgeList[m];
      TBaseCell *cell = boundedge->GetNeighbour(0);
      // get basis dimension and FE space data of the current cell
      // Here it is assumed that only one local FE space is used for the hole mesh!
      FE2D CurrentElement_FEID = FESpace2D->GetFE2D(0, cell);

      BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
      // calculate all needed values and derivatives of the basis functions
      BaseFunct2D BaseFunct = BaseFuncts[CurrentElement_FEID];
      if(TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim() > 1)
      {
        ErrMsg("The method "
            << "TFEFunction2D::GetL2BoundaryErrors"
            << "is not designed for vector valued functions");
        OutPut("No errors were computed\n");
        return;
      }

      int BaseVectDim = 1;
      int joint_id = boundedge->get_index_in_neighbour(cell);
      // get a quadrature formula good enough for the argument of the integral
      int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement_FEID);
      QuadFormula1D LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree((TDatabase::ParamDB->INPUT_QUAD_RULE < 2*fe_degree)? 2*fe_degree : TDatabase::ParamDB->INPUT_QUAD_RULE);
      std::vector<double> quadWeights, quadPoints;
      BoundaryAssembling2D::get_quadrature_formula_data(quadPoints, quadWeights, LineQuadFormula);
      // compute values of all basis functions and their first partial derivatives at all quadrature points
      std::vector< std::vector<double> > uorig, u_dx_orig ,u_dy_orig;
      BoundaryAssembling2D::get_original_values(CurrentElement_FEID, joint_id, cell, quadPoints, BaseVectDim, uorig, u_dx_orig, u_dy_orig, LineQuadFormula);

      int N_ = N_BaseFunct[CurrentElement_FEID];

      double FEFunctValues[MaxN_BaseFunctions2D];
      int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()];
      for(int l = 0; l < N_; l++)
      {
        FEFunctValues[l] = Values[DOF[l]];
      }

      double summing_boundary_error_l2_on_edge = 0;
      double edge_length = boundedge->get_length();
      double reference_edge_length = 2; // [-1,1] is the reference edge here

      for(size_t j = 0; j < quadPoints.size(); j++)
      {
        double value = 0;
        for(int l = 0; l < N_; l++)
        {
          value += FEFunctValues[l] * uorig[j][l]; // compute u_h|T = \sum \alpha_i \Phi_i
        }
        auto comp = boundedge->GetBoundComp();
        int comp_ID = comp->GetID();
        double t0, t1;
        boundedge->GetParameters(t0, t1);
        double t = t0 + 0.5 * (t1-t0) * (quadPoints[j]+1);
        double ExactVal, Diff;
        Exact(comp_ID, t, ExactVal);
        Diff = ExactVal - value;

        if(rescale_by_h_E)
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * 1/edge_length * Diff * Diff;
        }
        else
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * Diff * Diff;
        }
      }
      boundary_error_l2_on_Component += summing_boundary_error_l2_on_edge;
    }
    final_boundary_error_l2[0] +=  boundary_error_l2_on_Component;
 }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** calculate errors to given function */
/** parallel with MPI, Sashi : 10.10.09 */
void TFEFunction2D::GetErrors(DoubleFunct2D *Exact, int N_Derivatives,
    MultiIndex2D *NeededDerivatives,
    int N_Errors, ErrorMethod *ErrorMeth, 
    CoeffFct2D Coeff, 
    TAuxParam2D *Aux,
    int n_fespaces, const TFESpace2D **fespaces,
    double *errors, bool,
    std::function<bool(const TBaseCell*, int)>funct) const
{
  bool SecondDer[n_fespaces];
  // initialize the array SecondDer as either all true or all false.
  {
    // find out if one of the indices in NeededDerivatives refers to a second
    // order derivative
    bool need_second_derivatives = false;
    for(int i = 0; i < N_Derivatives; ++i)
    {
      if(NeededDerivatives[i] == D20 || NeededDerivatives[i] == D11 ||
          NeededDerivatives[i] == D02)
        need_second_derivatives = true;
    }
    for(int i = 0; i < n_fespaces; i++)
      SecondDer[i] = need_second_derivatives;
  }

  double *aux;
  double *Param[MaxN_QuadPoints_2D];
  int N_Parameters = Aux->GetN_Parameters();
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(int j = 0; j < MaxN_QuadPoints_2D; j++)
    Param[j] = aux + j*N_Parameters;

  double *Derivatives[MaxN_QuadPoints_2D];
  aux = new double [MaxN_QuadPoints_2D*N_Derivatives];
  for(int j = 0; j < MaxN_QuadPoints_2D; j++)
    Derivatives[j] = aux + j*N_Derivatives;

  double *ExactVal[MaxN_QuadPoints_2D];
  aux = new double [MaxN_QuadPoints_2D * 4];
  for(int j = 0; j < MaxN_QuadPoints_2D; j++)
    ExactVal[j] = aux + j*4;

  // 20 <= number of term
  double *AuxArray[MaxN_QuadPoints_2D];
  aux = new double [MaxN_QuadPoints_2D*20];
  for(int j = 0; j < MaxN_QuadPoints_2D; j++)
    AuxArray[j] = aux + j*20;

  BaseFunct2D *BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  int *N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // for L_infty-error there is one extra entry
  memset(errors, 0, (N_Errors + 1) * SizeOfDouble);

  // ########################################################################
  // loop over all cells
  // ########################################################################
  TCollection *Coll = fespaces[0]->GetCollection();  // all spaces use same Coll
  int N_Cells = Coll->GetN_Cells();
  double x_coord_max_error, y_coord_max_error;

  for(int i = 0; i < N_Cells; i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    if ( funct(cell,-4711) )
      continue;

#ifdef _MPI
    int ID, rank;
    MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
    ID  = cell->GetSubDomainNo();

    if(rank!=ID) // halo cells errors will not be calculated
    {
      continue; 
    }
#endif

    // get basis dimension and FE space data of cell i
    FE2D CurrentElement_FEID = FESpace2D->GetFE2D(i, cell);
    //CurrentElement_FEID is an identifier to some local FE Space (i.e., basis, functionals, polynomial degree,...)
    BaseFunct2D BaseFunct = BaseFuncts[CurrentElement_FEID];
    if(TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetBaseVectDim() > 1)
    {
      ErrMsg("For vector valued basis functions, you should use "
          << "TFEFunction2D::GetErrorsForVectorValuedFunction instead of "
          << "TFEFunction2D::GetErrors");
      OutPut("No errors were computed\n");
      return;
    }

    RefTrans2D RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(CurrentElement_FEID);
    BF2DRefElements ref_element = TFEDatabase2D::GetRefElementFromFE2D(
        CurrentElement_FEID);
    // get quadrature formula id, be carefull: TDatabase::ParamDB->INPUT_QUAD_RULE is initialized with 0. It might be the case that the resulting Quadrature rule is not accurate enough for the exact solution

    //QuadFormula2D qf_id = TFEDatabase2D::GetQFFromDegree(
    //    TDatabase::ParamDB->INPUT_QUAD_RULE, ref_element);
    int fe_degree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement_FEID);
    QuadFormula2D qf_id = TFEDatabase2D::GetQFFromDegree((TDatabase::ParamDB->INPUT_QUAD_RULE < 2*fe_degree)? 2*fe_degree : TDatabase::ParamDB->INPUT_QUAD_RULE, ref_element);

    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);

    TQuadFormula2D *quad_formula = TFEDatabase2D::GetQuadFormula2D(qf_id);
    int N_QuadPoints;
    const double *weights, *xi, *eta;
    quad_formula->GetFormulaData(N_QuadPoints, weights, xi, eta);

    // get quadrature coordinates on original cell (also AbsDetjk is filled)
    double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];
    // Compute quadrature points on original element (X, Y) according to RefTrans
    TFEDatabase2D::GetOrigFromRef(RefTrans, N_QuadPoints, xi, eta, X, Y, AbsDetjk);
    double hK = cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    int Used[N_FEs2D];
    memset(Used, 0, N_FEs2D*SizeOfInt);
    for(int j = 0; j < n_fespaces; j++)
    {
      CurrentElement_FEID = fespaces[j]->GetFE2D(i, cell);
      Used[CurrentElement_FEID] = 1;
    }

    FE2D LocalUsedElements[N_FEs2D];
    int N_LocalUsedElements = 0;
    memset(LocalUsedElements, 0, SizeOfInt*N_FEs2D);
    int j = 0;
    for(int k = 0; k < N_FEs2D; k++)
      if(Used[k])
      {
        LocalUsedElements[j] = (FE2D)k;
        j++;
      }
    N_LocalUsedElements = j;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    // Compute transformation of basis functions (and their derivatives) to original cell and the evaluation for all quadrature points 
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
        Coll, cell, SecondDer, N_QuadPoints, xi, eta, weights, X, Y, AbsDetjk);

    if(N_Parameters>0)
      Aux->GetParameters(N_QuadPoints, Coll, cell, i, xi, eta, X, Y, Param);

    // calculate all needed derivatives of this FE function
    CurrentElement_FEID = FESpace2D->GetFE2D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement_FEID];
    int N_ = N_BaseFunct[CurrentElement_FEID];

    double FEFunctValues[MaxN_BaseFunctions2D];
    int *DOF = FESpace2D->GetGlobalDOF(i);
    for(int l = 0; l < N_; l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }

    for(int k = 0; k < N_Derivatives; k++)
    { // create a pointer to the QuadPoint-values of the basis on original element and its derivatives
      double **OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
          NeededDerivatives[k]);
      for(int j = 0; j < N_QuadPoints; j++)
      {
        double *Orig = OrigFEValues[j];
        double value = 0;

        for(int l = 0; l < N_; l++)
        {
          value += FEFunctValues[l] * Orig[l];  // compute u_h|T = \sum \alpha_i \Phi_i
        }                                         // endfor l
        // Derivatives[j][0] contains values of the FEFunction at QuadPoints
        // Derivatives[j][1] contains values of x derivatives of the FEFunction at QuadPoints
        // Derivatives[j][2] contains values of y derivatives of the FEFunction at QuadPoints↲
        Derivatives[j][k] = value;
      }                                           // endfor j
    }                                             // endfor k

    /* Compute L^\inf-errors of the function in all quadrature points*/ 
    // This computation is the reason why errors has to have one 'entry' more
    for(int j = 0; j < N_QuadPoints; j++)
    {
      Exact(X[j], Y[j], ExactVal[j]);
      //OutPut("Pointwise_Error in: " << X[j] << " " << Y[j] << " is " << fabs(*ExactVal[j] - Derivatives[j][0]) << endl);
      if(fabs(*ExactVal[j] - Derivatives[j][0]) > errors[N_Errors]) // *ExactVal[j] is the same as ExactVal[j][0]
      {
        errors[N_Errors] = fabs(*ExactVal[j] - Derivatives[j][0]);
        x_coord_max_error = X[j];
        y_coord_max_error = Y[j];
      }
    }

    if(Coeff)
      Coeff(N_QuadPoints, X, Y, Param, AuxArray);
    // special rule for L1 error and P1 FE
#ifdef __2D__
    if ((ErrorMeth == &L1Error) && (LocalUsedElements[0] == C_P1_2D_T_A))
    {
      // compute coordinates of vertices
      for (int j = 0; j < 3; j++)
      {
        X[j] = cell->GetVertex(j)->GetX();
        Y[j] = cell->GetVertex(j)->GetY();
      }
      // set flag
      hK = -4711;
    }
#endif
    double LocError[N_Errors];
    ErrorMeth(N_QuadPoints, {{X, Y}}, AbsDetjk, weights, hK, Derivatives, 
        ExactVal, AuxArray, LocError);

#ifdef __2D__
    if (!(ErrorMeth == &conv_diff_l2_h1_linf_error<2>))
    {
      for(int j = 0; j < N_Errors; j++)
        errors[j] += LocError[j];
    }
    else
    {
      for(int j = 0; j < N_Errors-1; j++)
        errors[j] += LocError[j];
      // L_infty error
      if(errors[N_Errors-1] <  LocError[N_Errors-1])
        errors[N_Errors-1] = LocError[N_Errors-1];
    }
#endif

  } // end for i, loop over cells


#ifndef _MPI // sqrt(errors[j]) in the main programm after collecting error from all subdomain
#ifdef __2D__
  if (!(ErrorMeth == &L1Error))
  {
    if (!(ErrorMeth == &conv_diff_l2_h1_linf_error<2>))
    {
      for(int j = 0; j < N_Errors; j++)
        errors[j] = sqrt(errors[j]);
    }
    else
    {
      for(int j = 0; j < N_Errors-1; j++)
        errors[j] = sqrt(errors[j]);
    }
  }
#endif
#ifdef __3D__
  for(int j = 0; j < N_Errors; j++)
    errors[j] = sqrt(errors[j]);
#endif
#endif

  bool detailed_information = false;
  if (detailed_information)
    OutPut("Maximal pointwise error appeared in: " <<  x_coord_max_error << " " <<  y_coord_max_error << " value: " << errors[N_Errors]  << endl);

  delete [] AuxArray[0];
  delete [] ExactVal[0];
  delete [] Derivatives[0];
  delete [] Param[0];

}                                                 // TFEFunction2D::GetErrors


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double TFEFunction2D::get_L2_norm_on_boundary(int boundary_component) const
{
  auto collection = FESpace2D->GetCollection();
  auto n_cells = collection->GetN_Cells();
  double l2_norm_on_boundary = 0.;
  for(int i_cell = 0; i_cell < n_cells; ++i_cell)
  {
    auto cell = collection->GetCell(i_cell);
    auto n_joints = cell->GetN_Joints();
    for(int i_joint = 0; i_joint < n_joints; ++i_joint)
    {
      auto joint = cell->GetJoint(i_joint);
      if(joint->GetType() == BoundaryEdge)
      {
        auto *boundary_joint = dynamic_cast<TBoundEdge*>(joint);
        auto component = boundary_joint->GetBoundComp();
        if(boundary_component < 0 || component->GetID() == boundary_component)
        {
          // found a joint which is on the desired boundary
          auto fe = FESpace2D->get_fe(i_cell);
          FE2D fe_id = FESpace2D->GetFE2D(0, cell); // why is this not in 'fe'?
          auto basis_functions = fe.GetBaseFunct2D();
          auto fe_degree = basis_functions->GetPolynomialDegree();
          auto n_local_basis_functions = basis_functions->GetDimension();
          auto LineQuadFormula_id = TFEDatabase2D::GetQFLineFromDegree(fe_degree);
          std::vector<double> quadWeights, quadPoints;
          BoundaryAssembling2D::get_quadrature_formula_data(
            quadPoints, quadWeights, LineQuadFormula_id);
          std::vector<std::vector<double>> uorig, u_dx_orig ,u_dy_orig;
          BoundaryAssembling2D::get_original_values(
            fe_id, i_joint, cell, quadPoints, basis_functions->GetBaseVectDim(), 
            uorig, u_dx_orig, u_dy_orig, LineQuadFormula_id);
          
          auto local_to_global_dof = FESpace2D->GetGlobalDOF(i_cell);
          double FEFunctValues[n_local_basis_functions];
          for(int l = 0; l < n_local_basis_functions; l++)
          {
            FEFunctValues[l] = Values[local_to_global_dof[l]];
          }
          double edge_length = boundary_joint->get_length();
          double reference_edge_length = 2; // [-1,1] is the reference edge
          
          for(auto j = 0u; j < quadPoints.size(); j++)
          {
            double value = 0;
            for(int l = 0; l < n_local_basis_functions; l++)
            { 
              value += FEFunctValues[l] * uorig[j][l]; // compute u_h|T = \sum \alpha_i \Phi_i
            }
            l2_norm_on_boundary += quadWeights[j] * (value * value) 
                                   * edge_length/reference_edge_length;
          }
        }
      }
    }
  }
  return std::sqrt(l2_norm_on_boundary);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double TFEFunction2D::get_L2_norm() const
{
  auto collection = FESpace2D->GetCollection();
  auto n_cells = collection->GetN_Cells();
  double l2_norm = 0.;
  for(int i_cell = 0; i_cell < n_cells; ++i_cell)
  {
    auto cell = collection->GetCell(i_cell);
    auto fe = FESpace2D->get_fe(i_cell);
    FE2D fe_id = FESpace2D->GetFE2D(0, cell); // why is this not in 'fe'?
    auto basis_functions = fe.GetBaseFunct2D();
    auto n_local_basis_functions = basis_functions->GetDimension();
    int N_QuadPoints;
    const double *weights, *xi, *eta;
    // get quadrature coordinates on original cell (also AbsDetjk is filled)
    double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
    double AbsDetjk[MaxN_QuadPoints_2D];
    bool SecondDer = false;

    // Compute transformation of basis functions (and their derivatives) to 
    // original cell and the evaluation for all quadrature points 
    TFEDatabase2D::GetOrig(1, &fe_id, collection, cell, &SecondDer, 
                           N_QuadPoints, xi, eta, weights, X, Y, AbsDetjk);

    // calculate all needed derivatives of this FE function
    double FEFunctValues[n_local_basis_functions];
    int *DOF = FESpace2D->GetGlobalDOF(i_cell);
    for(int l = 0; l < n_local_basis_functions; l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }

    // create a pointer to the QuadPoint-values of the basis on original element and its derivatives
    double **OrigFEValues = TFEDatabase2D::GetOrigElementValues(
        basis_functions->GetID(), D00);
    for(int j = 0; j < N_QuadPoints; j++)
    {
      double *Orig = OrigFEValues[j];
      double value = 0;

      for(int l = 0; l < n_local_basis_functions; l++)
      {
        value += FEFunctValues[l] * Orig[l];
      }
      l2_norm += value * value * weights[j] * AbsDetjk[j];
    }
  }
  return l2_norm;
}



//==========================================================================

/** determine the value of function and its first derivatives at
  the given point */
void TFEFunction2D::FindGradient(double x, double y, double *values) const
{
  int i,j, N_Cells;
  double xi, eta;
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

      delete[] uorig;
      delete[] uxorig;
      delete[] uyorig;
      delete[] uref;
      delete[] uxiref;
      delete[] uetaref;

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
void TFEFunction2D::FindGradientLocal(const TBaseCell *cell, int cell_no, 
                                      double x, double y, double *values) const
{
  int j,k;
  double xi, eta, eps = 1e-20;
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
  int BaseVectDim = bf->GetBaseVectDim();

  // allocate arrays
  uorig = new double[N_BaseFunct*BaseVectDim];
  uxorig = new double[N_BaseFunct*BaseVectDim];
  uyorig = new double[N_BaseFunct*BaseVectDim];

  uref = new double[N_BaseFunct*BaseVectDim];
  uxiref = new double[N_BaseFunct*BaseVectDim];
  uetaref = new double[N_BaseFunct*BaseVectDim];

  // get values and derivatives of basis functions on the
  // reference mesh cell
  bf->GetDerivatives(D00, xi, eta, uref);
  bf->GetDerivatives(D10, xi, eta, uxiref);
  bf->GetDerivatives(D01, xi, eta, uetaref);

  // compute values on the original mesh cell
  TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
      uref, uxiref, uetaref,
      uorig, uxorig, uyorig);

  // for vector valued basis functions (Raviart-Thomas elements), some signs
  // must be changed
  if(BaseVectDim>1)
  {
    for (int i = 0; i < BaseVectDim; i++)
    {
      for(j=0;j<N_BaseFunct;j++)
      {
        int edge = FE_Obj->GetFEDesc2D()->GetJointOfThisDOF(j);
        if(edge!=-1) // edge ==-1 means inner dof
        {
          uorig[j+i*N_BaseFunct] *= cell->GetNormalOrientation(edge);
          uxorig[j+i*N_BaseFunct] *= cell->GetNormalOrientation(edge);
          uyorig[j+i*N_BaseFunct] *= cell->GetNormalOrientation(edge);
        }
      }
    }
  }

  Numbers = GlobalNumbers + BeginIndex[cell_no];
  for(int i = 0; i < BaseVectDim; i++)
  {
    u = 0;
    ux = 0;
    uy = 0;
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

    values[    3*i] = u;
    values[1 + 3*i] = ux;
    values[2 + 3*i] = uy;
  }

  // set memory free
  delete [] uorig;
  delete [] uxorig;
  delete [] uyorig;
  delete [] uref;
  delete [] uxiref;
  delete [] uetaref;
}


/** determine the value of function at the given point which lies within the
  cell *cell. This also works for vector valued basis functions as are used 
  for Raviart-Thomas elements.
 */
void TFEFunction2D::FindValueLocal(const TBaseCell *cell, int cell_no, double x,
                                   double y, double *values) const
{
  int i, j, k;
  double xi, eta;
  TCollection *Coll;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  int N_BaseFunct;
  int edge;
  int BaseVectDim;
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
  BaseVectDim = bf->GetBaseVectDim();

  // allocate arrays
  uorig = new double[N_BaseFunct*BaseVectDim];
  uxorig = new double[N_BaseFunct*BaseVectDim];
  uyorig = new double[N_BaseFunct*BaseVectDim];

  uref = new double[N_BaseFunct*BaseVectDim];
  uxiref = new double[N_BaseFunct*BaseVectDim];
  uetaref = new double[N_BaseFunct*BaseVectDim];

  // get values and derivatives of basis functions on the
  // reference mesh cell
  bf->GetDerivatives(D00, xi, eta, uref);
  bf->GetDerivatives(D10, xi, eta, uxiref);
  bf->GetDerivatives(D01, xi, eta, uetaref);

  // compute values on the original mesh cell
  TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
      uref, uxiref, uetaref,
      uorig, uxorig, uyorig);

  // for vector valued basis functions (Raviart-Thomas elements), some signs
  // must be changed
  if(BaseVectDim>1)
  {
    for (i=0; i<BaseVectDim; i++)
    {
      for(j=0;j<N_BaseFunct;j++)
      {
        edge = FE_Obj->GetFEDesc2D()->GetJointOfThisDOF(j);
        if(edge!=-1) // edge ==-1 means inner dof
          uorig[j+i*N_BaseFunct] *= cell->GetNormalOrientation(edge);
      }
    }
  }

  Numbers = GlobalNumbers + BeginIndex[cell_no];
  for (i=0; i<BaseVectDim; i++)
  {
    u = 0;
    for(j=0;j<N_BaseFunct;j++)
    {
      k = Numbers[j];
      val = Values[k];
      u += uorig[j+i*N_BaseFunct]*val;
    }

    values[i] = u;
  }

  // set memory free
  delete [] uorig;
  delete [] uxorig;
  delete [] uyorig;
  delete [] uref;
  delete [] uxiref;
  delete [] uetaref;
}


/** calculate the interpolation of an exact function */
void TFEFunction2D::Interpolate(DoubleFunct2D *Exact)
{
  int i,j;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  const double *xi, *eta;
  int *DOF;
  TRefTrans2D *rt;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D];
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double FctVal[4];
  int PolynomialDegree, ApproxOrder;
  QuadFormula2D QuadFormula = BaryCenterTria; //to avoid uninit warning
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Edges = 0;
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
        break;
      case QuadBilinear:
        rt = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)rt)->SetCell(cell);
        break;
      case QuadIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(QuadIsoparametric);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        break;
      case TriaAffin:
        rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)rt)->SetCell(cell);
        break;
      case TriaIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(TriaIsoparametric);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        break;
      default:
        cout << "unknown reftrans id: " << RefTrans << endl;
        break;
    }
    TFEDatabase2D::GetOrigFromRef(RefTrans, N_Points, xi, eta,
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
void TFEFunction2D::GetMeshCellParams(DoubleFunct2D* Exact, int N_Derivatives, 
    MultiIndex2D* NeededDerivatives,
    int N_Errors, ErrorMethod* ErrorMeth,
    CoeffFct2D Coeff, TAuxParam2D* Aux,
    int n_fespaces,
    const TFESpace2D** fespaces,
    double* errors, double* parameters)
{
  int i,j,k,l, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_;
  int Used[N_FEs2D], *N_BaseFunct;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  const double *weights;
  const double *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Param[MaxN_QuadPoints_2D], *aux;
  double *Derivatives[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF;
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
    SecondDer[i] = false;

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

    ErrorMeth(N_Points, {{X, Y}}, AbsDetjk, weights, hK, Derivatives,
        ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
    {
      errors[j] += LocError[j];
      parameters[i + j *N_Cells] = LocError[j];
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
  int i,j;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  const double *xi, *eta;
  int *DOF;
  RefTrans2D F_K = TriaAffin; //avoid uninit warning
  TRefTrans2D *rt;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D];
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double FctVal[4];
  int PolynomialDegree, ApproxOrder;
  QuadFormula2D QuadFormula = BaryCenterTria; //to avoid uninit warning
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  int N_Edges = 0;
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
        cout << "unknown reftrans id: " << endl;
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
        break;
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


/*
   Compute errors for one vector valued function. This is needed if you use for 
   example Raviart-Thomas elements and can't compute the error separetely for each
   component. This happens if you want to compute the L2-error of the divergence.

   Output is double *errors with length 3:
   errors[0] -- L2-Error for this function
   errors[1] -- L2-Error of divergence of this function 
   errors[2] -- H1-semi-Norm-Error of this function
 */
void TFEFunction2D::GetErrorsForVectorValuedFunction(
    DoubleFunct2D * const * const Exact, 
    ErrorMethod * const ErrorMeth, 
    double * const errors) const
{
  // set all errors to zero at first
  memset(errors,0,3*sizeof(double));
  TCollection *coll = FESpace2D->GetCollection(); 
  const int n_cells = coll->GetN_Cells();

  for(int i = 0; i < n_cells; i++)
  {
    TBaseCell *cell = coll->GetCell(i);
    FE2D element_id = FESpace2D->GetFE2D(i, cell);
    // number of basis functions
    int n_base_functs = TFEDatabase2D::GetN_BaseFunctFromFE2D(element_id);
    // id of basis function
    BaseFunct2D bf_id = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(element_id);
    TBaseFunct2D* bf = TFEDatabase2D::GetBaseFunct2D(bf_id);
    int * DOF = FESpace2D->GetGlobalDOF(i);
    RefTrans2D ref_trans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(element_id);
    TFEDatabase2D::SetCellForRefTrans(cell, ref_trans);
    BF2DRefElements ref_element = TFEDatabase2D::GetRefElementFromFE2D(
        element_id);
    // get quadrature formula id
    QuadFormula2D quad_formula_id = 
      TFEDatabase2D::GetQFFromDegree(TDatabase::ParamDB->INPUT_QUAD_RULE,
          ref_element);

    TQuadFormula2D *qf = TFEDatabase2D::GetQuadFormula2D(quad_formula_id);
    int n_points; // number of quadrature points in one cell
    const double *xi,*eta; // coordinates of quadrature points in reference cell
    const double *weights; // quadrature weights
    qf->GetFormulaData(n_points, weights, xi, eta);
    // modulus of determinant of reference transformation for each quadrature
    // point
    double AbsDetjk[n_points];
    // coordinates of quadrature points in original cell
    double X[n_points];
    double Y[n_points];
    // get quadrature coordinates on original cell (also AbsDetjk is filled)
    TFEDatabase2D::GetOrigFromRef(ref_trans, n_points, xi,eta, X, Y, AbsDetjk);

    // will store values of exact solution at all quadrature points in one cell
    double *ExactVal[n_points];
    for(int j = 0; j < n_points; j++)
    {
      ExactVal[j] = new double[8]; // 4 values for each of the two components 
      Exact[0](X[j], Y[j], ExactVal[j]); // x-component
      Exact[1](X[j], Y[j], ExactVal[j]+4); // y-component
    }
    // will store values of this FEFunction and its first derivatives at all 
    // quadrature points
    // the six means values and 2 first derivatives for both components
    double *Derivatives[n_points];
    for(int j = 0; j < n_points; j++)
    {
      Derivatives[j] = new double[6];  
    }
    // Get the function values of this FE-function at the local dofs.
    // some local dofs get a negative sign according to global orientation of 
    // the normals
    double FEFunctValues[n_base_functs];
    for(int l = 0; l < n_base_functs; l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
      // revert sign at some inner edges
      // edge number where a local DOF is on
      int edge = TFEDatabase2D::GetFE2D(element_id)
        ->GetFEDesc2D()->GetJointOfThisDOF(l);
      if(edge != -1) // edge==-1 means inner dof
      {
        FEFunctValues[l] *= cell->GetNormalOrientation(edge);
      }
    }

    double uref[n_base_functs*2]; // 2 because vector valued basis functions
    double uxiref[n_base_functs*2]; // 2 because vector valued basis functions
    double uetaref[n_base_functs*2]; // 2 because vector valued basis functions
    double AllOrigValues[3][n_points][2*n_base_functs];
    for(int k = 0; k < n_points; k++)
    {
      bf->GetDerivatives(D00, xi[k], eta[k], uref);
      bf->GetDerivatives(D10, xi[k], eta[k], uxiref);
      bf->GetDerivatives(D01, xi[k], eta[k], uetaref);
      // Piola transform
      TFEDatabase2D::GetOrigValues(ref_trans, xi[k], eta[k], bf, coll, 
          (TGridCell *) cell, uref, uxiref, uetaref,
          AllOrigValues[0][k], AllOrigValues[1][k],
          AllOrigValues[2][k]);
    }

    // loop over all needed derivatives
    for(int k = 0; k < 3; k++)
    {
      // loop over all quadrature points
      for(int j = 0; j < n_points; j++) 
      {
        double value_x = 0;
        double value_y = 0;
        // loop over all basis functions
        for(int l = 0; l < n_base_functs; l++)
        {
          value_x += FEFunctValues[l] * AllOrigValues[k][j][l];
          value_y += FEFunctValues[l] * AllOrigValues[k][j][l+n_base_functs];
        }
        Derivatives[j][k] = value_x;
        Derivatives[j][k+3] = value_y;
      }
    }
    double hK=1; // we set it one here since it is not needed in 
    // ErrorMeth=L2H1Errors
    double LocError[3]; // L^2 error in value, divergence and first derivative
    ErrorMeth(n_points, {{X, Y}}, AbsDetjk, weights, hK, Derivatives, 
        ExactVal, nullptr, LocError);
    for(int j = 0; j < 3; j++) 
    {
      errors[j] += LocError[j];
    }
    // delete everything which was created with "new" within this loop
    // otherwise one would get (many) memory leaks
    for (int j = 0; j < n_points; j++)
    {
      delete [] ExactVal[j]; ExactVal[j] = nullptr;
      delete [] Derivatives[j]; Derivatives[j] = nullptr;
    }
  } // end loop over all cells

  for(int j = 0; j < 3; j++)  
    errors[j] = sqrt(errors[j]);
}


/** write the solution into a data file - written by Sashi **/
void TFEFunction2D::WriteSol(std::string directory, std::string basename)
{
  int i, N_Joints, N_Cells;
  static int img=0;
  char Dquot;

  TCollection *Coll;
  TBaseCell *cell;

  Dquot = 34; //  see ASCII Chart
  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  i=0;
  cell =  Coll->GetCell(i);
  N_Joints = cell->GetN_Joints();
  const char* BaseName = basename.c_str();
  const char* output_directory = directory.c_str();

  std::ostringstream os;
  os << " ";

  OutPut("Writing solution into "<< output_directory << "/" << BaseName << img
      << ".Sol MooNMD file"<< endl);

  os.seekp(std::ios::beg);
  os << output_directory << "/" << BaseName << img << ".Sol" << ends;
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
void TFEFunction2D::ReadSol(std::string BaseName)
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
  int i,j, N_Cells, N_Edges = 0;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, *DOF;
  int *IntIndex;

  const double *xi, *eta;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D];
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double values[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TNodalFunctional2D *nf;
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

    RefElement = Element->GetBaseFunct2D()->GetRefElement();
    switch(RefElement)
    {
      case BFUnitSquare:
        N_Edges = 4;
        break;

      case BFUnitTriangle:
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

    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
    TFEDatabase2D::GetOrigFromRef(RefTrans, N_Points, xi, eta, X, Y, AbsDetjk);

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

/** ************************************************************************* */
void TFEFunction2D::add(AnalyticFunction f)
{
  int N_Points;
  const double *xi, *eta;
  // begin code

  TCollection *Coll = FESpace2D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  int N_DOFs = FESpace2D->GetN_DegreesOfFreedom();
  std::vector<int> IntIndex(N_DOFs, 0);

  for(int i = 0; i < N_Cells; i++)
  {
    const TBaseCell * cell = Coll->GetCell(i);
    const TFE2D& Element = FESpace2D->get_fe(i);
    auto nf = Element.GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);
    int N_LocalDOFs = Element.GetN_DOF();
    
    if(Element.GetBaseFunct2D()->GetBaseVectDim() != 1)
      ErrThrow("TFEFunction2D::add not implemented for vector-valued basis "
               "functions");

    BF2DRefElements RefElement = Element.GetBaseFunct2D()->GetRefElement();
    int N_Edges = (RefElement == BFUnitSquare) ? 4 : 3; // else BFUnitTriangle
    RefTrans2D RefTrans = Element.GetRefTransID();
    bool IsIsoparametric = false;
    if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(int j = 0; j < N_Edges; j++)
      {
        const TJoint * joint = cell->GetJoint(j);
        JointType jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          auto bdjoint = static_cast<const TBoundEdge *>(joint);
          BoundTypes bdtype = bdjoint->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == InterfaceJoint)
        {
          auto ifjoint = static_cast<const TInterfaceJoint*>(joint);
          BoundTypes bdtype = ifjoint->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
          IsIsoparametric = true;
      }
    }

    if(IsIsoparametric)
    {
      RefTrans = (RefElement == BFUnitSquare) ? QuadIsoparametric 
                                              : TriaIsoparametric;
    }

    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
    std::vector<double> X(N_Points, 0.);
    std::vector<double> Y(N_Points, 0.);
    std::vector<double> AbsDetjk(N_Points, 0.);
    TFEDatabase2D::GetOrigFromRef(RefTrans, N_Points, xi, eta, X.data(),
                                  Y.data(), AbsDetjk.data());
    std::vector<double> PointValues(N_Points, 0.);

    for(int j = 0; j < N_Points; j++)
    {
      PointValues[j] = f(cell, i, {{X[j], Y[j]}});
    }
    std::vector<double> FunctionalValues(N_LocalDOFs, 0.);
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues.data(),
        FunctionalValues.data());

    int * DOF = FESpace2D->GetGlobalDOF(i);

    for(int j = 0; j < N_LocalDOFs; j++)
    {
      if(IntIndex[DOF[j]] == 0)
        Values[DOF[j]] += FunctionalValues[j];
      IntIndex[DOF[j]] ++;
    }
  }

  for(int i = 0; i < N_DOFs; i++)
  {
    if(IntIndex[i] == 0)
    {     
      ErrThrow("TFEFunction2D::add: unable to set dof ", i);
    }
  }
}

/** ************************************************************************* */
void TFEFunction2D::compute_integral_and_measure(double& integral, 
    double& measure) const
{
  TCollection *coll = FESpace2D->GetCollection(); 

  integral = 0.0; // variable to store integral value of this TFEFunction2D
  measure = 0.0; // variable to store the measure of the domain

  // loop over all cells, find out integral value of this FEFunction2D and the`
  // measure of its domain
  const int n_cells = coll->GetN_Cells();
  for(int i = 0; i < n_cells; i++)
  {
    TBaseCell *cell = coll->GetCell(i); // current cell
    FE2D feID = FESpace2D->GetFE2D(i, cell); // id of finite element

    // calculate values on original element (i.e. prepare reference 
    // transformation)
    bool SecondDer = false;
    // quadrature weights and points in reference cell
    const double *weights, *xi, *eta;
    // ugly, we need to change GetOrig!!
    double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D]; // quadrature points
    double AbsDetjk[MaxN_QuadPoints_2D]; // determinant of transformation
    int n_points = 0;
    TFEDatabase2D::GetOrig(1, &feID, coll, cell, &SecondDer,
        n_points, xi, eta, weights, X, Y, AbsDetjk);

    // finite element on the current cell
    TFE2D *fe = TFEDatabase2D::GetFE2D(feID);
    const int n_loc_dof = fe->GetN_DOF(); // number of local dofs
    int * DOF = FESpace2D->GetGlobalDOF(i);

    // id of the local basis functions
    BaseFunct2D base_fc_id = fe->GetBaseFunct2D()->GetID();
    // transformed values of basis functions
    double **orig_values = TFEDatabase2D::GetOrigElementValues(base_fc_id, D00);
    // local integration (loop over all quadrature points)
    for(int j = 0; j < n_points; j++)
    {
      // local transformed values on this quadrature point
      double * orig = orig_values[j]; 
      double value = 0; // value of this TFEFunction2D at this quadrature point
      for(int l = 0; l < n_loc_dof; l++)
      {
        // entry in the vector of this TFEFunction2D times basis function
        value += Values[DOF[l]] * orig[l];
      } // endfor l

      const double w = weights[j] * AbsDetjk[j];
      integral += w * value;
      measure += w;
    } // endfor j
  }
}

/** ************************************************************************* */
double TFEFunction2D::compute_mean() const
{
  double integral, measure;
  this->compute_integral_and_measure(integral, measure);

  return integral/measure;
}

/** ************************************************************************* */
void TFEFunction2D::project_into_L20(double a)
{
  // compute current integral and measure of the domain:
  double integral, measure;
  this->compute_integral_and_measure(integral, measure);
  double new_mean = (integral - a)/measure;

  //OutPut("TFEFunction2D::project_into_L20 computed mean value before " << 
  //       integral/measure << endl);

  // vector of the same length as this TFEFunction2D. It represents a function
  // which has the constant value 'mean' for all nodal functionals. The last 
  // step in this projection will be to substract this vector from the vector of
  // this TFEFunction2D
  // for standard P_k or Q_k finite elements this is a constant function
  double *interpol = new double[Length];

  TCollection *coll = FESpace2D->GetCollection();
  const int n_cells = coll->GetN_Cells();
  for(int i = 0; i < n_cells; i++)
  {
    TBaseCell *cell = coll->GetCell(i); // current cell
    FE2D feID = FESpace2D->GetFE2D(i, cell); // id of finite element
    // finite element on the current cell
    TFE2D *fe = TFEDatabase2D::GetFE2D(feID);
    const int n_loc_dof = fe->GetN_DOF(); // number of local dofs
    int * DOF = FESpace2D->GetGlobalDOF(i);
    TNodalFunctional2D *nf = TFEDatabase2D::GetNodalFunctional2DFromFE2D(feID);
    int n_points; // number of evaluation points to compute nodal functionals
    const double *xi, *eta; //coordinates of evaluation points in reference cell
    nf->GetPointsForAll(n_points, xi, eta);
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
  // substract it from the vector of this TFEFunction2D
  for(int i = 0; i < Length; i++)
    Values[i] -= interpol[i];

  delete [] interpol;
}


TFEFunction2D& TFEFunction2D::operator*=(double alpha)
{
  int N_Active = FESpace2D->GetActiveBound();
  for (int i=0; i<N_Active; i++)
  {
    Values[i] *= alpha;
  }
  return *this;
}

TFEFunction2D & TFEFunction2D::operator+=(const TFEFunction2D & rhs)
{
  if(FESpace2D != rhs.FESpace2D)
  {
    OutPut("ERROR: TFEFunction2D::operator+=() The two arguments "
        << "have different fe spaces. Exiting" << endl);
    exit(1);
  }
  if(Length != rhs.Length)
  {
    OutPut("ERROR: TFEFunction2D::operator+=() The two arguments "
        << "have different lengths. Exiting" << endl);
    exit(1);
  }
  if(Values == rhs.Values)
  {
    OutPut("ERROR: TFEFunction2D::operator+=() The two arguments "
        << "have the same solution vector. This operation would be "
        << "equivalent to a multiplication by 2! Exiting" << endl);
    exit(1);
  }
  int N_Active = FESpace2D->GetActiveBound();
  for (int i=0; i<N_Active; i++)
  {
    Values[i] += rhs.Values[i];
  }
  return *this;
}

TFEFunction2D & TFEFunction2D::operator=(const TFEFunction2D & other)
{
  this->Name        = other.Name;
  this->Description = other.Description;
  this->FESpace2D   = other.FESpace2D;
  this->Values      = other.Values;
  this->Length      = other.Length;

  return *this;
}

void TFEFunction2D::MinMax(double & min, double & max) const
{
  double val;
  max = -1e100, min = 1e100;

  for(int i = 0; i<Length; i++)
  {
    val = Values[i];
    if(val>max) max = val;
    if(val<min) min = val;
  }
}

void TFEFunction2D::PrintMinMax(std::string name) const
{
  double min, max;
  this->MinMax(min, max);
  if( min <= max )
  {
    if(name.empty())
      Output::print<1>(this->Name, " min ", min, ", max ", max);
    else
      Output::print<1>(name, " min ", min, ", max ", max);
  }
  else
  {
    Output::print<1>("WARNING: TFEFunction2D::MinMax was not successful!");
  }
}


/** set Dirichlet values according to boundary conditions */
void TFEFunction2D::SetDirichletBC(BoundCondFunct2D *BoundaryCondition,
    BoundValueFunct2D *BoundaryValue)
{
  TCollection *Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc_Obj;
  TJoint *joint;
  TBoundEdge *boundedge;
  TIsoBoundEdge *isoboundedge;
  TNodalFunctional2D *nf;
  const TBoundComp2D* BoundComp;
  BoundCond Cond0, Cond1;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D];
  double t0, t1, s, t;
  int comp;
  int i, l, m;
  int N_Cells, N_Joints, N_EdgeDOF;
  int *BeginIndex, *GlobalNumbers;
  int *DOF;
  int *EdgeDOF, N_EdgePoints;
  const double *EdgePoints;
  double eps=1e-4;

  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    FEDesc_Obj = Element->GetFEDesc2D();

    nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

    DOF = GlobalNumbers+BeginIndex[i];

    N_Joints = cell->GetN_Edges();
    for(m=0;m<N_Joints;m++)
    {
      joint = cell->GetJoint(m);
      if(joint->GetType() == BoundaryEdge ||
          joint->GetType() == IsoBoundEdge)
      {
        if(joint->GetType() == BoundaryEdge)
        {
          boundedge = (TBoundEdge *)joint;
          BoundComp = boundedge->GetBoundComp();
          boundedge->GetParameters(t0, t1);
        }
        else
        {
          isoboundedge = (TIsoBoundEdge *)joint;
          BoundComp = isoboundedge->GetBoundComp();
          isoboundedge->GetParameters(t0, t1);
        }
        // get id of the boundary component
        comp=BoundComp->GetID();
        // get type of the boundary condition at the beginning
        // and at the end of the current edge
        if (t0 < t1)
        {
          BoundaryCondition(comp, t0+eps, Cond0);
          BoundaryCondition(comp, t1-eps, Cond1);
        }
        else
        {
          BoundaryCondition(comp, t0-eps, Cond0);
          BoundaryCondition(comp, t1+eps, Cond1);
        }
        // only one boundary condition per edge allowed
        if(Cond0 == Cond1)
        {
          if(Cond0 == DIRICHLET)
          {
            // read boundary values for each quadrature point
            for(l=0;l<N_EdgePoints;l++)
            {
              s = EdgePoints[l];
              t = 0.5*(t0*(1-s) + t1*(1+s));
              BoundaryValue(comp, t, PointValues[l]);
            } // endfor l
            // compute boundary values for each dof on the 
            // boundary edge with the nodal functionals
            nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
                FunctionalValues);
            EdgeDOF = FEDesc_Obj->GetJointDOF(m);
            N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
            // save boundary values of each dof on the boundary
            // edge in the rhs
            for(l=0;l<N_EdgeDOF;l++)
            {
              Values[DOF[EdgeDOF[l]]] = FunctionalValues[l];
              // cout << i <<  setw(25) << Values[DOF[EdgeDOF[l]]]<< endl;
            }
          } // endif Cond0==DIRICHLET
        } // endif (Cond0==Cond1)
        else
        {
          OutPut("different boundary condition on one edge ");
          OutPut("are not allowed!" << endl);
          exit(4711);
        }
      } // endif (boundary joint)
    } // endfor m<N_Joints
  } // endfor i<N_Cells
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
  int N_Cells, N_Joints;
  int *BeginIndex, *GlobalNumbers, *DOF, N_BF;

  double Mult, r_axial, val, mass, volume;
  const double *weights, *xi, *eta;
  double values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double AbsDetjk[MaxN_QuadPoints_2D], X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];

  TJoint *joint;
  TBaseCell *cell;
  TCollection *coll;
  JointType jointtype;
  BoundTypes bdtype;
  RefTrans2D RefTrans;
  bool IsIsoparametric;
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
    IsIsoparametric = false;
    if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Joints;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)  IsIsoparametric = true;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
          IsIsoparametric = true;

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


